// SPDX-License-Identifier: MIT
#include "dsp/biquad.hpp"
#include "dsp/svf.hpp"
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

using namespace shores_dsp;
constexpr std::size_t rate = 48000, frames = rate * 20, block = 256;
constexpr int trials = 9;
constexpr double tau = 6.2831853071795864769;
struct Control { double hz, q; };

// Compiler barrier, no runtime instruction on the tested GCC/Clang builds.
// Needed only in the apply-only loop to prevent dead-store elimination.
template<class T> static void observe(const T &value) {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(&value) : "memory");
#else
    (void)value;
#endif
}
struct SvfEngine {
    using Filter = Svf;
    using Coefficients = SvfCoefficients;
    static constexpr const char *name = "svf";
    static Coefficients design(Control c) { return design_svf(rate, {SvfType::LowPass, c.hz, c.q}); }
    static double checksum(Coefficients c) { return c.a1 + c.a2 + c.a3 + c.m0 + c.m1 + c.m2; }
};
struct BiquadEngine {
    using Filter = Biquad;
    using Coefficients = BiquadCoefficients;
    static constexpr const char *name = "biquad";
    static Coefficients design(Control c) { return design_biquad(rate, {BiquadType::LowPass, c.hz, c.q}); }
    static double checksum(Coefficients c) { return c.b0 + c.b1 + c.b2 + c.a1 + c.a2; }
};
static double elapsed(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start).count();
}
static void row(const char *engine, const char *mode, std::size_t stages, int trial,
                const char *unit, double ns, double checksum, double peak = 0, bool valid = true) {
    valid = valid && std::isfinite(checksum);
    std::printf("%s,%s,%zu,%d,%s,%.6f,%.9f,%.9f,%s\n", engine, mode, stages, trial, unit,
                valid ? ns : std::numeric_limits<double>::quiet_NaN(), checksum, peak,
                valid ? "finite" : "nonfinite");
    std::fflush(stdout);
}

template<class Engine> static void benchmark(const std::vector<float> &input,
                                            const std::vector<Control> &controls, int trial_count) {
    using Filter = typename Engine::Filter;
    using Coefficients = typename Engine::Coefficients;
    const auto initial = Engine::design({1200, 0.8});
    std::vector<Coefficients> events(frames);
    for (std::size_t i = 0; i < frames; ++i) events[i] = Engine::design(controls[i]);
    std::vector<Coefficients> block_events((frames + block - 1) / block);
    for (std::size_t b = 0; b < block_events.size(); ++b) {
        const double t = double(b * block) / rate;
        block_events[b] = Engine::design({200 * std::pow(40., 0.5 + 0.5 * std::sin(tau * 2 * t)),
                                         0.7 + 1.3 * (0.5 + 0.5 * std::cos(tau * 1.462 * t))});
    }
    std::vector<float> work(frames);
    // All event design is excluded from prepared processing modes. Both engines
    // use identical inputs, controls, event timing, stages and block boundaries.
    for (int mode = 0; mode < 4; ++mode) {
        if (mode == 3 && !std::is_same_v<Engine, SvfEngine>) continue;
        for (std::size_t stages : {1u, 8u}) {
            if (mode >= 2 && stages != 1) continue;
            const char *label = mode == 0 ? "fixed" : mode == 1 ? "prepared_block" :
                                mode == 2 ? "prepared_sample_matched" : "prepared_stream_svf";
            for (int trial = 0; trial < trial_count; ++trial) {
                work = input;
                std::vector<Filter> filters(stages, Filter(initial));
                const auto start = std::chrono::steady_clock::now();
                for (std::size_t at = 0; at < frames; at += block) {
                    const auto count = std::min(block, frames - at);
                    for (auto &filter : filters) {
                        if (mode == 1) filter.set_coefficients(block_events[at / block]);
                        if (mode <= 1) {
                            filter.process(work.data() + at, work.data() + at, count);
                        } else if (mode == 2) {
                            for (std::size_t j = at; j < at + count; ++j) {
                                filter.set_coefficients(events[j]);
                                filter.process(work.data() + j, work.data() + j, 1);
                            }
                        } else {
                            if constexpr (std::is_same_v<Engine, SvfEngine>)
                                filter.process_modulated(work.data() + at, work.data() + at,
                                                         events.data() + at, count);
                        }
                    }
                }
                const double ns = elapsed(start) / frames;
                double checksum = 0;
                for (std::size_t i = 0; i < frames; i += 997) checksum += work[i];
                double peak = 0;
                bool valid = true;
                for (float x : work) {
                    valid = valid && std::isfinite(x);
                    if (std::isfinite(x)) peak = std::max(peak, std::abs(double(x)));
                }
                row(Engine::name, label, stages, trial, "sample", ns, checksum, peak, valid);
            }
        }
    }
    // Cost of coefficient design + storing prepared results. Controls and output
    // storage already allocated. No audio processing in these measurements.
    for (int trial = 0; trial < trial_count; ++trial) {
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < frames; ++i) events[i] = Engine::design(controls[i]);
        const double ns = elapsed(start) / frames;
        double checksum = 0;
        for (std::size_t i = 0; i < frames; i += 997) checksum += Engine::checksum(events[i]);
        row(Engine::name, "design_and_store", 1, trial, "event", ns, checksum);
    }
    // Apply-only includes coefficient reads, filter stores, loop, and compiler
    // barriers. Cycle through 4096 prepared sets to bound working-set size.
    for (int trial = 0; trial < trial_count; ++trial) {
        Filter filter(initial);
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < frames; ++i) {
            filter.set_coefficients(events[i % 4096]);
            observe(filter);
        }
        const double ns = elapsed(start) / frames;
        row(Engine::name, "apply_only", 1, trial, "event", ns, Engine::checksum(filter.coefficients()));
    }
}

int main(int argc, char **argv) {
    const bool stress = argc == 2 && std::string(argv[1]) == "--stress";
    if (argc > 1 && !stress) { std::fputs("usage: dsp_svf_benchmark [--stress]\n", stderr); return 1; }
    std::vector<float> input(frames);
    std::vector<Control> controls(frames);
    std::uint32_t rng = 0x835f1234u;
    for (std::size_t i = 0; i < frames; ++i) {
        rng = 1664525u * rng + 1013904223u;
        input[i] = float(0.02 * (double(rng) / 4294967296.0 - 0.5));
        const double phase = tau * 6000 * double(i) / rate;
        // 6 kHz cutoff/resonance controls for per-sample paths. Block paths
        // use separate 2 Hz cutoff/1.462 Hz Q controls prepared above.
        controls[i] = stress ? Control{200 * std::pow(40., 0.5 + 0.5 * std::sin(phase)),
                                      0.7 + 1.3 * (0.5 + 0.5 * std::cos(phase * 0.731))} :
                              Control{800 * std::pow(2., 0.5 + 0.5 * std::sin(phase)),
                                      0.5 + 0.3 * (0.5 + 0.5 * std::cos(phase * 0.731))};
    }
    std::puts("engine,mode,stages,trial,unit,ns_per_unit,checksum,max_abs_output,status");
    benchmark<BiquadEngine>(input, controls, stress ? 1 : trials);
    benchmark<SvfEngine>(input, controls, stress ? 1 : trials);
}
