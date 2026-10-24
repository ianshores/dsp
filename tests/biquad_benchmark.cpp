// SPDX-License-Identifier: MIT
#include "dsp/biquad.hpp"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

int main() {
    constexpr std::size_t sample_rate = 48000, frames = sample_rate * 20, block = 256;
    const auto coefficients = shores_dsp::design_biquad(
        sample_rate, {shores_dsp::BiquadType::Peak, 2400, 0.8, 3, 1});
    const auto alternate = shores_dsp::design_biquad(
        sample_rate, {shores_dsp::BiquadType::Peak, 2400, 0.8, -3, 1});
    std::vector<float> input(frames), work(frames);
    std::uint32_t rng = 0x835f1234u;
    for (auto &x : input) {
        rng = 1664525u * rng + 1013904223u;
        x = float(double(rng) / 4294967296.0 - 0.5);
    }
    std::puts("mode,stages,median_ns_per_sample,percent_of_one_core_at_48khz,checksum");
    for (bool automated : {false, true}) for (std::size_t stages : {1u, 2u, 8u}) {
        std::vector<shores_dsp::Biquad> filters(stages,
                                                shores_dsp::Biquad(coefficients));
        std::vector<double> results;
        double checksum = 0;
        for (int trial = 0; trial < 9; ++trial) {
            work = input;
            for (auto &f : filters) f.reset();
            const auto start = std::chrono::steady_clock::now();
            for (std::size_t offset = 0; offset < frames; offset += block) {
                const auto n = std::min(block, frames - offset);
                for (auto &f : filters) {
                    if (automated)
                        f.ramp_to((offset / block) % 2 == 0 ? alternate : coefficients, n);
                    f.process(work.data() + offset, work.data() + offset, n);
                }
            }
            const auto end = std::chrono::steady_clock::now();
            for (std::size_t i = 0; i < frames; i += 997)
                checksum += work[i];
            const auto elapsed = std::chrono::duration<double, std::nano>(end - start).count();
            results.push_back(elapsed / double(frames));
        }
        std::sort(results.begin(), results.end());
        const double ns = results[results.size() / 2];
        std::printf("%s,%zu,%.3f,%.5f,%.6f\n", automated ? "ramp" : "fixed", stages, ns,
                    ns * double(sample_rate) * 1.e-7, checksum);
    }
}
