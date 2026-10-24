// SPDX-License-Identifier: MIT
#include "dsp/svf.hpp"
#include <array>
#include <complex>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

using namespace shores_dsp;
constexpr long double pi = 3.1415926535897932384626433832795L;
constexpr std::array<SvfType, 4> modes{
    SvfType::LowPass, SvfType::HighPass, SvfType::BandPass, SvfType::Notch};

static void require(bool condition, const char *why) {
    if (!condition) throw std::runtime_error(why);
}
static void near(double actual, double expected, double tolerance, const char *why) {
    require(std::isfinite(actual) && std::isfinite(expected) &&
            std::abs(actual - expected) <= tolerance, why);
}
static std::vector<float> signal(std::size_t n) {
    std::vector<float> x(n);
    std::uint32_t rng = 0x104f29u;
    for (std::size_t i = 0; i < n; ++i) {
        rng = 1664525u * rng + 1013904223u;
        x[i] = float(0.1 * (double(rng) / 4294967296.0 - 0.5) +
                     0.1 * std::sin(double(i) * 0.173) + 0.025);
    }
    x[0] += 0.5f;
    return x;
}

// Independent static reference: W3C/RBJ sine/cosine coefficients and long-double
// direct form I, not the production tan-based two-integrator realization.
struct ReferenceDF1 {
    long double b0{}, b1{}, b2{}, a1{}, a2{};
    long double x1{}, x2{}, y1{}, y2{};
    ReferenceDF1(double rate, SvfParameters p) {
        const auto hz = std::clamp(p.frequency_hz, 10.0, 0.49 * rate);
        const auto q = std::clamp(p.q, 0.1, 18.0);
        const auto w = 2 * pi * hz / rate;
        const auto c = std::cos(w), alpha = std::sin(w) / (2 * q);
        const auto a0 = 1 + alpha;
        a1 = -2 * c / a0; a2 = (1 - alpha) / a0;
        switch (p.type) {
        case SvfType::LowPass: b0 = (1 - c) / 2; b1 = 1 - c; b2 = b0; break;
        case SvfType::HighPass: b0 = (1 + c) / 2; b1 = -(1 + c); b2 = b0; break;
        case SvfType::BandPass: b0 = alpha; b1 = 0; b2 = -alpha; break;
        case SvfType::Notch: b0 = 1; b1 = -2 * c; b2 = 1; break;
        }
        b0 /= a0; b1 /= a0; b2 /= a0;
    }
    long double tick(float x) {
        const auto y = b0 * x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
        x2 = x1; x1 = x; y2 = y1; y1 = y;
        return y;
    }
    std::complex<long double> response(double hz, double rate) const {
        const auto z = std::polar(1.L, -2 * pi * hz / rate);
        return (b0 + b1 * z + b2 * z * z) / (1.L + a1 * z + a2 * z * z);
    }
};

// Independent time-varying reference: solve the original two simultaneous
// nodal equations by Gaussian elimination in long double at EVERY sample.
// It never uses prepared production coefficients or the optimized recurrence.
struct ReferenceNodal {
    long double state1 = 0, state2 = 0;
    long double tick(float x, double rate, SvfParameters p) {
        const long double g = std::tan(pi * std::clamp(p.frequency_hz, 10.0, rate * 0.49) / rate);
        const long double k = 1.L / std::clamp(p.q, 0.1, 18.0);
        const long double a = 1 + g * k, b = g, c = -g, d = 1;
        const long double rhs1 = state1 + g * x, rhs2 = state2;
        const long double eliminate = c / a;
        const long double low = (rhs2 - eliminate * rhs1) / (d - eliminate * b);
        const long double band = (rhs1 - b * low) / a;
        state1 = 2 * band - state1; state2 = 2 * low - state2;
        switch (p.type) {
        case SvfType::LowPass: return low;
        case SvfType::HighPass: return x - k * band - low;
        case SvfType::BandPass: return k * band;
        case SvfType::Notch: return x - k * band;
        }
        throw std::runtime_error("bad test mode");
    }
};

static double static_checks() {
    const auto input = signal(8192);
    std::vector<float> output(input.size());
    double worst = 0;
    for (double rate : {1000., 8000., 44100., 48000., 96000., 192000., 768000.})
        for (auto mode : modes)
            for (double hz : {10., 100., 1000., 0.49 * rate})
                for (double q : {0.1, 0.7071067811865476, 2., 18.}) {
                    const SvfParameters p{mode, hz, q};
                    Svf filter(design_svf(rate, p));
                    ReferenceDF1 reference(rate, p);
                    filter.process(input.data(), output.data(), input.size());
                    for (std::size_t i = 0; i < input.size(); ++i) {
                        const double expected = double(reference.tick(input[i]));
                        const double error = std::abs(output[i] - expected);
                        worst = std::max(worst, error);
                        near(output[i], expected, 3.e-6 * (1 + std::abs(expected)), "static DF1 mismatch");
                    }
                }

    // Measure the actual processed impulse spectrum, rather than evaluating
    // the production coefficients themselves. Complex comparison checks phase
    // as well as magnitude and covers DC, Nyquist and center normalization.
    constexpr double rate = 48000, cutoff = 1000;
    std::vector<float> impulse(65536, 0.f), response(impulse.size());
    impulse[0] = 1.f;
    for (auto mode : modes) for (double q : {0.1, 0.7071067811865476, 18.}) {
        Svf filter(design_svf(rate, {mode, cutoff, q}));
        filter.process(impulse.data(), response.data(), impulse.size());
        ReferenceDF1 reference(rate, {mode, cutoff, q});
        for (double hz : {0., 100., cutoff, 10000., 24000.}) {
            std::complex<long double> measured{};
            const auto step = std::polar(1.L, -2 * pi * hz / rate);
            std::complex<long double> phase{1, 0};
            for (float h : response) { measured += static_cast<long double>(h) * phase; phase *= step; }
            const auto expected = reference.response(hz, rate);
            require(std::abs(measured - expected) < 3.e-6L * (1 + std::abs(expected)),
                    "processed impulse response/phase mismatch");
            if (hz == cutoff && mode == SvfType::BandPass)
                near(double(std::abs(measured)), 1., 3.e-6, "BP must be unity at center for every Q");
            if (hz == cutoff && mode == SvfType::Notch)
                require(std::abs(measured) < 3.e-6, "notch center rejection");
        }
    }
    return worst;
}

static void contracts() {
    for (double bad : {-std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::quiet_NaN()}) {
        for (int field = 0; field < 3; ++field) {
            bool rejected = false;
            try {
                design_svf(field == 0 ? bad : 48000,
                           {SvfType::LowPass, field == 1 ? bad : 1000, field == 2 ? bad : 1});
            } catch (const std::invalid_argument &) { rejected = true; }
            require(rejected, "nonfinite control accepted");
        }
    }
    for (double bad_rate : {0., 999., 768001., 1.e9}) {
        bool rejected = false;
        try { design_svf(bad_rate, {}); } catch (const std::invalid_argument &) { rejected = true; }
        require(rejected, "invalid sample rate accepted");
    }
    bool rejected = false;
    try { design_svf(48000, {static_cast<SvfType>(42), 1000, 1}); }
    catch (const std::invalid_argument &) { rejected = true; }
    require(rejected, "unknown mode accepted");
    const auto same = [](const SvfCoefficients &a, const SvfCoefficients &b) {
        return a.a1 == b.a1 && a.a2 == b.a2 && a.a3 == b.a3 &&
               a.m0 == b.m0 && a.m1 == b.m1 && a.m2 == b.m2;
    };
    require(same(design_svf(48000, {SvfType::BandPass, -1.e20, -1.e20}),
                 design_svf(48000, {SvfType::BandPass, 10, 0.1})), "lower clamp");
    require(same(design_svf(48000, {SvfType::HighPass, 1.e20, 1.e20}),
                 design_svf(48000, {SvfType::HighPass, 23520, 18})), "upper clamp");
    const auto input = signal(8192);
    std::vector<float> output(input.size());
    Svf wire;
    wire.process(input.data(), output.data(), input.size());
    require(input == output, "default wire differs");
    for (auto mode : modes) {
        const auto c = design_svf(48000, {mode, 1000, 18});
        Svf full(c), chunked(c), samplewise(c), in_place(c);
        full.process(input.data(), output.data(), input.size());
        auto inplace = input;
        in_place.process(inplace.data(), inplace.data(), inplace.size());
        require(inplace == output, "fixed in-place mismatch");
        std::vector<float> chunks(input.size()), singles(input.size());
        std::size_t at = 0, iteration = 0;
        constexpr std::array<std::size_t, 6> blocks{1, 7, 63, 256, 17, 1024};
        while (at < input.size()) {
            const auto count = std::min(blocks[iteration++ % blocks.size()], input.size() - at);
            chunked.process(nullptr, nullptr, 0);
            chunked.process_modulated(nullptr, nullptr, nullptr, 0);
            chunked.process(input.data() + at, chunks.data() + at, count);
            at += count;
        }
        for (std::size_t i = 0; i < input.size(); ++i)
            samplewise.process(input.data() + i, singles.data() + i, 1);
        require(chunks == output && singles == output, "variable fixed block mismatch");
        full.reset();
        std::fill(chunks.begin(), chunks.end(), 0.f);
        full.process(chunks.data(), chunks.data(), chunks.size());
        require(std::all_of(chunks.begin(), chunks.end(), [](float x) { return x == 0; }),
                "reset/silence must be exact zero");
        require(same(full.coefficients(), c), "reset changed coefficients");
        full.process(input.data(), chunks.data(), chunks.size());
        require(chunks == output, "reset did not reproduce fresh output");
        // Ring down high-Q extremes, with block cleanup exercised repeatedly.
        for (double rate : {1000., 48000., 192000., 768000.})
            for (double hz : {10., rate * 0.49}) {
                Svf tail(design_svf(rate, {mode, hz, 18}));
                float impulse = 1, first{};
                tail.process(&impulse, &first, 1);
                double energy = 0;
                for (int block = 0; block < 4096; ++block) {
                    std::array<float, 256> zero{}, y{};
                    tail.process(zero.data(), y.data(), y.size());
                    for (float value : y) { require(std::isfinite(value), "nonfinite silence tail"); energy += value * value; }
                }
                require(std::isfinite(energy) && energy < 100, "unbounded static ringdown");
            }
    }
}

static double modulation_checks() {
    constexpr std::size_t n = 48000;
    const auto input = signal(n);
    double worst = 0, peak = 0;
    // 0.5 Hz, 80 Hz, and 6 kHz modulation of BOTH log cutoff and Q, followed
    // by alternating near-limit values, actual sample-rate changes, and mode
    // changes without reset. All coefficients are prepared before processing.
    for (auto mode : modes) for (int scenario = 0; scenario < 5; ++scenario) {
        std::vector<SvfCoefficients> events(n);
        std::vector<SvfParameters> parameters(n);
        std::vector<double> rates(n, 48000);
        std::vector<float> output(n), chunked(n), single(n);
        for (std::size_t i = 0; i < n; ++i) {
            const double mod_rate = scenario == 0 ? 0.5 : scenario == 1 ? 80 : 6000;
            const double phase = 2 * double(pi) * mod_rate * double(i) / 48000;
            auto &p = parameters[i];
            p = {mode, 100. * std::pow(180., 0.5 + 0.5 * std::sin(phase)),
                 0.1 + 17.9 * (0.5 + 0.5 * std::cos(phase * 0.731))};
            if (scenario >= 3) {
                rates[i] = scenario == 4 ? (i < n / 3 ? 48000 : i < 2 * n / 3 ? 96000 : 8000) : 48000;
                p.frequency_hz = (i % 2 == 0) ? 10 : rates[i] * 0.49;
                p.q = (i % 3 == 0) ? 0.1 : 18;
                if (scenario == 4) p.type = modes[(i / 127) % modes.size()];
            }
            events[i] = design_svf(rates[i], p);
        }
        Svf full, chunks, per_sample, inplace;
        full.process_modulated(input.data(), output.data(), events.data(), n);
        auto in_place = input;
        inplace.process_modulated(in_place.data(), in_place.data(), events.data(), n);
        require(in_place == output, "modulated in-place mismatch");
        std::size_t at = 0, iteration = 0;
        constexpr std::array<std::size_t, 6> blocks{1, 17, 256, 7, 1024, 63};
        while (at < n) {
            const auto count = std::min(blocks[iteration++ % blocks.size()], n - at);
            chunks.process_modulated(nullptr, nullptr, nullptr, 0);
            chunks.process_modulated(input.data() + at, chunked.data() + at, events.data() + at, count);
            at += count;
        }
        for (std::size_t i = 0; i < n; ++i) {
            per_sample.set_coefficients(events[i]);
            per_sample.process(input.data() + i, single.data() + i, 1);
        }
        require(output == chunked && output == single, "modulated event timing or block mismatch");
        ReferenceNodal reference;
        for (std::size_t i = 0; i < n; ++i) {
            const double expected = double(reference.tick(input[i], rates[i], parameters[i]));
            worst = std::max(worst, std::abs(output[i] - expected));
            peak = std::max(peak, std::abs(double(output[i])));
            near(output[i], expected, 3.e-6 * (1 + std::abs(expected)), "time-varying nodal mismatch");
        }
        // process() after modulated processing must retain the last event.
        float next{}, next_expected{};
        full.process(input.data(), &next, 1);
        per_sample.process(input.data(), &next_expected, 1);
        near(next, next_expected, 1.e-7, "last modulation event not retained");
        // Silence under aggressive controls must not create state from zero.
        std::vector<float> silence(n, 0.f);
        full.reset();
        full.process_modulated(silence.data(), silence.data(), events.data(), n);
        require(std::all_of(silence.begin(), silence.end(), [](float x) { return x == 0; }),
                "modulated silence generated signal");
    }
    require(peak < 100, "stress modulation exceeded test amplitude bound");
    std::cout << "modulation peak (test input only): " << peak << '\n';
    return worst;
}

int main() {
    try {
        contracts();
        const double static_error = static_checks();
        const double modulation_error = modulation_checks();
        std::cout << "PASS: four modes; independent long-double DF1/complex response and nodal modulation; "
                     "limits, Q, rates, silence/tails, reset, variable blocks, in-place, event timing; "
                     "max static abs error=" << static_error << "; max modulation abs error=" << modulation_error << '\n';
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
