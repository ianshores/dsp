// SPDX-License-Identifier: MIT
#include "dsp/biquad.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using shores_dsp::Biquad;
using shores_dsp::BiquadCoefficients;
using shores_dsp::BiquadParameters;
using shores_dsp::BiquadType;

static void require(bool condition, const char *reason) {
    if (!condition) throw std::runtime_error(reason);
}
static double gain_db(const BiquadCoefficients &c, double hz, double sample_rate) {
    const double angle = -6.2831853071795864769 * hz / sample_rate;
    const std::complex<double> z(std::cos(angle), std::sin(angle));
    const auto z2 = z * z;
    const auto numerator = double(c.b0) + double(c.b1) * z + double(c.b2) * z2;
    const auto denominator = 1.0 + double(c.a1) * z + double(c.a2) * z2;
    return 20.0 * std::log10(std::abs(numerator / denominator));
}
static BiquadCoefficients design(BiquadType type, double hz = 1000,
                                  double q = 0.7071067811865476, double db = 0,
                                  double slope = 1) {
    return shores_dsp::design_biquad(48000, {type, hz, q, db, slope});
}
static void near(double actual, double expected, double tolerance, const char *reason) {
    if (std::abs(actual - expected) > tolerance) throw std::runtime_error(reason);
}

int main() {
    try {
        const auto lp = design(BiquadType::LowPass);
        near(gain_db(lp, 1000, 48000), -3.0103, 0.02, "low-pass cutoff");
        require(gain_db(lp, 100, 48000) > -0.1, "low-pass bass passband");
        require(gain_db(lp, 10000, 48000) < -40, "low-pass rejection");
        const auto hp = design(BiquadType::HighPass);
        near(gain_db(hp, 1000, 48000), -3.0103, 0.02, "high-pass cutoff");
        require(gain_db(hp, 100, 48000) < -40, "high-pass rejection");
        require(gain_db(hp, 10000, 48000) > -0.1, "high-pass treble passband");
        const auto bp = design(BiquadType::BandPass);
        near(gain_db(bp, 1000, 48000), 0, 0.01, "band-pass center gain");
        require(gain_db(bp, 100, 48000) < -15, "band-pass bass rejection");
        require(gain_db(bp, 10000, 48000) < -15, "band-pass treble rejection");
        require(gain_db(design(BiquadType::Notch), 1000, 48000) < -80,
                "notch center rejection");
        const auto peak = design(BiquadType::Peak, 1000, 1, 12);
        near(gain_db(peak, 1000, 48000), 12, 0.01, "peak center boost");
        const auto cut = design(BiquadType::Peak, 1000, 1, -12);
        for (double hz : {20., 300., 1000., 5000., 18000.})
            near(gain_db(peak, hz, 48000) + gain_db(cut, hz, 48000), 0,
                 0.001, "matching peak boost/cut must cancel");
        const auto low_shelf = design(BiquadType::LowShelf, 1000, 1, 6);
        near(gain_db(low_shelf, 20, 48000), 6, 0.05, "low shelf bass gain");
        near(gain_db(low_shelf, 18000, 48000), 0, 0.05, "low shelf treble gain");
        const auto high_shelf = design(BiquadType::HighShelf, 1000, 1, -6);
        near(gain_db(high_shelf, 20, 48000), 0, 0.05, "high shelf bass gain");
        near(gain_db(high_shelf, 18000, 48000), -6, 0.05, "high shelf treble gain");
        const auto ap = design(BiquadType::AllPass);
        for (double hz : {20., 100., 1000., 5000., 18000.})
            near(gain_db(ap, hz, 48000), 0, 0.001, "all-pass magnitude");

        // Check the recurrence against an independent double-precision DF1.
        Biquad filter(peak);
        std::vector<float> input(8192), output(input.size());
        for (std::size_t i = 0; i < input.size(); ++i)
            input[i] = float(0.2 * std::sin(0.031 * i) + 0.13 * std::sin(0.183 * i));
        for (std::size_t offset = 0; offset < input.size(); offset += 256)
            filter.process(input.data() + offset, output.data() + offset, 256);
        double x1 = 0, x2 = 0, y1 = 0, y2 = 0, max_error = 0;
        for (std::size_t i = 0; i < input.size(); ++i) {
            const double x = input[i];
            const double y = peak.b0 * x + peak.b1 * x1 + peak.b2 * x2 -
                             peak.a1 * y1 - peak.a2 * y2;
            max_error = std::max(max_error, std::abs(y - output[i]));
            x2 = x1; x1 = x; y2 = y1; y1 = y;
        }
        require(max_error < 3.e-5, "biquad differs from independent DF1");
        filter.reset();
        auto in_place = input;
        for (std::size_t offset = 0; offset < input.size(); offset += 256)
            filter.process(in_place.data() + offset, in_place.data() + offset, 256);
        for (std::size_t i = 0; i < input.size(); ++i)
            require(in_place[i] == output[i], "in-place output differs");

        filter.reset();
        filter.set_coefficients(lp);
        filter.ramp_to(hp, 500);
        filter.process(input.data(), output.data(), 183);
        require(filter.ramping(), "coefficient ramp ended early");
        filter.process(input.data() + 183, output.data() + 183, 317);
        require(!filter.ramping(), "coefficient ramp did not end");
        const auto settled = filter.coefficients();
        require(settled.b0 == hp.b0 && settled.a2 == hp.a2, "ramp target not exact");
        filter.process(input.data() + 500, output.data() + 500, input.size() - 500);
        for (float sample : output) require(std::isfinite(sample), "ramp produced nonfinite output");

        // Sweep safe design limits, including high sample rates and Q.
        for (double sample_rate : {8000., 44100., 48000., 96000., 192000.})
            for (auto type : {BiquadType::LowPass, BiquadType::HighPass,
                              BiquadType::BandPass, BiquadType::Notch,
                              BiquadType::Peak, BiquadType::LowShelf,
                              BiquadType::HighShelf, BiquadType::AllPass})
                for (double frequency : {10., 100., 1000., 10000., 100000.})
                    for (double q : {0.1, 0.707, 18.}) {
                        const auto c = shores_dsp::design_biquad(
                            sample_rate, {type, frequency, q, 36., 1.});
                        for (double v : {c.b0, c.b1, c.b2, c.a1, c.a2})
                            require(std::isfinite(v), "nonfinite coefficient");
                        require(1. + c.a1 + c.a2 > 0. &&
                                1. - c.a1 + c.a2 > 0. &&
                                1. - c.a2 > 0. && 1. + c.a2 > 0.,
                                "unstable denominator");
                    }
        for (double bad : {std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::infinity()}) {
            try {
                shores_dsp::design_biquad(48000, {BiquadType::Peak, bad});
                throw std::runtime_error("nonfinite frequency accepted");
            } catch (const std::invalid_argument &) {}
        }
        try {
            shores_dsp::design_biquad(1.e9, {BiquadType::LowPass, 1000});
            throw std::runtime_error("unsupported sample rate accepted");
        } catch (const std::invalid_argument &) {}
        std::cout << "PASS: eight response modes, independent DF1, ramps, in-place, "
                     "safe parameter sweep; max DF1 error " << max_error << '\n';
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
