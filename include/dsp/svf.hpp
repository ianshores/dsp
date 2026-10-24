// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Ian Shores
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace shores_dsp {

enum class SvfType { LowPass, HighPass, BandPass, Notch };

struct SvfParameters {
    SvfType type = SvfType::LowPass;
    double frequency_hz = 1000.0;
    double q = 0.7071067811865476;
};

// Prepared trapezoidal state-variable coefficients. Use design_svf() to create
// these; arbitrary edits are outside the processing contract. Default is a
// silent-state wire. All coefficients and integrator states are double; I/O is
// float. The mix selects LP, HP, unity-center BP, or notch without a mode branch.
struct SvfCoefficients {
    double a1 = 1.0, a2 = 0.0, a3 = 0.0;
    double m0 = 1.0, m1 = 0.0, m2 = 0.0;
};

// Call outside the audio callback, including when preparing modulation events.
// Cutoff/center frequency: 10 Hz..0.49*rate, Q: 0.1..18 (Butterworth Q=1/sqrt(2)).
// Finite out-of-range controls are clamped. Nonfinite controls, unknown modes,
// and rates outside 1000..768000 Hz throw. Q is the analog prototype Q, not a
// normalized resonance knob. BP is multiplied by 1/Q for unity center gain.
// Derived from Andrew Simper's optimized trapezoidal SVF nodal equations:
// https://cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf (pp. 2-7).
inline SvfCoefficients design_svf(double sample_rate, const SvfParameters &parameters) {
    if (!std::isfinite(sample_rate) || sample_rate < 1000.0 || sample_rate > 768000.0 ||
        !std::isfinite(parameters.frequency_hz) || !std::isfinite(parameters.q))
        throw std::invalid_argument("SVF parameters must be finite; sample rate 1000..768000 Hz");
    const double frequency = std::clamp(parameters.frequency_hz, 10.0, 0.49 * sample_rate);
    const double q = std::clamp(parameters.q, 0.1, 18.0);
    constexpr double pi = 3.1415926535897932384626433832795;
    const double g = std::tan(pi * frequency / sample_rate);
    const double k = 1.0 / q;
    const double a1 = 1.0 / (1.0 + g * (g + k));
    SvfCoefficients c{a1, g * a1, g * g * a1, 0.0, 0.0, 0.0};
    switch (parameters.type) {
    case SvfType::LowPass:  c.m2 = 1.0; break;
    case SvfType::HighPass: c.m0 = 1.0; c.m1 = -k; c.m2 = -1.0; break;
    case SvfType::BandPass: c.m1 = k; break;
    case SvfType::Notch:    c.m0 = 1.0; c.m1 = -k; break;
    default: throw std::invalid_argument("unknown SVF type");
    }
    return c;
}

// One instance per channel. Designed coefficients may be shared across channels.
// Caller owns control/audio thread publication and synchronization. Applying
// prepared coefficients preserves integrator state; reset() clears it without
// changing coefficients. On a sample-rate change, prepare new coefficients and
// apply at an event boundary; reset/crossfade if the application needs it.
class Svf {
    SvfCoefficients coefficients_{};
    double ic1_ = 0.0, ic2_ = 0.0;

    static float tick(float input, const SvfCoefficients &c,
                      double &ic1, double &ic2) noexcept {
        const double x = input;
        const double v3 = x - ic2;
        const double band = c.a1 * ic1 + c.a2 * v3;
        const double low = ic2 + c.a2 * ic1 + c.a3 * v3;
        ic1 = 2.0 * band - ic1;
        ic2 = 2.0 * low - ic2;
        return float(c.m0 * x + c.m1 * band + c.m2 * low);
    }
    void finish_block(double ic1, double ic2) noexcept {
        // Clear inaudible state well above double denormals, as in Biquad.
        ic1_ = std::abs(ic1) < 1.e-30 ? 0.0 : ic1;
        ic2_ = std::abs(ic2) < 1.e-30 ? 0.0 : ic2;
    }

  public:
    Svf() noexcept = default;
    explicit Svf(SvfCoefficients coefficients) noexcept : coefficients_(coefficients) {}
    SvfCoefficients coefficients() const noexcept { return coefficients_; }
    void reset() noexcept { ic1_ = ic2_ = 0.0; }
    void set_coefficients(SvfCoefficients coefficients) noexcept { coefficients_ = coefficients; }

    // Both processing paths assume finite input and valid prepared coefficients.
    // Exact in-place buffers supported; partial overlap is not. Zero count is a
    // no-op, allowing null pointers. No allocation, locks, I/O, trig, or throws.
    void process(const float *input, float *output, std::size_t count) noexcept {
        if (count == 0) return;
        double ic1 = ic1_, ic2 = ic2_;
        const auto c = coefficients_;
        for (std::size_t i = 0; i < count; ++i)
            output[i] = tick(input[i], c, ic1, ic2);
        finish_block(ic1, ic2);
    }

    // One independently designed coefficient set per sample; applied BEFORE
    // that sample. The final set remains current for subsequent process().
    // Prepare/smooth frequency and Q off-thread, then design each event; do not
    // interpolate these internal coefficients arbitrarily. Event storage must
    // be caller-owned, immutable during processing, and disjoint from audio.
    // This finite testable update contract is not an artifact-free guarantee:
    // abrupt control/mode changes can click; smoothing/crossfades belong to the
    // caller. No hidden ramp, sample-rate conversion, or concurrency mechanism.
    void process_modulated(const float *input, float *output,
                           const SvfCoefficients *events, std::size_t count) noexcept {
        if (count == 0) return;
        double ic1 = ic1_, ic2 = ic2_;
        for (std::size_t i = 0; i < count; ++i)
            output[i] = tick(input[i], events[i], ic1, ic2);
        coefficients_ = events[count - 1];
        finish_block(ic1, ic2);
    }
};

} // namespace shores_dsp
