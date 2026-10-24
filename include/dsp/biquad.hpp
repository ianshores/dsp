// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Ian Shores
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace shores_dsp {

enum class BiquadType {
    LowPass,
    HighPass,
    BandPass,
    Notch,
    Peak,
    LowShelf,
    HighShelf,
    AllPass
};

struct BiquadParameters {
    BiquadType type = BiquadType::LowPass;
    double frequency_hz = 1000.0;
    double q = 0.7071067811865476;
    double gain_db = 0.0;
    double shelf_slope = 1.0;
};

// H(z) = (b0 + b1 z^-1 + b2 z^-2) / (1 + a1 z^-1 + a2 z^-2).
// Band-pass is normalized to 0 dB at its center frequency. Gain is used only
// by Peak and the shelves. Shelf slope replaces Q for shelf modes.
struct BiquadCoefficients {
    double b0 = 1.0, b1 = 0.0, b2 = 0.0, a1 = 0.0, a2 = 0.0;
};

// Design coefficients off the audio thread. The safe user ranges are 10 Hz to
// 0.49 * sample rate, Q 0.1..18, gain -36..36 dB, shelf slope 0.1..1.
// Finite out-of-range values are clamped; NaN/infinity and invalid rates throw.
inline BiquadCoefficients design_biquad(double sample_rate,
                                         const BiquadParameters &parameters) {
    if (!std::isfinite(sample_rate) || sample_rate < 1000.0 || sample_rate > 768000.0 ||
        !std::isfinite(parameters.frequency_hz) || !std::isfinite(parameters.q) ||
        !std::isfinite(parameters.gain_db) || !std::isfinite(parameters.shelf_slope))
        throw std::invalid_argument("biquad parameters must be finite; sample rate 1000..768000 Hz");

    const double frequency = std::clamp(parameters.frequency_hz, 10.0, sample_rate * 0.49);
    const double q = std::clamp(parameters.q, 0.1, 18.0);
    const double gain = std::clamp(parameters.gain_db, -36.0, 36.0);
    const double slope = std::clamp(parameters.shelf_slope, 0.1, 1.0);
    constexpr double tau = 6.283185307179586476925286766559;
    const double omega = tau * frequency / sample_rate;
    const double c = std::cos(omega), s = std::sin(omega);
    const double alpha = s / (2.0 * q);
    const bool uses_gain = parameters.type == BiquadType::Peak ||
                           parameters.type == BiquadType::LowShelf ||
                           parameters.type == BiquadType::HighShelf;
    const double A = uses_gain ? std::pow(10.0, gain / 40.0) : 1.0;
    double b0 = 1.0, b1 = 0.0, b2 = 0.0, a0 = 1.0, a1 = 0.0, a2 = 0.0;

    switch (parameters.type) {
    case BiquadType::LowPass:
        b0 = (1.0 - c) * 0.5; b1 = 1.0 - c; b2 = b0;
        a0 = 1.0 + alpha; a1 = -2.0 * c; a2 = 1.0 - alpha;
        break;
    case BiquadType::HighPass:
        b0 = (1.0 + c) * 0.5; b1 = -(1.0 + c); b2 = b0;
        a0 = 1.0 + alpha; a1 = -2.0 * c; a2 = 1.0 - alpha;
        break;
    case BiquadType::BandPass: // constant 0 dB peak gain
        b0 = alpha; b1 = 0.0; b2 = -alpha;
        a0 = 1.0 + alpha; a1 = -2.0 * c; a2 = 1.0 - alpha;
        break;
    case BiquadType::Notch:
        b0 = 1.0; b1 = -2.0 * c; b2 = 1.0;
        a0 = 1.0 + alpha; a1 = -2.0 * c; a2 = 1.0 - alpha;
        break;
    case BiquadType::Peak:
        b0 = 1.0 + alpha * A; b1 = -2.0 * c; b2 = 1.0 - alpha * A;
        a0 = 1.0 + alpha / A; a1 = -2.0 * c; a2 = 1.0 - alpha / A;
        break;
    case BiquadType::LowShelf:
    case BiquadType::HighShelf: {
        const double shelf_alpha = s * 0.5 *
            std::sqrt((A + 1.0 / A) * (1.0 / slope - 1.0) + 2.0);
        const double two_root_a_alpha = 2.0 * std::sqrt(A) * shelf_alpha;
        if (parameters.type == BiquadType::LowShelf) {
            b0 = A * ((A + 1.0) - (A - 1.0) * c + two_root_a_alpha);
            b1 = 2.0 * A * ((A - 1.0) - (A + 1.0) * c);
            b2 = A * ((A + 1.0) - (A - 1.0) * c - two_root_a_alpha);
            a0 = (A + 1.0) + (A - 1.0) * c + two_root_a_alpha;
            a1 = -2.0 * ((A - 1.0) + (A + 1.0) * c);
            a2 = (A + 1.0) + (A - 1.0) * c - two_root_a_alpha;
        } else {
            b0 = A * ((A + 1.0) + (A - 1.0) * c + two_root_a_alpha);
            b1 = -2.0 * A * ((A - 1.0) + (A + 1.0) * c);
            b2 = A * ((A + 1.0) + (A - 1.0) * c - two_root_a_alpha);
            a0 = (A + 1.0) - (A - 1.0) * c + two_root_a_alpha;
            a1 = 2.0 * ((A - 1.0) - (A + 1.0) * c);
            a2 = (A + 1.0) - (A - 1.0) * c - two_root_a_alpha;
        }
        break;
    }
    case BiquadType::AllPass:
        b0 = 1.0 - alpha; b1 = -2.0 * c; b2 = 1.0 + alpha;
        a0 = 1.0 + alpha; a1 = -2.0 * c; a2 = 1.0 - alpha;
        break;
    default:
        throw std::invalid_argument("unknown biquad type");
    }
    const double scale = 1.0 / a0;
    return {b0 * scale, b1 * scale, b2 * scale, a1 * scale, a2 * scale};
}

// One channel and two state values per instance. Share designed coefficients
// across channels, but give each channel its own Biquad. The normal block path
// has five multiplies per sample and no per-sample mode or ramp branch.
class Biquad {
    BiquadCoefficients coefficients_{}, target_{}, increment_{};
    double z1_ = 0.0, z2_ = 0.0;
    std::size_t ramp_remaining_ = 0;

    static double clean(double value) noexcept {
        // Clear inaudible state well above the double denormal range.
        return std::abs(value) < 1.e-30 ? 0.0 : value;
    }
    void finish_block(double z1, double z2) noexcept {
        z1_ = clean(z1);
        z2_ = clean(z2);
    }

  public:
    Biquad() noexcept = default;
    explicit Biquad(BiquadCoefficients coefficients) noexcept
        : coefficients_(coefficients), target_(coefficients) {}

    BiquadCoefficients coefficients() const noexcept { return coefficients_; }
    bool ramping() const noexcept { return ramp_remaining_ != 0; }
    void reset() noexcept { z1_ = z2_ = 0.0; }

    // Call from the audio thread at a known event boundary; the caller owns
    // publication/synchronization of coefficients from its control thread.
    void set_coefficients(BiquadCoefficients coefficients) noexcept {
        coefficients_ = target_ = coefficients;
        increment_ = {};
        ramp_remaining_ = 0;
    }
    void ramp_to(BiquadCoefficients target, std::size_t samples) noexcept {
        if (samples == 0) { set_coefficients(target); return; }
        target_ = target;
        const double inverse = 1.0 / double(samples);
        increment_ = {(target.b0 - coefficients_.b0) * inverse,
                      (target.b1 - coefficients_.b1) * inverse,
                      (target.b2 - coefficients_.b2) * inverse,
                      (target.a1 - coefficients_.a1) * inverse,
                      (target.a2 - coefficients_.a2) * inverse};
        ramp_remaining_ = samples;
    }

    // Exact in-place operation is supported; partially overlapping buffers
    // are not. Input is assumed finite. No allocations, locks or I/O.
    void process(const float *input, float *output, std::size_t count) noexcept {
        double z1 = z1_, z2 = z2_;
        while (count != 0) {
            if (ramp_remaining_ == 0) {
                const auto c = coefficients_;
                for (std::size_t i = 0; i < count; ++i) {
                    const double x = input[i];
                    const double y = c.b0 * x + z1;
                    z1 = c.b1 * x - c.a1 * y + z2;
                    z2 = c.b2 * x - c.a2 * y;
                    output[i] = float(y);
                }
                break;
            }
            const auto n = std::min(count, ramp_remaining_);
            for (std::size_t i = 0; i < n; ++i) {
                coefficients_.b0 += increment_.b0;
                coefficients_.b1 += increment_.b1;
                coefficients_.b2 += increment_.b2;
                coefficients_.a1 += increment_.a1;
                coefficients_.a2 += increment_.a2;
                const double x = input[i];
                const double y = coefficients_.b0 * x + z1;
                z1 = coefficients_.b1 * x - coefficients_.a1 * y + z2;
                z2 = coefficients_.b2 * x - coefficients_.a2 * y;
                output[i] = float(y);
            }
            ramp_remaining_ -= n;
            if (ramp_remaining_ == 0) coefficients_ = target_;
            input += n;
            output += n;
            count -= n;
        }
        finish_block(z1, z2);
    }
};

} // namespace shores_dsp
