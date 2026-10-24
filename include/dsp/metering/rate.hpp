// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Ian Shores
#pragma once
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace shores_dsp { namespace metering_detail {
inline void rate(double fs) {
    if (!std::isfinite(fs) || fs < 1000 || fs > 768000)
        throw std::invalid_argument("metering rate must be 1000..768000 Hz");
}

// Exact rational representation of the supplied binary64 rate. Bounds keep
// all integer products below uint64 limits; no elapsed-stream counter is used.
struct RateRatio {
    std::uint64_t numerator, denominator;
    explicit RateRatio(double fs) {
        rate(fs);
        static_assert(std::numeric_limits<double>::radix == 2 &&
                      std::numeric_limits<double>::digits == 53, "binary64 double required");
        int exponent = 0;
        const double fraction = std::frexp(fs, &exponent);
        numerator = std::uint64_t(std::ldexp(fraction, 53));
        denominator = std::uint64_t{1} << (53 - exponent);
    }
    std::size_t block_frames() const noexcept {
        // round_half_up(fs * 2/5), exactly, not a binary approximation to .4.
        const auto d = 5 * denominator;
        return std::size_t((2 * numerator + d / 2) / d);
    }
    std::size_t factor_for_384khz() const noexcept {
        const auto n = 384000 * denominator;
        return std::size_t((n + numerator - 1) / numerator);
    }
};

// Successive differences of round_half_up(j*fs/10), using exact integer
// remainder arithmetic. Bounded phase error <=0.5 input frame forever, including
// fractional/unusual rates. Intervals may differ by one; no accumulating drift.
class HopClock {
    std::uint64_t denominator_ = 1, remainder_ = 0, phase_ = 0;
    std::size_t whole_ = 0;
  public:
    explicit HopClock(double fs = 48000) {
        const RateRatio r(fs);
        denominator_ = 10 * r.denominator;
        whole_ = std::size_t(r.numerator / denominator_);
        remainder_ = r.numerator % denominator_;
        reset();
    }
    void reset() noexcept { phase_ = denominator_ / 2; }
    bool integral() const noexcept { return remainder_ == 0; }
    std::size_t floor_frames() const noexcept { return whole_; }
    std::size_t next() noexcept {
        phase_ += remainder_;
        const bool carry = phase_ >= denominator_;
        if (carry) phase_ -= denominator_;
        return whole_ + std::size_t(carry);
    }
};

inline std::size_t window_frames(double fs, double seconds) {
    rate(fs);
    if (!std::isfinite(seconds) || seconds < 0 || seconds > 10)
        throw std::invalid_argument("window time must be 0..10 seconds");
    double rounded = 0;
    if (seconds == 0.4) rounded = double(RateRatio(fs).block_frames());
    else {
        const double product = fs * seconds;
        const double lower = std::floor(product);
        const double fraction = product - lower;
        // When rounded multiplication lands exactly on a half, fma recovers
        // the sign of its error, preserving nearest/half-up for binary inputs.
        const bool up = fraction > 0.5 || (fraction == 0.5 && std::fma(fs, seconds, -product) >= 0);
        rounded = lower + double(up);
    }
    if (rounded < 1 || rounded > 2304000)
        throw std::invalid_argument("window must contain 1..2304000 samples");
    return std::size_t(rounded);
}
} } // namespace shores_dsp::metering_detail
