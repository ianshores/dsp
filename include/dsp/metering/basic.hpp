// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Ian Shores
#pragma once

#include "dsp/biquad.hpp"
#include <cstdint>
#include <limits>
#include <vector>

namespace shores_dsp {

// The existing general IIR remains the sole biquad implementation.
using BiquadFilter = Biquad;

namespace metering_detail {
inline void rate(double fs) {
    if (!std::isfinite(fs) || fs < 1000 || fs > 768000)
        throw std::invalid_argument("metering rate must be 1000..768000 Hz");
}
inline void seconds(double t, double maximum = 60) {
    if (!std::isfinite(t) || t < 0 || t > maximum)
        throw std::invalid_argument("metering time outside supported range");
}
// Nonfinite/out-of-domain input is replaced by zero, preserving stream time.
// Invalidity is sticky until reset; callers must inspect valid(). No exceptions
// or logging in callback paths. No class here provides inter-thread publication.
class Status {
  protected:
    bool invalid_ = false;
    double sample(float x) noexcept {
        if (!std::isfinite(x)) { invalid_ = true; return 0; }
        return x;
    }
    double energy(double x) noexcept {
        // Generous finite domain keeps bounded sums far below double overflow.
        if (!std::isfinite(x) || x < 0 || x > 1.e100) { invalid_ = true; return 0; }
        return x;
    }
  public:
    bool valid() const noexcept { return !invalid_; }
};
class Sum {
    double value_ = 0, correction_ = 0;
  public:
    void add(double x) noexcept {
        const double t = value_ + x;
        correction_ += std::abs(value_) >= std::abs(x) ? (value_ - t) + x : (x - t) + value_;
        value_ = t;
    }
    double value() const noexcept { return value_ + correction_; }
    void reset() noexcept { value_ = correction_ = 0; }
};
} // namespace metering_detail

class PeakDetector : public metering_detail::Status {
    double peak_ = 0;
  public:
    void reset() noexcept { peak_ = 0; invalid_ = false; }
    double value() const noexcept { return peak_; }
    double process(float x) noexcept { peak_ = std::max(peak_, std::abs(sample(x))); return peak_; }
    double process(const float *input, std::size_t n) noexcept {
        for (std::size_t i = 0; i < n; ++i) process(input[i]);
        return peak_;
    }
};

// Rectified-amplitude one-pole follower, not an RMS or EBU loudness time window.
// Attack/release are 1/e time constants in seconds; zero means instantaneous.
class EnvelopeFollower : public metering_detail::Status {
    double attack_ = 0, release_ = 0, value_ = 0;
  public:
    EnvelopeFollower(double fs = 48000, double attack_seconds = 0.01,
                     double release_seconds = 0.3) { configure(fs, attack_seconds, release_seconds); }
    // Setup/control thread only; changing ballistics preserves envelope state.
    void configure(double fs, double attack_seconds, double release_seconds) {
        metering_detail::rate(fs); metering_detail::seconds(attack_seconds);
        metering_detail::seconds(release_seconds);
        attack_ = attack_seconds == 0 ? 0 : std::exp(-1 / (fs * attack_seconds));
        release_ = release_seconds == 0 ? 0 : std::exp(-1 / (fs * release_seconds));
    }
    void reset() noexcept { value_ = 0; invalid_ = false; }
    double value() const noexcept { return value_; }
    double process(float x) noexcept {
        const double magnitude = std::abs(sample(x));
        const double c = magnitude > value_ ? attack_ : release_;
        value_ = magnitude + c * (value_ - magnitude);
        if (value_ < 1.e-30) value_ = 0;
        return value_;
    }
    void process(const float *input, double *output, std::size_t n) noexcept {
        for (std::size_t i = 0; i < n; ++i) output[i] = process(input[i]);
    }
};

// Input is a detected nonnegative linear amplitude, not signed audio. A new
// value >= the held value refreshes the hold. It holds for N FOLLOWING samples,
// then decays exponentially by decay_db_per_second; zero decay holds forever.
class PeakHold : public metering_detail::Status {
    std::uint64_t hold_ = 0, remaining_ = 0;
    double decay_ = 1, value_ = 0;
  public:
    PeakHold(double fs = 48000, double hold_seconds = 1, double decay_db_per_second = 20) {
        metering_detail::rate(fs); metering_detail::seconds(hold_seconds);
        if (!std::isfinite(decay_db_per_second) || decay_db_per_second < 0 || decay_db_per_second > 1000)
            throw std::invalid_argument("peak decay must be 0..1000 dB/second");
        hold_ = std::uint64_t(std::llround(fs * hold_seconds));
        decay_ = std::exp(-decay_db_per_second * std::log(10.) / (20 * fs));
    }
    void reset() noexcept { value_ = 0; remaining_ = 0; invalid_ = false; }
    double value() const noexcept { return value_; }
    double process(double detected_peak) noexcept {
        if (!std::isfinite(detected_peak) || detected_peak < 0 || detected_peak > 1.e50) {
            invalid_ = true; detected_peak = 0;
        }
        if (detected_peak >= value_) { value_ = detected_peak; remaining_ = hold_; }
        else if (remaining_ > 0) --remaining_;
        else value_ = std::max(detected_peak, value_ * decay_);
        if (value_ < 1.e-30) value_ = 0;
        return value_;
    }
};

// Running arithmetic mean of squared float input. No growing sum or storage.
// At capacity the prefix freezes and overflowed() becomes sticky; reset to
// start another measurement. The default 2^32 samples is ~24.9 hours at 48 kHz.
class MeanSquareAccumulator : public metering_detail::Status {
    std::uint64_t count_ = 0, capacity_;
    double mean_ = 0;
    bool overflow_ = false;
  public:
    explicit MeanSquareAccumulator(std::uint64_t capacity = std::uint64_t{1} << 32) : capacity_(capacity) {
        if (capacity == 0 || capacity > (std::uint64_t{1} << 53))
            throw std::invalid_argument("mean-square sample capacity must be 1..2^53");
    }
    void reset() noexcept { count_ = 0; mean_ = 0; overflow_ = invalid_ = false; }
    std::uint64_t count() const noexcept { return count_; }
    bool overflowed() const noexcept { return overflow_; }
    double mean_square() const noexcept { return mean_; }
    double rms() const noexcept { return std::sqrt(mean_); }
    bool process(float input) noexcept {
        const double x = sample(input);
        if (count_ == capacity_) { overflow_ = true; return false; }
        ++count_; mean_ += (x * x - mean_) / double(count_);
        return true;
    }
    void process(const float *input, std::size_t n) noexcept {
        for (std::size_t i = 0; i < n; ++i) process(input[i]);
    }
};

// Moving rectangular integration of NONNEGATIVE ENERGY, not audio amplitude.
// Setup allocates N doubles, with N=round(fs*seconds), 1..1,048,576 and time
// <=10 seconds. No callback allocation. Each push/reset is O(1); no ring clear
// on reset. Partial startup mean divides by available count; ready() identifies
// a complete window. Zero-count block calls throughout this API are no-ops.
class SlidingWindow : public metering_detail::Status {
    std::vector<double> ring_;
    std::size_t head_ = 0, filled_ = 0;
    metering_detail::Sum sum_;
    static std::size_t frames(double fs, double seconds) {
        metering_detail::rate(fs); metering_detail::seconds(seconds, 10);
        const auto n = std::llround(fs * seconds);
        if (n < 1 || n > 1048576) throw std::invalid_argument("window must contain 1..1048576 samples");
        return std::size_t(n);
    }
  public:
    SlidingWindow(double fs = 48000, double seconds = 0.4) : ring_(frames(fs, seconds), 0.) {}
    void reset() noexcept { head_ = filled_ = 0; sum_.reset(); invalid_ = false; }
    std::size_t capacity() const noexcept { return ring_.size(); }
    std::size_t count() const noexcept { return filled_; }
    bool ready() const noexcept { return filled_ == ring_.size(); }
    double mean() const noexcept { return filled_ ? std::max(0., sum_.value()) / double(filled_) : 0; }
    double process(double input_energy) noexcept {
        const double e = energy(input_energy);
        if (filled_ == ring_.size()) sum_.add(-ring_[head_]); else ++filled_;
        ring_[head_] = e; sum_.add(e);
        if (++head_ == ring_.size()) head_ = 0;
        return mean();
    }
    void process(const double *energy_values, double *means, std::size_t n) noexcept {
        for (std::size_t i = 0; i < n; ++i) means[i] = process(energy_values[i]);
    }
};

} // namespace shores_dsp
