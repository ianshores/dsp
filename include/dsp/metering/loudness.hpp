// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Ian Shores
#pragma once
#include "dsp/metering/basic.hpp"
#include <array>
#include <initializer_list>

namespace shores_dsp {

inline double energy_to_lufs(double weighted_mean_square) noexcept {
    return weighted_mean_square > 0 && std::isfinite(weighted_mean_square)
        ? -0.691 + 10 * std::log10(weighted_mean_square)
        : -std::numeric_limits<double>::infinity();
}

struct KWeightingCoefficients {
    BiquadCoefficients shelf{1.53512485958697, -2.69169618940638, 1.19839281085285,
                             -1.69065929318241, 0.73248077421585};
    BiquadCoefficients highpass{1., -2., 1., -1.99004745483398, 0.99007225036621};
};

// Inverse bilinear polynomial of a published section in u=(1-z^-1)/(1+z^-1).
// Recover its natural frequency in physical Hz, then prewarp at the new rate.
// This preserves prototype damping, numerator shape and endpoint gains (also
// the published RLB gain), rather than substituting a generic normalized shelf.
inline KWeightingCoefficients design_k_weighting(double fs) {
    metering_detail::rate(fs);
    if (fs < 4000) throw std::invalid_argument("K weighting needs 4000..768000 Hz: 1682 Hz shelf below Nyquist");
    KWeightingCoefficients result;
    if (fs == 48000) return result; // exact decimal tables, no recomputation
    const auto retune = [fs](BiquadCoefficients c) {
        constexpr long double pi = 3.1415926535897932384626433832795L;
        const long double a0 = 1.L + c.a1 + c.a2, a1 = 2.L * (1.L - c.a2), a2 = 1.L - c.a1 + c.a2;
        const long double g = std::sqrt(a0 / a2);
        const long double hz = std::atan(g) * 48000 / pi;
        const long double scale = g / std::tan(pi * hz / fs);
        const long double d0 = a0 + a1 * scale + a2 * scale * scale;
        const long double b0 = static_cast<long double>(c.b0) + c.b1 + c.b2;
        const long double b1 = 2.L * (static_cast<long double>(c.b0) - c.b2);
        const long double b2 = static_cast<long double>(c.b0) - c.b1 + c.b2;
        return BiquadCoefficients{
            double((b0 + b1 * scale + b2 * scale * scale) / d0),
            double(2 * (b0 - b2 * scale * scale) / d0),
            double((b0 - b1 * scale + b2 * scale * scale) / d0),
            double(2 * (a0 - a2 * scale * scale) / d0),
            double((a0 - a1 * scale + a2 * scale * scale) / d0)};
    };
    result.shelf = retune(result.shelf); result.highpass = retune(result.highpass);
    return result;
}

// Existing Biquad twice; all design runs outside callbacks. 48 kHz retains the
// published coefficients and processing path. configure() resets the epoch;
// failed validation leaves coefficients, states, rate and status unchanged.
class KWeightingFilter : public metering_detail::Status {
    BiquadFilter shelf_{{1.53512485958697, -2.69169618940638, 1.19839281085285,
                       -1.69065929318241, 0.73248077421585}};
    BiquadFilter highpass_{{1., -2., 1., -1.99004745483398, 0.99007225036621}};
    double rate_ = 48000;
  public:
    explicit KWeightingFilter(double fs = 48000) { configure(fs); }
    void configure(double fs) {
        const auto c = design_k_weighting(fs);
        shelf_.set_coefficients(c.shelf); highpass_.set_coefficients(c.highpass);
        rate_ = fs; reset();
    }
    double sample_rate() const noexcept { return rate_; }
    void reset() noexcept { shelf_.reset(); highpass_.reset(); invalid_ = false; }
    float process(float input) noexcept {
        float x = float(sample(input));
        shelf_.process(&x, &x, 1);
        if (!std::isfinite(x)) { reset(); invalid_ = true; x = 0; }
        highpass_.process(&x, &x, 1);
        if (!std::isfinite(x)) { reset(); invalid_ = true; x = 0; }
        return x;
    }
    void process(const float *input, float *output, std::size_t n) noexcept {
        for (std::size_t i = 0; i < n; ++i) output[i] = process(input[i]);
    }
};

enum class ChannelRole { Left, Right, Centre, LeftSurround, RightSurround, Lfe };

// Explicit channel order: caller passes one role per input channel. Annex 1
// roles only; no implicit mapping from count, no object/height/7.1 support.
// Energies are summed, not waveforms. LFE gets zero weight, not +10 dB.
class ChannelWeighting : public metering_detail::Status {
    std::array<double, 6> weights_{};
    std::size_t count_ = 0;
  public:
    explicit ChannelWeighting(std::initializer_list<ChannelRole> roles = {ChannelRole::Centre}) {
        if (roles.size() == 0 || roles.size() > 6) throw std::invalid_argument("need 1..6 channel roles");
        std::array<bool, 6> used{};
        for (auto role : roles) {
            const int i = int(role);
            if (i < 0 || i > 5 || used[std::size_t(i)]) throw std::invalid_argument("unknown or duplicate channel role");
            used[std::size_t(i)] = true;
            weights_[count_++] = role == ChannelRole::Lfe ? 0 :
                (role == ChannelRole::LeftSurround || role == ChannelRole::RightSurround ? 1.41 : 1.);
        }
    }
    std::size_t channels() const noexcept { return count_; }
    void reset() noexcept { invalid_ = false; }
    double process(const double *energies) noexcept {
        double sum = 0;
        for (std::size_t i = 0; i < count_; ++i) sum += weights_[i] * energy(energies[i]);
        return sum;
    }
};

struct IntegratedLoudness {
    double lufs = -std::numeric_limits<double>::infinity();
    double relative_threshold_lufs = -std::numeric_limits<double>::infinity();
    std::size_t stored_blocks = 0, absolute_blocks = 0, included_blocks = 0;
    bool complete = true; // false on invalid input or exhausted history capacity
};

// Full nearest-sample 400 ms blocks, starts round_half_up(j*fs/10). The exact
// rational clock bounds start error by half a sample, with no cumulative drift.
// Fractional hops/overlaps vary by one frame: exact 75% is not representable at
// arbitrary rates. Integral-hop rates retain the original four-hop arithmetic.
// Other rates use a fixed block-sized energy ring. No partial first/final block.
// Setup allocates capacity doubles (default 36,000, max 1,000,000). Audio work
// and reset are O(1) per frame; the history never grows or gets decimated.
class LoudnessGate : public metering_detail::Status {
    std::vector<double> blocks_;
    std::array<double, 4> hops_{};
    std::size_t stored_ = 0, hop_head_ = 0, hop_count_ = 0, in_hop_ = 0;
    metering_detail::Sum hop_sum_;
    metering_detail::HopClock clock_;
    std::vector<double> ring_;
    metering_detail::Sum window_sum_;
    std::size_t ring_head_ = 0, ring_filled_ = 0, until_block_ = 19200, block_frames_ = 19200;
    double rate_ = 48000;
    bool overflow_ = false;
    static std::size_t capacity(std::size_t n) {
        if (n == 0 || n > 1000000) throw std::invalid_argument("gate capacity must be 1..1000000 blocks");
        return n;
    }
  public:
    explicit LoudnessGate(std::size_t block_capacity = 36000, double fs = 48000)
        : blocks_(capacity(block_capacity), 0.) {
        configure(fs);
    }
    // Setup only; preserves configured history capacity, starts a new epoch.
    void configure(double fs) {
        metering_detail::HopClock clock(fs);
        const auto width = metering_detail::RateRatio(fs).block_frames();
        std::vector<double> ring(clock.integral() ? 0 : width, 0.);
        clock_ = clock; ring_ = std::move(ring); block_frames_ = width; rate_ = fs;
        reset();
    }
    double sample_rate() const noexcept { return rate_; }
    std::size_t block_frames() const noexcept { return block_frames_; }
    std::size_t hop_floor_frames() const noexcept { return clock_.floor_frames(); }
    bool integral_hops() const noexcept { return clock_.integral(); }
    void reset() noexcept {
        stored_ = hop_head_ = hop_count_ = in_hop_ = 0;
        hop_sum_.reset(); overflow_ = invalid_ = false;
        clock_.reset(); ring_head_ = ring_filled_ = 0; until_block_ = block_frames_; window_sum_.reset();
    }
    std::size_t stored_blocks() const noexcept { return stored_; }
    std::size_t block_capacity() const noexcept { return blocks_.size(); }
    bool overflowed() const noexcept { return overflow_; }
    // Returns true whenever a complete gating block boundary is reached,
    // including after capacity is exhausted. The stored prefix then freezes.
    bool process(double weighted_energy) noexcept {
        const double e = energy(weighted_energy);
        if (!clock_.integral()) {
            if (ring_filled_ == ring_.size()) window_sum_.add(-ring_[ring_head_]); else ++ring_filled_;
            ring_[ring_head_] = e; window_sum_.add(e);
            if (++ring_head_ == ring_.size()) ring_head_ = 0;
            if (--until_block_ != 0) return false;
            until_block_ = clock_.next();
            if (stored_ == blocks_.size()) { overflow_ = true; return true; }
            blocks_[stored_++] = std::max(0., window_sum_.value()) / double(block_frames_);
            return true;
        }
        hop_sum_.add(e);
        if (++in_hop_ != clock_.floor_frames()) return false;
        hops_[hop_head_] = hop_sum_.value() / double(clock_.floor_frames());
        hop_head_ = (hop_head_ + 1) % 4; hop_sum_.reset(); in_hop_ = 0;
        if (hop_count_ < 4) ++hop_count_;
        if (hop_count_ < 4) return false;
        if (stored_ == blocks_.size()) { overflow_ = true; return true; }
        blocks_[stored_++] = ((hops_[0] + hops_[1]) + (hops_[2] + hops_[3])) * 0.25;
        return true;
    }
    void process(const double *weighted_energy, std::size_t n) noexcept {
        for (std::size_t i = 0; i < n; ++i) process(weighted_energy[i]);
    }
    // Query at a caller-controlled analysis boundary, not in an audio callback:
    // two passes O(stored blocks), at most 2,000,000 iterations. No allocation.
    // No concurrent processing/query. Invalid/overflow results remain diagnostic
    // prefix readings only; complete=false prohibits a full-program claim.
    IntegratedLoudness result() const noexcept {
        IntegratedLoudness r;
        r.stored_blocks = stored_; r.complete = valid() && !overflow_;
        const double absolute = std::pow(10., (-70. + 0.691) / 10.);
        metering_detail::Sum first;
        for (std::size_t i = 0; i < stored_; ++i)
            if (blocks_[i] > absolute) { first.add(blocks_[i]); ++r.absolute_blocks; }
        if (r.absolute_blocks == 0) return r;
        const double relative = first.value() / double(r.absolute_blocks) * 0.1;
        r.relative_threshold_lufs = energy_to_lufs(relative);
        const double threshold = std::max(absolute, relative);
        metering_detail::Sum second;
        for (std::size_t i = 0; i < stored_; ++i)
            if (blocks_[i] > threshold) { second.add(blocks_[i]); ++r.included_blocks; }
        if (r.included_blocks) r.lufs = energy_to_lufs(second.value() / double(r.included_blocks));
        return r;
    }
};

} // namespace shores_dsp
