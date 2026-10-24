// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Ian Shores
#pragma once
#include "dsp/metering/basic.hpp"
#include <array>

namespace shores_dsp {

// 48 kHz: published 4x FIR. Other rates: setup-designed 64-tap/phase FIR.
// ITU-R BS.1770-5 Annex 2's published interpolation coefficients (numeric
// algorithm data, independently implemented here). Double MAC/state and peak.
// At 48 kHz: delay 5.875 frames, support 12. Other rates: delay 32, support 64.
// Generated Blackman-windowed sinc has unity DC per phase. Factor is
// max(2, ceil(384000/fs)); validation band min(20 kHz, .45*fs). Setup may allocate.
// No limiter or silent clipping to 0 dBFS. This is an estimate; fixture scope
// and finite-band accuracy must be read alongside docs/metering.md.
class TruePeakOversampler : public metering_detail::Status {
    inline static constexpr std::array<std::array<double, 4>, 12> fir_{{
        {{ .001708984375, -.0291748046875, -.0189208984375, -.00830078125 }},
        {{ .010986328125, .029296875, .0330810546875, .014892578125 }},
        {{ -.0196533203125, -.0517578125, -.0582275390625, -.026611328125 }},
        {{ .033203125, .089111328125, .1015625, .047607421875 }},
        {{ -.0594482421875, -.16650390625, -.2003173828125, -.102294921875 }},
        {{ .1373291015625, .465087890625, .77978515625, .97216796875 }},
        {{ .97216796875, .77978515625, .465087890625, .1373291015625 }},
        {{ -.102294921875, -.2003173828125, -.16650390625, -.0594482421875 }},
        {{ .047607421875, .1015625, .089111328125, .033203125 }},
        {{ -.026611328125, -.0582275390625, -.0517578125, -.0196533203125 }},
        {{ .014892578125, .0330810546875, .029296875, .010986328125 }},
        {{ -.00830078125, -.0189208984375, -.0291748046875, .001708984375 }}
    }};
    std::array<double, 64> history_{};
    std::vector<double> coefficients_;
    double rate_ = 48000;
    std::size_t factor_ = 4;
    std::size_t head_ = 0;
    double peak_ = 0;
  public:
    explicit TruePeakOversampler(double fs = 48000) { configure(fs); }
    // Setup only. Successful changes reset the epoch; failed setup preserves it.
    void configure(double fs) {
        const metering_detail::RateRatio ratio(fs);
        const auto factor = fs == 48000 ? 4 : std::max<std::size_t>(2, ratio.factor_for_384khz());
        std::vector<double> coefficients(fs == 48000 ? 0 : factor * 64);
        constexpr long double pi = 3.1415926535897932384626433832795L;
        if (fs != 48000) for (std::size_t p = 0; p < factor; ++p) {
            long double sum = 0;
            for (std::size_t k = 0; k < 64; ++k) {
                const long double t = static_cast<long double>(k) - 32 + static_cast<long double>(p) / factor;
                const auto sinc = t == 0 ? 1.L : std::sin(pi*t)/(pi*t);
                const auto window = .42L + .5L*std::cos(pi*t/32) + .08L*std::cos(2*pi*t/32);
                const auto c = sinc * window;
                coefficients[k*factor+p] = double(c); sum += c;
            }
            for (std::size_t k = 0; k < 64; ++k) coefficients[k*factor+p] /= double(sum);
        }
        coefficients_ = std::move(coefficients); rate_ = fs; factor_ = factor; reset();
    }
    double sample_rate() const noexcept { return rate_; }
    std::size_t factor() const noexcept { return factor_; }
    std::size_t taps_per_phase() const noexcept { return rate_ == 48000 ? 12 : 64; }
    double latency_frames() const noexcept { return rate_ == 48000 ? 5.875 : 32.; }
    std::size_t flush_frames() const noexcept { return taps_per_phase() - 1; }
    // Engineering validation band; no universal error guarantee above/below it.
    double validation_bandwidth_hz() const noexcept { return std::min(20000., .45*rate_); }
    void reset() noexcept { history_.fill(0); head_ = 0; peak_ = 0; invalid_ = false; }
    double value() const noexcept { return peak_; }
    double dbtp() const noexcept {
        return peak_ > 0 ? 20 * std::log10(peak_) : -std::numeric_limits<double>::infinity();
    }
    // Returns the maximum of the newly interpolated points, not the
    // lifetime maximum. value()/dbtp() retain the maximum since reset.
    double process(float input) noexcept {
        history_[head_] = sample(input);
        double current = 0;
        if (rate_ == 48000) {
            std::array<double, 4> out{};
            std::size_t at = head_;
            for (std::size_t tap = 0; tap < 12; ++tap) {
                for (std::size_t phase = 0; phase < 4; ++phase)
                    out[phase] += fir_[tap][phase] * history_[at];
                at = at == 0 ? 11 : at - 1;
            }
            if (++head_ == 12) head_ = 0;
            for (double x : out) current = std::max(current, std::abs(x));
        } else {
            std::array<double, 384> out;
            std::fill_n(out.begin(), factor_, 0.);
            std::size_t at = head_;
            for (std::size_t tap = 0; tap < 64; ++tap) {
                const auto c = coefficients_.data() + tap * factor_;
                for (std::size_t phase = 0; phase < factor_; ++phase)
                    out[phase] += c[phase] * history_[at];
                at = at == 0 ? 63 : at - 1;
            }
            if (++head_ == 64) head_ = 0;
            for (std::size_t p = 0; p < factor_; ++p) current = std::max(current, std::abs(out[p]));
        }
        peak_ = std::max(peak_, current); return current;
    }
    double process(const float *input, std::size_t n) noexcept {
        for (std::size_t i = 0; i < n; ++i) process(input[i]);
        return peak_;
    }
    // File/end-of-stream only: evaluate all delayed tail contributions under
    // zero extension. Advances flush_frames() frames; not for block boundaries.
    double flush() noexcept { for (std::size_t i = 0; i < flush_frames(); ++i) process(0); return peak_; }
};

} // namespace shores_dsp
