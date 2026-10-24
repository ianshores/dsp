// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Ian Shores
#pragma once
#include "dsp/metering/basic.hpp"
#include <array>

namespace shores_dsp {

// 48 kHz input only, 4x polyphase FIR, 12 taps/phase, 48 total coefficients.
// ITU-R BS.1770-5 Annex 2's published interpolation coefficients (numeric
// algorithm data, independently implemented here). Double MAC/state and peak.
// Group delay: 23.5/4 = 5.875 input samples; causal support spans 12 samples.
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
    std::array<double, 12> history_{};
    std::size_t head_ = 0;
    double peak_ = 0;
  public:
    explicit TruePeakOversampler(double fs = 48000) {
        if (fs != 48000) throw std::invalid_argument("true peak currently supports exactly 48000 Hz");
    }
    void reset() noexcept { history_.fill(0); head_ = 0; peak_ = 0; invalid_ = false; }
    double value() const noexcept { return peak_; }
    double dbtp() const noexcept {
        return peak_ > 0 ? 20 * std::log10(peak_) : -std::numeric_limits<double>::infinity();
    }
    // Returns the maximum of the four newly interpolated points, not the
    // lifetime maximum. value()/dbtp() retain the maximum since reset.
    double process(float input) noexcept {
        history_[head_] = sample(input);
        std::array<double, 4> out{};
        std::size_t at = head_;
        for (std::size_t tap = 0; tap < 12; ++tap) {
            for (std::size_t phase = 0; phase < 4; ++phase)
                out[phase] += fir_[tap][phase] * history_[at];
            at = at == 0 ? 11 : at - 1;
        }
        if (++head_ == 12) head_ = 0;
        double current = 0;
        for (double x : out) current = std::max(current, std::abs(x));
        peak_ = std::max(peak_, current); return current;
    }
    double process(const float *input, std::size_t n) noexcept {
        for (std::size_t i = 0; i < n; ++i) process(input[i]);
        return peak_;
    }
    // File/end-of-stream only: evaluate all delayed tail contributions under
    // zero extension. Advances time by 11 frames; not for block boundaries.
    double flush() noexcept { for (int i = 0; i < 11; ++i) process(0); return peak_; }
};

} // namespace shores_dsp
