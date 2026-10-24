# Metering primitives (local development)

Header-only C++17 primitives in `dsp/metering/`. `BiquadFilter` is an alias of
the existing `Biquad`, retaining its coefficient/state and synchronization API.
No host, UI, adapter, Heavy or runtime dependency is introduced.

## Basic stage

Include `dsp/metering/basic.hpp`. All processing/reset functions are `noexcept`,
allocate no memory and perform no locks or I/O. Work is constant per sample;
block calls are linear in the caller's frame count. Caller owns each instance
and all synchronization; no concurrent read/write is supported. Null pointers
are allowed only for zero-length blocks. Envelope output buffers are double;
window input/output are double energy/means. Exact same-type in-place window
buffers work; partial overlap is unsupported.

| Primitive | Semantics and defaults |
| --- | --- |
| `PeakDetector` | Maximum absolute float sample since reset; values above 1 allowed. |
| `EnvelopeFollower` | Rectified-amplitude one-pole; 10 ms attack/300 ms release at 48 kHz. Times are 1/e constants; zero is instantaneous. Not loudness smoothing. |
| `PeakHold` | Nonnegative detected amplitude; 1 s hold, 20 dB/s exponential decay at 48 kHz. A value >= held value refreshes the hold for the following N samples. Zero decay holds indefinitely. |
| `MeanSquareAccumulator` | Running mean of squared float audio; double mean, uint64 count. `rms()` takes its square root. No fixed window. |
| `SlidingWindow` | Rectangular nonnegative-energy mean, 400 ms at 48 kHz. `ready()` requires a full window; startup mean uses available frames. |

Rates for basic timing are 1000–768000 Hz. Ballistic times are 0–60 s;
hold duration is rounded to the nearest sample; decay is 0–1000 dB/s.
Configuration/constructors run outside callbacks and reject invalid controls.
Changing envelope ballistics preserves state. Envelope/hold tails below 1e-30
are cleared to avoid denormals.

Window setup allocates nearest-sample/half-up `rate*seconds` doubles, 1–2,304,000
frames and at most 10 s (~17.58 MiB maximum). Three seconds fits at every
supported rate through 768 kHz; longer windows remain capacity-limited. For
400 ms specifically, use exact rational `round_half_up(2*rate/5)`, not an
approximate binary `.4` multiplication. `sample_rate()` and `duration_seconds()`
report the configured rate and actual quantized duration. `configure(rate,time)`
allocates at setup and starts a fresh epoch; invalid setup leaves state intact. Push and logical reset are O(1); stale storage
is overwritten before it can contribute. A compensated double sum reduces
cancellation as old values leave, with normal binary64 rounding limitations.
Energy input domain is 0–1e100, avoiding overflow of any supported sum.

Mean-square capacity defaults to 2^32 samples (~24.9 h at 48 kHz), configurable
1–2^53. When full, extra input sets `overflowed()` and freezes the prefix until
reset; no overwrite/decimation/implicit restart. `valid()` indicates input
validity; check overflow separately. No count wrap or growing allocation.

Nonfinite audio or out-of-domain magnitude/energy becomes zero for that sample
and sets sticky `valid()==false`; elapsed time still advances. PeakHold accepts
0–1e50. Reset clears invalid/overflow status. Never describe an invalid or
overflowed measurement as a complete valid result. The raw reused Biquad keeps
its existing finite-input contract; it does not inherit these sanitation rules.

`tests/metering_basic.cpp` checks independent analytic attack/release and decay,
reference sums/windows, block equivalence, silence/reset, extreme-value removal,
nonfinite input and setup/capacity boundaries.

## Loudness stage

Include `dsp/metering/loudness.hpp` or the aggregate `dsp/metering.hpp`.
Rates are continuous finite binary64 values, including fractional/unusual rates:
**K weighting 4000–768000 Hz; gate/windows/true peak/basic timing 1000–768000 Hz**.
Values outside a component's range, NaN and infinity throw `invalid_argument`
during setup, without clamping. The 4 kHz K minimum keeps the 1682 Hz shelf
below Nyquist with a prewarp-conditioning margin; it does not promise a
20 kHz measurement bandwidth at low rates. This is numeric support, not a
certification of every rate. Official fixture evidence remains native 48 kHz.

`KWeightingFilter` uses the existing Biquad twice, with the exact decimal
coefficient sets in ITU-R BS.1770-5 Annex 1 Tables 1/2 **at exactly 48 kHz**.
For other rates `design_k_weighting(rate)` inverse-transforms each published
section into its analog bilinear polynomial, recovers natural frequency in Hz,
then prewarps that frequency at the new rate. With `u=(1-z^-1)/(1+z^-1)`,
the denominator polynomial is `(1+a1+a2) + 2(1-a2)u + (1-a1+a2)u²`.
Recover `g=sqrt(A0/A2)`, `f0=48000*atan(g)/pi`, then substitute
`u_old = [g/tan(pi*f0/rate)]*u_new` in numerator and denominator.
Normalize the constant digital denominator to one. This preserves damping,
numerator shape and endpoint gains, including the published RLB Nyquist gain
`4/(1+1.99004745483398+.99007225036621)`; it is deliberately not the alternative
raw `[1,-2,1]` numerator at every rate (which changes its endpoint gain).

This agrees with the physical-Hz shelf/RLB parameters in
[Brecht De Man's reference implementation](https://github.com/BrechtDeMan/loudness.py/blob/master/loudness.py)
to the tested numerical tolerance, after that explicit RLB gain convention.
The independent tests use its published frequency/Q/gain constants, not the
production inverse polynomial. No reference code or runtime dependency is
copied. Shelf: 1681.9744509555319 Hz, Q 0.7071752369554193, high gain
3.99984385397 dB, middle exponent 0.499666774155; RLB: 38.13547087613982 Hz,
Q 0.5003270373253953. This is a two-section bilinear design, not an exact
frequency-invariant resampling of the entire 48 kHz digital response; warping
is material near Nyquist at low rates. No generic RBJ shelf or SVF substitution.
`configure(rate)` designs outside callbacks and resets the epoch on success;
invalid setup preserves states/rate/status. Each channel has separate
state; there is no SVF substitution. Float intermediate output reflects the
existing biquad API; the coefficients and states remain double. An intermediate
float overflow clears the cascade and marks input invalid. Input NaN/infinity
becomes zero, with sticky invalidity as in the basics. Processing is bounded per
sample and in-place capable; reset clears states/status but retains coefficients.

`ChannelWeighting` requires an explicit ordered list of unique roles, maximum
six: Left, Right, Centre, LeftSurround, RightSurround, Lfe. Front energy weights
are 1, surrounds 1.41, LFE 0. It sums squared/mean-square energies, never audio
samples, and never infers dual-mono or layout from a channel count. Default is
one Centre channel. Subsets/permutations are allowed; caller owns the mapping.
No height/object/scene or Annex 3 layouts are supported. Check its sticky status
even for ignored LFE input. Channel energies must each lie in 0–1e100.

`LoudnessGate` takes **one already K-weighted, channel-summed energy per frame**.
A block has `N=round_half_up(2*rate/5)` frames, as BS.1770 specifies a
400 ms duration rounded to the nearest sample. Block j starts at
`S_j=round_half_up(j*rate/10)` and ends after frame `S_j+N`. Half ties go
upwards. The standard's 75% overlap is exactly representable only at some
rates; at others hop/overlap counts vary by one frame. Independently rounding
starts bounds phase error to half a frame rather than repeating a rounded hop
and accumulating drift. This scheduling choice extends the standard's nearest-
sample duration rule; the standard does not prescribe half-tie direction or an
exact fractional-start scheduler.

`rate.hpp` decomposes the supplied binary64 rate into an exact integer ratio.
An integer quotient/remainder phase accumulator generates successive differences
of S_j, without an ever-growing elapsed-frame counter or callback dependence.
At 48 kHz the first end is 19200 and subsequent hops are exactly 4800; the old
four-hop summation arithmetic is retained at all integral-hop rates. Fractional-
hop rates use an N-frame energy ring (at most 307200 doubles, 2.344 MiB), with
compensated add/remove sums. Partial first/final blocks are not used. Processing
and reset are O(1) per frame. There is no pre-gate attack/release follower.
`result()` recomputes the absolute and relative gate from all stored blocks,
using strictly greater energy tests at −70 LUFS and 10 LU below the absolute-
gated mean. Levels use `−0.691 + 10 log10(weighted mean square)`. Silence and
no eligible blocks return negative infinity; `stored_blocks==0` distinguishes
startup from a completed silent measurement.

Gate setup allocates a fixed history array and the fractional-hop ring if needed.
`configure(rate)` preserves history capacity and resets the epoch; invalid setup
preserves the previous configuration/state. `block_frames()`, `hop_floor_frames()`
and `integral_hops()` expose the timing choices. At the default history capacity
the fractional path needs less than 2.62 MiB total. The history is a fixed array of doubles. Default 36,000 blocks covers
3600.3 s at 48 kHz; maximum is 1,000,000 blocks (~8 MiB, 100000.3 s).
At the next completed block beyond capacity, the prefix freezes and sticky
`overflowed()` becomes true. Stream ingestion continues with bounded work;
`result().complete` is false. Invalid input also makes `complete` false.
The diagnostic prefix result is not a valid full-program result. Reset logically
clears history in O(1). There is no dropped-oldest history or histogram binning.

`result()` is an **analysis-thread/boundary query**, not a callback operation:
two bounded passes through stored blocks (at most 2,000,000 iterations), with
no allocation. Caller must stop/synchronize ingestion before querying; concurrent
reads/writes are unsupported. The library does not install locks or a publication
mechanism. `complete` indicates valid/exhaustion-free history, not EBU Mode
certification or readiness before the first full block.

`SlidingWindow(48000, .4)` and `(48000, 3)` provide ungated momentary and
short-term energy means for composition. Only report full-window readings when
`ready()` is true. Basic EnvelopeFollower is independent ballistics, not an
extra smoothing stage on these loudness windows.

## True peak stage

At **exactly 48 kHz**, `TruePeakOversampler` retains
the published BS.1770-5 Annex 2 four-phase FIR: **4x**, 12 taps/phase, 48 total
coefficients, double state/MACs, no callback allocations. Group delay is **5.875 input
samples (122.396 microseconds at 48 kHz)**; finite causal support is 12 input
samples. Floating point removes the fixed-point headroom shifts; no limiter or
clamp to full scale is applied. Use one instance for each measured channel,
including LFE when reporting that channel's raw true peak; loudness LFE weights
do not apply to this separate amplitude detector.

At every other supported rate the setup designer builds a **64 taps/phase**
Blackman-windowed sinc interpolator, **F=max(2,ceil(384000/rate))** phases.
The ceiling uses the exact rate ratio; 48 kHz alone retains its 4x table. Examples:
44.1/88.2/96/176.4/192 kHz use 9/5/4/3/2 phases. Output grids are at least
384 kHz, with at least 2x even above 192 kHz; this is an engineering accuracy
choice, not a claim that output-grid rate alone establishes conformance.

For tap k=0..63 and phase p=0..F-1, `t=k-32+p/F` and
`h=sinc(t)*(.42+.5*cos(pi*t/32)+.08*cos(2*pi*t/32))`, normalized to unity DC
independently per phase. Input-normalized bandwidth is rate-independent;
physical engineering validation band is **min(20 kHz, .45*rate)**. There is no
extra 20 kHz brick-wall input filter: callers must provide appropriately
bandlimited audio. No accuracy guarantee for higher-band signals, arbitrary
transients or untested rates follows. The wider transition/ripple near Nyquist
is a finite-FIR limit. This is an original generated alternative, not an ITU
published coefficient table for these rates.

Generated-path nominal group delay is **32 input frames** (`32/rate` seconds),
64-frame finite support, **63 zero frames** to flush. Per-phase normalization
and finite sampling mean this is nominal rather than an exact linear-phase
guarantee across every interpolated phase. Maximum F=384 at 1 kHz: 24576
coefficient doubles (192 KiB), 64 history doubles and a bounded 384-double
processing scratch array. Setup allocates/designs coefficients; process is
64*F double MACs plus F comparisons per input, with no allocation/trig/locks/I/O.
At 48 kHz it remains 12*4 MACs. Factor changes discontinuously at ceil thresholds
and at the preserved 48 kHz special case; configure is an epoch reset, not smooth
sample-rate modulation. Low rates cost appreciably more per input frame.

`configure(rate)` performs setup, resets measurement/signal state on success,
and preserves the old epoch if validation/allocation fails. Getters expose
`sample_rate()`, `factor()`, `taps_per_phase()`, `latency_frames()`,
`flush_frames()`, and `validation_bandwidth_hz()`. No lookahead buffer is hidden.

`process(sample)` returns the maximum of F new interpolated points.
`value()`/`dbtp()` retain the peak since reset. `flush()` feeds `flush_frames()` zero frames
to include all delayed end-of-file contributions under zero extension. Use it
once at the stream end, never at ordinary callback boundaries. It advances time;
it is not a lookahead or sample-preserving adapter. Silence is −infinity dBTP.

This finite FIR estimates bandlimited reconstruction, with pass-band ripple and
finite interpolation-grid under-read. The analytic 4x sinusoidal worst-grid
under-read at normalized frequency 0.5 is about 0.688 dB before FIR error; this
is not a universal bound on arbitrary transients. A test grid of 480 Hz–20 kHz
tones/phases measured −0.292675..+0.220440 dB relative to analytic peak. That grid
does not establish a full-band tolerance. The new generated path has independent
convolution, analytic-tone, dense multisine and finite shifted-sinc tests across
a rate grid; see [rate validation](validation/metering-rates-zeus/README.md).
These are generated references, not official resampled fixtures. No infinite-sinc,
universal ±0.3 dB or blanket R128 accuracy claim is made.

## Reference validation and conformance boundary

Primary references checked 24 October 2026:

- [ITU BS.1770-5, currently in force](https://www.itu.int/rec/R-REC-BS.1770), Annexes 1/2 and numerical tables.
- [EBU Tech 3341 v4, November 2023](https://tech.ebu.ch/docs/tech/tech3341.pdf) and [R128 v5, November 2023](https://tech.ebu.ch/docs/r/r128.pdf).
- [Current ITU BS.2217-2 compliance report](https://www.itu.int/pub/R-REP-BS.2217-2-2016) and [official fixtures](https://www.itu.int/oth/R1102000001/en).

29 downloaded official ITU loudness fixtures were tested at 48 kHz: 12 stereo
frequency/level tones, 12 explicit 5.1 channel checks including LFE, two summing
tests, one mono sweep, absolute gate and relative gate. All applicable integrated
readings pass ±0.1 LU, with LFE-only returning −infinity. Source files remain
outside this repository at `/tmp/audio-metering-reference/itu-fixtures/` and
are not redistributed. Results/hashes, not audio, are in the validation evidence.

Independently generated EBU Table 1 cases 1–6 pass ±0.1 LU; basic stationary
M/S readings are included for cases 1/2. Generated true-peak tone cases 15–19
pass +0.2/−0.4 dB. The FIR is independently checked against zero-stuffed
long-double direct convolution, including delayed tail. Generated material is
**not** described as official EBU audio. The EBU fixture download returned HTTP
403; authentic programme cases 7/8, official transient true-peak cases 20–23,
full M/S minimum-requirements sequences, advanced-layout and remaining ITU
programme fixtures have not passed this validation because they were not run.
The [EBU terms](https://tech.ebu.ch/files/live/sites/tech/files/shared/testmaterial/use%20of%20EBU%20AUDIO%20test%20sequences.pdf)
limit their sequences to internal R&D and prohibit redistribution. None is
bundled, copied into MIT-licensed assets, or used for an audible demo.

This is a set of reusable primitives, **not an EBU Mode compliant meter**.
LRA, start/pause controls, UI/scales and adapters are outside the requested table.
No blanket ITU/EBU certification follows from the scoped passing tests.

## Build, reproduce and compose

```sh
uv tool run --from cmake==4.4.4 cmake -S . -B build-metering-gcc \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ \
  -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_EXTENSIONS=OFF \
  -DCMAKE_CXX_FLAGS='-Wall -Wextra -Werror' \
  -DDSP_BUILD_BENCHMARKS=ON -DDSP_BUILD_EXAMPLES=ON \
  -DDSP_METERING_ITU_FIXTURES=/tmp/audio-metering-reference/itu-fixtures
uv tool run --from cmake==4.4.4 cmake --build build-metering-gcc --parallel 2
uv tool run --from cmake==4.4.4 ctest --test-dir build-metering-gcc --output-on-failure
./build-metering-gcc/dsp_metering_example
./build-metering-gcc/dsp_metering_benchmark
```

Omit `DSP_METERING_ITU_FIXTURES` if licensed external files are unavailable;
CTest then has 11 entries instead of 12, and must not be reported as testing
the official fixtures. Python is used only by the optional fixture test, never
by library processing. The fixture driver expects 48 kHz PCM16LE stdin and
explicit channels 1/2/5/6 with documented role order; it is test-only I/O.

`examples/metering.cpp` composes the exact requested primitives on synthetic
stereo tone data without hardware/hosts. Its generation, printing, end-of-stream
flush and integrated query are outside a real host callback contract. Host
applications supply audio buffers and handle synchronization separately.

## Arbitrary-rate reproduction

`tests/metering_rates.cpp` covers standard and unusual/fractional rates, numerical
bounds, independently parameterized K responses and energy, one million exact
clock hops per rate, raw-window gate scans, channel/LFE mapping, setup changes,
invalid/silence/reset/block invariance and true-peak impulse delay/tail.
`metering_realtime` also instruments the fractional ring/generated FIR paths and
maximum 384-phase processing. Setup/reconfigure allocation is deliberately outside
that observation interval. No new runtime dependency.

[Rate-specific final evidence](validation/metering-rates-zeus/README.md) and
[measured setup/reconfigure/process costs](performance/metering-rates-zeus/README.md)
extend, rather than replace, the historical 48 kHz evidence. The available 29
official ITU files remain native/unmodified 48 kHz; the EBU HTTP 403 and unrun
official programme/transient/M/S/advanced-layout gaps remain open.
