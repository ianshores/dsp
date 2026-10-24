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

Window setup allocates `round(rate*seconds)` doubles, 1–1,048,576 frames and
at most 10 s (~8 MiB maximum). Push and logical reset are O(1); stale storage
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
Standards-related filtering and gating currently support **48,000 Hz only**.
Other rates throw during setup; no unvalidated resampling or approximate shelf
design is substituted. Basic primitives retain their wider timing ranges.

`KWeightingFilter` uses the existing Biquad twice, with the exact decimal
coefficient sets in ITU-R BS.1770-5 Annex 1 Tables 1/2. Each channel has separate
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
Four completed 100 ms hops produce the first full 400 ms block at frame 19200;
subsequent blocks end every 4800 frames, giving 75% overlap. Partial first/final
blocks are not used and callbacks do not define gate boundaries. Processing
and reset are O(1) per frame. There is no pre-gate attack/release follower.
`result()` recomputes the absolute and relative gate from all stored blocks,
using strictly greater energy tests at −70 LUFS and 10 LU below the absolute-
gated mean. Levels use `−0.691 + 10 log10(weighted mean square)`. Silence and
no eligible blocks return negative infinity; `stored_blocks==0` distinguishes
startup from a completed silent measurement.

Gate setup allocates a fixed array of doubles. Default 36,000 blocks covers
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

`TruePeakOversampler` currently accepts 48 kHz only. It independently implements
the published BS.1770-5 Annex 2 four-phase FIR: **4x**, 12 taps/phase, 48 total
coefficients, double state/MACs, no allocations. Group delay is **5.875 input
samples (122.396 microseconds at 48 kHz)**; finite causal support is 12 input
samples. Floating point removes the fixed-point headroom shifts; no limiter or
clamp to full scale is applied. Use one instance for each measured channel,
including LFE when reporting that channel's raw true peak; loudness LFE weights
do not apply to this separate amplitude detector.

`process(sample)` returns the maximum of four new interpolated points.
`value()`/`dbtp()` retain the peak since reset. `flush()` feeds 11 zero frames
to include all delayed end-of-file contributions under zero extension. Use it
once at the stream end, never at ordinary callback boundaries. It advances time;
it is not a lookahead or sample-preserving adapter. Silence is −infinity dBTP.

This finite FIR estimates bandlimited reconstruction, with pass-band ripple and
finite interpolation-grid under-read. The analytic 4x sinusoidal worst-grid
under-read at normalized frequency 0.5 is about 0.688 dB before FIR error; this
is not a universal bound on arbitrary transients. A test grid of 480 Hz–20 kHz
tones/phases measured −0.292675..+0.220440 dB relative to analytic peak. That grid
does not establish a full-band tolerance. No 44.1/96/192 kHz, 8x/16x factor,
infinite-sinc, or universal ±0.3 dB/R128 accuracy claim is made.

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
