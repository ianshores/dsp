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
