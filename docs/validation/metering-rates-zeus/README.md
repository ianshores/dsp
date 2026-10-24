# Arbitrary-rate metering validation — Zeus, 24 October 2026

Base: `feature/metering-primitives` at published
`59fe52d48ee94a587130899f4f78876c87389628`. The next single DSP commit is local;
its exact ID is recorded in the project receipt outside this repository. Existing
48 kHz K tables, FIR table/recurrence, gate arithmetic and prior regressions are
preserved. No Convoy/Foundry pin, plugin adapter, deployment or runtime dependency.

## Final matrix

| Build | Settings | CTest |
| --- | --- | --- |
| GCC 16.2.1 Release | C++17, `-O3 -DNDEBUG -Wall -Wextra -Werror` | 13/13 pass |
| Clang 22.1.8 Release | same | 13/13 pass |
| Clang ASan+UBSan | C++17, `-O1 -g -Wall -Wextra -Werror -fsanitize=address,undefined -fno-omit-frame-pointer` | 13/13 pass |

[CMake build logs](build-gcc.log), [GCC CTest](ctest-gcc.log),
[Clang CTest](ctest-clang.log), [sanitizers](ctest-sanitize.log).
CMake 4.4.4, `--parallel 2`, Release benchmarks/examples enabled; sanitizers have
examples enabled and benchmarks disabled. No sanitizer findings. Thirteen entries
include all seven original DSP regressions, five metering entries and optional
external official fixtures. Default suite without fixtures has twelve entries.

## Bounds, math and capacity

Finite continuous binary64 rates: K 4000–768000 Hz, gate/windows/true peak/basic
timing 1000–768000 Hz. Invalid setup throws without clamping. K's lower bound
keeps the 1682 Hz shelf below Nyquist with conditioning margin; low rates cannot
provide a 20 kHz band. Only numerical/rate-grid support is claimed.

K uses the existing biquad twice. Recover each published analog bilinear
polynomial, its natural frequency in Hz and damping, then prewarp at the new
rate, preserving numerator shape and endpoint gains. Exact 48 kHz bypasses
redesign. In particular, the RLB published Nyquist gain is retained rather than
silently changing gain with rate through an unscaled `[1,-2,1]` numerator.
Independent math uses De Man's frequency/Q/gain constants and explicit endpoint
gain, not the production inverse-polynomial implementation. See
[API/reference equations](../../metering.md).

400 ms blocks use exact `N=round_half_up(2*fs/5)`. Starts use
`round_half_up(j*fs/10)`; ends are start+N. Integer quotient/remainder scheduling
uses the exact rational representation of fs, with no accumulating rounded-hop
or global elapsed-counter overflow. Half ties go upward. Fractional overlap
necessarily varies by a frame; this is the documented scheduling choice where
BS.1770's duration/75% overlap cannot both be exact integer samples. Integral
hops retain four-sum arithmetic; fractional hops use a fixed N-frame energy ring.

SlidingWindow cap is now **2304000 doubles**, so three seconds fits through
768 kHz. Maximum setup storage ~17.58 MiB/window; at most ten seconds subject to
that cap. Fractional gate ring adds at most 307200 doubles (~2.344 MiB), plus
fixed configured block history (36000 default, 1000000 maximum). Overflow still
freezes history and makes the result incomplete. Configuration resets the epoch;
invalid setup preserves it. Setup may allocate; callback ingestion/reset do not.

New true peak at rates other than exactly 48 kHz: F=max(2,ceil(384000/fs)),
64 Blackman-windowed sinc taps/phase, per-phase unity-DC normalization, double
MAC/state, nominal delay 32 input frames, end flush 63 zero frames. Exact 48 retains
4x/12 taps, delay 5.875 and flush 11. Engineering validation band min(20 kHz,.45fs),
not an automatic input bandlimit or universal error guarantee. Max factor 384,
24576 coefficient doubles, bounded 384-double scratch. No callbacks design FIRs.

## Independent results

[Full numerical output](rates-numeric.txt), original
[K/gate/generated loudness cases](loudness-numeric.txt),
[48 kHz FIR regressions](true_peak-numeric.txt), [allocation observation](realtime-numeric.txt).

Rates include 4/8/22.05/44.1/48/88.2/96/176.4/192/768 kHz, 4000.25,
4001.25, 44117.3, 48000.25, 96000.123, 349525.5 Hz; true peak additionally
1 kHz, 1000.25, 192000.125 and the immediate binary64 neighbours of 48 kHz.
All applicable components accept these without a rate whitelist.

- Independent K complex frequency responses (10 Hz through 20 kHz where below
  Nyquist): maximum absolute difference **2.91726e-9 dB**. Independent filtered
  997 Hz steady-tone energy: maximum difference **3.11699e-8 dB** on the grid.
  Old long-double DF1 48 kHz test remains unchanged apart from obsolete rate
  rejection cases.
- One million clock hops per rate: independent absolute-time long-double check
  bounds start phase to **0.50001 sample** (tolerance includes reference rounding).
  Explicit exact rational block half-ties and adjacent representable rates pass.
  Raw-frame independent gate scans differ by at most **3.55271e-15 LU**.
  Variable block ingestion is identical, with tested first/last boundaries,
  silence/reset, invalid status, capacity freeze and setup epoch changes.
- Generated FIR analytic tones across phase/frequency grid: error
  **−0.0698461..+0.0187400 dB**. Dense analytic multisine reference maximum
  absolute error **0.0363158 dB**. Finite shifted-sinc pulse maximum absolute
  peak error **0.0357592 dB**. These tests enforce 0.3 dB regression limits on
  their generated cases only. Pulse truncation is explicit; it is not an
  infinite-sinc or official transient test.
- Independent long-double direct generated-kernel convolution, full tail:
  maximum difference **2.22045e-16**. Exact 48 kHz direct zero-stuffed published
  FIR test remains intact. Impulse timing, delay getters, tail drain, variable
  blocks, invalid controls, reset/silence and rate reconfiguration pass.
- Explicit channel order, surround 1.41/LFE 0 semantics remain unchanged across
  rates. Fractional envelope ballistics match analytic smoothing. Three-second
  windows allocate successfully at the maximum rate; oversized windows reject.
- Allocation instrumentation observes zero heap allocations during processing
  and resets, including fractional gate/window, generated FIR and max 384-phase
  processing. Setup is excluded. Replacement new/delete use a compiler-guarded noinline macro/fallback to keep
  instrumentation boundaries and avoid GCC's mismatched-new/delete inlining
  false positive. The benchmark GNU observation barrier is likewise guarded,
  with a no-op fallback. GCC/Clang were rebuilt and all three matrices rerun
  after both portability fixes; no MSVC build is claimed. No warning suppression. This is not a scheduler/concurrency
  or hard real-time proof; source inspection also finds no callback locks/I/O.

## Official fixtures and remaining gaps

[All 29 per-file results and hashes](itu-fixtures-gcc.json) come from **unmodified
native 48 kHz official ITU loudness WAVs**, external to Git. Maximum absolute
deviation remains **0.0484443932 LU**, passing applicable ±0.1 LU; LFE-only −infinity.
The same 29-file adapter passes in all three matrices.
[Original download provenance](../metering-zeus/itu-download-provenance.json).

Generated EBU-described loudness 1–6 and true-peak tones 15–19 remain independently
synthesized tests, not official fixtures. No official files were resampled or
relabelled to certify other rates. EBU retrieval remains HTTP 403; official
programme 7/8, transient true-peak 20–23, full M/S minimum sequences, advanced
layouts and remaining ITU programme items are **unrun**. No blanket ITU/EBU Mode
or arbitrary-rate conformance claim; no human listening acceptance.

## Reproduce and compose

Use the CMake commands in [the main API guide](../../metering.md), with external
`DSP_METERING_ITU_FIXTURES=/tmp/audio-metering-reference/itu-fixtures` when available.
This machine uses `uv tool run --from cmake==4.4.4 cmake`/`ctest` because cmake
is not on its default PATH. C++17/extensions-off and warnings flags are recorded
in [environment metadata](../../performance/metering-rates-zeus/environment.txt).
[Source digest manifest](source-sha256.txt) identifies tested code/test files.

`./build-metering-gcc/dsp_metering_example 48000.25` composes five seconds of
synthetic stereo headlessly. [Recorded runs](composition-examples.txt) also cover
44.1/48/96/192/768 kHz with full M/S/I windows and complete history. No audio
hardware/UI is opened. This numerical example is neither audition material nor
a human sound acceptance. [Measured performance](../../performance/metering-rates-zeus/README.md)
includes setup/reconfiguration and full-window processing separately.
