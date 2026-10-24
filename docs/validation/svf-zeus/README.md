# SVF validation on Zeus — 24 October 2026

The local feature branch adds a dependency-free C++17 header-only state-variable
filter while preserving the existing biquad and FFT sources and tests. The API
separates control-thread design from audio-thread application of prepared events.
It supports LP/HP/unity-center BP/notch, float I/O, double state/coefficients,
independent channel state, exact in-place processing, and variable blocks.

## Results

| Build | Warnings | CTest |
| --- | --- | --- |
| GCC 16.2.1, Release, C++17 | `-Wall -Wextra -Werror` | 7/7 pass |
| Clang 22.1.8, Release, C++17 | `-Wall -Wextra -Werror` | 7/7 pass |
| Clang 22.1.8, `-O1 -g`, C++17, ASan + UBSan | `-Wall -Wextra -Werror` | 7/7 pass |

The suite includes the existing biquad and all five FFT variants. Both Release
builds also compiled the benchmark and audition generator with warnings as
errors. No sanitizer finding was reported. Raw logs:
[GCC](gcc-ctest.txt), [Clang](clang-ctest.txt), [ASan/UBSan](sanitizers-ctest.txt).
The [environment record](environment.json) captures CPU, OS, compiler versions,
flags, sample rate, channel/stage counts and callback size for measurements.

`tests/svf.cpp` uses an independent long-double W3C/RBJ direct form I reference
for static output, and measures the actual float impulse response's complex
spectrum against that reference. It covers all four modes at sample rates
1000, 8000, 44100, 48000, 96000, 192000 and 768000 Hz; cutoffs at 10 Hz,
100 Hz, 1000 Hz and 0.49 × rate (with specified clamping); and Q 0.1,
1/√2, 2 and 18. Response measurements include DC, center, pass/stop bands,
Nyquist, notch rejection, phase and band-pass unity-center normalization.

Time-varying checks use a separate long-double Gaussian-elimination solver for
the original two nodal equations. It recomputes its equations from cutoff/Q;
it never reads the production coefficients or uses the optimized recurrence.
Tests exercise simultaneous log cutoff and Q motion at 0.5, 80 and 6000 Hz
cutoff modulation (Q modulation at 0.731 × those rates), alternating cutoff/Q
limits, real sample-rate transitions 48→96→8 kHz, and mode changes without
reset. Each trajectory covers 48,000 samples with deterministic noise, tone,
DC and an initial transient. Slow modulation covers one second, not a long
listening or endurance trial.

All three builds measured the same maxima: **2.97363e-8 static absolute error**,
**1.19048e-7 modulation absolute error**, and **3.21476 maximum stress output**
for the test signal. Assertions use `3e-6 × (1 + abs(reference))`; peak output
must also stay below the explicit test bound of 100. These are tested-signal
results, not universal gain bounds. Detailed SVF outputs:
[GCC](gcc-svf-detail.txt), [Clang](clang-svf-detail.txt),
[sanitizers](sanitizers-svf-detail.txt).

Additional checks cover invalid modes/rates/NaN/infinity, exact finite-control
clamping, default passthrough, reset reproducibility, fixed and modulated silence,
long high-Q impulse tails near frequency/rate limits, in-place operation, zero
counts with null pointers, and varied blocks (including 1, 7, 17, 63, 256, 1024).
Per-sample `set_coefficients()` agrees exactly with the prepared event stream
for these signals; the final event remains active afterward.

## Design references and limits

Checked 24 October 2026:

- [Andrew Simper's optimized trapezoidal SVF derivation](https://cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf), nodal equations and equivalent-integrator state updates, pp. 2–7.
- [W3C Audio EQ Cookbook](https://www.w3.org/TR/audio-eq-cookbook/), LP/HP/notch and constant-0-dB band-pass conventions plus the independent static reference.
- [JUCE TPT filter behavior](https://docs.juce.com/master/classjuce_1_1dsp_1_1StateVariableTPTFilter.html), suitability for fast modulation, possible need for smoothing, and differing raw band-pass normalization. No JUCE implementation was copied or linked.

LP/HP retain prototype Q; multiplying the internal band signal by `1/Q`
normalizes BP. Caller-side smoothing and synchronization remain necessary as
appropriate. Arbitrary coefficient edits, nonfinite input, partially overlapping
buffers, extreme signal amplitudes, concurrent mutation, every possible control
trajectory, and perceptual quality are not covered by the tests. Finite float
input near full numeric range can overflow a resonant float output. Block state
cleanup can discard tails below 1e-30. No allocation, lock, I/O or design occurs
in either processing loop, as verified by source inspection; this was not a
thread-race, allocator-interception or hard real-time scheduling experiment.

The [benchmark report](../../performance/svf-zeus/README.md) records a reproduced
nonfinite biquad trajectory and distinguishes it from usable timing runs.
[Audition material](audition.md) is provided for Ian's listening decision.
No final radio-engine choice, host adapter, public code push, plugin release
or deployment is part of this change.

## Reproduce builds

```sh
uv tool run --from cmake==4.4.4 cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_EXTENSIONS=OFF \
  -DCMAKE_CXX_FLAGS='-Wall -Wextra -Werror' \
  -DDSP_BUILD_BENCHMARKS=ON -DDSP_BUILD_EXAMPLES=ON
uv tool run --from cmake==4.4.4 cmake --build build --parallel 2
uv tool run --from cmake==4.4.4 ctest --test-dir build --output-on-failure
```

For Clang Release, use `-B build-clang -DCMAKE_CXX_COMPILER=clang++` with the
same options. For the sanitizer build, use `-B build-sanitize`, Clang,
`-DCMAKE_BUILD_TYPE=Debug`, and replace the flags with
`-O1 -g -Wall -Wextra -Werror -fsanitize=address,undefined -fno-omit-frame-pointer`;
also set `-DCMAKE_EXE_LINKER_FLAGS='-fsanitize=address,undefined'`. Each build used
`--parallel 2`. Sanitizer runtime defaults were retained.
