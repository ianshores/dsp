# DSP library

MIT-licensed, header-only C++17 building blocks for real-time audio. The library
contains a multimode biquad and Convoy's custom real FFT. Host adapters, plugin
packaging and deployment belong in [Foundry](https://github.com/ianshores/foundry).

## Biquad filter

`dsp/biquad.hpp` provides one filter instance per audio channel. Modes are
low-pass, high-pass, band-pass, notch, peak, low shelf, high shelf, and all-pass.
The designer-facing parameters are frequency in Hz, Q (resonance or width), and
gain in dB for peaks and shelves. Shelves use a separate slope parameter instead
of Q. Band-pass has 0 dB gain at its center frequency.

```cpp
#include <dsp/biquad.hpp>

using namespace shores_dsp;
const auto high_pass = design_biquad(48000, {BiquadType::HighPass, 300});
const auto low_pass = design_biquad(48000, {BiquadType::LowPass, 4000});
Biquad high(high_pass), low(low_pass); // one pair per channel

// In the audio callback, with caller-owned buffers:
high.process(input, output, frame_count);
low.process(output, output, frame_count); // exact in-place processing is supported
```

Design coefficients on a control thread, then apply them at an audio event
boundary with `set_coefficients()` or `ramp_to(coefficients, sample_count)`.
The latter interpolates coefficients for a defined number of samples; a host
can crossfade whole filters for a particularly large mode change. The caller
controls communication between threads. The fixed-setting processing path uses
five multiplies and two state values per sample, with no allocation or lock.
Double-precision coefficients and state preserve stability at low cutoffs and
high sample rates while input and output remain float. Inaudible state is
cleared at block boundaries to avoid denormal slowdowns.

Coefficient design follows the [W3C Audio EQ Cookbook](https://www.w3.org/TR/audio-eq-cookbook/).
Finite controls are limited to 10 Hz–0.49 × sample rate, Q 0.1–18, gain
−36–36 dB, and shelf slope 0.1–1. Sample rates must be 1–768 kHz. Invalid
rates and nonfinite controls are rejected. `tests/biquad.cpp` checks the
frequency-response landmarks,
independent double-precision processing, automation, in-place processing, and
the safe parameter range.

**Measured on Intel Core Ultra 9 185H, Linux, GCC 16.2.1**, `-O3` C++17,
48 kHz, 256-sample blocks, noise input, nine 20-second-audio trials per run:
one fixed peak filter took about 2.2–2.4 ns/sample (roughly 0.011% of one
core's real-time budget), and eight cascaded filters took about 17–18 ns/sample.
These are DSP-only measurements on this machine, not plugin or host CPU readings.
Run `dsp_biquad_benchmark` to measure fixed and per-block ramped operation on
your machine. [Method and raw baseline results](docs/performance/biquad-baseline.md).

## Real FFT

Include `dsp/fft.hpp` and construct `shores_dsp::RealFFT` during setup, outside the
audio callback. Sizes must be powers of two of at least 32 samples. Forward and
inverse processing allocate no memory or acquire locks. Plans own their scratch
buffers, so use one plan per concurrent processing stream.

```cpp
#include <dsp/fft.hpp>

shores_dsp::RealFFT fft(1024);
float time[1024] = {};
float real[513], imag[513], recovered[1024];
fft.forward(time, real, imag);
fft.inverse(real, imag, recovered);
```

The spectrum has `N/2 + 1` bins. This FFT stores logical bin `k` at
`fft.slot(k)`, a digit-reversed layout used to avoid reorder work in repeated
convolution. DC and Nyquist have zero imaginary component. The inverse is
normalized so a forward/inverse round trip has unit gain. The half-input
forward method reads only `N/2` time samples and treats the rest as zero.
Packed variants use an eight-real/eight-imaginary group layout in an `N + 16`
float buffer; see the header contract for details. These layouts are explicit
parts of the API so consumers do not accidentally treat the output as ordinary
natural-order spectrum data.

The implementation uses runtime AVX on supported GNU/Clang x86 builds, with
SSE and scalar fallbacks. Define `DSP_DISABLE_AVX` or `DSP_FFT_SCALAR` before
including the header to force a fallback for testing. `DSP_FFT_FMA` requires
an FMA compiler target and compatible CPU.

With CMake, add this repository as a subdirectory and link
`shores_dsp::shores_dsp`. Installed packages offer the same target via
`find_package(ShoresDSP CONFIG REQUIRED)`.

## Verify

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

The tests compare transforms against an independent double-precision FFT and
selected direct DFT bins, and check round trips over impulse, DC, Nyquist,
tone, noise and decaying signals from 32 through 16384 samples. Convoy's
convolution tests and benchmarks provide integration and performance checks.
