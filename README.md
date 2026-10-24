# DSP library

MIT-licensed, header-only C++17 building blocks for real-time audio. The library
contains a state-variable filter, a multimode biquad, and Convoy's custom real FFT. Host adapters, plugin
packaging and deployment belong in [Foundry](https://github.com/ianshores/foundry).

## State-variable filter

`dsp/svf.hpp` provides `SvfType`, `SvfParameters`, `SvfCoefficients`,
`design_svf(rate, parameters)`, and one `Svf` per channel. Modes are low-pass,
high-pass, band-pass and notch. Cutoff/center frequency and Q use the same
conventions as the corresponding biquad modes. Band-pass is scaled by `1/Q`
for **0 dB center gain**, including at high Q. Low/high-pass have the usual
−3.01 dB corner at Q = 1/√2; higher Q can boost around the cutoff.

```cpp
#include <dsp/svf.hpp>
using namespace shores_dsp;
// Setup/control thread:
const auto c = design_svf(48000, {SvfType::LowPass, 3600, 0.7071067811865476});
Svf filter(c); // independent state per channel
// Audio thread, after safe publication by the caller:
filter.process(input, output, frame_count);
filter.set_coefficients(prepared_next_event); // preserves integrator state
```

This is an optimized trapezoidal two-integrator design based on
[Simper's nodal derivation](https://cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf).
Its state and coefficients are double, with float input/output. Processing is
`noexcept` and performs no allocation, lock, I/O, coefficient design or trig.
Exact in-place blocks are supported; partially overlapping buffers are not.
Default construction passes input through. `reset()` clears state and retains
the current coefficients; inaudible state is cleared at block boundaries.

For sample-rate control motion, prepare an immutable array off the audio thread:
`process_modulated(input, output, events, frame_count)` applies `events[i]`
**before** sample `i` and retains the last set for later fixed processing.
Smooth cutoff/Q before calling `design_svf()` for each desired event. The caller
owns event memory, delivery deadlines, and synchronization. `set_coefficients()`
supports prepared events at less frequent boundaries. Design costs and the memory
bandwidth of sample-rate event arrays are separate from fixed processing cost.
Do not arbitrarily interpolate the internal SVF coefficients.

Finite frequency/Q controls clamp to 10 Hz–0.49 × sample rate and Q 0.1–18.
Rates outside 1–768 kHz, unknown modes, and nonfinite controls are rejected during
design. Both paths require finite audio and coefficients returned by the designer
(or the default wire). On a sample-rate change, prepare new coefficients and apply
at a known boundary; reset/crossfade if needed for the application.

Rapid modulation is explicitly tested, including against an independent nodal
solver, but this is not a promise of stability or inaudibility for every possible
signal/control sequence. Abrupt cutoff/Q or mode changes can click. As the
[JUCE behavioral documentation](https://docs.juce.com/master/classjuce_1_1dsp_1_1StateVariableTPTFilter.html)
also explains, TPT filters can still need smoothing. No JUCE code or dependency
is used. Biquad remains useful for EQ/shelves and has its own state and API.

See the [Zeus validation report](docs/validation/svf-zeus/README.md),
[matched performance measurements](docs/performance/svf-zeus/README.md), and
[audition instructions](docs/validation/svf-zeus/audition.md).
The SVF is a DSP primitive; host/PlugData/Heavy adapters and packaged plugins
require separate implementation and verification.

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
Configure with `-DDSP_BUILD_BENCHMARKS=ON` and run `dsp_biquad_benchmark` to
measure fixed and per-block ramped operation on your machine.
[Method and raw baseline results](docs/performance/biquad-baseline.md).

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
