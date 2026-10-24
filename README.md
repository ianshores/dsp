# DSP library

MIT-licensed, header-only C++17 building blocks for real-time audio. The first
module is Convoy's custom real FFT, extracted without changing its transform
algorithm. Filters and other reusable modules can join it here. Host adapters,
plugin packaging and deployment belong in [Foundry](https://github.com/ianshores/foundry).

## Real FFT

Include `dsp/fft.hpp` and construct `dsp::RealFFT` during setup, outside the
audio callback. Sizes must be powers of two of at least 32 samples. Forward and
inverse processing allocate no memory or acquire locks. Plans own their scratch
buffers, so use one plan per concurrent processing stream.

```cpp
#include <dsp/fft.hpp>

dsp::RealFFT fft(1024);
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
