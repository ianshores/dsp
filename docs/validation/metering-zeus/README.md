# Metering validation — Zeus, 24 October 2026

This milestone adds exactly the requested primitive set (with BiquadFilter
reusing Biquad) on local branch `feature/metering-primitives`. Base is completed
SVF commit `5e563031240f1ee7cd6c99f511f8541612fe58eb`; that branch is preserved.
The first basic-primitives commit is `ee27afc69af43e8aa50993b8324949d724fa9b2c`.
The final local receipt records subsequent commits. No code/docs push, release,
main merge, deployment or dependency-pin promotion belongs to this milestone.

## Tests

| Build | Language / warnings | Final CTest |
| --- | --- | --- |
| GCC 16.2.1, Release | C++17, Wall/Wextra/Werror | 12/12 pass |
| Clang 22.1.8, Release | C++17, Wall/Wextra/Werror | 12/12 pass |
| Clang 22.1.8, O1/debug, ASan + UBSan | C++17, Wall/Wextra/Werror | 12/12 pass |

Raw logs: [GCC](gcc-ctest.txt), [Clang](clang-ctest.txt),
[sanitizers](sanitizers-ctest.txt). Builds use CMake 4.4.4 and `--parallel 2`.
All seven original SVF/biquad/FFT entries pass, alongside four metering entries
and the optional official-fixture harness. Original DSP headers/tests/benchmarks
remain unchanged. Without external fixtures, the default suite has 11 entries.
GCC/Clang also compile the new example and benchmark with Werror.

- Basic tests use independent analytic 1/e attack/release, exact hold timing,
  dB decay, reference energy sums and brute-force moving means. Covered: silence,
  logical reset, variable blocks, invalid setup/nonfinite inputs, accumulator
  capacity freeze and large-peak removal from the moving window.
- K weighting is compared to an independent long-double direct-form-I cascade:
  maximum observed absolute error **6.1386e-8**. Explicit roles, LFE exclusion,
  permutations, bad roles, input sanitation and intermediate float overflow are
  checked. Variable blocks and exact in-place output agree after reset.
- Gate tests independently scan raw frame energies into overlapping full
  blocks and recompute both gates. Covered: partial blocks, exact first/hop
  boundaries, silence/reset, invalid input, fixed-capacity freeze, variable block
  invariance and complete/incomplete result flags.
- The allocation observer found **zero heap allocations** while processing and
  resetting all primitives over 96,000 frames. Setup is excluded. This is not a
  proof of thread safety, scheduling bounds or absence of every allocator type;
  relevant sources are also inspected for no dynamic allocation/locks/I/O.
- True peak agrees with independent long-double, explicitly zero-stuffed
  full-rate convolution (observed maximum difference **0** for the test vector).
  End-of-stream flush, phase ordering, reset, silence, invalid input, block
  invariance, >0 dBFS values and rate rejection are covered.

## Official loudness fixtures

29 unmodified WAVs from the [official ITU BS.2217 service](https://www.itu.int/oth/R1102000001/en)
were supplied externally to the test adapter, at 48 kHz PCM16. The harness uses
the reference expected readings in the [current BS.2217-2 report](https://www.itu.int/pub/R-REP-BS.2217-2-2016).
All 29 pass ±0.1 LU; LFE-only cases return −infinity. Largest absolute deviation
is **0.0484443932 LU**, on the absolute gate file. This covers mono frequency
sweep, 12 stereo frequency/level tones, 12 six-channel checks, two summing files
and both gate files. Advanced channel layouts and programme conformance items
were not downloaded or tested.

[Per-file results and WAV hashes](itu-fixtures-gcc.json) were captured with the
GCC standalone test adapter; the same 29-file harness passes in each CMake matrix
build. [Download URLs/archive hashes](itu-download-provenance.json) preserve
provenance. WAV/ZIP assets remain outside Git under
`/tmp/audio-metering-reference/itu-fixtures/`; no third-party audio is bundled.
The fixture harness does not download anything during CTest.

## Synthetic loudness and true-peak checks

Independently generated EBU Table 1 loudness cases 1–6 pass their ±0.1 LU
expectations. Integrated results respectively: **−22.9933, −32.9933, −23.0139,
−23.0139, −22.9787, −23.0163 LUFS**. Stationary M/S checks are also included
for cases 1/2. Gate energy/reference tests are separate from this composition.

Generated true-peak cases 15–19 yield **−6.21706, −5.97591, −6.31328,
−6.02649, +3.02907 dBTP**, within +0.2/−0.4 dB of the specified −6/+3
nominals. A 480 Hz–20 kHz finite tone/phase grid measured **−0.292675..+0.220440
dB** against analytic tone peaks. These are generated vectors, not official
EBU fixture files or a global reconstruction-error bound.

EBU fixture retrieval returned HTTP 403 at the public download attempt. Official
authentic programme cases 7/8, true-peak transient cases 20–23, and the full
M/S alignment sequences were **not run**, so their tolerances are unverified.
The [EBU terms](https://tech.ebu.ch/files/live/sites/tech/files/shared/testmaterial/use%20of%20EBU%20AUDIO%20test%20sequences.pdf)
restrict fixture use to internal R&D and prohibit redistribution. No access
restriction was bypassed, and none of that audio was adopted as project content.

## Scope and limits

See [API/reference/capacity documentation](../../metering.md) for exact supported
rates, role mapping, nonfinite/overflow policy and threading contracts.
Standards-related primitives currently support 48 kHz and six Annex 1 roles
only; basic timing supports 1–768 kHz within its frame/time caps. Gate history
is exact stored block energies, not a histogram, and fails its completeness flag
on exhaustion. Integrated querying is stopped-boundary work outside callbacks.
True peak uses published 4x/48-coefficient FIR, 5.875 input-sample group delay,
and requires the documented tail flush for file ends.

This does **not** establish full ITU/EBU Mode compliance. LRA, UI/scales,
start/pause operation, adapters and packaged plugins are outside the requested
table. No human audio acceptance, plugin-host test or subjective preference is
claimed. Remaining reference gaps are explicit above.

[Headless composition output](composition-example.txt), [performance method and
measurements](../../performance/metering-zeus/README.md), [environment](environment.json),
and [exact source digests](source-sha256.txt) make the milestone reviewable.
