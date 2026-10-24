# Metering rate performance — Zeus, 24 October 2026

Measured timings, not deadline guarantees or general claims about faster DSP.
[Environment and exact compiler flags](environment.txt): Intel Core i9-9900K
@3.60 GHz, 8 cores/16 threads, Linux 7.1.9-arch1-2 x86-64, GCC 16.2.1.
C++17 Release `-O3 -DNDEBUG -Wall -Wextra -Werror`; no native targeting,
fast-math, forced FMA, affinity or CPU clock changes. Background/system load
was not measured. Benchmarks run sequentially, without this task's tests/builds
competing during their timing intervals. Results include timer/loop/checksum
cost and remain sensitive to cache, frequency and system activity.

## Sequential unchanged-workload 48 kHz comparison

[Before raw](before-48k.csv), [before summary](before-48k-summary.csv),
[after raw](after-48k.csv), [after summary](after-48k-summary.csv).
Baseline headers and unchanged `tests/metering_benchmark.cpp` came from exactly
published 59fe52d48ee94a587130899f4f78876c87389628 via git archive into
`/tmp/dsp-metering-before-rates`. Baseline compiled with the same GCC/flags.
The after executable uses the local new headers and the identical harness.
Baseline ran first, then after; no randomized/interleaved trial order.

Each row is the median of nine 20-second-audio trials, 48 kHz, blocks of 256,
deterministic MT19937(seed 41) noise in [-.5,.5]. Mono primitive tests, and
explicit stereo/5.1 compositions. Preparation/input generation excluded;
setup measured separately, process includes per-block checksum observation.
No end flush/analysis query included in callback process timings. Composition
includes all requested primitives and K uses the existing biquad twice.

| Processing | Channels | Before ns/frame | After ns/frame | Change |
| --- | ---: | ---: | ---: | ---: |
| K weighting | 1 | 9.273935 | 8.282214 | −10.69% |
| Gate | 1 | 4.974302 | 4.929012 | −0.91% |
| Published 4x true peak | 1 | 13.215925 | 14.227760 | +7.66% |
| Composition | 2 | 65.549436 | 69.743484 | +6.40% |
| Composition | 6 | **165.608529** | **186.324449** | **+12.51%** |

The six-channel composition regression is real in this saved measurement:
about 0.795% → 0.894% of one core at 48 kHz under this harness. Numerical
preservation of the 48 kHz coefficients/output paths is **not a performance
guarantee**. Dispatch/state/storage changes and measurement variability can
influence timings; this sequential experiment does not isolate their causes.
No claim follows that K became inherently faster, or that all hosts see these
percentages. All raw setup/update/query rows are retained in the summaries.

## Rate-aware workload

[Harness](../../../tests/metering_rates_benchmark.cpp), [raw 9 trials](rates.csv),
[median/min/max summary](rates-summary.csv), [summary script](summarize.py).
Mono, blocks of 256, deterministic seed 41 noise. Per-rate frame count is
`max(8192,min(384000,2*fs))`: normally two seconds of audio, floor of 8192 frames
for short low-rate timings, cap of 384000 frames at high rates. Exact frame counts
are present in every CSV row. No hidden universal 20-second workload claim.
Gate history is 100 blocks (not default 36000); this workload never overflows it.
Energy conversion x*x is included in gate/window timing. K writes a caller
256-float block; true peak retains a scalar peak. Per-block checksums and guarded GNU compiler observation barriers are inside
processing timings. Anti-optimization/timing claims apply only to the tested
GCC/Clang builds; other compilers use the existing no-op fallback. MSVC was not
built or timed. Query/flush excluded.

Window3s is filled to capacity with .01 energy **outside both timers**, so
process measures steady add/remove work, not just partial startup. Other
primitives start reset. Setup measures one constructor/allocation/design,
excluding destruction. Reconfiguration is separately timed over 32 alternating
rates fs and fs+.25 (fs−.25 at the upper bound), including epoch reset,
allocation/design/replacement and old-storage destruction. In particular the
48 kHz true-peak reconfiguration row mixes the published path with generated
48000.25 Hz; it is not the cost of repeated exact48 kHz setup. Gate update rows
likewise mix integral/fractional storage. Tiny constructors approach clock floor.

| Rate Hz | K ns/frame | TP ns/frame | Gate ns/frame | Full 3s window ns/frame |
| ---: | ---: | ---: | ---: | ---: |
| 1000 | — | 5494.541992 | 4.586548 | 9.376831 |
| 4000 | 7.824829 | 1338.522339 | 4.750000 | 9.531006 |
| 44100 | 9.204921 | 282.890941 | 5.102177 | 10.000998 |
| 48000 | 8.589063 | 14.177271 | 4.748083 | 9.565781 |
| 48000.25 | 8.393865 | 250.854719 | 8.593906 | 9.577000 |
| 88200 | 8.267874 | 218.385181 | 4.787137 | 9.580051 |
| 96000 | 7.991130 | 183.873057 | 4.739396 | 9.576432 |
| 176400 | 8.059141 | 189.644889 | 4.837171 | 9.604255 |
| 192000 | 9.501734 | 99.134612 | 4.799172 | 9.650385 |
| 768000 | 8.225956 | 99.735305 | 4.807471 | 9.651445 |

Generated TP has 64 taps/phase; F=384/96/9/8/5/4/3/2 as determined by
max(2,ceil(384000/fs)). Exactly 48 kHz remains 12 taps/phase, 4x. That filter-length/cost
choice causes a substantial cost discontinuity around 48 kHz; generated paths
must not be described as as-cheap-as the published short FIR. Odd/small SIMD
loop lengths also affect the non-monotonic timings. At 768kHz TP process here
is ~7.66% of one core per channel; input 1000 Hz costs ~0.55% despite its large
per-frame count. No subjective sound superiority or general accuracy claim.

| Rate Hz | TP setup ns/instance | TP reconfigure ns/event | 3s window setup ns | 3s window reconfigure ns |
| ---: | ---: | ---: | ---: | ---: |
| 1000 | 8538708.000000 | 8782499.750000 | 198.000000 | 366.281250 |
| 4000 | 2134927.000000 | 2220236.906250 | 1530.000000 | 1613.062500 |
| 44100 | 208191.000000 | 217987.218750 | 273167.000000 | 34209.531250 |
| 48000 | 128.000000 | 91788.062500 | 294656.000000 | 31225.625000 |
| 48000.25 | 192857.000000 | 190638.906250 | 222909.000000 | 46381.218750 |
| 96000 | 92448.000000 | 93108.468750 | 496309.000000 | 61992.093750 |
| 192000 | 46355.000000 | 45093.687500 | 269241.000000 | 96790.093750 |
| 768000 | 46539.000000 | 46875.656250 | 2683492.000000 | 784514.281250 |

The maximum 3 s window allocates 2304000 doubles (~17.58 MiB) during setup;
process/reset never allocates. K design costs, all rate-specific gate setup/
reconfigure rows and trial spreads are in rates-summary.csv. Generated-FIR
setup includes long-double trig and may take milliseconds at low rates; it
belongs outside callbacks. State changes intentionally start a new epoch.

## Reproduce

Configure Release with `DSP_BUILD_BENCHMARKS=ON`, build with `--parallel 2`.
Run `dsp_metering_benchmark` before/after under the same environment. Run
`dsp_metering_rates_benchmark > rates.csv`, then
`python3 summarize.py rates.csv > rates-summary.csv` from this evidence folder.
For the original schema use `../metering-zeus/summarize.py`.
Source/header digests are in the [validation manifest](../../validation/metering-rates-zeus/source-sha256.txt).
No plugin/device/GUI/host, wall-clock deadline or full concurrency measurement.
