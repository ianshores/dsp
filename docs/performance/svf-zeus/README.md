# SVF and biquad on Zeus — 24 October 2026

On an **Intel Core i9-9900K @ 3.60 GHz**, Linux 7.1.9-arch1-2 x86-64,
**GCC 16.2.1**, `-O3 -DNDEBUG -std=c++17 -Wall -Wextra -Werror`, the fixed
SVF was slower than the existing biquad in this workload. No native CPU target,
fast-math or forced FMA was used. These are local DSP-loop measurements, not
plugin CPU readings or worst-case callback deadlines.

All processing rows use **48 kHz, one channel, 256-sample callbacks**, 20 seconds
(960,000 frames) of deterministic noise per trial, nine trials per row, with one
or eight cascaded filters. A stage count of eight still means one channel, not
eight channels. Timing includes processing and the updates named in each row;
input preparation, allocation, design and checksums are excluded from processing
timings. Coefficient arrays are prepared beforehand. Every output sample is
checked for finiteness outside the timed region.

| Path | Stages | Biquad ns/sample | SVF ns/sample |
| --- | ---: | ---: | ---: |
| Fixed, 1200 Hz LP, Q 0.8 | 1 | 3.466 | 4.686 |
| Fixed, same filter | 8 | 27.250 | 36.913 |
| Prepared event every block | 1 | 3.476 | 4.687 |
| Prepared event every block | 8 | 27.282 | 36.986 |
| Prepared event every sample, matched API calls | 1 | 6.266 | 7.246 |
| Prepared event stream, SVF-specific API | 1 | — | 6.008 |

Each value is the median of nine trials; [raw trials](raw.csv) and
[summary with min/max](summary.csv) retain the variation. At 48 kHz the one-stage
fixed rows equal about **0.0166% (biquad)** and **0.0225% (SVF)** of one core's
sample budget (`ns/sample × 48000 × 1e-7`). This conversion excludes all host,
UI, device, event transport and other callback work.

Both engines use LP with identical cutoff/Q at each event and the same input.
Per-block controls sweep cutoff 200–8000 Hz at 2 Hz and Q 0.7–2 at 1.462 Hz;
updates are immediate, with no ramp. Matched per-sample paths call
`set_coefficients()` and `process(..., 1)` on both engines. Their 6 kHz cutoff
motion spans 800–1600 Hz; Q varies 0.5–0.8 at 4386 Hz. This modest per-sample
range was chosen after the aggressive workload below failed for the biquad.
The additional SVF stream row uses `process_modulated()` in 256-sample blocks,
so its cleanup/call overhead differs from the matched row. Coefficient-stream
reads and their memory bandwidth are part of these timings. It does not measure
the cost of delivering live parameters from another thread.

## Design and update costs

| Operation | Biquad ns/event | SVF ns/event |
| --- | ---: | ---: |
| Design and store prepared coefficients | 13.599 | 13.215 |
| Apply prepared coefficients only | 3.046 | 0.889 |

Same CPU/compiler/flags as above; median of nine trials of 960,000 events.
Design-and-store includes control reads, validation/clamping, coefficient math,
and stores into preallocated arrays. Apply-only cycles through 4096 prepared
sets and includes reads, filter stores, loop and compiler barriers that prevent
dead-store elimination. It performs no audio processing. SVF and biquad structs
have different sizes and update internals; these are complete API costs, not
isolated arithmetic latency. Neither operation's timing is a thread-safety or
deadline guarantee.

## Aggressive modulation finding

An initial timing attempt used a much wider 6 kHz cutoff sweep, 200–8000 Hz,
with Q 0.7–2 at 4386 Hz and immediate updates every sample. The biquad produced
nonfinite audio; the first benchmark harness threw on its nonfinite checksum
and aborted before flushing CSV output. No useful timing was retained from that
attempt. The harness was corrected to check every output sample and write an
explicit `nonfinite` result with `nan` timing instead of aborting.

[stress.csv](stress.csv) reproduces that workload with one trial per row. Its
biquad `prepared_sample_matched` row is nonfinite; both SVF sample-update paths
remain finite (tested peak about 0.0162 for input bounded by 0.01). Invalid
biquad timing is deliberately excluded. Static and per-block paths remain
finite in this run. This illustrates a limit of this particular direct biquad
coefficient-swap trajectory, not a proof about all biquads, all modulation or
perceptual quality. The original biquad and its tests were left unchanged.

## Reproduce

Use CMake 4.4.4 through the isolated tool environment on Zeus (there was no
`cmake` in PATH initially):

```sh
uv tool run --from cmake==4.4.4 cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_EXTENSIONS=OFF \
  -DCMAKE_CXX_FLAGS='-Wall -Wextra -Werror' -DDSP_BUILD_BENCHMARKS=ON
uv tool run --from cmake==4.4.4 cmake --build build --parallel 2
./build/dsp_svf_benchmark > /tmp/svf-trials.csv
python docs/performance/svf-zeus/summarize.py /tmp/svf-trials.csv
./build/dsp_svf_benchmark --stress > /tmp/svf-stress.csv
./build/dsp_biquad_benchmark
```

The unchanged existing peak-EQ benchmark was also rerun on Zeus:
[existing-biquad.csv](existing-biquad.csv). It measures a different filter and
gain-ramp workload, so it is contextual evidence rather than the matched LP
comparison above. The old Core Ultra 9 185H baseline was not reused as Zeus data.

There was no CPU affinity, clock/governor control or desktop interruption; other
load was not measured. The engines ran sequentially. No universal speed ranking
or artifact-free claim follows from these measurements. Exact environment and
source hashes are in the [validation evidence](../../validation/svf-zeus/README.md).
