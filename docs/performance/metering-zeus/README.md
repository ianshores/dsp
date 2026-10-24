# Metering measurements on Zeus — 24 October 2026

**CPU:** Intel Core i9-9900K @ 3.60 GHz, 8 cores/16 threads.
**OS:** Linux 7.1.9-arch1-2 x86-64. **Compiler:** GCC 16.2.1.
**Flags:** `-Wall -Wextra -Werror -O3 -DNDEBUG -std=c++17`.
No native CPU target, fast-math or forced FMA. CMake 4.4.4 via isolated uv tool.

Each process row is the median of **nine trials**, each processing **20 seconds
(960,000 frames), 48 kHz, blocks of 256**. Mono except explicit six-channel
weighting and stereo/5.1 compositions. Random input is deterministic (seed 41,
uniform −0.5..0.5); squared-energy input for energy primitives is prepared
outside timing. No host, audio device or UI is involved. Per-block return-value
checksum accumulation and callback-equivalent loops are inside process timing.
Input generation, allocations, setup, result printing and destructor work are
outside process timing. Compiler barriers keep state observable without runtime
instructions on the tested compiler. No CPU affinity/clock changes or desktop
interruption; background load was not measured. Engines/rows ran sequentially.

| Primitive/path | Channels | Median ns/frame |
| --- | ---: | ---: |
| Absolute sample maximum | 1 | 1.027696 |
| Attack/release envelope | 1 | 3.657695 |
| Peak hold/decay | 1 | 1.439886 |
| Running mean square | 1 | 4.759391 |
| 400 ms sliding energy window | 1 | 9.466938 |
| Existing biquad, HP 300 Hz | 1 | 3.724448 |
| Exact K-weighting cascade | 1 | 8.201239 |
| Explicit 5.1 energy weighting | 6 | 3.286471 |
| Gate ingestion (query excluded) | 1 | 4.914388 |
| True peak, 4x/48-coefficient FIR | 1 | 11.828198 |
| Full requested-primitive composition | 2 | 63.588873 |
| Full requested-primitive composition | 6 | 164.496977 |

The six-channel weighting row includes assembling six identical energy inputs
per frame; its number is per frame, not per channel. The composition includes
sample peak, envelope, hold, mean square, true peak and K weighting independently
for each input channel, then channel energy sum, M/S windows and gate ingestion.
It feeds identical noise with alternating signs into each channel. This is a
representative composition workload, not an additional public meter class.
K weighting exercises the reused biquad twice. LFE contributes to per-channel
detectors and has zero integrated-loudness energy weight.

At 48 kHz, composition cost is about **0.3052% (stereo)** or **0.7896% (5.1)**
of one core (`ns/frame × 48000 × 1e-7`). This conversion excludes host/device,
UI, cross-thread transport and integrated query overhead and does not prove
worst-case callback deadlines or performance on other machines.

## Separate setup, update and query costs

| Operation | Unit | Median ns |
| --- | --- | ---: |
| Construct 400 ms window, allocate/initialize 19,200 doubles | instance | 4940 |
| Construct gate, allocate/initialize 36,000 doubles | instance | 8649 |
| Construct stereo composition | instance | 43144 |
| Construct 5.1 composition | instance | 44169 |
| Configure attack/release ballistics (two exp calls) | event | 11.735302 |
| Integrated query after 20 s (197 stored blocks) | query | 566 |

All medians are nine trials on the same environment. Setup is one construction
per trial, including allocation/default initialization and timer overhead;
destruction is excluded. Trivial-constructor setup readings are near the clock
overhead floor and should not be interpreted as precise instruction latency;
all are retained in raw data. Allocator/cache warmup and background activity
affect setup variation. Composition allocates M/S windows and a gate history
and contains six state slots, even for its stereo run. BiQuad coefficient
design is outside the setup timing for its fixed prepared-coefficient row.
K-weighting coefficients and true-peak FIR are fixed published constants.

Envelope updates are nine 960,000-event trials, including prepared energy reads,
two parameter additions, validation, two exponentials, stores and compiler
barriers. Query is one complete two-pass scan per trial, at a stopped analysis
boundary. It measures 197 blocks only, not worst-case million-block capacity.
Queries must not be called concurrently with processing or placed in callbacks.

## Evidence and reproduction

[Raw per-trial results](raw.csv), [median/min/max summary](summary.csv),
[environment record](../../validation/metering-zeus/environment.json), and
[source hashes](../../validation/metering-zeus/source-sha256.txt).

```sh
# Configure/build as documented in docs/metering.md, with benchmarks ON.
./build-metering-gcc/dsp_metering_benchmark > /tmp/metering-trials.csv
python docs/performance/metering-zeus/summarize.py /tmp/metering-trials.csv
```

Processing is bounded and source-inspected for no allocation/locks/I/O; a
separate test also observes zero heap allocations during 96,000 processing/reset
frames. These timing trials are not a hard real-time, concurrent publication,
compliance or subjective-quality test.
