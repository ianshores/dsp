# Audible comparison

Four generated WAV files are available locally in
`/home/ian/dev/code/dsp/build/audition-zeus/`:

1. `01-dry-synthetic.wav` — original harmonic/vowel-like and noise/fricative bursts.
2. `02-biquad.wav` — the existing biquad with prepared coefficient changes.
3. `03-svf.wav` — the new SVF with the same cutoff/Q at each event.
4. `04-left-biquad-right-svf.wav` — biquad on the left, SVF on the right for channel soloing.

Open files 02 and 03 in a player of your choice and compare from the same
timestamp at a comfortable level; file 01 gives the source. File 04 supports
soloing either channel. No playback or REAPER session was opened by the agent.
The material is synthetic, not intelligible recorded radio speech, and is not a
blind or loudness-matched perceptual evaluation.

| Time | Scene |
| --- | --- |
| 0–4 s | Static radio band: HP 300 Hz + LP 3.6 kHz, Q 1/√2 |
| 4–8 s | Slow LP sweep: 0.5 Hz cutoff motion, 300–8000 Hz; Q 0.8–4 |
| 8–12 s | Rapid BP motion: 80 Hz, 300–6000 Hz; Q 0.7–8, unity center BP |
| 12–16 s | Audio-rate LP motion: 2 kHz, 600–8000 Hz; Q 0.7–4 |

Both engines receive identical independently designed per-sample control events.
The biquad uses immediate coefficient swaps, not its coefficient-ramp API.
Scenes reset state at their boundaries; boundaries and fast changes can click.
No extra smoothing, limiter, saturation, impulse response or third-party asset
is included. This comparison does not select the final radio filter for Ian.

All files are **16 seconds, 48 kHz, PCM16**. They share one gain of about
**0.193106** to bring the overall raw peak of **4.14281** to 0.8 full scale.
Individual files/scenes were not normalized, preserving relative level behavior.
The raw A/B difference RMS is 0.136907 over the full clip; that is a numerical
description, not an audible-quality score. Decoded-file verification found zero
clipped samples; headers, lengths, peaks, per-scene RMS and SHA256 digests are
retained in [audition-files.json](audition-files.json). These measurements do not
claim the agent listened to the audio.

Reproduce without audio hardware or host software:

```sh
uv tool run --from cmake==4.4.4 cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_EXTENSIONS=OFF \
  -DDSP_BUILD_EXAMPLES=ON
uv tool run --from cmake==4.4.4 cmake --build build --parallel 2 --target dsp_svf_audition
./build/dsp_svf_audition /tmp/svf-audition-new
```

The generator refuses to overwrite existing files. Generated WAVs remain in
the ignored local build directory; source and digest evidence are committed on
the local feature branch. No public code push or audio publication occurred.
