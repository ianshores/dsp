"""Validate locally supplied, unmodified official ITU BS.2217 loudness WAVs.

No downloads, fixture generation, redistribution, or runtime dependencies.
Usage: python tests/metering_itu_fixtures.py DRIVER FIXTURE_DIRECTORY
"""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import wave

driver, directory = sys.argv[1:]
results = []
for path in sorted(Path(directory).glob("*.wav")):
    name = path.name
    if not name.startswith("1770-2_Comp_"):
        continue
    expected = (None if "ChannelCheckLFE" in name else
                -69.5 if "AbsGateTest" in name else
                -10.0 if "RelGateTest" in name else
                -18.0 if "FrequencySweep" in name else
                -23.0 if "23LKFS" in name else
                -24.0 if "24LKFS" in name else "unknown")
    if expected == "unknown":
        raise SystemExit("Unmapped fixture: " + name)
    with wave.open(str(path)) as source:
        channels = source.getnchannels()
        assert source.getframerate() == 48000 and source.getsampwidth() == 2
        assert channels in (1, 2, 5, 6)
        frames = source.getnframes()
        pcm = source.readframes(frames)
    completed = subprocess.run([driver, str(channels)], input=pcm,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    got = json.loads(completed.stdout)
    actual = got["integrated_lufs"]
    passed = (actual is None if expected is None else
              actual is not None and math.isfinite(actual) and abs(actual - expected) <= 0.1)
    assert got["complete"]
    assert passed, (name, expected, got)
    results.append({"file": name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    "channels": channels, "frames": frames, "sample_rate": 48000,
                    "expected_lufs": expected, "tolerance_lu": 0.1,
                    "actual_lufs": actual, "passed": passed})
assert len(results) == 29, "Expected exactly the 29 supported Annex 1 fixture files"
print(json.dumps({"source": "https://www.itu.int/oth/R1102000001/en",
                  "scope": "29 supplied BS.2217 loudness fixtures only; excludes larger layouts, programme conformance items and true-peak fixtures",
                  "passed": len(results), "results": results}, indent=2))
