# Biquad baseline — 2026-10-24

Source: `tests/biquad_benchmark.cpp`; raw results: [biquad-baseline.csv](biquad-baseline.csv).

Machine: Intel Core Ultra 9 185H, Linux x86-64, GCC 16.2.1. Built with
`g++ -O3 -std=c++17 -Wall -Wextra -Werror` (no `-march=native`, fast-math, or
forced FMA). Each CSV row is the median of nine timed trials processing
20 seconds of 48 kHz noise in 256-sample blocks. The `fixed` mode processes
one, two, or eight cascaded 2.4 kHz/+3 dB peak filters. The `ramp` mode
alternates between +3 and -3 dB peak coefficients every block and ramps
throughout the following block. Input preparation, coefficient design, and
checksum calculation are outside the timed region. Filter processing and
per-block ramp scheduling are inside it.

Across three runs, one fixed filter took 2.273–2.371 ns/sample, equivalent to
about 0.011% of one core's 48 kHz real-time budget. Eight cascaded filters
took 17.732–18.053 ns/sample, about 0.085–0.087% of that budget. Results measure
the library processing loop, not plugin, host, UI, or audio-device overhead.
Real-time deadlines and performance on other CPUs require separate testing.
