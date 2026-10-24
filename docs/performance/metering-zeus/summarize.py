"""Summarize the metering benchmark's raw per-trial CSV, with no dependencies."""
import collections
import csv
import statistics
import sys

groups = collections.defaultdict(list)
with open(sys.argv[1], newline="") as source:
    for row in csv.DictReader(source):
        groups[row["primitive"], row["mode"], row["channels"], row["unit"]].append(float(row["ns_per_unit"]))
writer = csv.writer(sys.stdout, lineterminator="\n")
writer.writerow(["primitive", "mode", "channels", "unit", "trials", "median_ns", "min_ns", "max_ns", "percent_one_core_48khz"])
for key, values in groups.items():
    median = statistics.median(values)
    writer.writerow([*key, len(values), f"{median:.6f}", f"{min(values):.6f}", f"{max(values):.6f}",
                     f"{median * 48000 * 1e-7:.6f}" if key[1] == "process" else ""])
