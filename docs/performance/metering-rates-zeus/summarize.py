"""Rate benchmark summary; Python standard library only."""
import collections
import csv
import statistics
import sys

groups = collections.defaultdict(list)
with open(sys.argv[1], newline="") as source:
    for row in csv.DictReader(source):
        key = tuple(row[k] for k in ("primitive", "mode", "sample_rate", "channels", "block", "frames"))
        groups[key].append(float(row["ns_per_unit"]))
writer = csv.writer(sys.stdout, lineterminator="\n")
writer.writerow(["primitive", "mode", "sample_rate", "channels", "block", "frames", "trials", "median_ns", "min_ns", "max_ns"])
for key, values in groups.items():
    writer.writerow([*key, len(values), f"{statistics.median(values):.6f}", f"{min(values):.6f}", f"{max(values):.6f}"])
