"""Summarize finite trial rows; failed workloads never receive a timing."""
import csv
import statistics
import sys
from collections import defaultdict

groups = defaultdict(list)
with open(sys.argv[1], newline="") as source:
    for row in csv.DictReader(source):
        if row["status"] != "finite":
            raise SystemExit("nonfinite workload: " + row["engine"] + "/" + row["mode"])
        key = (row["engine"], row["mode"], row["stages"], row["unit"])
        groups[key].append(float(row["ns_per_unit"]))
writer = csv.writer(sys.stdout, lineterminator="\n")
writer.writerow(["engine", "mode", "stages", "unit", "trials", "median_ns", "min_ns", "max_ns", "percent_one_core_48khz"])
for (engine, mode, stages, unit), values in groups.items():
    median = statistics.median(values)
    writer.writerow([engine, mode, stages, unit, len(values), f"{median:.6f}",
                     f"{min(values):.6f}", f"{max(values):.6f}",
                     f"{median * 48000 * 1e-7:.6f}" if unit == "sample" else ""])
