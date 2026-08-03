#!/usr/bin/env python3
# Copyright 2026 PyReMatching contributors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Renders the M1 exit artifact CSV as the tables the M1 go/no-go read is made from.

    python benchmarks/two_phase/summarize_m1_artifact.py \
        benchmarks/two_phase/results/m1_exit_artifact.csv [--error-rate 0.003]

Every derived quantity here is read straight out of the CSV, which the C++ harness fills from
`pm::two_phase::summarize`. Nothing is re-derived: that is the whole point of having one
implementation of the summary table.
"""

import argparse
import csv
import sys
from collections import defaultdict


def load(path):
    with open(path) as handle:
        return list(csv.DictReader(handle))


def fmt_us(ns):
    return f"{float(ns) / 1000:.1f}"


def sweep_table(rows, error_rate):
    """Residual, fallback rate and cost, per (d, T). The core of the exit artifact."""
    print(f"\n## Horizon sweep at p = {error_rate}\n")
    print(
        f"{'d':>3} {'T/edge':>7} {'q':>6} {'residual':>9} {'maxres':>7} {'density':>8} "
        f"{'tree/shot':>10} {'blossom%':>9} {'phase1':>8} {'harvest':>8} {'no-trunc':>9} "
        f"{'exact':>8} {'P1/exact':>9}"
    )
    for row in rows:
        if row["error_rate"] != error_rate:
            continue
        phase1 = float(row["phase1_only_ns"])
        exact = float(row["exact_mean_ns"])
        print(
            f"{row['distance']:>3} {float(row['horizon_multiple']):>7.2f} "
            f"{float(row['q_fallback_rate']):>6.3f} {float(row['mean_residual_size']):>9.2f} "
            f"{row['max_residual_size']:>7} {float(row['mean_residual_density']):>8.4f} "
            f"{float(row['committed_pairs_tree_per_shot']):>10.2f} "
            f"{float(row['exposed_root_blossom_rate']) * 100:>9.1f} "
            f"{fmt_us(phase1):>8} {fmt_us(row['harvest_ns']):>8} "
            f"{fmt_us(row['untruncated_pipeline_ns']):>9} {fmt_us(exact):>8} "
            f"{phase1 / exact if exact else 0:>9.3f}"
        )


def dual_table(rows, error_rate):
    """max_S y_S against the cluster weight-diameter, both in units of T."""
    print(f"\n## max_S y_S vs cluster weight-diameter at p = {error_rate}\n")
    print(f"{'d':>3} {'T/edge':>7} {'max y_S / T':>12} {'mean max y_S / T':>17} {'diameter / T':>13} {'clusters':>9} {'skipped':>8}")
    for row in rows:
        if row["error_rate"] != error_rate:
            continue
        horizon = float(row["horizon"])
        print(
            f"{row['distance']:>3} {float(row['horizon_multiple']):>7.2f} "
            f"{float(row['max_region_dual']) / horizon:>12.2f} "
            f"{float(row['mean_max_region_dual']) / horizon:>17.2f} "
            f"{float(row['max_cluster_weight_diameter']) / horizon:>13.2f} "
            f"{row['clusters_measured']:>9} {row['clusters_skipped_for_budget']:>8}"
        )


def residual_histogram(rows, error_rate, distance):
    print(f"\n## Residual size distribution at d = {distance}, p = {error_rate}\n")
    for row in rows:
        if row["error_rate"] != error_rate or row["distance"] != str(distance):
            continue
        bins = [int(v) for v in row["residual_hist"].split()]
        shots = float(row["shots"])
        head = " ".join(f"{i}:{100 * b / shots:.0f}%" for i, b in enumerate(bins[:8]) if b)
        overflow = bins[-1]
        print(f"  T/edge={float(row['horizon_multiple']):>5.2f}  {head}" + (f"  >={len(bins) - 1}:{overflow}" if overflow else ""))


def scaling_table(rows):
    """Does the residual stay bounded as d grows at fixed T? That is the locality claim."""
    print("\n## Residual vs distance, at fixed T (the locality read)\n")
    by_multiple = defaultdict(list)
    for row in rows:
        by_multiple[float(row["horizon_multiple"])].append(row)
    for multiple in sorted(by_multiple):
        for error_rate in sorted({r["error_rate"] for r in by_multiple[multiple]}):
            selected = [r for r in by_multiple[multiple] if r["error_rate"] == error_rate]
            selected.sort(key=lambda r: int(r["distance"]))
            residuals = " ".join(f"d{r['distance']}:{float(r['mean_residual_size']):.2f}" for r in selected)
            print(f"  T/edge={multiple:>5.2f}  p={error_rate:<7} {residuals}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", nargs="?", default="benchmarks/two_phase/results/m1_exit_artifact.csv")
    parser.add_argument("--error-rate", default="0.003")
    parser.add_argument("--distance", type=int, default=21)
    args = parser.parse_args()

    rows = load(args.csv)
    if not rows:
        print("no rows", file=sys.stderr)
        return 1

    sweep_table(rows, args.error_rate)
    dual_table(rows, args.error_rate)
    residual_histogram(rows, args.error_rate, args.distance)
    scaling_table(rows)
    print("\nAll times in microseconds per shot.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
