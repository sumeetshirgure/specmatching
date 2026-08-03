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

"""Renders the M2 exit artifact CSV as the tables the M2 exit read is made from.

    python benchmarks/two_phase/summarize_m2_artifact.py \
        benchmarks/two_phase/results/m2_exit_artifact.csv

The CSV is long-format — `section,d,p,T,key,value` — so one file holds the ball-table cost, the
`speedup_vs_m1` profile and the level-3 LER comparison. Every derived quantity is read straight out
of it; the C++ harness fills them from `pm::two_phase::summarize_ball`, and nothing is re-derived
here.
"""

import argparse
import csv
from collections import defaultdict


def load(path):
    table = defaultdict(dict)
    with open(path) as handle:
        for row in csv.DictReader(handle):
            key = (row["section"], row["d"], row["p"], row["T"])
            table[key][row["key"]] = float(row["value"])
    return table


def rows_of(table, section):
    for (sec, d, p, t), values in table.items():
        if sec == section:
            yield int(d), float(p), float(t), values


def us(ns):
    """Nanoseconds per shot -> microseconds per shot. Every wall time here is a whole-shot total."""
    return ns / 1000


def us_per_round(ns, d):
    """Nanoseconds per shot -> microseconds per round of syndrome extraction.

    The harness generates `rotated_memory_x` with `CircuitGenParameters(distance, distance, ...)`,
    so a distance-d shot is d rounds: d = 13 decodes 13 rounds. Only the steady-state per-round cost
    is reported this way -- the front end (ball intersect + H build) is a once-per-shot cost, so
    dividing it by d would not mean anything.
    """
    return ns / 1000 / d


def ball_table(table):
    print("\n## Ball table cost (the local-memory budget)\n")
    print("`compile_s` is a one-time build cost, not a per-shot one.\n")
    print(f"{'d':>4} {'nodes':>8} {'entries':>10} {'MiB':>8} {'mean|B|':>9} {'max|B|':>8} "
          f"{'words':>7} {'ambig':>7} {'compile_s':>10}")
    for d, _, _, values in sorted(rows_of(table, "ball")):
        print(f"{d:>4} {values['nodes']:>8.0f} {values['entries']:>10.0f} "
              f"{values['bytes_total'] / 2**20:>8.2f} {values['mean_ball_size']:>9.1f} "
              f"{values['max_ball_size']:>8.0f} {values['mean_ball_word_len']:>7.1f} "
              f"{values['ambiguous_mask_pairs']:>7.0f} {values['compile_wall_seconds']:>10.2f}")


def net_of_frontend_ns(values):
    """Mean total ns per shot with the M2-only front end netted out.

    `intersect` (the ball table lookups) and `h_build` (assembling H) are the two costs M1 does not
    pay at all, so they are what separates "is the pipeline faster today" from "how much of the win
    is the front end eating".
    """
    total = values["mean_total_ns"]
    return total - total * values["frac_intersect"] - total * values["frac_h_build"]


def speedup_table(table):
    """Both readings of `speedup_vs_m1`, plus the un-truncated baseline.

    `speedup` is the measured one — whole-shot M2 against whole-shot M1. `net_spd` divides out the
    ball-intersect and H-build time first; it is a ceiling, not a measurement, since it is what the
    speedup would be if the front end were free. `net/rd` is that same netted total spread over the
    d rounds of a distance-d shot — the number that has to clear the 1 us per-round budget on
    hardware. The front end itself is a once-per-shot cost, so it is never divided by d.

    `stock_us` is stock sparse blossom run to completion on `G` — no horizon, no residual. It is
    context, not a like-for-like comparison: both truncated columns stop at `T` and hand back a
    residual that M3-M6 still have to clean up, so `spd_stk = stock_us / m2_us` is an upper bound on
    the end-to-end win, not the win itself. `m1_us` is the comparison that holds the job fixed.
    Older artifacts have no `mean_stock_ns` row and print `-` in both columns.
    """
    print("\n## Phase 1 on H vs Phase 1 on G (`speedup_vs_m1`), us per shot\n")
    print(f"{'d':>4} {'p':>8} {'T':>5} {'m2_us':>8} {'m1_us':>8} {'stock_us':>9} {'speedup':>8} "
          f"{'spd_stk':>8} {'net_spd':>8} "
          f"{'net/rd':>8} {'h_nodes':>8} {'h_edges':>8} {'deg':>6} "
          f"{'isect':>6} {'hbld':>6} {'mwpm':>6} {'blsm':>6} {'hrvst':>6} {'unattr':>7}")
    for d, p, t, values in sorted(rows_of(table, "profile")):
        split = [values[k] for k in
                 ("frac_intersect", "frac_h_build", "frac_mwpm_build", "frac_blossom_on_h", "frac_harvest")]
        net = net_of_frontend_ns(values)
        reference = values["mean_g_reference_ns"]
        total = values["mean_total_ns"]
        stock = values.get("mean_stock_ns")
        stock_us = f"{us(stock):>9.2f}" if stock else f"{'-':>9}"
        stock_speedup = f"{stock / total:>8.2f}" if stock and total else f"{'-':>8}"
        print(f"{d:>4} {p:>8.4f} {t:>5.2f} {us(total):>8.2f} "
              f"{us(reference):>8.2f} {stock_us} {values['speedup_vs_m1']:>8.2f} "
              f"{stock_speedup} "
              f"{(reference / net if net else 0):>8.2f} "
              f"{us_per_round(net, d):>8.3f} "
              f"{values['mean_h_nodes']:>8.1f} {values['mean_h_edges']:>8.1f} "
              f"{values['mean_degree']:>6.2f} " + " ".join(f"{x:>6.3f}" for x in split) +
              f" {1 - sum(split):>7.3f}")


def frontend_cost_table(table):
    """What the M2-only front end costs in absolute terms; the speedups are in `speedup_table`."""
    print("\n## Cost of the ball-intersect + H-build front end, us per shot\n")
    print(f"{'d':>4} {'p':>8} {'T':>5} {'m1_us':>8} {'m2_us':>8} {'isect_us':>9} {'hbld_us':>8} "
          f"{'net_us':>8} {'front%':>7}")
    for d, p, t, values in sorted(rows_of(table, "profile")):
        total = values["mean_total_ns"]
        print(f"{d:>4} {p:>8.4f} {t:>5.2f} {us(values['mean_g_reference_ns']):>8.2f} "
              f"{us(total):>8.2f} "
              f"{us(total * values['frac_intersect']):>9.2f} "
              f"{us(total * values['frac_h_build']):>8.2f} "
              f"{us(net_of_frontend_ns(values)):>8.2f} "
              f"{values['frac_intersect'] + values['frac_h_build']:>7.3f}")


def tie_table(table):
    """§M2.6 divergence rates. These are ties between optima, not errors — see the M2 notes."""
    print("\n## §M2.6 tie rates (H resolving a degeneracy differently from G)\n")
    print(f"{'d':>4} {'p':>8} {'T':>5} {'residual':>9} {'pairing':>8}")
    for d, p, t, values in sorted(rows_of(table, "profile")):
        print(f"{d:>4} {p:>8.4f} {t:>5.2f} {values['residual_tie_rate']:>9.4f} "
              f"{values['pairing_tie_rate']:>8.4f}")


def ler_table(table):
    print("\n## Level 3: end-to-end logical error rate at T = 2\n")
    print(f"{'d':>4} {'p':>8} {'shots':>8} {'exact':>11} {'M1+cleanup':>12} {'M2+cleanup':>12} {'M2==M1':>7}")
    for d, p, _, values in sorted(rows_of(table, "ler")):
        same = "yes" if values["ler_m1"] == values["ler_m2"] else "NO"
        print(f"{d:>4} {p:>8.4f} {values['shots']:>8.0f} {values['ler_exact']:>11.3e} "
              f"{values['ler_m1']:>12.3e} {values['ler_m2']:>12.3e} {same:>7}")

    print("\n## Level 3: slope of ln(LER) vs d (proportional to `d_eff/d`)\n")
    print(f"{'p':>8} {'exact':>20} {'M1':>20} {'M2':>20} {'M2-M1':>10}")
    for _, p, _, values in sorted(rows_of(table, "ler_fit")):
        delta = values["slope_m2"] - values["slope_m1"]
        print(f"{p:>8.4f} "
              f"{values['slope_exact']:>12.4f}+-{values['stderr_exact']:<7.4f} "
              f"{values['slope_m1']:>12.4f}+-{values['stderr_m1']:<7.4f} "
              f"{values['slope_m2']:>12.4f}+-{values['stderr_m2']:<7.4f} "
              f"{delta:>10.4f}")


def exit_read(table):
    """The M2 exit checkpoint's read: stop and profile if `speedup_vs_m1 < 1.5` at the lowest p."""
    profile = list(rows_of(table, "profile"))
    if not profile:
        return
    lowest_p = min(p for _, p, _, _ in profile)
    at_lowest = [(d, t, v["speedup_vs_m1"]) for d, p, t, v in profile if p == lowest_p]
    best = max(speedup for _, _, speedup in at_lowest)
    print(f"\n## Exit read\n")
    print(f"Lowest p in scope: {lowest_p}. Best `speedup_vs_m1` there: {best:.2f}.")
    if best < 1.5:
        print("=> Below the 1.5 threshold. The M2 exit checkpoint says: stop and profile before")
        print("   starting M3. The front end is not the bottleneck and compression will not rescue it.")
    else:
        print("=> At or above the 1.5 threshold; M3 may proceed.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path")
    args = parser.parse_args()

    table = load(args.csv_path)
    ball_table(table)
    speedup_table(table)
    frontend_cost_table(table)
    tie_table(table)
    ler_table(table)
    exit_read(table)
    print("\nAll wall times are microseconds per shot, except `net/rd`: a distance-d shot is d "
          "rounds, so net/rd divides the front-end-free total by d.")


if __name__ == "__main__":
    main()
