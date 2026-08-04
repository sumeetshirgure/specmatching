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
timing profile and the level-3 LER comparison. Every measured quantity is read straight out of it;
the C++ harness fills them from `pm::two_phase::summarize_ball`, and nothing is re-measured here.
The one derived number is the speedup: `stock_us / (m2_us * (blsm + hrvst))`, i.e. stock sparse
blossom against the M2 critical path, with everything outside blossom and harvest taken as
pipelined out.
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


def critical_ns(values):
    """Mean ns per shot on the critical path: blossom-on-H plus harvest, and nothing else.

    The other four buckets — ball intersect, H build, MWPM build and the unattributed remainder —
    are taken to be pipelined out: they either run ahead of the syndrome they serve or overlap the
    previous shot's matching, so they do not sit between a syndrome arriving and a correction
    leaving. What is left is blossom on H plus the harvest that reads the pairing back out.
    """
    return values["mean_total_ns"] * (values["frac_blossom_on_h"] + values["frac_harvest"])


def speedup_table(table):
    """The one speedup: stock sparse blossom against the M2 critical path.

    `speedup = stock_us / crit_us`, where `crit_us = m2_us * (blsm + hrvst)` — see `critical_ns`.
    `stock_us` is stock sparse blossom run to completion on `G`: no horizon, no residual. That is
    not a like-for-like job, since M2 stops at `T` and hands back a residual for M3-M6 to clean up,
    so this is an upper bound on the end-to-end win rather than the win itself. `m1_us` (Phase 1 on
    `G`, truncated the same way) and the fraction split are printed as context for it.

    Artifacts with no `mean_stock_ns` row print `-` for both `stock_us` and `speedup`.
    """
    print("\n## Stock sparse blossom vs the M2 critical path (blossom + harvest), us per shot\n")
    print(f"{'d':>4} {'p':>8} {'T':>5} {'m2_us':>8} {'m1_us':>8} {'stock_us':>9} {'crit_us':>8} "
          f"{'speedup':>8} {'h_nodes':>8} {'h_edges':>8} {'deg':>6} "
          f"{'isect':>6} {'hbld':>6} {'mwpm':>6} {'blsm':>6} {'hrvst':>6} {'unattr':>7}")
    for d, p, t, values in sorted(rows_of(table, "profile")):
        split = [values[k] for k in
                 ("frac_intersect", "frac_h_build", "frac_mwpm_build", "frac_blossom_on_h", "frac_harvest")]
        critical = critical_ns(values)
        total = values["mean_total_ns"]
        stock = values.get("mean_stock_ns")
        stock_us = f"{us(stock):>9.2f}" if stock else f"{'-':>9}"
        speedup = f"{stock / critical:>8.2f}" if stock and critical else f"{'-':>8}"
        print(f"{d:>4} {p:>8.4f} {t:>5.2f} {us(total):>8.2f} "
              f"{us(values['mean_g_reference_ns']):>8.2f} {stock_us} {us(critical):>8.2f} "
              f"{speedup} "
              f"{values['mean_h_nodes']:>8.1f} {values['mean_h_edges']:>8.1f} "
              f"{values['mean_degree']:>6.2f} " + " ".join(f"{x:>6.3f}" for x in split) +
              f" {1 - sum(split):>7.3f}")


def pipeline_split_table(table):
    """The critical path and the pipelined-out work in absolute terms, split by bucket.

    `pipe_us` is everything `critical_ns` drops — intersect, H build, MWPM build and the
    unattributed remainder — so `crit_us + pipe_us` is `m2_us`. It is not free work, it is work
    assumed to be hidden; if the pipelining does not hold, it comes back onto the critical path.
    """
    print("\n## Critical path vs pipelined-out work, us per shot\n")
    print(f"{'d':>4} {'p':>8} {'T':>5} {'m2_us':>8} {'crit_us':>8} {'pipe_us':>8} "
          f"{'isect_us':>9} {'hbld_us':>8} {'mwpm_us':>8} {'unattr_us':>10} {'crit%':>7}")
    for d, p, t, values in sorted(rows_of(table, "profile")):
        total = values["mean_total_ns"]
        critical = critical_ns(values)
        unattributed = 1 - sum(values[k] for k in
                               ("frac_intersect", "frac_h_build", "frac_mwpm_build",
                                "frac_blossom_on_h", "frac_harvest"))
        print(f"{d:>4} {p:>8.4f} {t:>5.2f} {us(total):>8.2f} {us(critical):>8.2f} "
              f"{us(total - critical):>8.2f} "
              f"{us(total * values['frac_intersect']):>9.2f} "
              f"{us(total * values['frac_h_build']):>8.2f} "
              f"{us(total * values['frac_mwpm_build']):>8.2f} "
              f"{us(total * unattributed):>10.2f} "
              f"{(critical / total if total else 0):>7.3f}")


def structural_table(table):
    """The hardware-budget counters: bytes moved and elements written, not nanoseconds.

    `isect`, `hbld` and `mwpm_rebuild` are the three stages meant to move off the CPU, so their
    wall time measures this laptop's DRAM latency rather than the target architecture. These are
    the machine-independent quantities behind them, and they feed a latency model directly.

    `isect_B/def` divides by `h_nodes` — the post-preamble defects, which is exactly the
    intersection loop's trip count — and is the per-PE local-memory budget the M2 exit checkpoint
    asks to be quoted. `other_B` is what the *other* build mode's traversal would have read on the
    same shots, derived from the tables without running it: that is the SCAN/BITSET crossover.

    Artifacts predating the counters print `-` throughout.
    """
    print("\n## Structural counters (bytes moved, elements written), per shot\n")
    print(f"{'d':>4} {'p':>8} {'T':>5} {'mode':>7} {'h_nodes':>8} {'h_edges':>8} "
          f"{'isect_B':>10} {'B/def':>9} {'scan_B':>10} {'hit_B':>9} {'other_B':>10} "
          f"{'hbld_edges':>10} {'maxdeg':>6} {'mwpm_init':>9} {'max_isect_B':>11}")
    for d, p, t, values in sorted(rows_of(table, "profile")):
        if "isect_bytes" not in values:
            print(f"{d:>4} {p:>8.4f} {t:>5.2f} " + " ".join(f"{'-':>{w}}" for w in
                  (7, 8, 8, 10, 9, 10, 9, 10, 10, 6, 9, 11)))
            continue
        mode = "bitset" if values.get("isect_mode_bitset", 0) else "scan"
        print(f"{d:>4} {p:>8.4f} {t:>5.2f} {mode:>7} "
              f"{values['mean_h_nodes']:>8.1f} {values['mean_h_edges']:>8.1f} "
              f"{values['isect_bytes']:>10.1f} {values['isect_bytes_per_defect']:>9.1f} "
              f"{values['isect_scan_bytes']:>10.1f} {values['isect_hit_bytes']:>9.1f} "
              f"{values['isect_bytes_other_mode']:>10.1f} "
              f"{values['hbld_edges_written']:>10.2f} {values['h_max_degree']:>6.0f} "
              f"{values['mwpm_init_elements']:>9.2f} {values['max_isect_bytes']:>11.0f}")


def structural_checks_table(table):
    """The three read-offs the counters exist for, as ratios rather than raw numbers.

    - `scan_B/def` against `words x 8` (BITSET) or `mean|B| x 8` (SCAN) says whether the traversal
      reads each defect's window exactly once. A ratio far from 1 means it does not.
    - `scan_B/def` at fixed `d` should barely move with `p`; `hit_B/def` should grow with it.
    - `edges/node` against the reported mean degree says which convention the edge count is in:
      `hbld_edges_written` counts each undirected pair once, so it should sit near `deg / 2` plus
      the boundary edges.
    """
    rows = [r for r in rows_of(table, "profile") if "isect_bytes" in r[3]]
    if not rows:
        return
    ball = {d: values for d, _, _, values in rows_of(table, "ball")}
    print("\n## Structural counter read-offs\n")
    print(f"{'d':>4} {'p':>8} {'T':>5} {'mode':>7} {'scan_B/def':>10} {'window_B':>9} {'ratio':>6} "
          f"{'hit_B/def':>10} {'edges/node':>10} {'deg':>6} {'2*edges/node/deg':>17}")
    for d, p, t, values in sorted(rows):
        bitset = bool(values.get("isect_mode_bitset", 0))
        table_row = ball.get(d, {})
        window = (table_row.get("mean_ball_word_len", 0) * 8 if bitset
                  else table_row.get("mean_ball_size", 0) * 8)
        scan_per_defect = values["isect_scan_bytes_per_defect"]
        nodes = values["mean_h_nodes"]
        edges_per_node = values["hbld_edges_written"] / nodes if nodes else 0
        degree = values["mean_degree"]
        print(f"{d:>4} {p:>8.4f} {t:>5.2f} {'bitset' if bitset else 'scan':>7} "
              f"{scan_per_defect:>10.1f} {window:>9.1f} "
              f"{(scan_per_defect / window if window else 0):>6.2f} "
              f"{values['isect_hit_bytes_per_defect']:>10.1f} {edges_per_node:>10.2f} "
              f"{degree:>6.2f} "
              f"{(2 * edges_per_node / degree if degree else 0):>17.2f}")
    print("\n`window_B` is the ball section's node-average window, and the ball section is compiled "
          "once per `d`,\nat the highest `p` in the sweep. The tables themselves are recompiled per "
          "`(d, p)` — `R` is quoted in\nmedian edge weights and the median edge weight moves with "
          "`p` — so a `ratio` that steps with `p` is that\nmismatch, not the traversal: at fixed "
          "`(d, p)` the scan component is identical across horizons.")


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
    """The M2 exit checkpoint's read, applied to the critical-path speedup at the lowest p.

    The threshold is the checkpoint's 1.5; the quantity it is read against is now
    `stock_us / crit_us` rather than `speedup_vs_m1`.
    """
    profile = list(rows_of(table, "profile"))
    if not profile:
        return
    lowest_p = min(p for _, p, _, _ in profile)
    at_lowest = [v.get("mean_stock_ns", 0) / critical_ns(v)
                 for _, p, _, v in profile if p == lowest_p and critical_ns(v)]
    print(f"\n## Exit read\n")
    if not any(at_lowest):
        print(f"Lowest p in scope: {lowest_p}. No `mean_stock_ns` in this artifact — no read.")
        return
    best = max(at_lowest)
    print(f"Lowest p in scope: {lowest_p}. Best `stock_us / crit_us` there: {best:.2f}.")
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
    pipeline_split_table(table)
    structural_table(table)
    structural_checks_table(table)
    tie_table(table)
    ler_table(table)
    exit_read(table)
    print("\nAll wall times are microseconds per shot. The only speedup reported is "
          "`stock_us / crit_us`, with `crit_us = m2_us * (blsm + hrvst)`: intersect, H build, "
          "MWPM build and the unattributed remainder are taken as pipelined out. The structural "
          "counters are bytes and element writes per shot, and are machine-independent: they are "
          "what the off-CPU stages move, whatever the memory system underneath.")


if __name__ == "__main__":
    main()
