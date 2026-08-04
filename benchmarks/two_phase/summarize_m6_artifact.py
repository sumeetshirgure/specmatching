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

"""Renders the M6 exit artifact CSV as the sections of §M6.4's benchmark list.

    python benchmarks/two_phase/summarize_m6_artifact.py \
        benchmarks/two_phase/results/m6_exit_artifact.csv

The CSV is long-format — `section,d,p,T,key,value`. The escalation and streaming budgets are the M3
artifact's job and are not repeated here.
"""

import argparse
import csv
from collections import defaultdict


def load(path):
    table = defaultdict(dict)
    heavy = defaultdict(dict)
    with open(path) as handle:
        for row in csv.DictReader(handle):
            key = (row["section"], int(row["d"]), float(row["p"]), float(row["T"]))
            if row["section"] == "heavy":
                name = row["key"]
                if "_" not in name.removeprefix("shot_"):
                    continue
                shot, _, field = name.removeprefix("shot_").partition("_")
                heavy[(key[1], key[2], key[3], int(shot))][field] = float(row["value"])
                continue
            table[key][row["key"]] = float(row["value"])
    return table, heavy


def points(table, section):
    return sorted(key for key in table if key[0] == section)


def report_latency(table):
    print("\n== latency vs stock exact ==")
    print(
        f"{'d':>4} {'p':>8} {'T':>5} {'mean':>10} {'p99':>10} {'p999':>10} {'p9999':>10} "
        f"{'max':>10} {'stock':>10} {'speedup':>8} {'contam':>9} {'gap':>7}"
    )
    for key in points(table, "latency"):
        row = table[key]
        print(
            f"{key[1]:>4} {key[2]:>8g} {key[3]:>5g} {row['mean_ns']:>10.1f} {row['p99_ns']:>10.1f} "
            f"{row['p999_ns']:>10.1f} {row['p9999_ns']:>10.1f} {row['max_ns']:>10.1f} "
            f"{row['stock_mean_ns']:>10.1f} {row['speedup_vs_stock']:>8.3f} "
            f"{100 * row['contaminated_shot_rate']:>8.4f}% {row['amortisation_gap']:>7.4f}"
        )
    print(
        "\n  `speedup` below 1 means the two-phase decoder is slower than stock exact on this"
        " machine.\n  `gap` is |amortised - measured| / measured, the M6 exit checkpoint's"
        " reconciliation."
    )


def report_frontend(table):
    print("\n== front-end comparison: stock exact vs M1 Phase 1 on G vs M2 Phase 1 on H ==")
    print(f"{'d':>4} {'p':>8} {'stock_ns':>11} {'m1_on_G_ns':>12} {'m2_on_H_ns':>12} {'m2 vs m1':>10}")
    for key in points(table, "frontend"):
        row = table[key]
        print(
            f"{key[1]:>4} {key[2]:>8g} {row['stock_mean_ns']:>11.1f} {row['m1_phase1_on_g_ns']:>12.1f} "
            f"{row['m2_phase1_on_h_ns']:>12.1f} {row['speedup_vs_m1']:>10.3f}"
        )
    print(
        "\n  Stock is a different job from the other two — it returns a full matching where Phase 1"
        "\n  returns a partial one — so it is a column, never a denominator of a Phase-1 ratio."
    )


def report_ball(table):
    print("\n== ball table cost and structural counters ==")
    print(
        f"{'d':>4} {'bytes':>12} {'paths_B':>12} {'mean|B|':>9} {'deg@2T':>8} "
        f"{'isect_B/def':>12} {'other_mode':>11} {'hbld_edges':>11} {'mwpm_init':>10}"
    )
    for key in points(table, "ball"):
        row = table[key]
        mode = "bitset" if row.get("mode_bitset") else "scan"
        print(
            f"{key[1]:>4} {row['bytes_total']:>12.0f} {row['bytes_paths']:>12.0f} "
            f"{row['mean_ball_size']:>9.1f} {row['mean_degree_at_2T']:>8.2f} "
            f"{row['isect_bytes_per_defect']:>12.1f} {row['isect_bytes_other_mode_per_defect']:>11.1f} "
            f"{row['hbld_edges_written']:>11.1f} {row['mwpm_init_elements']:>10.1f}  ({mode})"
        )
    print(
        "\n  `isect_B/def` is the local memory a per-defect processing element has to stream, which"
        " is the\n  headline number for the local-memory hardware argument. `other_mode` is what the"
        " other build\n  mode would have read on the same shots, derived from the tables rather than"
        " measured."
    )


def report_ler(table):
    print("\n== effective distance ==")
    print(f"{'p':>8} {'d_points':>9} {'slope_two_phase':>18} {'slope_exact':>18} {'disagreements':>14}")
    for key in points(table, "ler"):
        if key[1] != 0:
            continue  # per-distance LER rows; the fit rows carry d = 0
        row = table[key]
        if not row.get("fit_points"):
            continue
        print(
            f"{key[2]:>8g} {row['fit_points']:>9.0f} "
            f"{row['slope_two_phase']:>11.4f}+-{row['slope_two_phase_stderr']:<6.4f} "
            f"{row['slope_exact']:>11.4f}+-{row['slope_exact_stderr']:<6.4f} "
            f"{row.get('prediction_disagreements', float('nan')):>14.0f}"
        )
    print(
        "\n  `disagreements` is the number of shots where the two decoders predicted differently."
        "\n  It is zero by construction — the output is exact MWPM on every shot — and it is the"
        " number to\n  read. The slope fit can only be a weaker statistical restatement of it, and"
        " at p <= 1e-3 it has\n  no signal at all: every LER in the campaign is zero at these shot"
        " counts."
    )


def report_heavy(heavy):
    if not heavy:
        return
    print("\n== heavy-shot drill-down (top 0.1% by defect count, plus every escalating shot) ==")
    print(f"{'d':>4} {'p':>8} {'shot':>7} {'defects':>8} {'phase1':>10} {'harvest':>9} {'escal':>10} {'total':>10} {'stock':>10}")
    rows = sorted(heavy.items())
    for (d, p, _t, shot), fields in rows[:40]:
        print(
            f"{d:>4} {p:>8g} {shot:>7} {fields.get('num_defects', 0):>8.0f} "
            f"{fields.get('phase1_ns', 0):>10.0f} {fields.get('harvest_ns', 0):>9.0f} "
            f"{fields.get('escalation_ns', 0):>10.0f} {fields.get('total_ns', 0):>10.0f} "
            f"{fields.get('exact_reference_ns', 0):>10.0f}"
            + ("  escalated" if fields.get("escalated") else "")
        )
    if len(rows) > 40:
        print(f"  ... {len(rows) - 40} more rows in the CSV")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path")
    args = parser.parse_args()
    table, heavy = load(args.csv_path)
    report_latency(table)
    report_frontend(table)
    report_ball(table)
    report_ler(table)
    report_heavy(heavy)


if __name__ == "__main__":
    main()
