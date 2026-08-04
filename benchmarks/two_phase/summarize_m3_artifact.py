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

"""Renders the M3 exit artifact CSV as the four measured lines of the M3 exit checkpoint.

    python benchmarks/two_phase/summarize_m3_artifact.py \
        benchmarks/two_phase/results/m3_exit_artifact.csv

The CSV is long-format — `section,d,rounds,p,T,key,value`. Nothing is re-derived here that the C++
harness already computed; this aggregates across the grid and flags the two things a reader of an
escalation table has to be told rather than left to notice:

  * a `q` at or below `1 / shots` is the campaign's resolution floor, not a measurement (§M3.0);
  * a percentile quoted without `contaminated_shot_rate` is a measurement of the scheduler (§M6.4).
"""

import argparse
import csv
from collections import defaultdict


def load(path):
    table = defaultdict(dict)
    with open(path) as handle:
        for row in csv.DictReader(handle):
            key = (row["section"], int(row["d"]), int(row["rounds"]), float(row["p"]), float(row["T"]))
            table[key][row["key"]] = float(row["value"])
    return table


def points(table, section):
    return sorted(key for key in table if key[0] == section)


def report_q(table):
    print("\n== escalation rate, with its denominator ==")
    print(f"{'d':>4} {'rounds':>7} {'p':>8} {'T':>5} {'shots':>10} {'escalated':>10} {'q':>12} {'note':<30}")
    floored = 0
    for key in points(table, "q"):
        row = table[key]
        shots = row["shots"]
        escalated = row["shots_escalated"]
        note = ""
        if escalated == 0:
            note = "below the 1/shots floor"
            floored += 1
        elif escalated < 30:
            plural = "event" if escalated == 1 else "events"
            note = f"only {int(escalated)} {plural}; +-{escalated ** 0.5 / shots:.1e}"
        print(
            f"{key[1]:>4} {key[2]:>7} {key[3]:>8g} {key[4]:>5g} {shots:>10.0f} {escalated:>10.0f} "
            f"{row['q']:>12.3e} {note:<30}"
        )
    if floored:
        print(
            f"\n  {floored} point(s) reported q = 0. That is 'none observed at this shot count', not"
            f" 'never happens': raise --shots before quoting it."
        )


def report_cost(table):
    print("\n== escalation cost, and §M3.2's estimate replaced by a measurement ==")
    print(
        f"{'d':>4} {'p':>8} {'T':>5} {'C_phase1':>11} {'C_escal':>11} {'stock|esc':>11} "
        f"{'stock|all':>11} {'ratio':>7} {'penalty':>10}"
    )
    for key in points(table, "cost"):
        row = table[key]
        on_escalated = row["mean_stock_ns_on_escalated"]
        all_shots = row["mean_stock_ns_all_shots"]
        ratio = on_escalated / all_shots if all_shots else 0
        print(
            f"{key[1]:>4} {key[3]:>8g} {key[4]:>5g} {row['c_phase1_ns']:>11.1f} "
            f"{row['c_escalation_ns']:>11.1f} {on_escalated:>11.1f} {all_shots:>11.1f} "
            f"{ratio:>7.2f} {100 * row['amortised_penalty']:>9.4f}%"
        )
    print(
        "\n  `ratio` is the measured stock cost on escalating shots over the mean-over-all-shots"
        " estimate\n  §M3.2 used. Above 1 is expected: escalating shots have surviving trees at T"
        " and are harder\n  than average, so the estimate was a lower bound."
    )
    gaps = [table[key]["amortisation_gap"] for key in points(table, "cost")]
    if gaps:
        print(f"  worst amortisation gap (|amortised - measured| / measured): {max(gaps):.4f}")


def report_latency(table):
    print("\n== end-to-end latency with escalation live ==")
    print(
        f"{'d':>4} {'p':>8} {'T':>5} {'mean':>10} {'p99':>10} {'p999':>10} {'p9999':>10} "
        f"{'max':>10} {'stock':>10} {'contam':>9}"
    )
    for key in points(table, "latency"):
        row = table[key]
        print(
            f"{key[1]:>4} {key[3]:>8g} {key[4]:>5g} {row['mean_ns']:>10.1f} {row['p99_ns']:>10.1f} "
            f"{row['p999_ns']:>10.1f} {row['p9999_ns']:>10.1f} {row['max_ns']:>10.1f} "
            f"{row['stock_mean_ns']:>10.1f} {100 * row['contaminated_shot_rate']:>8.4f}%"
        )
    print(
        "\n  The escalation spike sits near p99.97 at q ~ 3e-4, so p99 and p999 do not show it."
        "\n  Compare p9999 and max against the mean, not against p99."
    )


def report_streaming(table):
    print("\n== streaming budget ==")
    print(
        f"{'d':>4} {'p':>8} {'T':>5} {'round_ns':>9} {'shot_budget':>12} {'worst_esc':>12} "
        f"{'worst_any':>12} {'overrun':>11} {'buffer':>7}"
    )
    for key in points(table, "streaming"):
        row = table[key]
        print(
            f"{key[1]:>4} {key[3]:>8g} {key[4]:>5g} {row['round_ns']:>9.0f} "
            f"{row['shot_budget_ns']:>12.1f} {row['worst_escalated_ns']:>12.1f} "
            f"{row['worst_any_ns']:>12.1f} {row['overrun_ns']:>11.1f} {row['buffer_depth_shots']:>7.0f}"
        )
    print(
        "\n  `buffer` is the input depth in shots that absorbs the worst observed shot at the"
        " quoted round\n  rate. q is extensive in spacetime volume (M1 result 3), so re-run this at"
        " the intended\n  (d, rounds) rather than scaling this table."
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path")
    args = parser.parse_args()
    table = load(args.csv_path)
    report_q(table)
    report_cost(table)
    report_latency(table)
    report_streaming(table)


if __name__ == "__main__":
    main()
