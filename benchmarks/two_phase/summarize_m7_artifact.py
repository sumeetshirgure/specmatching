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

"""Renders the M7 exit artifact CSV as the two lines of the §M7.8 checkpoint.

    python benchmarks/two_phase/summarize_m7_artifact.py \
        benchmarks/two_phase/results/m7_exit_artifact.csv

The CSV is long-format — `section,d,p,T,key,value`. Nothing is re-derived that the C++ harness
already computed; this aggregates across the grid and states the things a reader of these tables
has to be told rather than left to infer:

  * the identity count is the gate. A non-zero entry voids every timing below it.
  * the speedup is `stock_ns(G) / crit_ns(H)` with `crit_ns = blossom + dual_scan + harvest`, both
    sides running *stock* blossom. It discounts ball intersect, H build and Mwpm(H) build entirely,
    so it is an upper bound on an end-to-end win rather than the win — which is why the end-to-end
    ratio is printed beside it and the stage split says how far apart the two are.
  * `escal_us` is `q_this * C_escalation` at the speedup campaign's own shot count. A `q_this` at or
    below `1 / shots` is that campaign's resolution floor rather than a measurement (§M3.0), so the
    escalation charge is flagged as a bound where it sits on one; §M3.2 caps what it can be worth at
    ~0.02% amortised either way.
"""

import argparse
import csv
from collections import defaultdict


def load(path):
    table = defaultdict(dict)
    with open(path) as handle:
        for row in csv.DictReader(handle):
            key = (row["section"], int(row["d"]), float(row["p"]), float(row["T"]))
            table[key][row["key"]] = float(row["value"])
    return table


def points(table, section):
    return sorted(key for key in table if key[0] == section)


def report_identity(table):
    print("\n== §M7.6 identity against stock exact decode on G ==")
    print(
        f"{'d':>4} {'p':>8} {'T':>5} {'shots':>9} {'wt_disagr':>10} {'obs_disagr':>11} "
        f"{'certified':>10} {'escalated':>10} {'no_pm':>7}"
    )
    total = 0
    for key in points(table, "identity"):
        row = table[key]
        total += row["weight_disagreements"] + row["obs_disagreements"]
        print(
            f"{key[1]:>4} {key[2]:>8g} {key[3]:>5g} {row['shots']:>9.0f} "
            f"{row['weight_disagreements']:>10.0f} {row['obs_disagreements']:>11.0f} "
            f"{row['certified']:>10.0f} {row['escalated']:>10.0f} {row['h_no_perfect_matching']:>7.0f}"
        )
    if total:
        print(f"\n  {total:.0f} DISAGREEMENT(S). This is a certificate bug; do not read the timings below.")
    else:
        print("\n  0 disagreements. Every shot, certified or escalated, is exact MWPM.")
    print(
        "  `no_pm` counts shots that escalated because H had no perfect matching, as against"
        " completing\n  with a dual over T. Only the second trigger is silent, and it is the one"
        " invariant 4 guards."
    )


def report_speedup(table):
    """The read: `stock_us / crit_us`, both sides stock blossom, one on `G` and one on `H`.

    `crit_us` is blossom-on-H plus the terminal dual scan plus the harvest, and nothing else — the
    same definition `summarize_m2_artifact.py::critical_ns` uses, with the certificate added because
    it sits between the solve and the emission and gates it. Ball intersect, H build and MWPM build
    are discounted as pipelined out. That is an assumption about a machine that does not exist yet,
    so the end-to-end number is the honest one for a CPU and this is the architectural one.
    """
    keys = points(table, "speedup")
    if not keys:
        return
    print("\n== the read: stock exact on G vs the M7 critical path on H (blossom + dual scan + harvest) ==")
    print(
        f"{'d':>4} {'p':>8} {'T':>5} {'shots':>8} {'stock_us':>9} {'crit_us':>9} {'e2e_us':>9} "
        f"{'speedup':>8} {'e2e':>8} {'escal_us':>10} {'q_this':>10} {'q note':<24}"
    )
    speedups = []
    e2e_speedups = []
    for key in keys:
        row = table[key]
        shots = row["shots"] or 1.0
        stock_us = row["stock_g_ns"] / shots / 1000.0
        crit_us = row["crit_ns"] / shots / 1000.0
        e2e_us = row["total_ns"] / shots / 1000.0
        escal_us = row["amortised_escalation_ns"] / 1000.0
        speedup = stock_us / crit_us if crit_us else float("nan")
        e2e = stock_us / e2e_us if e2e_us else float("nan")
        speedups.append(speedup)
        e2e_speedups.append(e2e)
        # `q_this` is never printed without `shots`: a zero below `1 / shots` is this campaign's
        # resolution floor, so `escal_us` on that row is a bound and not a measured charge.
        q_this = row.get("q_this", 0.0)
        events = q_this * shots
        if events == 0:
            note = f"floor: q < {1.0 / shots:.1e}"
        elif events < 30:
            note = f"{events:.0f} events; +-{events ** 0.5 / shots:.1e}"
        else:
            note = ""
        print(
            f"{key[1]:>4} {key[2]:>8g} {key[3]:>5g} {shots:>8.0f} {stock_us:>9.3f} {crit_us:>9.3f} "
            f"{e2e_us:>9.3f} {speedup:>8.2f} {e2e:>8.2f} {escal_us:>10.4f} {q_this:>10.3e} {note:<24}"
        )
    print(
        f"\n  stock_us / crit_us: min {min(speedups):.2f}, mean {sum(speedups) / len(speedups):.2f},"
        f" max {max(speedups):.2f}"
    )
    print(
        f"  end-to-end:         min {min(e2e_speedups):.2f},"
        f" mean {sum(e2e_speedups) / len(e2e_speedups):.2f}, max {max(e2e_speedups):.2f}"
    )
    print(
        "  Both sides run stock blossom — the same solver, one on G and one on the ball graph — so"
        " the\n  ratio is H's own win and nothing else. It is an upper bound on the end-to-end win,"
        " not the win:\n  it charges the shot only the stages assumed to sit between a syndrome"
        " arriving and a correction\n  leaving, and compares against a stock decode that does the"
        " whole job in one place. `escal_us` is\n  `q_this * C_escalation`, the Phase-2 fallback"
        " amortised over all shots — outside the Phase-1\n  critical path, so it is beside it, not"
        " in it. Where `q note` says `floor`, this campaign saw no\n  escalation at all and"
        " `escal_us` is a bound rather than a charge; §M3.2 caps it at ~0.02%\n  amortised either"
        " way, which is why it is reported and not fitted."
    )

    print("\n  where the shot goes, as a share of the whole two-phase decode:")
    print(
        f"{'d':>4} {'p':>8} {'T':>5} {'isect%':>8} {'hbuild%':>8} {'mwpm%':>8} {'blossom%':>9} "
        f"{'dscan%':>8} {'hrvst%':>8} {'charged%':>9}"
    )
    for key in keys:
        row = table[key]
        total = row["total_ns"] or 1.0
        charged = row["blossom_on_h_ns"] + row["dual_scan_ns"] + row["harvest_ns"]
        print(
            f"{key[1]:>4} {key[2]:>8g} {key[3]:>5g} {100 * row['intersect_ns'] / total:>7.2f}% "
            f"{100 * row['h_build_ns'] / total:>7.2f}% {100 * row['mwpm_build_ns'] / total:>7.2f}% "
            f"{100 * row['blossom_on_h_ns'] / total:>8.2f}% {100 * row['dual_scan_ns'] / total:>7.3f}% "
            f"{100 * row['harvest_ns'] / total:>7.2f}% {100 * charged / total:>8.2f}%"
        )
    print(
        "  `charged%` is exactly what `speedup` bills and `isect%` + `hbuild%` + `mwpm%` is what it"
        " discounts\n  as pipelined out — which is the whole distance between `speedup` and `e2e`."
        " The lower `charged%`\n  is, the more of the read rests on the pipelining assumption rather"
        " than on this machine."
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path")
    args = parser.parse_args()
    table = load(args.csv_path)
    report_identity(table)
    report_speedup(table)


if __name__ == "__main__":
    main()
