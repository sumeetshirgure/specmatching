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

"""Renders the M2.9 exit artifact CSV as the four answers the M2.9 exit checkpoint asks for.

    python benchmarks/two_phase/summarize_m29_artifact.py \
        benchmarks/two_phase/results/m29_exit_artifact.csv

The CSV is long-format — `section,d,p,T,key,value` — so one file holds the §M2.9.6 measurements,
the enumeration A/B, the exposed-root-blossom depth histogram and the §M2.9.5 E1/E2 experiments.
Nothing is re-derived here that the C++ harness already computed; this only aggregates across the
grid, because the decisions §M2.9.4 and §M2.9.5 turn on are grid-wide rather than per point.
"""

import argparse
import csv
import statistics
from collections import Counter, defaultdict


def load(path):
    table = defaultdict(dict)
    histogram = Counter()
    with open(path) as handle:
        for row in csv.DictReader(handle):
            value = float(row["value"])
            if row["section"] == "depth" and row["key"].startswith("exposed_depth_"):
                histogram[int(row["key"].rsplit("_", 1)[1])] += value
                continue
            table[(row["section"], row["d"], row["p"], row["T"])][row["key"]] = value
    return table, histogram


def column(table, section, key):
    return [row[key] for (name, *_), row in table.items() if name == section and key in row]


def span(values, fmt="{:.3f}"):
    if not values:
        return "n/a"
    return "{}..{} (median {})".format(
        fmt.format(min(values)), fmt.format(max(values)), fmt.format(statistics.median(values))
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path")
    options = parser.parse_args()
    table, histogram = load(options.csv_path)

    print("=== M2.9.6 measurements ===")
    print("1. largest_tree_size          max over grid : {:.0f}".format(max(column(table, "measure", "max_largest_tree_size"))))
    print("2. blossom nesting depth      max, all      : {:.0f}".format(max(column(table, "measure", "max_blossom_nesting_depth"))))
    print("   blossom members            max, all      : {:.0f}".format(max(column(table, "measure", "max_blossom_members"))))
    print("   exposed-root depth         max           : {:.0f}".format(max(column(table, "measure", "max_exposed_blossom_depth"))))
    print("   exposed-root members       max           : {:.0f}".format(max(column(table, "measure", "max_exposed_blossom_members"))))
    print("3. blossom formations / shot                : " + span(column(table, "measure", "mean_blossom_formations"), "{:.2f}"))
    print("4. solve dependent depth, mean              : " + span(column(table, "measure", "mean_solve_dependent_depth"), "{:.1f}"))
    print("   solve events, mean                       : " + span(column(table, "measure", "mean_solve_events"), "{:.1f}"))
    print("   harvest dependent depth, mean            : " + span(column(table, "measure", "mean_harvest_dependent_depth"), "{:.1f}"))
    print()
    print("   harvest share of the critical path       : " + span(column(table, "measure", "harvest_share_of_critical_path")))
    print("   (the M2.9 exit checkpoint's replacement for the section's 'roughly a third' estimate)")

    print()
    print("=== exposed-root-blossom nesting depth, the distribution that decides M2.9.4 ===")
    total = sum(histogram.values())
    for depth in sorted(histogram):
        count = histogram[depth]
        if count:
            print("  depth {:<2} {:8.0f} shots  {:7.2%}".format(depth, count, count / total))
    shallow = sum(count for depth, count in histogram.items() if depth <= 3)
    print("  total {:.0f} shots with an exposed root blossom; depth <= 3 in {:.2%}".format(total, shallow / total))

    print()
    print("=== M2.9.6 A/B: M1.3's enumeration vs M2.9.1's ===")
    print("  harvest speedup                           : " + span(column(table, "ab", "harvest_speedup"), "{:.2f}"))
    print("  harvest share of the shot, M1.3           : " + span(column(table, "ab", "frac_harvest_legacy")))
    print("  harvest share of the shot, M2.9           : " + span(column(table, "ab", "frac_harvest_direct")))
    print("  stage split of the M2.9 harvest:")
    for label, key in [
        ("enumerate  ", "frac_harvest_enumerate"),
        ("base descent", "frac_harvest_base_descent"),
        ("shatter    ", "frac_harvest_shatter"),
        ("reduce     ", "frac_harvest_reduce"),
    ]:
        print("    {}                          : {}".format(label, span(column(table, "ab", key))))

    print()
    print("=== M2.9.5 E1 / E2 ===")
    matches = sum(column(table, "split", "splittable_matches"))
    e1 = sum(column(table, "split", "e1_obs_differed"))
    interior = sum(column(table, "split", "interior_obs_nonzero"))
    e2 = sum(column(table, "split", "e2_weight_differed"))
    cycles = sum(column(table, "split", "cycles_examined"))
    nonzero = sum(column(table, "split", "cycles_with_nonzero_obs"))
    if matches:
        print("  splittable matches                        : {:.0f}".format(matches))
        print("  E1  obs moved under a rotated split       : {:.0f} ({:.4%})".format(e1, e1 / matches))
        print("  E1' interior pairing flips an observable  : {:.0f} ({:.4%})".format(interior, interior / matches))
        print("  E2  weight moved                          : {:.0f} ({:.4%})".format(e2, e2 / matches))
    if cycles:
        print("  blossom cycles with non-zero obs mask     : {:.0f}/{:.0f} ({:.5%})".format(nonzero, cycles, nonzero / cycles))
    print()
    print("  M2.9.5 is open only if E1 and E2 both report zero.")


if __name__ == "__main__":
    main()
