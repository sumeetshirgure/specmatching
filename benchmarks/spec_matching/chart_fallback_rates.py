#!/usr/bin/env python3
# Copyright 2026 SpecMatching contributors
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

"""Tabulates the fraction of shots escalated over the requested `(d, p, T)` settings.

    python benchmarks/spec_matching/chart_fallback_rates.py --in RESULTS_DIR [MORE_DIRS ...]
        --cells d=17,p=1e-3,T=1.5 [d=21,p=1e-3,T=2 ...]
        [--format text|markdown] [--out FILE]

An input is a `sparse_graph_stats` output directory holding `summary.csv`. The fraction is that
row's `shots_escalated / shots`, and both counts stand beside it. Rows are grouped by `p`, then by
`d`, then by `T`, with each group's label written on its first row.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


def fmt_p(value):
    """`1e-3` style, falling back to `%g` for anything that is not a round mantissa."""
    if value == 0:
        return "0"
    exponent = math.floor(math.log10(abs(value)))
    mantissa = value / (10.0**exponent)
    if abs(mantissa - round(mantissa)) < 1e-9:
        mantissa = int(round(mantissa))
        return f"1e{exponent}" if mantissa == 1 else f"{mantissa}e{exponent}"
    return f"{value:g}"


def cell_key(d, p, t):
    """The canonical `(d, p, T)` key.

    `summary.csv` writes its floats at six significant digits. Rounding a value given on the command
    line to that same precision is what lets `1e-3` and `0.001` name one cell.
    """
    return (int(float(d)), float(f"{float(p):.6g}"), float(f"{float(t):.6g}"))


def title_of(key):
    d, p, t = key
    return f"d = {d}, p = {fmt_p(p)}, T = {t:g}"


def parse_cell(text):
    """`d=17,p=1e-3,T=1.5` into a canonical key."""
    fields = {}
    for token in text.split(","):
        if not token.strip():
            continue
        name, separator, value = token.partition("=")
        if not separator:
            raise SystemExit(f"error: '{token}' in '{text}' is not name=value")
        fields[name.strip()] = value.strip()
    missing = [name for name in ("d", "p", "T") if name not in fields]
    if missing:
        raise SystemExit(f"error: '{text}' has no {', '.join(missing)}")
    unknown = sorted(set(fields) - {"d", "p", "T"})
    if unknown:
        raise SystemExit(f"error: '{text}' has unrecognised {', '.join(unknown)}")
    try:
        return cell_key(fields["d"], fields["p"], fields["T"])
    except ValueError as error:
        raise SystemExit(f"error: '{text}': {error}")


def load_rows(in_dirs):
    """Every `summary.csv` row, keyed `(d, p, T)`, carrying its two shot counts."""
    rows = {}
    for directory in in_dirs:
        path = Path(directory) / "summary.csv"
        if not path.is_file():
            raise SystemExit(f"error: {path} does not exist")
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            for name in ("d", "p", "T", "shots", "shots_escalated"):
                if name not in (reader.fieldnames or ()):
                    raise SystemExit(f"error: {path} has no '{name}' column")
            for row in reader:
                key = cell_key(row["d"], row["p"], row["T"])
                rows[key] = (int(row["shots"]), int(row["shots_escalated"]))
    return rows


HEADERS = ("p", "d", "T", "shots", "escalated", "fraction escalated")

#: The group labels read down the left, so they are left aligned; the counts read against each
#: other, so they are right aligned.
ALIGN = ("<", "<", "<", ">", ">", ">")


def table_cells(keys, rows):
    """One display row per key, in `p`, `d`, `T` order, with a repeated group label left blank."""
    body = []
    previous_p = None
    previous_d = None
    for key in sorted(keys, key=lambda key: (key[1], key[0], key[2])):
        d, p, t = key
        shots, escalated = rows[key]
        fraction = escalated / shots if shots else float("nan")
        body.append(
            (
                fmt_p(p) if p != previous_p else "",
                f"{d}" if (p, d) != (previous_p, previous_d) else "",
                f"{t:g}",
                f"{shots:,}",
                f"{escalated:,}",
                f"{fraction:.6g}",
            )
        )
        previous_p, previous_d = p, d
    return body


def render_text(body):
    widths = [max([len(HEADERS[i])] + [len(row[i]) for row in body]) for i in range(len(HEADERS))]

    def line(cells):
        return "  ".join(f"{cells[i]:{ALIGN[i]}{widths[i]}}" for i in range(len(HEADERS))).rstrip()

    lines = [line(HEADERS), "  ".join("-" * widths[i] for i in range(len(HEADERS)))]
    lines.extend(line(row) for row in body)
    return "\n".join(lines) + "\n"


def render_markdown(body):
    lines = ["| " + " | ".join(HEADERS) + " |"]
    lines.append("| " + " | ".join("---" for _ in HEADERS) + " |")
    for row in body:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Tabulate the fraction of shots escalated, grouped by p, then d, then T.",
    )
    parser.add_argument(
        "--in",
        dest="in_dirs",
        nargs="+",
        required=True,
        metavar="DIR",
        help="one or more sparse_graph_stats output dirs, each holding a summary.csv",
    )
    parser.add_argument(
        "--cells",
        nargs="+",
        required=True,
        metavar="d=..,p=..,T=..",
        help="the (d, p, T) settings to tabulate",
    )
    parser.add_argument("--format", choices=["text", "markdown"], default="text")
    parser.add_argument("--out", default=None, metavar="FILE", help="default: stdout")
    args = parser.parse_args(argv)

    rows = load_rows(args.in_dirs)
    keys = []
    for text in args.cells:
        key = parse_cell(text)
        if key not in rows:
            raise SystemExit(f"error: no summary.csv row for {title_of(key)} under {', '.join(args.in_dirs)}")
        if key not in keys:
            keys.append(key)

    body = table_cells(keys, rows)
    table = render_markdown(body) if args.format == "markdown" else render_text(body)

    if args.out:
        path = Path(args.out)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(table)
        print(f"wrote {path}")
    else:
        sys.stdout.write(table)
    return 0


if __name__ == "__main__":
    sys.exit(main())
