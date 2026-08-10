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

"""Draws the per-shot latency distributions `two_phase_latency_profiler` logs.

    python benchmarks/two_phase/plot_latency_histograms.py \
        benchmarks/two_phase/results/latency

One figure per log file, two histograms in each: stock exact decode on the original graph `G`
against the M7 front end on the sparsified graph `H`. Both sides are stock sparse blossom, timed on
the same shots, so the gap between the two distributions is `H`'s own win and nothing else.

The sparsified series is the three stages that are on the critical path once the front end is
running,

    sparse = blossom_ns + dscan_ns + hrvst_ns  (+ escal_stock_ns on an escalating shot)

i.e. the solve on `H`, §M7.7's terminal dual scan and the harvest, plus — on the shots the
certificate rejects — the full Phase-2 stock re-decode on `G`. The escalation tail is therefore in
the histogram rather than amortised away beside it. Ball intersect, `H` build and `Mwpm(H)` build are
**not** in it: §M2's critical-path read discounts them as pipelined out. The profiler logs their sum
as `excluded_ns`, so `--include-excluded` puts them back and draws the honest CPU number instead; the
caption says which of the two is on the page.

What the figure states rather than leaves to be inferred:

  * the means, direct-labelled, because the speedup is a ratio of means and the reader should be able
    to see the two numbers it is a ratio of;
  * how many shots escalated, since those are the whole right tail of the sparsified series;
  * how many shots the plot clips off the right edge, and how many contaminated shots were dropped
    before anything was computed. Neither is ever silent.

No percentiles are reported. The one place a percentile survives is `--x-max-percentile`, which
chooses where the x axis stops — a view control, not a statistic, and the caption states as a count
how many shots fell past it.

The same numbers also go to a text table — `latency_speedups.txt`, written into the same directory as
the figures — so the run is readable without opening an image, and diffable between runs. `--table`
prints that table to the terminal too.
"""

import argparse
import csv
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

# Categorical slots 1 and 2, light and dark steps. Two series, fixed order: `G` is always slot 1 and
# the sparsified graph always slot 2, at every grid point, so a reader moving between figures never
# has to re-read the legend.
THEMES = {
    "light": {
        "surface": "#fcfcfb",
        "text_primary": "#0b0b0b",
        "text_secondary": "#52514e",
        "grid": "#dcdcd8",
        "series": ("#2a78d6", "#eb6834"),
    },
    "dark": {
        "surface": "#1a1a19",
        "text_primary": "#ffffff",
        "text_secondary": "#c3c2b7",
        "grid": "#3a3a37",
        "series": ("#3987e5", "#d95926"),
    },
}

STOCK_LABEL = "stock decode on G"
SPARSE_LABEL = "stock decode on sparsified H"


class Log:
    """One `(d, p, T, mode)` log file: its `#` metadata and its per-shot rows."""

    def __init__(self, path, meta, rows):
        self.path = path
        self.meta = meta
        self.rows = rows

    @property
    def title(self):
        return (
            f"d = {self.meta.get('d', '?')},  p = {self.meta.get('p', '?')},  "
            f"T = {self.meta.get('T', '?')},  isect = {self.meta.get('mode', '?')}"
        )

    @property
    def stem(self):
        return os.path.splitext(os.path.basename(self.path))[0]


def load(path):
    meta = {}
    rows = []
    with open(path) as handle:
        lines = []
        for line in handle:
            if line.startswith("#"):
                key, _, value = line[1:].strip().partition("=")
                meta[key.strip()] = value.strip()
            else:
                lines.append(line)
    for row in csv.DictReader(lines):
        rows.append({key: int(value) for key, value in row.items()})
    if not rows:
        raise ValueError(f"{path}: no shot rows")
    return Log(path, meta, rows)


def collect(paths):
    """Expands directories to the `latency_*.csv` files in them, and sorts for a stable figure order."""
    files = []
    for path in paths:
        if os.path.isdir(path):
            found = sorted(
                os.path.join(path, name)
                for name in os.listdir(path)
                if name.startswith("latency_") and name.endswith(".csv")
            )
            if not found:
                print(f"  note: no latency_*.csv under {path}")
            files.extend(found)
        else:
            files.append(path)
    return files


def series_of(log, include_contaminated, include_excluded):
    """The two series, in microseconds, plus the counts the caption has to state.

    The sparsified side is summed here, from the stage columns rather than from a pre-summed one, so
    that what is charged is visible in this function: solve on `H`, terminal dual scan, harvest, and
    on an escalating shot the Phase-2 re-decode on `G`. `escal_stock_ns` is zero on every shot that
    did not escalate, so no branch is needed to say "escalated shots pay both".

    Contaminated shots — the thread was descheduled mid-shot — are dropped from **both** series
    before anything is computed, so the two stay paired shot for shot. The count is returned rather
    than swallowed.
    """
    stock = []
    sparse = []
    escalated = 0
    dropped = 0
    for row in log.rows:
        if row["contaminated"] and not include_contaminated:
            dropped += 1
            continue
        stock.append(row["stock_g_ns"] / 1000.0)
        value = row["blossom_ns"] + row["dscan_ns"] + row["hrvst_ns"] + row["escal_stock_ns"]
        if include_excluded:
            value += row["excluded_ns"]
        sparse.append(value / 1000.0)
        escalated += row["escalated"]
    return stock, sparse, escalated, dropped


def percentile(values, fraction):
    """Only ever used to choose where the x axis stops; nothing reported is a percentile."""
    if not values:
        return float("nan")
    ordered = sorted(values)
    index = min(int(fraction * (len(ordered) - 1)), len(ordered) - 1)
    return ordered[index]


def mean_of(values):
    return sum(values) / len(values) if values else float("nan")


def stats_of(values):
    return {
        "n": len(values),
        "mean": mean_of(values),
        "max": max(values) if values else float("nan"),
    }


def bin_edges(values, x_max, bins, log_x):
    """Shared edges for both series — two histograms on different bins are not comparable."""
    low = min(v for v in values if v > 0) if log_x else min(values)
    low = max(low, 1e-3)
    high = max(x_max, low * (1.0 + 1e-6))
    if log_x:
        return [10 ** (math.log10(low) + i * (math.log10(high) - math.log10(low)) / bins) for i in range(bins + 1)]
    return [low + i * (high - low) / bins for i in range(bins + 1)]


def draw(log, args, theme):
    stock, sparse, escalated, dropped = series_of(log, args.include_contaminated, args.include_excluded)
    if not stock:
        print(f"  {log.stem}: every shot was contaminated; nothing to plot")
        return None

    both = stock + sparse
    # The axis has to reach past both means with room for their labels, or the escalation tail pushes
    # a mean off the right edge and the figure loses the number it is built around.
    x_max = max(percentile(both, args.x_max_percentile / 100.0), 1.15 * max(mean_of(stock), mean_of(sparse)))
    clipped = sum(1 for value in both if value > x_max)
    edges = bin_edges(both, x_max, args.bins, args.log_x)

    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    fig.patch.set_facecolor(theme["surface"])
    ax.set_facecolor(theme["surface"])

    for values, colour, label in (
        (stock, theme["series"][0], STOCK_LABEL),
        (sparse, theme["series"][1], SPARSE_LABEL),
    ):
        # Filled at low alpha so the overlap is legible either way round, with a 2px edge that keeps
        # each series readable where the fills stack.
        ax.hist(values, bins=edges, color=colour, alpha=0.45, label=label, zorder=2)
        ax.hist(values, bins=edges, color=colour, histtype="step", linewidth=2.0, zorder=3)

    # Selective direct labels: the two means, and nothing else. A number on every bar is noise.
    # Stacked at different heights and always set to the right of their line, because the two means
    # sit close together at every grid point and would otherwise collide or run off the axis.
    #
    # The mean of a distribution with an escalation tail does not sit at the mode, so the line will
    # often be to the right of the peak. That is the point of drawing it: the ratio in the caption is
    # a ratio of these two lines, not of the bumps under them.
    for values, colour, height in ((stock, theme["series"][0], 0.97), (sparse, theme["series"][1], 0.88)):
        mean = mean_of(values)
        ax.axvline(mean, color=colour, linewidth=1.5, linestyle=(0, (4, 3)), zorder=4)
        ax.annotate(
            f"mean {mean:.2f} us",
            xy=(mean, height),
            xycoords=("data", "axes fraction"),
            ha="left",
            va="top",
            fontsize=9,
            # Text stays in ink; the dashed line beside it is what carries the series identity.
            color=theme["text_primary"],
            xytext=(5, 0),
            textcoords="offset points",
            zorder=5,
        )

    if args.log_x:
        ax.set_xscale("log")
    ax.set_xlim(edges[0], edges[-1])
    ax.set_xlabel("per-shot latency (microseconds)", color=theme["text_secondary"], fontsize=10)
    ax.set_ylabel("shots", color=theme["text_secondary"], fontsize=10)
    if args.log_y:
        ax.set_yscale("log")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
    ax.grid(axis="y", color=theme["grid"], linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(theme["grid"])
    ax.tick_params(colors=theme["text_secondary"], labelsize=9)

    ax.set_title(log.title, color=theme["text_primary"], fontsize=13, loc="left", pad=14)
    legend = ax.legend(frameon=False, loc="upper right", fontsize=10)
    for text in legend.get_texts():
        text.set_color(theme["text_primary"])

    stock_stats = stats_of(stock)
    sparse_stats = stats_of(sparse)
    speedup = stock_stats["mean"] / sparse_stats["mean"] if sparse_stats["mean"] else float("nan")
    excluded_note = (
        "sparsified = blossom + dscan + hrvst + intersect + H build + Mwpm(H) build,"
        " plus the Phase-2 re-decode on escalating shots"
        if args.include_excluded
        else "sparsified = blossom + dscan + hrvst, plus the Phase-2 re-decode on escalating shots"
    )
    caption = (
        f"{stock_stats['n']:,} shots  ·  mean ratio {speedup:.2f}x  ·  "
        f"max {stock_stats['max']:.2f} vs {sparse_stats['max']:.2f} us  ·  "
        f"{escalated:,} escalated  ·  {dropped:,} contaminated dropped  ·  "
        f"{clipped:,} beyond the right edge\n{excluded_note}  ·  "
        f"timer {log.meta.get('timer_backend', '?')}"
        f"{'' if log.meta.get('timer_thread_scoped') == '1' else ' (WALL CLOCK)'}"
    )
    fig.text(0.012, 0.005, caption, color=theme["text_secondary"], fontsize=8.5, va="bottom")
    fig.tight_layout(rect=(0, 0.075, 1, 1))

    written = []
    for extension in args.formats:
        out_path = os.path.join(args.out_dir or os.path.dirname(log.path) or ".", f"{log.stem}.{extension}")
        fig.savefig(out_path, dpi=args.dpi, facecolor=theme["surface"])
        written.append(out_path)
    plt.close(fig)
    return {
        "log": log,
        "stock": stock_stats,
        "sparse": sparse_stats,
        "escalated": escalated,
        "dropped": dropped,
        "clipped": clipped,
        "written": written,
    }


def table_lines(results, args):
    """The table view: everything the figures show, as text, so the plots are never the only read.

    One list of lines, so the same table can go to the terminal and to the file beside the figures
    without the two drifting apart.
    """
    lines = [
        f"{'d':>4} {'p':>8} {'T':>5} {'mode':>7} {'shots':>8} "
        f"{'stock mean':>11} {'sparse mean':>12} {'stock max':>10} {'sparse max':>11} "
        f"{'speedup':>9} {'escal':>7} {'contam':>7}"
    ]
    for result in results:
        meta = result["log"].meta
        stock = result["stock"]
        sparse = result["sparse"]
        ratio = stock["mean"] / sparse["mean"] if sparse["mean"] else float("nan")
        lines.append(
            f"{meta.get('d', '?'):>4} {meta.get('p', '?'):>8} {meta.get('T', '?'):>5} "
            f"{meta.get('mode', '?'):>7} {stock['n']:>8,} "
            f"{stock['mean']:>11.3f} {sparse['mean']:>12.3f} "
            f"{stock['max']:>10.3f} {sparse['max']:>11.3f} {ratio:>9.2f} "
            f"{result['escalated']:>7,} {result['dropped']:>7,}"
        )
    lines.append("")
    lines.append(
        "  All times in microseconds, over uncontaminated shots. `speedup` is stock mean over"
        " sparsified\n  mean, both taken over every uncontaminated shot — the escalating ones"
        " included, each charged its\n  Phase-2 re-decode on G on top of its time on H. Dropping"
        " them would price the front end at a\n  rate no deployment gets, so the escalation tail is"
        " inside this ratio rather than beside it."
    )
    lines.append(
        "  sparsified = blossom + dscan + hrvst + intersect + H build + Mwpm(H) build, plus the"
        " Phase-2\n  re-decode on escalating shots."
        if args.include_excluded
        else "  sparsified = blossom + dscan + hrvst, plus the Phase-2 re-decode on escalating"
        " shots. Ball\n  intersect, H build and Mwpm(H) build are not charged (--include-excluded"
        " adds them back)."
    )
    if args.include_contaminated:
        lines.append("  Contaminated shots were KEPT in both series (--include-contaminated).")
    return lines


def write_tables(results, args):
    """Writes the table as text beside the figures it summarises, one file per output directory.

    The figures land in `--out-dir` or beside each log, so the table follows the same rule: a reader
    who has the pictures in a directory has the numbers for exactly those pictures there too.
    """
    by_dir = {}
    for result in results:
        out_dir = args.out_dir or os.path.dirname(result["log"].path) or "."
        by_dir.setdefault(out_dir, []).append(result)

    written = []
    for out_dir, group in sorted(by_dir.items()):
        out_path = os.path.join(out_dir, args.table_name)
        timers = sorted({result["log"].meta.get("timer_backend", "?") for result in group})
        header = [
            "Two-phase per-shot latency: stock decode on G vs stock decode on sparsified H.",
            f"{len(group)} log file(s); timer backend {', '.join(timers)}.",
            "",
        ]
        if any(result["log"].meta.get("timer_thread_scoped") != "1" for result in group):
            header.insert(2, "WARNING: at least one log was timed on a WALL CLOCK, not a thread-scoped counter.")
        with open(out_path, "w") as handle:
            handle.write("\n".join(header + table_lines(group, args)) + "\n")
        written.append(out_path)
    return written


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "paths",
        nargs="+",
        metavar="PATH",
        help="log files, or directories to take every latency_*.csv from",
    )
    parser.add_argument("--out-dir", help="where the figures go; default is beside each log file")
    parser.add_argument("--formats", default="png", help="comma-separated: png,pdf,svg")
    parser.add_argument("--bins", type=int, default=80)
    parser.add_argument("--dpi", type=int, default=160)
    parser.add_argument(
        "--x-max-percentile",
        type=float,
        default=99.9,
        metavar="Q",
        help="where the x axis stops, as a percentile of both series pooled — a view control, not a"
        " reported statistic; the axis is always extended past both means, and the caption says how"
        " many shots fell beyond the edge (default 99.9)",
    )
    parser.add_argument("--log-x", action="store_true", help="log latency axis; the escalation tail is long")
    parser.add_argument("--log-y", action="store_true", help="log count axis, which is where a thin tail shows up")
    parser.add_argument(
        "--include-excluded",
        action="store_true",
        help="add `excluded_ns` (intersect + H build + Mwpm(H) build) back into the sparsified"
        " series, i.e. charge the stages §M2 discounts as pipelined out rather than taking the"
        " critical-path read",
    )
    parser.add_argument(
        "--include-contaminated",
        action="store_true",
        help="keep shots the scheduler interfered with; they are dropped from both series by default",
    )
    parser.add_argument("--theme", choices=sorted(THEMES), default="light")
    parser.add_argument("--table", action="store_true", help="print the summary table to stdout as well")
    parser.add_argument(
        "--table-name",
        default="latency_speedups.txt",
        metavar="NAME",
        help="filename of the text table written beside the figures (default latency_speedups.txt)",
    )
    args = parser.parse_args()
    args.formats = [item.strip() for item in args.formats.split(",") if item.strip()]

    files = collect(args.paths)
    if not files:
        print("  nothing to plot.")
        return 1
    if args.out_dir:
        os.makedirs(args.out_dir, exist_ok=True)

    theme = THEMES[args.theme]
    results = []
    for path in files:
        try:
            log = load(path)
        except (OSError, ValueError, KeyError) as error:
            print(f"  skipping {path}: {error}")
            continue
        # v1 logged one pre-summed `sparse_ns` column and no stage split, so it cannot be re-read
        # under the definition above. Named rather than guessed at.
        if log.meta.get("schema") not in (None, "two_phase_latency_v2"):
            print(f"  skipping {path}: schema {log.meta['schema']}; this script reads two_phase_latency_v2")
            continue
        result = draw(log, args, theme)
        if result is None:
            continue
        results.append(result)
        for out_path in result["written"]:
            print(f"  wrote {out_path}")

    if not results:
        print("  nothing plotted.")
        return 1
    for out_path in write_tables(results, args):
        print(f"  wrote {out_path}")
    if args.table:
        print()
        print("\n".join(table_lines(results, args)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
