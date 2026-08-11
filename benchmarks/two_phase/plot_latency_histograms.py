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

A million-shot campaign is read a chunk at a time and never held as Python objects: see `load` and
`series_of` for what that costs and what it buys.
"""

import argparse
import itertools
import math
import os
import sys

import matplotlib
import numpy as np

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

# Rows are parsed in blocks of this many, so the transient cost of a read is set by the block and not
# by the length of the campaign. Big enough that the per-block overhead is lost in the parse, small
# enough that a block of text is a few megabytes rather than the whole file.
CHUNK_ROWS = 100_000

# The only columns any series needs. Everything else in the row — `shot`, `defects`, `certified`,
# `total_ns` — is never read, so it is never parsed either.
NEEDED = (
    "contaminated",
    "escalated",
    "stock_g_ns",
    "blossom_ns",
    "dscan_ns",
    "hrvst_ns",
    "escal_stock_ns",
    "excluded_ns",
)
AT = {name: position for position, name in enumerate(NEEDED)}


class Log:
    """One `(d, p, T, mode)` log file: its `#` metadata and its column header.

    Deliberately not its rows. The rows stay on disk and are streamed by `series_of`; a `Log` is
    small enough that `main` can hold one per file for the whole run, which is what the summary table
    at the end needs.
    """

    def __init__(self, path, meta, columns):
        self.path = path
        self.meta = meta
        self.columns = columns

    @property
    def title(self):
        return (
            f"d = {self.meta.get('d', '?')},  p = {self.meta.get('p', '?')},  "
            f"T = {self.meta.get('T', '?')},  isect = {self.meta.get('mode', '?')}"
        )

    @property
    def stem(self):
        return os.path.splitext(os.path.basename(self.path))[0]


def read_header(handle):
    """Consumes the `#` metadata and the column header, leaving `handle` on the first shot row."""
    meta = {}
    for line in handle:
        if line.startswith("#"):
            key, _, value = line[1:].strip().partition("=")
            meta[key.strip()] = value.strip()
        elif line.strip():
            return meta, [name.strip() for name in line.strip().split(",")]
    return meta, []


def load(path):
    """Reads the header. The shot rows are left on disk for `series_of` to stream.

    A million-shot log is ~50 MB of text, and materialising it the obvious way — a list of lines,
    then a list of per-shot dicts — costs upwards of a gigabyte for that one file, because a dict of
    twelve boxed ints per row is nearly two orders of magnitude wider than the two numbers actually
    wanted from it. Holding one of those per log file, which is what a summary table over a campaign
    implies, is what runs the machine out of memory rather than any single figure.

    So a `Log` is its header and nothing else, and the rows are read once, in blocks, straight into
    two `float64` arrays. Whether the columns needed are even present is settled here, before a
    figure is started, so a malformed log is reported by name like any other bad file.
    """
    with open(path) as handle:
        meta, columns = read_header(handle)
        if not columns:
            raise ValueError(f"{path}: no column header")
        if next(handle, None) is None:
            raise ValueError(f"{path}: no shot rows")
    # v1 logged one pre-summed `sparse_ns` column and no stage split, so it cannot be re-read under
    # the definition in `series_of`. Named rather than reported as a pile of absent columns.
    if meta.get("schema") not in (None, "two_phase_latency_v2"):
        raise ValueError(f"schema {meta['schema']}; this script reads two_phase_latency_v2")
    missing = [name for name in NEEDED if name not in columns]
    if missing:
        raise ValueError(f"{path}: header is missing {', '.join(missing)}")
    return Log(path, meta, columns)


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

    The file is walked once, `CHUNK_ROWS` rows at a time, and each block is dropped as soon as the
    two columns it contributes have been appended. The filter and the stage sum are done on the block
    while it is still integer nanoseconds, so a contaminated shot is never converted at all. What
    survives the read is two `float64` arrays — 8 MB each per million kept shots, against roughly a
    hundred times that for the same rows as Python objects.
    """
    usecols = tuple(log.columns.index(name) for name in NEEDED)
    stock_chunks = []
    sparse_chunks = []
    escalated = 0
    dropped = 0
    with open(log.path) as handle:
        read_header(handle)
        while True:
            text = list(itertools.islice(handle, CHUNK_ROWS))
            if not text:
                break
            block = np.loadtxt(text, delimiter=",", usecols=usecols, dtype=np.int64, ndmin=2)
            del text
            if not block.size:
                continue
            if include_contaminated:
                kept = block
            else:
                clean = block[:, AT["contaminated"]] == 0
                dropped += int(block.shape[0] - np.count_nonzero(clean))
                kept = block[clean]
            del block
            if not kept.shape[0]:
                continue
            escalated += int(kept[:, AT["escalated"]].sum())
            value = (
                kept[:, AT["blossom_ns"]]
                + kept[:, AT["dscan_ns"]]
                + kept[:, AT["hrvst_ns"]]
                + kept[:, AT["escal_stock_ns"]]
            )
            if include_excluded:
                value += kept[:, AT["excluded_ns"]]
            stock_chunks.append(kept[:, AT["stock_g_ns"]] / 1000.0)
            sparse_chunks.append(value / 1000.0)
    empty = np.empty(0, dtype=np.float64)
    stock = np.concatenate(stock_chunks) if stock_chunks else empty
    stock_chunks.clear()
    sparse = np.concatenate(sparse_chunks) if sparse_chunks else empty
    sparse_chunks.clear()
    return stock, sparse, escalated, dropped


def percentile(values, fraction):
    """Only ever used to choose where the x axis stops; nothing reported is a percentile.

    A partial sort, not a full one: the axis needs one order statistic, and selecting it costs a
    linear pass over one scratch copy instead of sorting a two-million-element pool.
    """
    if values.size == 0:
        return float("nan")
    index = min(int(fraction * (values.size - 1)), values.size - 1)
    return float(np.partition(values, index)[index])


def mean_of(values):
    return float(values.mean()) if values.size else float("nan")


def stats_of(values):
    return {
        "n": int(values.size),
        "mean": mean_of(values),
        "max": float(values.max()) if values.size else float("nan"),
    }


def bin_edges(stock, sparse, x_max, bins, log_x):
    """Shared edges for both series — two histograms on different bins are not comparable.

    Takes the two series rather than a pooled copy of them: the low edge is a minimum, and a minimum
    over a union is the smaller of the two minima.
    """
    if log_x:
        lows = [float(side[side > 0].min()) for side in (stock, sparse) if np.any(side > 0)]
    else:
        lows = [float(side.min()) for side in (stock, sparse) if side.size]
    low = max(min(lows), 1e-3) if lows else 1e-3
    high = max(x_max, low * (1.0 + 1e-6))
    if log_x:
        return np.logspace(math.log10(low), math.log10(high), bins + 1)
    return np.linspace(low, high, bins + 1)


def draw(log, args, theme):
    stock, sparse, escalated, dropped = series_of(log, args.include_contaminated, args.include_excluded)
    if not stock.size:
        print(f"  {log.stem}: every shot was contaminated; nothing to plot")
        return None

    # The axis has to reach past both means with room for their labels, or the escalation tail pushes
    # a mean off the right edge and the figure loses the number it is built around.
    #
    # The pooled copy exists only for that one order statistic and is released before anything is
    # drawn; the count past the edge is two counts summed, which needs no pool at all.
    pooled = np.concatenate((stock, sparse))
    x_max = max(percentile(pooled, args.x_max_percentile / 100.0), 1.15 * max(mean_of(stock), mean_of(sparse)))
    del pooled
    clipped = int(np.count_nonzero(stock > x_max) + np.count_nonzero(sparse > x_max))
    edges = bin_edges(stock, sparse, x_max, args.bins, args.log_x)

    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    fig.patch.set_facecolor(theme["surface"])
    ax.set_facecolor(theme["surface"])

    for values, colour, label in (
        (stock, theme["series"][0], STOCK_LABEL),
        (sparse, theme["series"][1], SPARSE_LABEL),
    ):
        # Binned once and drawn twice. `ax.hist` would re-bin the whole series for each of the two
        # passes and keep a copy of it inside the axes; past a few hundred thousand shots that is the
        # difference between a figure and a swap storm, and the counts are identical either way.
        #
        # Filled at low alpha so the overlap is legible either way round, with a 2px edge that keeps
        # each series readable where the fills stack.
        counts, _ = np.histogram(values, bins=edges)
        ax.stairs(counts, edges, color=colour, alpha=0.45, fill=True, label=label, zorder=2)
        ax.stairs(counts, edges, color=colour, linewidth=2.0, zorder=3)

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
    parser.add_argument(
        "--bins",
        type=int,
        default=200,
        help="bins per series, shared between the two so they stay comparable; the default suits a"
        " campaign of ~1e6 shots, where a coarser binning flattens the escalation tail into the bulk"
        " (default 200)",
    )
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
        # `draw` is inside the same guard as `load`, because with the rows streamed rather than
        # pre-parsed a malformed row is first seen while the figure is being built. One bad file
        # still costs the run one file.
        try:
            result = draw(load(path), args, theme)
        except (OSError, ValueError, KeyError) as error:
            print(f"  skipping {path}: {error}")
            continue
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
