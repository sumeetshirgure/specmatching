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

"""Draws the per-shot latency distributions `profiler_driver` logs.

    python benchmarks/spec_matching/plot_latency_histograms.py \
        benchmarks/spec_matching/results/latency

One figure per log file, two histograms in each: stock exact decode on the original graph `G`, and
the decoder system that runs `G` and the M7 front end on the sparsified graph `H` at once. Every side
is stock sparse blossom, timed on the same shots, so the gap between the distributions is `H`'s own
win and nothing else.

The serial `sparse` machine described next is still computed and still reported in the text table —
it is simply not drawn: the figure is the two distributions the speedup is a ratio of, and the third
curve shared its bulk with the system series exactly.

The sparsified series is the three stages that are on the critical path once the front end is
running,

    sparse = blossom_ns + dscan_ns + hrvst_ns  (+ escal_stock_ns on an escalating shot)

i.e. the solve on `H`, §M7.7's terminal dual scan and the harvest, plus — on the shots the
certificate rejects — the full Phase-2 stock re-decode on `G`, *after* them. That is one serial
machine: it starts `G` only once `H` has been rejected. Ball intersect, `H` build and `Mwpm(H)` build
are **not** in it: §M2's critical-path read discounts them as pipelined out. The profiler logs their
sum as `excluded_ns`, so `--include-excluded` puts them back and draws the honest CPU number instead;
the caption says which of the two is on the page.

The system series is the machine anyone would actually deploy: `H` and `G` solved **concurrently**,
the shot ending the moment a usable matching exists.

    system = stock_g_ns                                         on an escalating shot
           = min(blossom_ns + dscan_ns + hrvst_ns, stock_g_ns)  otherwise

An escalating shot is charged the `G` decode alone — that decode was already running, and `H`'s work
bought nothing — so escalation costs a lost race rather than a second decode stacked on the first.
Every speedup here is `mean(stock) / mean(system)`: the decoder system without graph sparsification
against the same system predicating on `H`. The serial `sparse` mean is reported beside it and is
never a denominator.

What the figure states rather than leaves to be inferred:

  * the means, direct-labelled, because the speedup is a ratio of means and the reader should be able
    to see the two numbers it is a ratio of;
  * how many shots escalated, since those are the shots the system is charged the `G` decode alone
    for, and the whole right tail of the serial series in the table;
  * how many shots the plot clips off the right edge, and how many contaminated shots were dropped
    before anything was computed. Neither is ever silent;
  * which timer produced the numbers, and — loudly, in the caption, in the table header and on
    stderr — whether it was a wall clock rather than a thread-scoped counter.

No percentiles are reported. The one place a percentile survives is `--x-max-percentile`, which
chooses where the x axis stops — a view control, not a statistic, and the caption states as a count
how many shots fell past it.

Logs written under the v3 schema also carry the **component structure** of `H` — its connected
components, their sizes, weighted diameters, degrees, `H` edge weights and boundary costs — measured
outside every timed window and charged to no latency number. Two things are drawn from it:

  * a component figure per grid point, from the profiler's companion `components_*.csv`: the five
    distributions as small multiples, plus the headline "how much of the defect set never needs the
    solver at all" numbers. `--no-components` skips it;
  * a second table beside the latency one, from the per-shot rows, so both tables describe the same
    uncontaminated shots.

A v2 log still draws its latency figure; it simply has no component table row and no component
figure. So does a v3 log written under `--no-component-stats`, which has the columns and zeros in
them — "not measured" is read off the header, never off a column of zeros.

A v4 log adds `k` — the profiler ran one independent experiment per `k`, and the solver saw only
components of size `> k`. `--k` and `--horizons` select which of them to draw, and every v4 log in a
run is additionally **overlaid** on one figure per `(d, p, mode)`: the solver+harvest distribution,
one curve per `(k, T)`, which is what comparing `k` values means here. This script does not sweep
`k` — the profiler does, and the human runs it. The size-`<= k` resolve is in none of these numbers
and never was: it is a serial pre-pass on the critical path, outside the profiler's measurement
scope by intent, which the caption states rather than leaves to be assumed.

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
import textwrap

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

# Categorical slots 1 to 3, light and dark steps. Fixed order: `G` is always slot 1, the sparsified
# graph always slot 2 — which the latency figure no longer draws, but the component figure still
# wears — and the concurrent system always slot 3, at every grid point, so a reader moving between
# figures never has to re-read the legend.
THEMES = {
    "light": {
        "surface": "#fcfcfb",
        "text_primary": "#0b0b0b",
        "text_secondary": "#52514e",
        "grid": "#dcdcd8",
        "series": ("#2a78d6", "#eb6834", "#2f8f5b"),
        # The `k` overlay's categorical slots, taken in `(k, T)` order so a run keeps its colour
        # between figures. Cycled if a figure carries more runs than there are slots.
        "k_series": ("#2a78d6", "#eb6834", "#2f8f5b", "#8a5cd6", "#c0a02c", "#4a4a46"),
    },
    "dark": {
        "surface": "#1a1a19",
        "text_primary": "#ffffff",
        "text_secondary": "#c3c2b7",
        "grid": "#3a3a37",
        "series": ("#3987e5", "#d95926", "#3f9e69"),
        "k_series": ("#3987e5", "#d95926", "#3f9e69", "#9a6ee0", "#cfae37", "#9a9a90"),
    },
}

STOCK_LABEL = "stock decode on G"
SYSTEM_LABEL = "system: H and G at once, first match wins"
SOLVE_HARVEST_LABEL = "solver + harvest, on the components of size > k"

# Rows are parsed in blocks of this many, so the transient cost of a read is set by the block and not
# by the length of the campaign. Big enough that the per-block overhead is lost in the parse, small
# enough that a block of text is a few megabytes rather than the whole file.
CHUNK_ROWS = 100_000

# Latency schemas this script reads. v1 logged one pre-summed `sparse_ns` column and no stage split,
# so it cannot be re-read under the definition in `series_of`; v2 is the stage split; v3 adds the
# component columns, which are optional here — a v2 log still draws, without the component panels.
# The `two_phase_*` spellings are the same schemas under the name the driver wrote before the
# rename to `spec_matching`. Logs already on disk are the record of a run that cannot be repeated
# cheaply, so they are read, not rejected.
LATENCY_SCHEMAS = (
    "spec_matching_latency_v2",
    "spec_matching_latency_v3",
    "spec_matching_latency_v4",
    "two_phase_latency_v2",
    "two_phase_latency_v3",
    "two_phase_latency_v4",
)
COMPONENT_SCHEMAS = ("spec_matching_components_v1", "two_phase_components_v1")

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

# v3's component columns: the structure of the sparsified graph the stages above ran on, measured
# outside every timed window. Read as one block or not at all — a log that has some of them and not
# others is malformed, not half-supported.
#
# `h_defects` is the denominator of every fraction taken from these: the nodes of `H`, which the
# component counts partition. Not `defects`, which is the raw detection-event count before the
# negative-weight preamble.
COMPONENT_NEEDED = (
    "h_defects",
    "components",
    "singleton_components",
    "pair_components",
    "nontrivial_components",
    "largest_component",
    "defects_in_trivial",
    "defects_to_solver",
)


class Log:
    """One `(d, p, T, mode)` log file: its `#` metadata and its column header.

    Deliberately not its rows. The rows stay on disk and are streamed by `series_of`; a `Log` is
    small enough that `main` can hold one per file for the whole run, which is what the summary table
    at the end needs.
    """

    def __init__(self, path, meta, columns, has_components=False):
        self.path = path
        self.meta = meta
        self.columns = columns
        # Whether the per-shot component columns are there *and* the run actually collected them.
        # A v3 log written under `--no-component-stats` has the columns and zeros in them, which
        # would otherwise be reported as "every component is trivial" rather than "not measured".
        self.has_components = has_components

    @property
    def title(self):
        k = self.meta.get("k")
        k_part = f",  k = {k}" if k is not None else ""
        return (
            f"d = {self.meta.get('d', '?')},  p = {self.meta.get('p', '?')},  "
            f"T = {self.meta.get('T', '?')}{k_part},  isect = {self.meta.get('mode', '?')}"
        )

    @property
    def k(self):
        """The run's `k`, or None on a pre-v4 log — which is not the same statement as `k = 0`."""
        try:
            return int(self.meta["k"])
        except (KeyError, ValueError):
            return None

    @property
    def horizon(self):
        try:
            return float(self.meta["T"])
        except (KeyError, ValueError):
            return None

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
    if meta.get("schema") not in (None,) + LATENCY_SCHEMAS:
        raise ValueError(f"schema {meta['schema']}; this script reads {' or '.join(LATENCY_SCHEMAS)}")
    missing = [name for name in NEEDED if name not in columns]
    if missing:
        raise ValueError(f"{path}: header is missing {', '.join(missing)}")
    present = [name for name in COMPONENT_NEEDED if name in columns]
    if present and len(present) != len(COMPONENT_NEEDED):
        absent = [name for name in COMPONENT_NEEDED if name not in columns]
        raise ValueError(f"{path}: has some component columns but not {', '.join(absent)}")
    has_components = len(present) == len(COMPONENT_NEEDED) and meta.get("component_stats") != "0"
    return Log(path, meta, columns, has_components)


def parse_number_list(text, cast):
    """`--k 0,2,4` / `--horizons 1.5,2.0` as a set, or None when the flag was not given.

    None and the empty set are deliberately different: "no filter" and "a filter nothing passes"
    are different requests, and an empty `--k ,` is a typo worth reporting rather than a silent
    run that plots nothing.
    """
    if text is None:
        return None
    values = set()
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            values.add(cast(item))
        except ValueError:
            raise ValueError(f"cannot read {item!r} as a {cast.__name__}") from None
    if not values:
        raise ValueError(f"no values in {text!r}")
    return values


def selected(log, wanted_k, wanted_T):
    """Whether a log passes `--k` / `--horizons`.

    A pre-v4 log has no `k` at all, and `--k` therefore excludes it: it was written before `k`
    existed, so claiming it is any particular `k` would be inventing the label.
    """
    if wanted_k is not None and log.k not in wanted_k:
        return False
    if wanted_T is not None and log.horizon not in wanted_T:
        return False
    return True


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
    """The three series, in microseconds, plus the counts the caption has to state.

    The sparsified side is summed here, from the stage columns rather than from a pre-summed one, so
    that what is charged is visible in this function: solve on `H`, terminal dual scan, harvest, and
    on an escalating shot the Phase-2 re-decode on `G` after them. `escal_stock_ns` is zero on every
    shot that did not escalate, so no branch is needed to say "escalated shots pay both".

    The system side is the concurrent machine, and is where the branch is real: an escalating shot is
    charged `stock_g_ns` and nothing else, because in that machine the `G` decode has been running
    since the start of the shot and the `H` work is simply discarded. Every other shot is charged
    whichever of the two finished first. Both are derived here rather than logged, from columns the
    profiler already writes, and the profiler states the same definition in each log's header.

    Contaminated shots — the thread was descheduled mid-shot — are dropped from **all three** series
    before anything is computed, so they stay paired shot for shot. The count is returned rather
    than swallowed.

    The file is walked once, `CHUNK_ROWS` rows at a time, and each block is dropped as soon as the
    two columns it contributes have been appended. The filter and the stage sum are done on the block
    while it is still integer nanoseconds, so a contaminated shot is never converted at all. What
    survives the read is two `float64` arrays — 8 MB each per million kept shots, against roughly a
    hundred times that for the same rows as Python objects.

    v3's component columns are reduced to running totals *inside* the same loop and never kept per
    shot: the component structure is reported as means and rates, and the shape of its distributions
    comes from the companion `components_*.csv`, which is already aggregated. So the per-shot memory
    cost of reading them is nothing at all.
    """
    names = NEEDED + (COMPONENT_NEEDED if log.has_components else ())
    at = {name: position for position, name in enumerate(names)}
    usecols = tuple(log.columns.index(name) for name in names)
    stock_chunks = []
    sparse_chunks = []
    system_chunks = []
    escalated = 0
    dropped = 0
    components = dict.fromkeys(COMPONENT_NEEDED, 0) if log.has_components else None
    if components is not None:
        components.update(shots=0, shots_with_defects=0, solver_empty=0, max_component=0)
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
                clean = block[:, at["contaminated"]] == 0
                dropped += int(block.shape[0] - np.count_nonzero(clean))
                kept = block[clean]
            del block
            if not kept.shape[0]:
                continue
            escalated += int(kept[:, at["escalated"]].sum())
            stock = kept[:, at["stock_g_ns"]]
            front_end = kept[:, at["blossom_ns"]] + kept[:, at["dscan_ns"]] + kept[:, at["hrvst_ns"]]
            value = front_end + kept[:, at["escal_stock_ns"]]
            if include_excluded:
                # The pipelined-out stages are charged to the front end itself, so they land in both
                # the serial series and the system's `H` side — in the concurrent machine they are
                # still work that has to happen before `H` can answer.
                front_end = front_end + kept[:, at["excluded_ns"]]
                value = value + kept[:, at["excluded_ns"]]
            system = np.where(
                kept[:, at["escalated"]] != 0,
                stock,
                np.minimum(front_end, stock),
            )
            stock_chunks.append(stock / 1000.0)
            sparse_chunks.append(value / 1000.0)
            system_chunks.append(system / 1000.0)
            if components is not None:
                accumulate_components(components, kept, at)
    empty = np.empty(0, dtype=np.float64)
    stock = np.concatenate(stock_chunks) if stock_chunks else empty
    stock_chunks.clear()
    sparse = np.concatenate(sparse_chunks) if sparse_chunks else empty
    sparse_chunks.clear()
    system = np.concatenate(system_chunks) if system_chunks else empty
    system_chunks.clear()
    return stock, sparse, system, escalated, dropped, components


def solve_harvest_series(log, include_contaminated):
    """`blossom_ns + hrvst_ns` per shot, in microseconds — the figure of merit `k` moves.

    Read on its own rather than carried out of `series_of`, so that the overlay holds one campaign's
    array at a time instead of one per log file for the whole run. The cost is a second streaming
    pass over the logs in a group, which is a read of two integer columns and no allocation beyond
    the block.

    It is the solve on the size-`> k` graph and the harvest, and nothing else: not the dual scan,
    not the Phase-2 re-decode on an escalating shot, and not the size-`<= k` resolve, which the
    profiler never timed — it is a serial pre-pass outside its measurement scope.
    """
    names = ("contaminated", "blossom_ns", "hrvst_ns")
    at = {name: position for position, name in enumerate(names)}
    usecols = tuple(log.columns.index(name) for name in names)
    chunks = []
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
                clean = block[:, at["contaminated"]] == 0
                dropped += int(block.shape[0] - np.count_nonzero(clean))
                kept = block[clean]
            del block
            if kept.shape[0]:
                chunks.append((kept[:, at["blossom_ns"]] + kept[:, at["hrvst_ns"]]) / 1000.0)
    values = np.concatenate(chunks) if chunks else np.empty(0, dtype=np.float64)
    chunks.clear()
    return values, dropped


def draw_k_overlay(logs, args, theme):
    """One figure per `(d, p, mode)`: the solver+harvest distribution, one curve per `(k, T)`.

    This is the comparison the `k` experiment is for, and it is an overlay rather than a sweep: the
    profiler ran each `k` as an independent experiment and this draws the runs it was given, in
    ascending `(k, T)` order so a curve keeps its colour and its legend position between figures.

    Curves are outlines rather than fills. Two filled histograms read as a comparison; five read as
    a stack, and the quantity being compared here is where each distribution sits, not how they
    overlap.
    """
    runs = []
    for log in logs:
        values, dropped = solve_harvest_series(log, args.include_contaminated)
        if values.size:
            runs.append({"log": log, "values": values, "dropped": dropped})
    if len(runs) < 2:
        # One curve is the per-log figure again, drawn worse. The overlay exists to compare runs.
        for run in runs:
            del run["values"]
        return []
    runs.sort(key=lambda run: (run["log"].k if run["log"].k is not None else -1, run["log"].horizon or 0.0))

    pooled_max = max(percentile(run["values"], args.x_max_percentile / 100.0) for run in runs)
    x_max = max(pooled_max, 1.15 * max(mean_of(run["values"]) for run in runs))
    lows = [
        float(run["values"][run["values"] > 0].min()) if args.log_x else float(run["values"].min())
        for run in runs
        if run["values"].size and (not args.log_x or np.any(run["values"] > 0))
    ]
    low = max(min(lows), 1e-3) if lows else 1e-3
    high = max(x_max, low * (1.0 + 1e-6))
    edges = (
        np.logspace(math.log10(low), math.log10(high), args.bins + 1)
        if args.log_x
        else np.linspace(low, high, args.bins + 1)
    )

    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    fig.patch.set_facecolor(theme["surface"])
    ax.set_facecolor(theme["surface"])
    palette = theme["k_series"]
    clipped = 0
    for index, run in enumerate(runs):
        values = run["values"]
        colour = palette[index % len(palette)]
        clipped += int(np.count_nonzero(values > x_max))
        counts, _ = np.histogram(values, bins=edges)
        label = f"k = {run['log'].k},  T = {run['log'].meta.get('T', '?')}  ·  mean {mean_of(values):.2f} us"
        ax.stairs(counts, edges, color=colour, linewidth=2.0, label=label, zorder=3)
        ax.axvline(mean_of(values), color=colour, linewidth=1.2, linestyle=(0, (4, 3)), zorder=2)
        del run["values"]

    if args.log_x:
        ax.set_xscale("log")
    ax.set_xlim(edges[0], edges[-1])
    ax.set_xlabel("solver + harvest, per shot (microseconds)", color=theme["text_secondary"], fontsize=10)
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

    first = runs[0]["log"]
    ax.set_title(
        f"d = {first.meta.get('d', '?')},  p = {first.meta.get('p', '?')},  "
        f"isect = {first.meta.get('mode', '?')}  ·  {SOLVE_HARVEST_LABEL}",
        color=theme["text_primary"],
        fontsize=13,
        loc="left",
        pad=14,
    )
    legend = ax.legend(frameon=False, loc="upper right", fontsize=9)
    for text in legend.get_texts():
        text.set_color(theme["text_primary"])

    caption = (
        f"{len(runs)} independent runs on the same shots, one per (k, T)  ·  "
        f"{clipped:,} shots beyond the right edge\n"
        "solver + harvest only: the resolve of the size <= k components runs before it, in series, on the"
        " critical path, and is in none of these numbers —"
        " it is outside the profiler's measurement scope by intent, not free and not concurrent."
    )
    caption = wrap_caption(caption)
    fig.text(0.012, 0.005, caption, color=theme["text_secondary"], fontsize=8.5, va="bottom")
    # The axes give up exactly as much of the figure as the caption turned out to need, so a longer
    # caption pushes the plot up rather than being drawn over it.
    fig.tight_layout(rect=(0, min(0.3, 0.02 + 0.033 * (caption.count("\n") + 1)), 1, 1))

    stem = f"bh_by_k_d{first.meta.get('d', 'NA')}_p{first.meta.get('p', 'NA')}_{first.meta.get('mode', 'NA')}"
    written = []
    for extension in args.formats:
        out_path = os.path.join(args.out_dir or os.path.dirname(first.path) or ".", f"{stem}.{extension}")
        fig.savefig(out_path, dpi=args.dpi, facecolor=theme["surface"])
        written.append(out_path)
    plt.close(fig)
    return written


def accumulate_components(totals, kept, at):
    """One block's component columns, folded into the running totals. O(1) memory in the campaign.

    `shots_with_defects` and `solver_empty` are counted apart from `shots` because the rate worth
    quoting — how often the solver and the harvest would not have run at all — has to be conditioned
    on the shot having had a defect. A shot with an empty syndrome skips the solve for a reason that
    has nothing to do with the prune, and pooling the two would report the corpus's zero-defect rate
    as if it were the prune's win.
    """
    for name in COMPONENT_NEEDED:
        totals[name] += int(kept[:, at[name]].sum())
    has_defects = kept[:, at["h_defects"]] > 0
    totals["shots"] += int(kept.shape[0])
    totals["shots_with_defects"] += int(np.count_nonzero(has_defects))
    totals["solver_empty"] += int(np.count_nonzero(has_defects & (kept[:, at["defects_to_solver"]] == 0)))
    totals["max_component"] = max(totals["max_component"], int(kept[:, at["largest_component"]].max()))


def component_file_for(log):
    """The `components_*.csv` the profiler writes beside a `latency_*.csv`, or None.

    The two files share everything after the prefix, so the pairing is a rename rather than a match
    on parsed metadata — a log renamed by hand loses its companion, which is the intended failure:
    guessing which of several component files belongs to it would be worse.
    """
    name = os.path.basename(log.path)
    if not name.startswith("latency_"):
        return None
    path = os.path.join(os.path.dirname(log.path), "components_" + name[len("latency_") :])
    return path if os.path.exists(path) else None


def read_component_file(path):
    """The aggregated component structure of one grid point: its scalars and its histograms.

    Long form on disk — `kind,key,value` — so one schema covers a table of scalars and five
    histograms of different lengths. Bins arrive in order but are placed by index rather than
    appended, so a file that is sorted differently, or that omits an empty trailing bin, still reads
    correctly.
    """
    scalars = {}
    raw = {}
    with open(path) as handle:
        meta, columns = read_header(handle)
        if meta.get("schema") not in COMPONENT_SCHEMAS:
            raise ValueError(f"{path}: schema {meta.get('schema')}; this script reads {' or '.join(COMPONENT_SCHEMAS)}")
        if columns != ["kind", "key", "value"]:
            raise ValueError(f"{path}: expected a kind,key,value header")
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            kind, _, rest = line.partition(",")
            key, _, value = rest.partition(",")
            if kind == "scalar":
                scalars[key] = float(value)
            else:
                raw.setdefault(kind, {})[int(key)] = float(value)
    histograms = {}
    for kind, bins in raw.items():
        counts = np.zeros(max(bins) + 1, dtype=np.float64)
        for index, count in bins.items():
            counts[index] = count
        histograms[kind] = counts
    return {"path": path, "meta": meta, "scalars": scalars, "hists": histograms}


# The five distributions the component file carries, in reading order, with the unit of their bins.
# `T` means the bin index is in units of the horizon divided by `bins_per_T`; `count` means the bin
# index *is* the value. The last bin of every one of them overflows, which the axis labels say.
COMPONENT_PANELS = (
    ("size_hist", "components", "component size (defects)", "count"),
    ("degree_hist", "H nodes", "degree in H", "count"),
    ("edge_weight_hist", "H edges", "edge weight $d_G(u,v)$  (units of T)", "T"),
    ("bcost_hist", "defects", "boundary cost  (units of T)", "T"),
    ("diameter_hist", "components", "weighted component diameter  (units of T)", "T"),
)


def draw_component_panel(ax, counts, unit, bins_per_T, label, colour, theme):
    """One histogram, drawn like the latency figure's: filled at low alpha under a 2px edge.

    The count axis is logarithmic because these distributions are the point: the bulk is one or two
    bins tall and everything the pruning experiment is about — the large components, the long edges —
    lives three or four decades down. On a linear axis that tail is a flat line at zero.
    """
    edges = np.arange(counts.size + 1, dtype=np.float64)
    if unit == "T":
        edges /= float(bins_per_T)
    ax.stairs(counts, edges, color=colour, alpha=0.45, fill=True, zorder=2)
    ax.stairs(counts, edges, color=colour, linewidth=2.0, zorder=3)
    ax.set_yscale("log")
    ax.set_xlim(edges[0], edges[-1])
    ax.set_xlabel(label, color=theme["text_secondary"], fontsize=9)
    ax.grid(axis="y", color=theme["grid"], linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(theme["grid"])
    ax.tick_params(colors=theme["text_secondary"], labelsize=8)


def draw_components(log, component, args, theme):
    """The component structure of `H` at one grid point: five distributions and the headline numbers.

    Small multiples of one hue, not a five-series chart: these are five different populations —
    components, nodes, edges — and colouring them apart would imply a comparison between quantities
    that share no axis. The one hue is the same one the latency figure gives the sparsified series,
    because every panel here describes that same graph.

    The sixth cell is deliberately not a chart. "How much of the defect set never needs the solver"
    is a handful of numbers, and a bar chart of six unrelated percentages is harder to read than the
    percentages.
    """
    hists = component["hists"]
    scalars = component["scalars"]
    bins_per_T = float(component["meta"].get("bins_per_T", 16))
    colour = theme["series"][1]

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.0))
    fig.patch.set_facecolor(theme["surface"])
    flat = axes.ravel()
    for ax in flat:
        ax.set_facecolor(theme["surface"])

    for ax, (kind, of_what, label, unit) in zip(flat, COMPONENT_PANELS):
        counts = hists.get(kind)
        if counts is None or not counts.sum():
            ax.text(
                0.5,
                0.5,
                f"no {kind}",
                transform=ax.transAxes,
                ha="center",
                va="center",
                color=theme["text_secondary"],
                fontsize=9,
            )
            ax.set_axis_off()
            continue
        draw_component_panel(ax, counts, unit, bins_per_T, label, colour, theme)
        ax.set_ylabel(of_what, color=theme["text_secondary"], fontsize=9)
        # One direct label per panel, and only where a threshold means something: the size axis is
        # where "trivial" is defined, and `H` holds no edge longer than `2T` by construction, so a
        # bar at the right edge of that panel is the overflow bin and not a longer edge.
        if kind == "size_hist":
            # The line sits at this run's own `k`, so the panel says which components the solver was
            # spared rather than restating a size-2 rule the run may not have used.
            k = log.k if log.k is not None else 2
            ax.axvline(k + 1.0, color=theme["text_secondary"], linewidth=1.0, linestyle=(0, (4, 3)), zorder=4)
            ax.annotate(
                f"size $\\leq$ {k}: no solver",
                xy=(k + 1.0, 0.94),
                xycoords=("data", "axes fraction"),
                ha="left",
                va="top",
                fontsize=8,
                color=theme["text_primary"],
                xytext=(4, 0),
                textcoords="offset points",
                zorder=5,
            )
        elif kind == "edge_weight_hist":
            ax.axvline(2.0, color=theme["text_secondary"], linewidth=1.0, linestyle=(0, (4, 3)), zorder=4)
            ax.annotate(
                "2T: H's own cutoff",
                xy=(2.0, 0.94),
                xycoords=("data", "axes fraction"),
                ha="right",
                va="top",
                fontsize=8,
                color=theme["text_primary"],
                xytext=(-4, 0),
                textcoords="offset points",
                zorder=5,
            )

    summary_ax = flat[len(COMPONENT_PANELS)]
    summary_ax.set_axis_off()
    # Every percentage here states the population it is over. The three defect shares are over the
    # same denominator — `H`'s nodes — and do not sum to 1: a component can be small enough to be
    # trivial and still ambiguous enough to go to the solver, and that gap is the interesting number.
    analysed = scalars.get("shots_with_component_stats", 0.0)
    with_defects = scalars.get("shots_with_component_defects", 0.0)
    lines = [
        f"what the prune removes at k = {log.k}" if log.k is not None else "what the prune would remove",
        "",
        "share of H's defects",
        f"  in a component of size $\\leq$ 2  {100 * scalars.get('frac_defects_in_trivial_components', 0):>6.1f}%",
        f"  resolved off the solver        {100 * scalars.get('frac_defects_committed_trivially', 0):>6.1f}%",
        f"  left to the solver anyway      {100 * scalars.get('frac_defects_to_solver', 0):>6.1f}%",
        "",
        "share of shots with a defect",
        f"  no solve, no harvest at all    {100 * scalars.get('solver_set_empty_rate', 0):>6.1f}%",
        "",
        f"defects resolved per shot        {scalars.get('mean_defects_resolved_small', 0):>6.2f}",
        f"components per shot              {scalars.get('mean_components', 0):>6.2f}",
        f"largest, mean over shots         {scalars.get('mean_largest_component_size', 0):>6.2f}",
        f"largest, over the campaign       {scalars.get('max_component_size', 0):>6.0f}",
        f"touching the boundary            {100 * scalars.get('boundary_touching_component_fraction', 0):>6.1f}%",
        "",
        f"{analysed:,.0f} shots analysed, {with_defects:,.0f} with a defect",
    ]
    # Sized so the block cannot reach the caption: seventeen lines at this size and spacing are
    # shorter than the cell, which a figure with a longer caption or a taller font would not be.
    summary_ax.text(
        0.0,
        1.0,
        "\n".join(lines),
        transform=summary_ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        family="monospace",
        color=theme["text_primary"],
        linespacing=1.4,
    )

    fig.suptitle(
        f"component structure of H  —  {log.title}",
        color=theme["text_primary"],
        fontsize=13,
        x=0.012,
        ha="left",
    )
    uncomputed = scalars.get("diameter_uncomputed_component_fraction", 0.0)
    clauses = [
        "Counts on a log axis; the last bin of every histogram is an overflow bin",
        f"weight bins are T/{bins_per_T:.0f}",
        "diameters are shortest paths confined to the component, so they are H-subgraph distances"
        " rather than distances in G",
    ]
    if uncomputed:
        clauses.append(f"{100 * uncomputed:.2f}% of components were too large to measure one")
    # Wrapped rather than left to run: at this width the clause list is longer than the figure, and
    # a caption that leaves the page is a caption that was not read.
    caption = textwrap.fill("  ·  ".join(clauses), width=170, break_long_words=False)
    caption += "\nStructural, measured outside every timed window: none of this is charged to any latency number."
    fig.text(0.012, 0.005, caption, color=theme["text_secondary"], fontsize=8.5, va="bottom")
    fig.tight_layout(rect=(0, 0.075, 1, 0.955))

    stem = os.path.splitext(os.path.basename(component["path"]))[0]
    written = []
    for extension in args.formats:
        out_path = os.path.join(args.out_dir or os.path.dirname(log.path) or ".", f"{stem}.{extension}")
        fig.savefig(out_path, dpi=args.dpi, facecolor=theme["surface"])
        written.append(out_path)
    plt.close(fig)
    return written


def wrap_caption(text, width=135):
    """Re-wraps a caption to the figure width, keeping the author's own line breaks.

    A caption that runs off the right edge is a caption that was not read, and these say things —
    which shots were dropped, which timer produced the numbers, what the speedup is a ratio of — that
    the figure is not honest without. The default width is what fits a 9-inch figure at 8.5 pt.
    """
    return "\n".join(
        textwrap.fill(line, width=width, break_long_words=False) if line else ""
        for line in text.split("\n")
    )


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


def bin_edges(sides, x_max, bins, log_x):
    """Shared edges for every series — histograms on different bins are not comparable.

    Takes the series rather than a pooled copy of them: the low edge is a minimum, and a minimum over
    a union is the smallest of the individual minima.
    """
    if log_x:
        lows = [float(side[side > 0].min()) for side in sides if np.any(side > 0)]
    else:
        lows = [float(side.min()) for side in sides if side.size]
    low = max(min(lows), 1e-3) if lows else 1e-3
    high = max(x_max, low * (1.0 + 1e-6))
    if log_x:
        return np.logspace(math.log10(low), math.log10(high), bins + 1)
    return np.linspace(low, high, bins + 1)


def draw(log, args, theme):
    stock, sparse, system, escalated, dropped, components = series_of(
        log, args.include_contaminated, args.include_excluded
    )
    if not stock.size:
        print(f"  {log.stem}: every shot was contaminated; nothing to plot")
        return None

    # The axis has to reach past every mean with room for its label, or the escalation tail pushes a
    # mean off the right edge and the figure loses the number it is built around.
    #
    # Only the two drawn series set it. The serial `sparse` series is still computed, for the table,
    # but an axis stretched to cover a curve nobody can see would only compress the ones they can.
    #
    # The pooled copy exists only for that one order statistic and is released before anything is
    # drawn; the count past the edge is two counts summed, which needs no pool at all.
    pooled = np.concatenate((stock, system))
    x_max = max(
        percentile(pooled, args.x_max_percentile / 100.0),
        1.15 * max(mean_of(stock), mean_of(system)),
    )
    del pooled
    clipped = int(np.count_nonzero(stock > x_max) + np.count_nonzero(system > x_max))
    edges = bin_edges((stock, system), x_max, args.bins, args.log_x)

    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    fig.patch.set_facecolor(theme["surface"])
    ax.set_facecolor(theme["surface"])

    # `stock` and `system` are the two the speedup is a ratio of, and they are the only two drawn.
    # The serial series is not: it shares its bulk with the system series exactly — the two differ on
    # escalating shots and nowhere else — so its curve was a second tracing of the same body, and the
    # escalation count in the caption says how many shots the two would have parted company on. Its
    # mean and max are in the text table, where they are read rather than compared by eye.
    for values, colour, label in (
        (stock, theme["series"][0], STOCK_LABEL),
        (system, theme["series"][2], SYSTEM_LABEL),
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
    # Stacked at different heights and always set to the right of their line, because the means sit
    # close together at every grid point and would otherwise collide or run off the axis.
    #
    # The mean of a distribution with an escalation tail does not sit at the mode, so the line will
    # often be to the right of the peak. That is the point of drawing it: the ratio in the caption is
    # a ratio of the stock line and the system line, not of the bumps under them.
    for values, colour, height in (
        (stock, theme["series"][0], 0.97),
        (system, theme["series"][2], 0.88),
    ):
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
    system_stats = stats_of(system)
    # The one speedup on the page: the decoder system without sparsification over the same system
    # predicating on H. The serial mean is in the table beside it, and is never a denominator.
    speedup = stock_stats["mean"] / system_stats["mean"] if system_stats["mean"] else float("nan")
    excluded_note = (
        "front end = blossom + dscan + hrvst + intersect + H build + Mwpm(H) build"
        if args.include_excluded
        else "front end = blossom + dscan + hrvst"
    )
    caption = (
        f"{stock_stats['n']:,} shots  ·  speedup {speedup:.2f}x = stock mean / system mean  ·  "
        f"max {stock_stats['max']:.2f} vs {system_stats['max']:.2f} us  ·  "
        f"{escalated:,} escalated  ·  {dropped:,} contaminated dropped  ·  "
        f"{clipped:,} beyond the right edge\n{excluded_note}; system = min(front end, G) on an"
        " ordinary shot and G alone on an escalating one, since that decode was already running  ·  "
        f"timer {log.meta.get('timer_backend', '?')}"
        f"{'' if log.meta.get('timer_thread_scoped') == '1' else ' — WALL CLOCK, timings include off-CPU time'}"
        + (
            ""
            if log.k is None
            else f"\nk = {log.k}: the solver saw only components of size > k. The size-<= k resolve"
            " runs before it, in series, and is in none of these numbers — outside the profiler's"
            " measurement scope by intent."
        )
    )
    caption = wrap_caption(caption)
    fig.text(0.012, 0.005, caption, color=theme["text_secondary"], fontsize=8.5, va="bottom")
    # As on the overlay: the caption states what the figure would not be honest without, so the axes
    # yield to it rather than the other way round.
    fig.tight_layout(rect=(0, min(0.3, 0.02 + 0.033 * (caption.count("\n") + 1)), 1, 1))

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
        "system": system_stats,
        "escalated": escalated,
        "dropped": dropped,
        "clipped": clipped,
        "components": components,
        "written": written,
    }


def table_lines(results, args):
    """The table view: everything the figures show, as text, so the plots are never the only read.

    One list of lines, so the same table can go to the terminal and to the file beside the figures
    without the two drifting apart.
    """
    lines = [
        f"{'d':>4} {'p':>8} {'T':>5} {'k':>3} {'mode':>7} {'shots':>8} "
        f"{'stock mean':>11} {'sparse mean':>12} {'system mean':>12} "
        f"{'stock max':>10} {'sparse max':>11} {'system max':>11} "
        f"{'speedup':>9} {'escal':>7} {'contam':>7}"
    ]
    for result in results:
        log = result["log"]
        meta = log.meta
        stock = result["stock"]
        sparse = result["sparse"]
        system = result["system"]
        ratio = stock["mean"] / system["mean"] if system["mean"] else float("nan")
        lines.append(
            f"{meta.get('d', '?'):>4} {meta.get('p', '?'):>8} {meta.get('T', '?'):>5} "
            f"{'-' if log.k is None else log.k:>3} "
            f"{meta.get('mode', '?'):>7} {stock['n']:>8,} "
            f"{stock['mean']:>11.3f} {sparse['mean']:>12.3f} {system['mean']:>12.3f} "
            f"{stock['max']:>10.3f} {sparse['max']:>11.3f} {system['max']:>11.3f} {ratio:>9.2f} "
            f"{result['escalated']:>7,} {result['dropped']:>7,}"
        )
    lines.append("")
    lines.append(
        "  All times in microseconds, over uncontaminated shots. `speedup` is stock mean over"
        " SYSTEM mean —\n  the decoder system without graph sparsification against the same system"
        " predicating on H — and\n  never over the serial `sparse` mean. Both are taken over every"
        " uncontaminated shot, the escalating\n  ones included: dropping them would price the front"
        " end at a rate no deployment gets, so the cost\n  of escalation is inside this ratio rather"
        " than beside it."
    )
    lines.append(
        "  sparse = front end, then the whole Phase-2 re-decode on G after it on an escalating shot."
        "\n  system = H and G solved concurrently, the shot ending at the first usable matching:"
        " stock_g on an\n  escalating shot (that decode was already running, so escalation costs a"
        " lost race and not a\n  second decode), min(front end, stock_g) on any other."
    )
    lines.append(
        "  front end = blossom + dscan + hrvst + intersect + H build + Mwpm(H) build."
        if args.include_excluded
        else "  front end = blossom + dscan + hrvst. Ball intersect, H build and Mwpm(H) build are"
        " not charged\n  (--include-excluded adds them back to both series)."
    )
    if any(result["log"].k is not None for result in results):
        lines.append(
            "  `k`: the solver saw only components of size > k; sizes 1..k were resolved exactly off"
            " it. That\n  resolve is a serial pre-pass on the critical path and is in none of these"
            " numbers — the profiler\n  scoped its measurement to the solver and the harvest, and"
            " left the resolve to be measured separately.\n  `-` is a log written before k existed."
        )
    if args.include_contaminated:
        lines.append("  Contaminated shots were KEPT in both series (--include-contaminated).")
    lines.extend(component_table_lines(results))
    return lines


def component_table_lines(results):
    """The component structure as text, over the same uncontaminated shots as the table above.

    Every column comes from the per-shot rows rather than from the companion file, so this table and
    the latency table describe exactly the same set of shots: the companion file is aggregated over
    *all* shots, contaminated ones included, and quoting the two side by side would be quoting two
    populations under one heading. The companion file's own numbers are on the component figure,
    where they are the whole subject and the difference is stated.
    """
    rows = [result for result in results if result.get("components")]
    if not rows:
        return []
    lines = [
        "",
        "Component structure of H — same uncontaminated shots as above. `triv` is the share of H's"
        " defects in a",
        "component of size <= 2, which is a fixed size class and not this run's k; `solver` is the"
        " share the",
        "resolver at this run's k had to hand over; `no solve` is the share of shots with at least"
        " one defect",
        "where it handed over nothing, so the Mwpm(H) build, the solve and the harvest did not run at"
        " all.",
        "Structural, and measured outside every timed window.",
        "",
        f"{'d':>4} {'p':>8} {'T':>5} {'k':>3} {'mode':>7} {'shots':>8} "
        f"{'comps':>7} {'single':>7} {'pairs':>7} {'nontriv':>8} "
        f"{'largest':>8} {'max':>5} {'triv':>7} {'solver':>7} {'no solve':>9}",
    ]
    for result in rows:
        log = result["log"]
        meta = log.meta
        totals = result["components"]
        shots = totals["shots"] or 1
        defects = totals["h_defects"] or 1
        with_defects = totals["shots_with_defects"] or 1
        lines.append(
            f"{meta.get('d', '?'):>4} {meta.get('p', '?'):>8} {meta.get('T', '?'):>5} "
            f"{'-' if log.k is None else log.k:>3} "
            f"{meta.get('mode', '?'):>7} {totals['shots']:>8,} "
            f"{totals['components'] / shots:>7.2f} "
            f"{totals['singleton_components'] / shots:>7.2f} "
            f"{totals['pair_components'] / shots:>7.2f} "
            f"{totals['nontrivial_components'] / shots:>8.2f} "
            f"{totals['largest_component'] / shots:>8.2f} "
            f"{totals['max_component']:>5,} "
            f"{100 * totals['defects_in_trivial'] / defects:>6.1f}% "
            f"{100 * totals['defects_to_solver'] / defects:>6.1f}% "
            f"{100 * totals['solver_empty'] / with_defects:>8.1f}%"
        )
    missing = len(results) - len(rows)
    if missing:
        lines.append(
            f"  {missing} log(s) carried no component columns and are absent from this table"
            " (schema v2, or --no-component-stats)."
        )
    return lines


def warn_about_wall_clock(results):
    """Says on stderr that some of these numbers came off a wall clock. Returns how many logs did.

    On stderr and as a banner, for the same reason the profiler does it: stdout here is a list of
    written files and a table, routinely redirected or piped, and this is the one property of a
    campaign that cannot be recovered from the figures afterwards. A wall-clock log does not look
    broken — it looks like a slower machine.
    """
    guilty = [result["log"] for result in results if result["log"].meta.get("timer_thread_scoped") != "1"]
    if not guilty:
        return 0
    backends = sorted({log.meta.get("timer_backend", "?") for log in guilty})
    print(
        "\n"
        "  ##########################################################################\n"
        f"  WARNING: {len(guilty)} of {len(results)} log(s) were timed with a WALL CLOCK"
        f" ({', '.join(backends)}).\n"
        "  ##########################################################################\n"
        "  Their timings include whatever time the decoding thread spent off the CPU,\n"
        "  so every mean, tail and speedup drawn from them is an upper bound with the\n"
        "  scheduler's interference inside it. Re-run the profiler on a host with a\n"
        "  thread-scoped clock before quoting these numbers.\n",
        file=sys.stderr,
    )
    for log in guilty:
        print(f"    {log.stem}: timer {log.meta.get('timer_backend', '?')}", file=sys.stderr)
    return len(guilty)


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
        probes = sorted({result["log"].meta.get("preemption_probe", "?") for result in group})
        header = [
            "Spec-matching per-shot latency: stock decode on G, against the same decoder predicating on"
            " the sparsified graph H.",
            f"{len(group)} log file(s); timer backend {', '.join(timers)};"
            f" contamination probe {', '.join(probes)}.",
            "",
        ]
        if any(result["log"].meta.get("timer_thread_scoped") != "1" for result in group):
            header.insert(
                2,
                "WARNING: at least one log was timed on a WALL CLOCK, not a thread-scoped counter."
                " Its timings include\n         time the decoding thread spent off the CPU, so every"
                " mean, tail and speedup below that\n         came from it is an upper bound with the"
                " scheduler's interference inside it.",
            )
        if any(result["log"].meta.get("preemption_probe", "?") == "none" for result in group):
            header.insert(
                2,
                "NOTE: at least one log came from a platform with no contamination probe, so its"
                " `contam` of 0 means\n      nothing was measured rather than that the machine was"
                " quiet.",
            )
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
    parser.add_argument(
        "--no-components",
        dest="components",
        action="store_false",
        help="skip the component-structure figure drawn from each log's companion components_*.csv;"
        " the component table stays, since it comes from the per-shot rows",
    )
    parser.add_argument(
        "--k",
        metavar="LIST",
        help="comma-separated k values to draw, e.g. 0,2,4; logs at any other k are skipped."
        " Selection only — this script does not sweep k, the profiler runs one experiment per k and"
        " writes one log each (default: draw every log given)",
    )
    parser.add_argument(
        "--horizons",
        "--T",
        dest="horizons",
        metavar="LIST",
        help="comma-separated T values to draw; logs at any other T are skipped (default: all)",
    )
    parser.add_argument(
        "--no-k-overlay",
        dest="k_overlay",
        action="store_false",
        help="skip the per-(d, p, mode) figure that overlays the solver+harvest distribution of every"
        " (k, T) run given, which is how k values are compared",
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
    try:
        wanted_k = parse_number_list(args.k, int)
        wanted_T = parse_number_list(args.horizons, float)
    except ValueError as error:
        print(f"  {error}")
        return 1

    files = collect(args.paths)
    if not files:
        print("  nothing to plot.")
        return 1
    if args.out_dir:
        os.makedirs(args.out_dir, exist_ok=True)

    theme = THEMES[args.theme]
    results = []
    overlay_groups = {}
    for path in files:
        # `draw` is inside the same guard as `load`, because with the rows streamed rather than
        # pre-parsed a malformed row is first seen while the figure is being built. One bad file
        # still costs the run one file.
        try:
            log = load(path)
            # Selection happens after the header is read and before anything is drawn, so a log is
            # skipped by what it says it is rather than by what its filename looks like.
            if not selected(log, wanted_k, wanted_T):
                continue
            result = draw(log, args, theme)
        except (OSError, ValueError, KeyError) as error:
            print(f"  skipping {path}: {error}")
            continue
        if result is None:
            continue
        if log.k is not None:
            key = (log.meta.get("d"), log.meta.get("p"), log.meta.get("mode"))
            overlay_groups.setdefault(key, []).append(log)
        results.append(result)
        for out_path in result["written"]:
            print(f"  wrote {out_path}")

        # The component figure is drawn from the companion file, so a missing or malformed one costs
        # the run that figure and nothing else: the latency figure beside it is already written.
        if args.components:
            component_path = component_file_for(log)
            if component_path is None:
                if log.has_components:
                    print(f"  note: no components_*.csv beside {log.stem}; component figure skipped")
                continue
            try:
                for out_path in draw_components(log, read_component_file(component_path), args, theme):
                    print(f"  wrote {out_path}")
            except (OSError, ValueError, KeyError) as error:
                print(f"  skipping {component_path}: {error}")

    if not results:
        print("  nothing plotted.")
        return 1

    # The `k` comparison, drawn last because it re-reads the logs of a group and holds one campaign's
    # series at a time rather than every campaign's at once.
    if args.k_overlay:
        for key in sorted(overlay_groups, key=lambda item: tuple(str(part) for part in item)):
            try:
                for out_path in draw_k_overlay(overlay_groups[key], args, theme):
                    print(f"  wrote {out_path}")
            except (OSError, ValueError, KeyError) as error:
                print(f"  skipping the k overlay for d={key[0]} p={key[1]} {key[2]}: {error}")

    warn_about_wall_clock(results)
    for out_path in write_tables(results, args):
        print(f"  wrote {out_path}")
    if args.table:
        print()
        print("\n".join(table_lines(results, args)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
