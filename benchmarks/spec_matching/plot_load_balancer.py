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

"""Draws the k-core critical-path latencies `load_balancer_profiler` logs.

    python benchmarks/spec_matching/plot_load_balancer.py path/to/run_dir

A run directory is `shots.csv` + `agg.json` + `run.log`, exactly as the profiler writes it. Three
kinds of figure come out of it, plus a text table.

## The two series, and why only two

Per shot the profiler logs three latencies, and `run.log` states the definitions:

    sparse_k = uf + balance + max_c (build_c + solve_c + extract_c) + combine
    fallback = stock exact decode on G, from the full syndrome, on its own core
    system   = escalated ? fallback : min(fallback, sparse_k)

The histograms draw **`fallback` and `sparse_k`** — the two inputs the model combines, and the pair
whose separation is the whole point of a `k` sweep. `system` is not drawn, for the same reason
`plot_latency_histograms.py` does not draw its serial series: it is a per-shot `min` of the two, so
its bulk is a second tracing of whichever side is faster, and its curve would sit on top of
`sparse_k` everywhere except the escalating tail. Its mean is on the figure as a reference rule, and
its mean, max and percentiles are in the table. Nothing about it is hidden; it is simply not a third
body of ink over the two it is derived from.

There is also a colour reason, and it is worth stating because it drove the layout. The corpus
palette's slot 2 (green) against slot 1 (orange) is ΔE 4.4 under protanopia — below the ΔE 6 floor,
so no amount of direct labelling makes a third filled series legible. Slot 0 (blue) against slot 1
(orange) is ΔE 24.7 and passes every check in both themes. So every figure here carries **at most
two hues**, and every comparison that needs more than two categories is a small multiple instead of
an overlay. That is why the `k` sweep is a row of panels rather than `|k|` curves on one axis.

## The figures

  * **one histogram per `(d, p, T, k)` cell** — `fallback` against `sparse_k`, means direct-labelled,
    the `system` mean as an ink reference rule. The caption carries the speedup, the escalation
    count, the mean LPT imbalance, how many shots fell past the right edge, and which timer produced
    the numbers;
  * **one `k`-sweep figure per `(d, p, T)`** — the same two histograms as small multiples, one panel
    per `k`, on shared bins and shared axes so the panels are comparable by eye. `fallback` is the
    internal control here: it is the same shots decoded the same way in every panel, so it should
    *not* move with `k`, and a panel where it does is telling you about the machine rather than about
    the balancer;
  * **one module figure per `(d, p, T)`** — where the critical path actually goes, one panel per `k`,
    the six timed regions as horizontal bars on a shared linear axis. Linear and not log, because a
    bar encodes length from zero and a log axis would inflate `combine` — twelve nanoseconds — into
    looking comparable with `extract_crit`. `combine`'s bar is a sliver and its direct label says so;
    that is the correct reading, not a defect of the scale.

The six module means sum to the `sparse_k` mean exactly, by the definition above. The table checks
that rather than assuming it, and says so when it fails.

## Overhead, and why the columns are not read-symmetric

Rows are raw, uncorrected tick deltas; `run.log` carries `timer_overhead_ns` and says the analysis is
where the subtraction happens. This is that analysis, and the subtraction is **not** a flat one.

A bracketed region costs about one clock read, so a column carries one `timer_overhead_ns` per
*interval summed into it* — and `sparse_k` sums six intervals (`uf`, `balance`, `build_crit`,
`solve_crit`, `extract_crit`, `combine`) while `fallback` is one. Subtracting a flat overhead from
both would leave five of them in `sparse_k` and bias the speedup against the sparse path.
`--subtract-overhead` does the per-column thing: `INTERVALS` below is the count per column, five
rather than six on a shot whose critical core truncated and therefore had no extraction. Off by
default, so the raw numbers are what you get unless you ask; the caption always says which is on the
page.

On this corpus's Linux `rdpmc` backend the correction is tens of nanoseconds and changes nothing. On
macOS the same profiler lands on `clock_gettime(CLOCK_THREAD_CPUTIME_ID)` at roughly 200 ns a read,
where five intervals is about a microsecond — the same order as the whole `sparse_k` at small `d`.
That is the case the flag exists for.

## What travels with every number

The profiler's `run.log` states two model assumptions and a caveat, and they are not optional
context: the `k` loads are solved **serially on one thread** and reported as a max, no cost is
modelled for moving a core's slice to it, and a load solved straight after the previous one sees a
warmer cache than a real core would. The caption repeats the caveat and the table repeats all three.
A figure from this script without them is a figure being read as something it is not.

A million-shot campaign is read a chunk at a time and never held as Python objects. The profiler
writes cells in `(d, p, T, k)` order with `k` innermost, so all the `k` of one `(d, p, T)` are
contiguous in the file and the reader holds one such group at a time — `|k|` cells, not the campaign.
"""

import argparse
import itertools
import json
import math
import os
import sys
import textwrap

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

# The corpus palette, taken from `plot_latency_histograms.py` so a reader moving between the two
# scripts' output sees one visual language.
#
# Only slots 0 and 1 are used for series, and that is a checked decision rather than a stylistic one:
# `validate_palette.js` puts slot 1 against slot 2 at ΔE 4.4 under protanopia (below the ΔE 6 floor,
# so secondary encoding cannot rescue it), and slot 0 against slot 1 at ΔE 24.7, passing every check
# against both surfaces. Anything here needing more than two categories is a small multiple.
THEMES = {
    "light": {
        "surface": "#fcfcfb",
        "text_primary": "#0b0b0b",
        "text_secondary": "#52514e",
        "grid": "#dcdcd8",
        # (fallback on G, the sparse k-core path).
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

FALLBACK_LABEL = "fallback: stock decode on G"
SPARSE_LABEL = "sparse_k: k-core critical path on H"

# Rows are parsed in blocks of this many, so the transient cost of a read is set by the block and not
# by the length of the campaign.
CHUNK_ROWS = 100_000

# The only columns any series needs. `shot`, `n_def`, `pred_*`, `crit_*` and the per-module ticks are
# never read from the CSV — the module figure takes its means from `agg.json`, which the profiler
# already accumulated in 128-bit integers.
#
# `ticks_system` is read even though it is recomputed here, so that the recomputation can be checked
# against what the profiler logged rather than trusted.
NEEDED = (
    "d",
    "p",
    "T",
    "k",
    "escalated",
    "imbalance",
    "ticks_extract_crit",
    "ticks_sparse_k",
    "ticks_fallback",
    "ticks_system",
)

# Timed intervals summed into each column, which is how many `timer_overhead_ns` it carries. See the
# module docstring: this is the whole reason the correction is per-column and not flat.
#
# `sparse_k` is six on an ordinary shot and five when the critical core truncated, since a truncated
# core runs no extraction and `ticks_extract_crit` is then exactly zero.
INTERVALS_FALLBACK = 1
INTERVALS_SPARSE_K = 6
INTERVALS_SPARSE_K_TRUNCATED_CRIT = 5
INTERVALS_MODULE = 1

# In the order the model composes them, which is the order the module panels read top to bottom.
MODULES = ("uf", "balance", "build_crit", "solve_crit", "extract_crit", "combine")


class Run:
    """One profiler output directory: its `run.log` metadata and its `agg.json` cells.

    Deliberately not its rows. The rows stay on disk and are streamed by `stream_groups`; a `Run` is
    small enough that `main` can hold one per directory for the whole invocation, which is what the
    summary table at the end needs.
    """

    def __init__(self, path, meta, agg):
        self.path = path
        self.meta = meta
        self.agg = agg
        # `agg.json` keyed the way the CSV rows are keyed, so a cell's aggregate can be found from
        # the four numbers the streaming reader already has.
        self.cells = {
            (float(cell["d"]), float(cell["p"]), float(cell["T"]), float(cell["k"])): cell
            for cell in agg.get("cells", ())
        }

    @property
    def shots_csv(self):
        return os.path.join(self.path, "shots.csv")

    @property
    def name(self):
        return os.path.basename(os.path.normpath(self.path)) or self.path

    @property
    def timer(self):
        return self.meta.get("hires_timer_name", "?")

    @property
    def thread_scoped(self):
        return self.meta.get("timer_is_thread_scoped") == "1"

    @property
    def overhead_ns(self):
        """`timer_overhead_ns=9.2114 (10^6 back-to-back ...)` — the number, not the parenthetical."""
        text = self.meta.get("timer_overhead_ns", "")
        try:
            return float(text.split()[0])
        except (IndexError, ValueError):
            return float("nan")

    @property
    def alpha(self):
        try:
            return float(self.meta.get("alpha", "nan"))
        except ValueError:
            return float("nan")


def read_run_log(path):
    """The top-level `key=value` lines of `run.log`.

    The indented `corpus` / `cell` lines below them are per-cell provenance the figures do not read —
    `agg.json` carries the same cells in a form that does not need parsing — so they are skipped by
    the leading-whitespace test rather than half-parsed. Values may themselves contain `=`, which is
    why this partitions on the first one instead of splitting.
    """
    meta = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line[:1].isspace():
                continue
            key, sep, value = line.partition("=")
            if sep:
                meta[key.strip()] = value.strip()
    return meta


def load(path):
    """Reads one run directory's metadata. The shot rows are left on disk for `stream_groups`.

    Everything that would make a figure wrong is settled here, before a figure is started: that the
    three files exist, that `shots.csv` has the columns the series are built from, and that
    `agg.json` parses. A malformed run is reported by name like any other bad input.
    """
    if os.path.isfile(path) and os.path.basename(path) == "shots.csv":
        path = os.path.dirname(path) or "."
    if not os.path.isdir(path):
        raise ValueError("not a directory (expected a load_balancer_profiler --out dir)")
    shots = os.path.join(path, "shots.csv")
    log_path = os.path.join(path, "run.log")
    agg_path = os.path.join(path, "agg.json")
    for needed in (shots, log_path, agg_path):
        if not os.path.isfile(needed):
            raise ValueError(f"{path}: no {os.path.basename(needed)}")
    with open(shots) as handle:
        header = handle.readline()
        if not header.strip():
            raise ValueError(f"{shots}: no column header")
        columns = [name.strip() for name in header.strip().split(",")]
        if handle.readline() == "":
            raise ValueError(f"{shots}: no shot rows")
    missing = [name for name in NEEDED if name not in columns]
    if missing:
        raise ValueError(f"{shots}: header is missing {', '.join(missing)}")
    meta = read_run_log(log_path)
    with open(agg_path) as handle:
        agg = json.load(handle)
    if agg.get("binary") not in (None, "load_balancer_profiler"):
        raise ValueError(f"{agg_path}: binary is {agg['binary']}; this script reads load_balancer_profiler")
    run = Run(path, meta, agg)
    run.columns = columns
    return run


def parse_number_list(text, cast):
    """`--k 3,5` / `--T 1.5,2` as a set, or None when the flag was not given."""
    if text is None:
        return None
    values = set()
    for item in text.split(","):
        item = item.strip()
        if item:
            values.add(cast(item))
    return values or None


class Cell:
    """One `(d, p, T, k)` cell's drawable series, accumulated a chunk at a time.

    Holds three `float64` arrays in microseconds and four scalars. It never holds the module columns
    or the balancer statistics: those are means, and `agg.json` already has them accumulated over the
    whole cell in 128-bit integers, which is a better number than anything re-summed here in floats.
    """

    def __init__(self, key, agg):
        self.key = key
        self.d, self.p, self.T, self.k = key
        self.agg = agg
        self._fallback = []
        self._sparse = []
        self._system = []
        self.escalated = 0
        self.imbalance_sum = 0.0
        self.shots = 0
        # Rows where the recomputed `system` disagreed with the column the profiler logged. Only
        # meaningful on the raw path; the corrected path recomputes from corrected inputs by design.
        self.system_mismatch = 0
        # Corrected values that came out below zero and were clamped. A negative latency is not a
        # thing, but silently clamping one is worse than saying how many.
        self.clamped = 0

    def add(self, fallback, sparse, system, escalated, imbalance):
        self._fallback.append(fallback)
        self._sparse.append(sparse)
        self._system.append(system)
        self.escalated += int(np.count_nonzero(escalated))
        self.imbalance_sum += float(imbalance.sum())
        self.shots += int(fallback.size)

    def finish(self):
        empty = np.empty(0, dtype=np.float64)
        self.fallback = np.concatenate(self._fallback) if self._fallback else empty
        self.sparse = np.concatenate(self._sparse) if self._sparse else empty
        self.system = np.concatenate(self._system) if self._system else empty
        self._fallback.clear()
        self._sparse.clear()
        self._system.clear()
        return self

    @property
    def imbalance(self):
        return self.imbalance_sum / self.shots if self.shots else float("nan")

    @property
    def title(self):
        return f"d = {self.d:g},  p = {self.p:g},  T = {self.T:g},  k = {self.k:g}"

    @property
    def slug(self):
        return f"d{self.d:g}_p{self.p:g}_T{self.T:g}_k{self.k:g}"


def group_slug(key):
    d, p, T = key
    return f"d{d:g}_p{p:g}_T{T:g}"


def group_title(key):
    d, p, T = key
    return f"d = {d:g},  p = {p:g},  T = {T:g}"


def stream_groups(run, args):
    """Yields `((d, p, T), [Cell, ...])` in file order, one `(d, p, T)` group at a time.

    The profiler writes cells in `(d, p, T, k)` order with `k` innermost, so every `k` of a group is
    contiguous and a group can be completed, handed to the figures and released before the next one
    starts. What is resident is `|k|` cells — a few megabytes at a hundred thousand shots — rather
    than the campaign.

    That ordering is a property of the profiler's loop nest, not an assumption about the file: a
    group is flushed when the `(d, p, T)` triple changes, so a file that interleaved them would
    produce a group per run of rows and the figures would be visibly wrong rather than silently so.
    Runs within a chunk are found by comparing adjacent key rows, which also handles a chunk boundary
    landing mid-cell without a special case.

    The correction, when asked for, is applied here — while the block is still a block — so a
    corrected campaign costs one extra vector op per chunk rather than a second pass.
    """
    at = {name: position for position, name in enumerate(NEEDED)}
    usecols = tuple(run.columns.index(name) for name in NEEDED)
    overhead = run.overhead_ns if args.subtract_overhead else 0.0
    if args.subtract_overhead and not math.isfinite(overhead):
        raise ValueError(f"{run.path}: --subtract-overhead but run.log has no usable timer_overhead_ns")

    current_key = None
    pending = {}

    def flush():
        if current_key is None or not pending:
            return None
        return current_key, [pending[k].finish() for k in sorted(pending)]

    with open(run.shots_csv) as handle:
        handle.readline()
        while True:
            text = list(itertools.islice(handle, CHUNK_ROWS))
            if not text:
                break
            block = np.loadtxt(text, delimiter=",", usecols=usecols, dtype=np.float64, ndmin=2)
            del text
            if not block.size:
                continue
            keys = block[:, (at["d"], at["p"], at["T"], at["k"])]
            change = np.any(keys[1:] != keys[:-1], axis=1)
            starts = np.concatenate(([0], np.flatnonzero(change) + 1))
            ends = np.concatenate((starts[1:], [keys.shape[0]]))
            for start, end in zip(starts, ends):
                rows = block[start:end]
                key = tuple(float(value) for value in keys[start])
                group = key[:3]
                if args.d is not None and key[0] not in args.d:
                    continue
                if args.T is not None and key[2] not in args.T:
                    continue
                if args.k is not None and key[3] not in args.k:
                    continue
                if group != current_key:
                    done = flush()
                    if done is not None:
                        yield done
                    current_key = group
                    pending = {}
                cell = pending.get(key[3])
                if cell is None:
                    cell = Cell(key, run.cells.get(key))
                    pending[key[3]] = cell

                escalated = rows[:, at["escalated"]] != 0
                fallback = rows[:, at["ticks_fallback"]]
                sparse = rows[:, at["ticks_sparse_k"]]
                if overhead:
                    # Five intervals rather than six when the critical core truncated: it ran no
                    # extraction, so `ticks_extract_crit` is exactly zero and the region was never
                    # opened. See INTERVALS_* above.
                    truncated_crit = rows[:, at["ticks_extract_crit"]] == 0
                    sparse_intervals = np.where(
                        truncated_crit, INTERVALS_SPARSE_K_TRUNCATED_CRIT, INTERVALS_SPARSE_K
                    )
                    fallback = fallback - INTERVALS_FALLBACK * overhead
                    sparse = sparse - sparse_intervals * overhead
                    below = int(np.count_nonzero(fallback < 0) + np.count_nonzero(sparse < 0))
                    if below:
                        cell.clamped += below
                        np.clip(fallback, 0.0, None, out=fallback)
                        np.clip(sparse, 0.0, None, out=sparse)
                # The definition `run.log` states, recomputed here rather than read, so that the
                # corrected path and the raw path produce it the same way.
                system = np.where(escalated, fallback, np.minimum(fallback, sparse))
                if not overhead:
                    cell.system_mismatch += int(np.count_nonzero(system != rows[:, at["ticks_system"]]))
                cell.add(
                    fallback / 1000.0,
                    sparse / 1000.0,
                    system / 1000.0,
                    escalated,
                    rows[:, at["imbalance"]],
                )
            del block
    done = flush()
    if done is not None:
        yield done


def wrap_caption(text, width=135):
    """Re-wraps a caption to the figure width, keeping the author's own line breaks.

    A caption that runs off the right edge is a caption that was not read, and these say things — the
    serial-emulation model, the warm-cache caveat, which shots fell past the edge — that the figure
    is not honest without.
    """
    return "\n".join(
        textwrap.fill(line, width=width, break_long_words=False) if line else ""
        for line in text.split("\n")
    )


def percentile(values, fraction):
    """Only ever used to choose where the x axis stops; nothing reported here is a percentile.

    A partial sort, not a full one: the axis needs one order statistic, and selecting it costs a
    linear pass over one scratch copy instead of sorting the whole pool. The percentiles that *are*
    reported come from `agg.json`, where the profiler computed them over the raw ticks.
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
    """Shared edges for every series and every panel — histograms on different bins are not
    comparable, and the `k` sweep exists precisely to be compared across panels.

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


def axis_limits(cells, args):
    """The x range and bins shared by every panel of a group, from the two drawn series only.

    `system` does not widen the axis: it is bounded above by `fallback` on every shot, so an axis
    that covers the two inputs covers it too.
    """
    sides = []
    for cell in cells:
        sides.extend((cell.fallback, cell.sparse))
    pooled = np.concatenate(sides) if sides else np.empty(0)
    x_max = max(
        percentile(pooled, args.x_max_percentile / 100.0),
        1.15 * max((mean_of(side) for side in sides if side.size), default=1.0),
    )
    del pooled
    return x_max, bin_edges(tuple(sides), x_max, args.bins, args.log_x)


def style_axes(ax, theme, log_x, log_y):
    """The recessive frame every panel wears: grid behind the data, two spines, muted ticks."""
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
    ax.grid(axis="y", color=theme["grid"], linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(theme["grid"])
    ax.tick_params(colors=theme["text_secondary"], labelsize=9)


def draw_histograms(ax, cell, edges, theme, log_x=False, label=True):
    """The two series into one axes, and the `system` mean as a reference rule.

    Binned once and drawn twice per series. `ax.hist` would re-bin the whole array for each of the
    two passes and keep a copy inside the axes; past a few hundred thousand shots that is the
    difference between a figure and a swap storm, and the counts are identical either way.

    Filled at low alpha with a 2px edge, so the overlap is legible either way round.
    """
    for values, colour, series_label in (
        (cell.fallback, theme["series"][0], FALLBACK_LABEL),
        (cell.sparse, theme["series"][1], SPARSE_LABEL),
    ):
        counts, _ = np.histogram(values, bins=edges)
        ax.stairs(counts, edges, color=colour, alpha=0.45, fill=True, label=series_label, zorder=2)
        ax.stairs(counts, edges, color=colour, linewidth=2.0, zorder=3)

    def side_of(value):
        """Which side of its rule a label goes, as `(ha, x offset)`.

        A label is set to the right of its line by default, and flipped to the left once the line is
        past ~60% of the axis — where "to the right" runs it into the legend box or off the edge
        entirely. The means sit close together at most grid points, so this is decided per label
        rather than once for the figure.
        """
        low, high = (math.log10(edges[0]), math.log10(edges[-1])) if log_x else (edges[0], edges[-1])
        at = (math.log10(value) if log_x else value) if value > 0 else low
        fraction = (at - low) / (high - low) if high > low else 0.0
        return ("right", -6) if fraction > 0.6 else ("left", 5)

    # Selective direct labels: the two means the figure is about, and the derived `system` mean as a
    # rule in text ink rather than a third series colour — it is a reference, not an identity.
    annotations = (
        (mean_of(cell.fallback), theme["series"][0], theme["text_primary"], 0.97, "mean {:,.2f} us"),
        (mean_of(cell.sparse), theme["series"][1], theme["text_primary"], 0.88, "mean {:,.2f} us"),
        (
            mean_of(cell.system),
            theme["text_secondary"],
            theme["text_secondary"],
            0.79,
            "system mean {:,.2f} us",
        ),
    )
    for index, (mean, rule_colour, text_colour, height, template) in enumerate(annotations):
        system = index == 2
        ax.axvline(
            mean,
            color=rule_colour,
            linewidth=1.2 if system else 1.5,
            linestyle=(0, (1, 2)) if system else (0, (4, 3)),
            zorder=4,
        )
        if not label:
            continue
        ha, offset = side_of(mean)
        ax.annotate(
            template.format(mean),
            xy=(mean, height),
            xycoords=("data", "axes fraction"),
            ha=ha,
            va="top",
            fontsize=9,
            color=text_colour,
            xytext=(offset, 0),
            textcoords="offset points",
            zorder=5,
        )


def model_caveat(run):
    """The two assumptions and the caveat that travel with every number, in one short line.

    Taken from `run.log` rather than restated here, so a run written under a changed model cannot be
    captioned with this script's memory of the old one.
    """
    parts = []
    if run.meta.get("k_core_critical_path_model", "").startswith("1"):
        parts.append(
            "k_core_critical_path_model: the k loads are solved SERIALLY on one thread and reported"
            " as a max; no cost is modelled for moving a core's slice to it"
        )
    if run.meta.get("warm_cache_caveat", "").startswith("1"):
        parts.append("a load solved straight after the previous one sees a warmer cache than a real core would")
    return "  ·  ".join(parts)


def timer_note(run, args):
    note = f"timer {run.timer}"
    if not run.thread_scoped:
        note += " — WALL CLOCK, timings include off-CPU time"
    if args.subtract_overhead:
        note += (
            f"  ·  overhead-corrected: {run.overhead_ns:.3f} ns x {INTERVALS_SPARSE_K} intervals off"
            f" sparse_k ({INTERVALS_SPARSE_K_TRUNCATED_CRIT} on a truncated critical core) and x"
            f" {INTERVALS_FALLBACK} off fallback"
        )
    else:
        note += "  ·  raw uncorrected ticks (--subtract-overhead applies the per-column correction)"
    return note


def save(fig, args, run, stem, theme):
    written = []
    out_dir = args.out_dir or run.path
    os.makedirs(out_dir, exist_ok=True)
    for extension in args.formats:
        out_path = os.path.join(out_dir, f"{stem}.{extension}")
        fig.savefig(out_path, dpi=args.dpi, facecolor=theme["surface"])
        written.append(out_path)
    plt.close(fig)
    return written


def draw_cell(run, cell, edges, args, theme):
    """One `(d, p, T, k)` cell: the two distributions, their means, and what the caption must state."""
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    fig.patch.set_facecolor(theme["surface"])
    ax.set_facecolor(theme["surface"])

    draw_histograms(ax, cell, edges, theme, log_x=args.log_x)
    style_axes(ax, theme, args.log_x, args.log_y)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_xlabel("per-shot latency (microseconds)", color=theme["text_secondary"], fontsize=10)
    ax.set_ylabel("shots", color=theme["text_secondary"], fontsize=10)
    ax.set_title(cell.title, color=theme["text_primary"], fontsize=13, loc="left", pad=14)
    legend = ax.legend(frameon=False, loc="upper right", fontsize=10)
    for text in legend.get_texts():
        text.set_color(theme["text_primary"])

    fallback = stats_of(cell.fallback)
    sparse = stats_of(cell.sparse)
    system = stats_of(cell.system)
    speedup = fallback["mean"] / system["mean"] if system["mean"] else float("nan")
    clipped = int(
        np.count_nonzero(cell.fallback > edges[-1]) + np.count_nonzero(cell.sparse > edges[-1])
    )
    caption = (
        f"{fallback['n']:,} shots  ·  speedup {speedup:.2f}x = fallback mean / system mean  ·  "
        f"max {fallback['max']:,.2f} vs {sparse['max']:,.2f} us  ·  "
        f"{cell.escalated:,} escalated  ·  mean LPT imbalance {cell.imbalance:.4f} (1.0 = perfect)  ·  "
        f"{clipped:,} beyond the right edge\n"
        "system = min(fallback, sparse_k) on an ordinary shot and fallback alone on an escalating"
        " one, since that decode was already running; it is not drawn because it is a per-shot min of"
        f" the two that are\n{timer_note(run, args)}"
    )
    caveat = model_caveat(run)
    if caveat:
        caption += f"\n{caveat}"
    if cell.clamped:
        caption += f"\n{cell.clamped:,} corrected values fell below zero and were clamped to it"
    if cell.system_mismatch:
        caption += (
            f"\nWARNING: {cell.system_mismatch:,} rows where the logged ticks_system disagrees with"
            " the definition recomputed from ticks_fallback and ticks_sparse_k"
        )
    caption = wrap_caption(caption)
    fig.text(0.012, 0.005, caption, color=theme["text_secondary"], fontsize=8.5, va="bottom")
    fig.tight_layout(rect=(0, min(0.34, 0.02 + 0.033 * (caption.count("\n") + 1)), 1, 1))
    return save(fig, args, run, f"lbp_{cell.slug}", theme), {
        "cell": cell,
        "fallback": fallback,
        "sparse": sparse,
        "system": system,
        "speedup": speedup,
        "clipped": clipped,
    }


def draw_k_sweep(run, key, cells, edges, args, theme):
    """The `k` sweep as small multiples: one panel per `k`, shared bins, shared axes.

    Small multiples rather than `|k|` overlaid curves because the palette will not carry more than
    two hues — see the module docstring. It is also the better read: the question is what the
    `sparse_k` body does as `k` grows, and a panel per `k` answers it without asking anyone to
    untangle four curves.

    `fallback` is repeated in every panel deliberately. It is the same shots decoded the same way at
    every `k`, so it is the control: it should sit still, and a panel where it has moved is telling
    you about the machine rather than about the balancer.
    """
    if len(cells) < 2:
        return None
    columns = len(cells)
    fig, axes = plt.subplots(
        1, columns, figsize=(min(4.6 * columns, 18.0), 5.2), sharex=True, sharey=True
    )
    axes = np.atleast_1d(axes)
    fig.patch.set_facecolor(theme["surface"])

    for ax, cell in zip(axes, cells):
        ax.set_facecolor(theme["surface"])
        draw_histograms(ax, cell, edges, theme, log_x=args.log_x, label=False)
        style_axes(ax, theme, args.log_x, args.log_y)
        ax.set_xlim(edges[0], edges[-1])
        speedup = mean_of(cell.fallback) / mean_of(cell.system) if mean_of(cell.system) else float("nan")
        ax.set_title(
            f"k = {cell.k:g}\nsparse_k mean {mean_of(cell.sparse):,.2f} us  ·  speedup {speedup:.2f}x",
            color=theme["text_primary"],
            fontsize=10,
            loc="left",
            pad=10,
        )
        ax.set_xlabel("per-shot latency (microseconds)", color=theme["text_secondary"], fontsize=9)
    axes[0].set_ylabel("shots", color=theme["text_secondary"], fontsize=10)

    fig.suptitle(
        f"{group_title(key)}  —  what k buys",
        color=theme["text_primary"],
        fontsize=13,
        x=0.012,
        y=0.985,
        ha="left",
    )
    # At the figure level rather than inside a panel: the identity is the same in every panel, and a
    # legend box inside one of them sits on that panel's data and on no other's.
    handles, labels = axes[0].get_legend_handles_labels()
    legend = fig.legend(
        handles, labels, frameon=False, fontsize=9, loc="upper right", bbox_to_anchor=(0.995, 1.0)
    )
    for text in legend.get_texts():
        text.set_color(theme["text_primary"])
    caption = (
        "One panel per k, shared bins and shared axes. The dashed rules are the two means; the dotted"
        " rule is the system mean. fallback is the control — the same shots decoded the same way in"
        " every panel, so it should not move with k.\n"
        f"{timer_note(run, args)}"
    )
    caveat = model_caveat(run)
    if caveat:
        caption += f"\n{caveat}"
    caption = wrap_caption(caption)
    fig.text(0.012, 0.005, caption, color=theme["text_secondary"], fontsize=8.5, va="bottom")
    fig.tight_layout(rect=(0, min(0.3, 0.02 + 0.031 * (caption.count("\n") + 1)), 1, 0.94))
    return save(fig, args, run, f"lbp_ksweep_{group_slug(key)}", theme)


def draw_modules(run, key, cells, args, theme):
    """Where the critical path goes: the six timed regions, one panel per `k`, on a shared axis.

    The means come from `agg.json`, where the profiler accumulated them in 128-bit integers over the
    whole cell, rather than being re-summed here in floats.

    Linear, not log. A bar encodes length from a zero baseline, and on a log axis `combine` — twelve
    nanoseconds against `extract_crit`'s tens of thousands — would grow a bar comparable to the ones
    that matter. Its sliver plus its direct label is the honest rendering of "this is nothing".

    One hue for every bar: the categories are identity, not magnitude, and the row label is what
    carries the identity. A ramp here would double-encode length as lightness.
    """
    usable = [cell for cell in cells if cell.agg and cell.agg.get("module_means")]
    if not usable:
        return None
    overhead = run.overhead_ns if args.subtract_overhead else 0.0
    columns = len(usable)
    fig, axes = plt.subplots(1, columns, figsize=(min(5.2 * columns, 18.0), 5.0), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    fig.patch.set_facecolor(theme["surface"])

    positions = np.arange(len(MODULES))
    series = []
    for cell in usable:
        means = cell.agg["module_means"]
        series.append([max(float(means.get(name, 0.0)) - INTERVALS_MODULE * overhead, 0.0) for name in MODULES])
    x_max = max(max(row) for row in series) if series else 1.0

    for ax, cell, values in zip(axes, usable, series):
        ax.set_facecolor(theme["surface"])
        # Top-to-bottom in the order the model composes the regions, which reads naturally only if
        # the y axis is inverted — matplotlib puts position 0 at the bottom.
        ax.barh(positions, values, height=0.62, color=theme["series"][1], zorder=2)
        ax.set_yticks(positions)
        ax.set_yticklabels(MODULES)
        ax.set_xlim(0, x_max * 1.32)
        # Five ticks, not matplotlib's default: at `d = 25` the labels run to six digits and the
        # default spacing butts them into one another.
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax.grid(axis="x", color=theme["grid"], linewidth=0.8, zorder=0)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(theme["grid"])
        ax.tick_params(colors=theme["text_secondary"], labelsize=9)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
        # A direct label on every bar, which is the one place that is right: six bars, three orders
        # of magnitude between them, and the smallest is the point.
        for position, value in zip(positions, values):
            ax.annotate(
                f"{value:,.0f}",
                xy=(value, position),
                ha="left",
                va="center",
                fontsize=8.5,
                color=theme["text_primary"],
                xytext=(4, 0),
                textcoords="offset points",
                zorder=4,
            )
        total = sum(values)
        ax.set_title(
            f"k = {cell.k:g}\nsum {total:,.0f} ns = sparse_k mean",
            color=theme["text_primary"],
            fontsize=10,
            loc="left",
            pad=10,
        )
        ax.set_xlabel("mean ticks (nanoseconds)", color=theme["text_secondary"], fontsize=9)

    # Once, not once per panel. The panels share a y axis, so inverting inside the loop inverts the
    # same axis `|k|` times and an even `|k|` silently leaves the modules bottom-up.
    axes[0].invert_yaxis()

    fig.suptitle(
        f"{group_title(key)}  —  where the critical path goes",
        color=theme["text_primary"],
        fontsize=13,
        x=0.012,
        y=0.985,
        ha="left",
    )
    caption = (
        "The six timed regions of the critical path, as means. They sum to the sparse_k mean by"
        " construction: sparse_k = uf + balance + max_c(build_c + solve_c + extract_c) + combine, and"
        " the build/solve/extract bars are the critical core's. uf, balance and combine are the"
        " serial parts — they do not shrink with k, and balance and combine grow with it.\n"
        f"{timer_note(run, args)}"
    )
    caption = wrap_caption(caption)
    fig.text(0.012, 0.005, caption, color=theme["text_secondary"], fontsize=8.5, va="bottom")
    fig.tight_layout(rect=(0, min(0.3, 0.02 + 0.033 * (caption.count("\n") + 1)), 1, 0.94))
    return save(fig, args, run, f"lbp_modules_{group_slug(key)}", theme)


def table_lines(run, results, args):
    """Everything the figures show, as text, so the plots are never the only read.

    One list of lines, so the same table can go to the terminal and to the file beside the figures
    without the two drifting apart.
    """
    lines = []
    lines.append(f"load_balancer_profiler — {run.path}")
    lines.append(f"  git {run.meta.get('git_hash', '?')}  ·  stim {run.meta.get('stim_version', '?')}")
    lines.append(
        f"  timer {run.timer}  ·  thread-scoped {'yes' if run.thread_scoped else 'NO — WALL CLOCK'}"
        f"  ·  overhead {run.overhead_ns:.3f} ns/read"
    )
    # `run.log` writes these as `value (parenthetical)`; the parenthetical is prose the table has no
    # room for, so it is trimmed here rather than wrapped.
    def first_word(key):
        return run.meta.get(key, "?").split()[0] if run.meta.get(key) else "?"

    lines.append(
        f"  seed {first_word('seed')}  ·  alpha {run.alpha:g}"
        f"  ·  shots/cell {first_word('shots_per_cell')}"
        f"  ·  warmup/cell {first_word('warmup_shots_per_cell')}"
    )
    lines.append(
        "  correction: "
        + (
            f"{run.overhead_ns:.3f} ns x per-column intervals subtracted"
            if args.subtract_overhead
            else "none — raw uncorrected ticks"
        )
    )
    for note in (run.meta.get("sparse_k"), run.meta.get("system")):
        if note:
            lines.append(f"  {note}")
    caveat = model_caveat(run)
    if caveat:
        for part in caveat.split("  ·  "):
            lines.append(f"  ! {part}")
    lines.append("")
    header = (
        f"{'d':>4} {'p':>8} {'T':>5} {'k':>3} {'shots':>9} {'escal':>7} {'imbal':>7} "
        f"{'fallback':>10} {'sparse_k':>10} {'system':>10} {'speedup':>8} "
        f"{'fb max':>10} {'sp max':>10} {'sys p99':>10} {'sys p99.9':>10}"
    )
    lines.append(header)
    lines.append("-" * len(header))
    for result in results:
        cell = result["cell"]
        agg = cell.agg or {}
        system_agg = agg.get("system", {})
        # The percentiles are the profiler's own, over raw ticks, so they are shown in raw
        # microseconds and are NOT moved by --subtract-overhead. Said in the footer rather than
        # silently mixed with corrected means.
        p99 = system_agg.get("p99")
        p999 = system_agg.get("p99_9")
        lines.append(
            f"{cell.d:>4.0f} {cell.p:>8g} {cell.T:>5g} {cell.k:>3.0f} {result['fallback']['n']:>9,} "
            f"{cell.escalated:>7,} {cell.imbalance:>7.4f} "
            f"{result['fallback']['mean']:>10,.2f} {result['sparse']['mean']:>10,.2f} "
            f"{result['system']['mean']:>10,.2f} {result['speedup']:>8.2f} "
            f"{result['fallback']['max']:>10,.2f} {result['sparse']['max']:>10,.2f} "
            f"{(p99 / 1000.0 if p99 is not None else float('nan')):>10,.2f} "
            f"{(p999 / 1000.0 if p999 is not None else float('nan')):>10,.2f}"
        )
    lines.append("")
    lines.append("  all latencies in microseconds; percentile columns are agg.json's, over RAW ticks")
    lines.append("  imbal = mean k * pred_max / pred_sum, the LPT prediction's imbalance; 1.0 = perfect")

    # The module split, and the identity it has to satisfy. Checked rather than assumed: if the six
    # regions stop summing to sparse_k the model in the caption is no longer the model in the data.
    lines.append("")
    module_header = f"{'d':>4} {'p':>8} {'T':>5} {'k':>3} " + " ".join(f"{name:>13}" for name in MODULES) + f" {'sum':>13} {'sparse_k':>13}"
    lines.append(module_header)
    lines.append("-" * len(module_header))
    for result in results:
        cell = result["cell"]
        agg = cell.agg or {}
        means = agg.get("module_means")
        if not means:
            continue
        values = [float(means.get(name, 0.0)) for name in MODULES]
        total = sum(values)
        sparse_mean = float(agg.get("sparse_k", {}).get("mean", float("nan")))
        flag = "" if not math.isfinite(sparse_mean) or abs(total - sparse_mean) <= 1.0 else "  <-- MISMATCH"
        lines.append(
            f"{cell.d:>4.0f} {cell.p:>8g} {cell.T:>5g} {cell.k:>3.0f} "
            + " ".join(f"{value:>13,.1f}" for value in values)
            + f" {total:>13,.1f} {sparse_mean:>13,.1f}{flag}"
        )
    lines.append("")
    lines.append("  module means in RAW nanoseconds, from agg.json's 128-bit accumulators")
    lines.append("  the six sum to the sparse_k mean by construction; MISMATCH means they no longer do")

    mismatches = sum(result["cell"].system_mismatch for result in results)
    if mismatches:
        lines.append("")
        lines.append(
            f"  WARNING: {mismatches:,} rows where the logged ticks_system disagrees with the"
            " definition recomputed from ticks_fallback, ticks_sparse_k and escalated"
        )
    clamped = sum(result["cell"].clamped for result in results)
    if clamped:
        lines.append(f"  {clamped:,} corrected values fell below zero and were clamped to it")
    return lines


def warn_about_wall_clock(run):
    if run.thread_scoped:
        return
    print(
        f"\n  WARNING: {run.path} was measured with a WALL CLOCK, backend"
        f" '{run.timer}'.\n  Every tick includes time the thread spent off the CPU and is an upper"
        " bound only.\n",
        file=sys.stderr,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "paths",
        nargs="+",
        help="load_balancer_profiler --out directories (shots.csv + agg.json + run.log)",
    )
    parser.add_argument("--out-dir", help="where the figures go; default is inside each run directory")
    parser.add_argument("--formats", default="png", help="comma-separated: png,pdf,svg")
    parser.add_argument("--dpi", type=int, default=160)
    parser.add_argument("--bins", type=int, default=90)
    parser.add_argument("--log-x", action="store_true", help="log latency axis; the escalation tail is long")
    parser.add_argument("--log-y", action="store_true", help="log count axis, which is where a thin tail shows up")
    parser.add_argument(
        "--x-max-percentile",
        type=float,
        default=99.5,
        help="where the latency axis stops; a view control, and the caption counts the shots past it",
    )
    parser.add_argument(
        "--subtract-overhead",
        action="store_true",
        help="subtract timer_overhead_ns x the interval count of each column (6 for sparse_k, 1 for"
        " fallback); off by default, and the caption always says which is on the page",
    )
    parser.add_argument("--d", help="only these distances, comma separated")
    parser.add_argument("--T", help="only these horizons, comma separated")
    parser.add_argument("--k", help="only these core counts, comma separated")
    parser.add_argument("--no-cells", action="store_true", help="skip the per-cell histograms")
    parser.add_argument("--no-ksweep", action="store_true", help="skip the k-sweep small multiples")
    parser.add_argument("--no-modules", action="store_true", help="skip the module breakdown")
    parser.add_argument("--theme", choices=sorted(THEMES), default="light")
    parser.add_argument("--table", action="store_true", help="print the summary table to stdout as well")
    args = parser.parse_args()

    args.formats = [item.strip() for item in args.formats.split(",") if item.strip()]
    args.d = parse_number_list(args.d, float)
    args.T = parse_number_list(args.T, float)
    args.k = parse_number_list(args.k, float)
    theme = THEMES[args.theme]

    failed = False
    for path in args.paths:
        try:
            run = load(path)
        except (ValueError, OSError, json.JSONDecodeError) as error:
            print(f"  {path}: {error}", file=sys.stderr)
            failed = True
            continue
        warn_about_wall_clock(run)
        print(f"{run.path}")
        results = []
        for key, cells in stream_groups(run, args):
            if not cells:
                continue
            x_max, edges = axis_limits(cells, args)
            for cell in cells:
                if cell.shots == 0:
                    continue
                if not args.no_cells:
                    written, result = draw_cell(run, cell, edges, args, theme)
                    for out_path in written:
                        print(f"  wrote {out_path}")
                else:
                    result = {
                        "cell": cell,
                        "fallback": stats_of(cell.fallback),
                        "sparse": stats_of(cell.sparse),
                        "system": stats_of(cell.system),
                        "clipped": 0,
                    }
                    result["speedup"] = (
                        result["fallback"]["mean"] / result["system"]["mean"] if result["system"]["mean"] else float("nan")
                    )
                results.append(result)
            if not args.no_ksweep:
                for out_path in draw_k_sweep(run, key, cells, edges, args, theme) or ():
                    print(f"  wrote {out_path}")
            if not args.no_modules:
                for out_path in draw_modules(run, key, cells, args, theme) or ():
                    print(f"  wrote {out_path}")
            # The group's arrays are released here, before the next one is read.
            for cell in cells:
                cell.fallback = cell.sparse = cell.system = np.empty(0)

        if not results:
            print("  no cells matched the selectors; nothing to plot", file=sys.stderr)
            failed = True
            continue
        lines = table_lines(run, results, args)
        out_dir = args.out_dir or run.path
        table_path = os.path.join(out_dir, "load_balancer_latency.txt")
        with open(table_path, "w") as handle:
            handle.write("\n".join(lines) + "\n")
        print(f"  wrote {table_path}")
        if args.table:
            print()
            print("\n".join(lines))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
