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

"""Draws the `summary.csv` / `hists.json` that `sparse_graph_stats` writes.

    python benchmarks/spec_matching/plot_sparse_graph_stats.py --in RESULTS_DIR [MORE_DIRS ...]
        [--out DIR] [--format png|pdf|svg] [--d 5,7,9] [--p 1e-3,5e-4] [--T 1.5,2]
        [--which all|escalation|structure|hists|joint] [--show]

Several `--in` dirs are concatenated (one campaign per `T`, or per machine); their rows may not
collide on `(d, p, T)`. Figures land in `--out` (default: the first `--in` dir + `/figures`)
alongside a `figures_index.md` that lists every one of them.

This script computes nothing the binary did not already compute, beyond the arithmetic needed to
draw: normalising a histogram, dividing a count by the shot count, and the Wilson interval on the
two counts of a `size x status` bin. It does not fit, extrapolate, smooth, or annotate with
conclusions, and it imports nothing from this repo -- it runs on a machine that has only the
artifacts.

Read `design/plot_sparse_graph_stats_design.md` for what each figure is meant to say, and
`design/branch_component_blossom_and_sparse_graph_stats.md` §3 for what the inputs mean.

Two notes on where the artifacts and the design doc differ, taken rather than papered over:

  * the doc says `benchmarks/two_phase/`; that directory is `benchmarks/spec_matching/` since the
    rename, and this file sits beside the binary that writes its inputs;
  * `hists.json` nests the cells under a `"cells"` key beside a top-level `"caps"` block. Both the
    nested and a bare `{"d=..,p=..,T=..": {...}}` mapping are accepted.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import textwrap
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

# The backend has to be chosen before `pyplot` is imported, and the only thing that changes it is
# `--show`. Peeking at `sys.argv` is what keeps a headless campaign machine from needing a display.
if "--show" not in sys.argv:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402  (must follow the backend selection above)
from matplotlib.lines import Line2D  # noqa: E402

# §6.3: the same artifacts must produce identical figures. `svg.hashsalt` fixes the element ids an
# SVG would otherwise seed from the process, and the two `metadata` blocks below drop the creation
# date out of PDF and SVG. Nothing else in the pipeline reads a clock.
matplotlib.rcParams["svg.hashsalt"] = "sparse_graph_stats"
matplotlib.rcParams["figure.max_open_warning"] = 0

#: The four histograms of §4.3, with the unit of one bin and what the bin counts are a fraction of.
HIST_KINDS = {
    "component_size_hist": ("component size", "count", "components"),
    "degree_hist": ("H degree (defect-defect)", "count", "H nodes"),
    "hop_diameter_hist": ("hop diameter", "count", "components"),
    "wint_diameter_hist": ("w_int diameter", "T_int/16", "components"),
}

#: §4.2's scalars, as `(design name, csv mean column, csv max column)`.
STRUCTURE_SCALARS = [
    ("n_defects", "mean_n_defects", "max_n_defects"),
    ("h_edges", "mean_h_edges", "max_h_edges"),
    ("h_boundary_edges", "mean_h_boundary_edges", "max_h_boundary_edges"),
    ("nodes_with_boundary_edge", "mean_nodes_with_boundary_edge", "max_nodes_with_boundary_edge"),
    ("num_components", "mean_num_components", "max_num_components"),
    ("largest_component_size", "mean_largest_component_size", "max_largest_component_size"),
]

#: §4.4's grid summary reads `P(TRUNCATED | size)` at these bins. `-1` is the overflow bin.
GRID_SIZE_BINS = [1, 2, 3, 4, -1]

MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]

#: §6.4, said once per figure rather than once per panel, where it lands on whatever the legend
#: did not want.
FLOOR_SENTENCE = "0 rates drawn at the resolution floor (hollow marker + arrow), never dropped"


# ---------------------------------------------------------------------------------------------
# Formatting. Every label and filename goes through these, so a `p` written `0.001` by the binary
# and one written `1e-3` by hand name the same figure.
# ---------------------------------------------------------------------------------------------


def fmt_p(value: float) -> str:
    """`1e-3` style, per §3. Falls back to `%g` for anything that is not a round mantissa."""
    if value == 0:
        return "0"
    exponent = math.floor(math.log10(abs(value)))
    mantissa = value / (10.0**exponent)
    if abs(mantissa - round(mantissa)) < 1e-9:
        mantissa = int(round(mantissa))
        return f"1e{exponent}" if mantissa == 1 else f"{mantissa}e{exponent}"
    return f"{value:g}"


def fmt_t(value: float) -> str:
    """`T` as given -- `2`, `1.5`."""
    return f"{value:g}"


def fmt_d(value) -> str:
    return f"{int(value)}"


def cell_key(d, p, t) -> tuple:
    """The canonical `(d, p, T)` key.

    `summary.csv` writes its floats at the default six significant digits and `hists.json` writes
    its key with `%g`, which is the same six. Rounding both to that precision is what lets a row and
    a cell that came from one run join without a tolerance.
    """
    return (int(d), float(f"{float(p):.6g}"), float(f"{float(t):.6g}"))


def label_of(**fixed) -> str:
    """`d=17_p=1e-3_T=2`, in that order, for filenames and index rows."""
    parts = []
    for name in ("d", "p", "T"):
        if name in fixed and fixed[name] is not None:
            value = fixed[name]
            formatted = {"d": fmt_d, "p": fmt_p, "T": fmt_t}[name](value)
            parts.append(f"{name}={formatted}")
    return "_".join(parts)


def wilson(successes, trials, z: float = 1.959963984540054):
    """The 95% Wilson interval, elementwise, matching the binary's own `wilson_interval`."""
    successes = np.asarray(successes, dtype=float)
    trials = np.asarray(trials, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        phat = np.where(trials > 0, successes / np.where(trials > 0, trials, 1.0), 0.0)
        denominator = 1.0 + z * z / np.where(trials > 0, trials, 1.0)
        centre = (phat + z * z / (2 * np.where(trials > 0, trials, 1.0))) / denominator
        half = (z / denominator) * np.sqrt(
            phat * (1 - phat) / np.where(trials > 0, trials, 1.0)
            + z * z / (4 * np.where(trials > 0, trials, 1.0) ** 2)
        )
    lo = np.where(successes == 0, 0.0, np.clip(centre - half, 0.0, 1.0))
    hi = np.where(successes == trials, 1.0, np.clip(centre + half, 0.0, 1.0))
    return np.where(trials > 0, lo, np.nan), np.where(trials > 0, hi, np.nan)


# ---------------------------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------------------------


@dataclass
class RunLog:
    """What §3 reads out of a `run.log`: the header, and the generator call per `(d, p)`."""

    path: Path
    git_hash: str = ""
    stim_version: str = ""
    calls: dict = field(default_factory=dict)
    present: bool = False


@dataclass
class Artifact:
    """One or more `sparse_graph_stats` output dirs, concatenated."""

    summary: pd.DataFrame
    cells: dict
    caps: dict
    logs: list
    source_dirs: list

    def row(self, d, p, t):
        key = cell_key(d, p, t)
        matches = self.summary[self.summary["_key"] == key]
        return None if matches.empty else matches.iloc[0]

    def cell(self, d, p, t):
        return self.cells.get(cell_key(d, p, t))

    def distances(self):
        return sorted(self.summary["d"].unique().tolist())

    def error_rates(self):
        return sorted(self.summary["p"].unique().tolist())

    def horizons(self):
        return sorted(self.summary["T"].unique().tolist())

    def has(self, column: str) -> bool:
        return column in self.summary.columns


def read_run_log(path: Path) -> RunLog:
    log = RunLog(path=path)
    if not path.is_file():
        return log
    log.present = True
    pending = None
    for line in path.read_text(errors="replace").splitlines():
        stripped = line.strip()
        if stripped.startswith("git_hash="):
            log.git_hash = stripped[len("git_hash=") :]
        elif stripped.startswith("stim_version="):
            log.stim_version = stripped[len("stim_version=") :]
        elif stripped.startswith("corpus=stim "):
            # `corpus=stim d=5 rounds=5 p=0.001` -- the labels the next generator_call belongs to.
            fields = {}
            for token in stripped.split():
                if "=" in token:
                    name, _, value = token.partition("=")
                    fields[name] = value
            try:
                pending = (int(fields["d"]), float(fields["p"]))
            except (KeyError, ValueError):
                pending = None
        elif stripped.startswith("generator_call=") and pending is not None:
            log.calls[cell_key(pending[0], pending[1], 0)[:2]] = stripped[len("generator_call=") :]
    return log


def load_artifacts(dirs) -> Artifact:
    """Reads every `--in` dir and concatenates them, asserting no `(d, p, T)` collision (§2)."""
    frames = []
    cells = {}
    caps = {}
    logs = []
    owner = {}
    for directory in dirs:
        directory = Path(directory)
        summary_path = directory / "summary.csv"
        if not summary_path.is_file():
            raise SystemExit(f"error: {summary_path} does not exist")
        frame = pd.read_csv(summary_path)
        for column in ("d", "p", "T", "shots"):
            if column not in frame.columns:
                raise SystemExit(f"error: {summary_path} has no '{column}' column")
        frame["_source"] = str(directory)
        frame["_key"] = [cell_key(d, p, t) for d, p, t in zip(frame["d"], frame["p"], frame["T"])]
        for key in frame["_key"]:
            if key in owner:
                d, p, t = key
                raise SystemExit(
                    f"error: (d={fmt_d(d)}, p={fmt_p(p)}, T={fmt_t(t)}) appears in both "
                    f"{owner[key]} and {directory}; --in dirs must not collide on (d, p, T)"
                )
            owner[key] = str(directory)
        frames.append(frame)

        hists_path = directory / "hists.json"
        if hists_path.is_file():
            document = json.loads(hists_path.read_text())
            # The binary nests the cells under "cells" beside a top-level "caps"; a bare mapping of
            # cell keys is accepted too, since that is what §3 of the design describes.
            raw_cells = document.get("cells", document)
            directory_caps = document.get("caps", {})
            for name, value in directory_caps.items():
                if name in caps and caps[name] != value:
                    raise SystemExit(
                        f"error: --in dirs disagree on {name} ({caps[name]} vs {value}); "
                        "histograms binned differently cannot share an axis"
                    )
                caps[name] = value
            for raw_key, cell in raw_cells.items():
                if not isinstance(cell, dict) or not raw_key.startswith("d="):
                    continue
                fields = dict(token.split("=", 1) for token in raw_key.split(","))
                key = cell_key(fields["d"], fields["p"], fields["T"])
                cell = dict(cell)
                cell["_source"] = str(directory)
                cells[key] = cell

        logs.append(read_run_log(directory / "run.log"))

    summary = pd.concat(frames, ignore_index=True)
    summary = summary.sort_values(["T", "d", "p"], kind="stable").reset_index(drop=True)
    return Artifact(
        summary=summary,
        cells=cells,
        caps=caps,
        logs=logs,
        source_dirs=[str(Path(directory)) for directory in dirs],
    )


def filter_artifact(artifact: Artifact, distances, rates, horizons) -> Artifact:
    """§6.5: filters only subset. Nothing is re-aggregated and no row is merged into another."""
    frame = artifact.summary
    if distances is not None:
        wanted = {int(value) for value in distances}
        frame = frame[frame["d"].astype(int).isin(wanted)]
    if rates is not None:
        wanted = {cell_key(0, value, 0)[1] for value in rates}
        frame = frame[[cell_key(0, value, 0)[1] in wanted for value in frame["p"]]]
    if horizons is not None:
        wanted = {cell_key(0, 0, value)[2] for value in horizons}
        frame = frame[[cell_key(0, 0, value)[2] in wanted for value in frame["T"]]]
    frame = frame.reset_index(drop=True)
    keys = set(frame["_key"])
    return Artifact(
        summary=frame,
        cells={key: cell for key, cell in artifact.cells.items() if key in keys},
        caps=artifact.caps,
        logs=artifact.logs,
        source_dirs=artifact.source_dirs,
    )


# ---------------------------------------------------------------------------------------------
# Figure furniture
# ---------------------------------------------------------------------------------------------


class Style:
    """Colour encodes `d`, marker encodes `T`, over the whole artifact rather than per figure.

    Fixing the assignment once means `d = 17` is the same colour on every page, which is the only
    way a reader can carry a series from one figure to the next.
    """

    def __init__(self, artifact: Artifact):
        # `tab10` is the colour-blind-acceptable palette §4 names.
        palette = plt.get_cmap("tab10")
        self.distances = artifact.distances()
        self.horizons = artifact.horizons()
        self.rates = artifact.error_rates()
        self._colour = {d: palette(i % 10) for i, d in enumerate(self.distances)}
        self._marker = {t: MARKERS[i % len(MARKERS)] for i, t in enumerate(self.horizons)}
        self._rate_colour = {p: palette(i % 10) for i, p in enumerate(self.rates)}
        self._horizon_colour = {t: palette(i % 10) for i, t in enumerate(self.horizons)}

    def colour(self, d):
        return self._colour[d]

    def marker(self, t):
        return self._marker[t]

    def rate_colour(self, p):
        return self._rate_colour[p]

    def horizon_colour(self, t):
        """For the two figures whose x axis is `d`, where colour is free and `T` is the series."""
        return self._horizon_colour[t]


def footer_text(artifact: Artifact, d=None, p=None) -> str:
    """§3: the git hash and the stim call, off `run.log`'s header. Says so when there is no log."""
    present = [log for log in artifact.logs if log.present]
    if not present:
        return "run.log missing"
    hashes = sorted({log.git_hash for log in present if log.git_hash})
    versions = sorted({log.stim_version for log in present if log.stim_version})
    parts = []
    parts.append("git=" + (",".join(h[:12] for h in hashes) if hashes else "unknown"))
    parts.append("stim=" + (",".join(versions) if versions else "unknown"))

    calls = {}
    for log in present:
        calls.update(log.calls)
    call = None
    if d is not None and p is not None:
        call = calls.get(cell_key(d, p, 0)[:2])
    if call is None:
        distinct = sorted(set(calls.values()))
        if len(distinct) == 1:
            call = distinct[0]
        elif distinct:
            call = f"{len(distinct)} stim generator calls over the (d, p) grid; see run.log"
    if call:
        if len(call) > 300:
            call = call[:297] + "... (full call in run.log)"
        parts.append(call)
    missing = len(artifact.logs) - len(present)
    if missing:
        parts.append(f"{missing} run.log missing")
    return "  ".join(parts)


def panel_grid(count: int, width: float = 4.6, height: float = 3.6, footer: bool = True):
    """A grid of at most three columns, sized so panels stay the same size as the count grows."""
    columns = min(count, 3)
    rows = math.ceil(count / columns)
    figure, axes = plt.subplots(
        rows,
        columns,
        figsize=(width * columns, height * rows + (0.7 if footer else 0.0)),
        squeeze=False,
    )
    flat = [axes[r][c] for r in range(rows) for c in range(columns)]
    for axis in flat[count:]:
        axis.set_visible(False)
    return figure, flat[:count]


def finish(figure, title: str, subtitle: str, footer: str):
    """Title, subtitle and footer. Titles name the quantity and the fixed variables, nothing else.

    The three are placed in inches off the edges rather than in figure fractions, so a one-row
    figure and a three-row one give them the same amount of room instead of the short one crushing
    its subtitle into its title.
    """
    height = figure.get_figheight()
    figure.suptitle(title, fontsize=12, y=1.0 - 0.28 / height)
    subtitle_lines = 0
    if subtitle:
        # An 8-point DejaVu line fits about 14 characters per inch of figure width; wrapping on that
        # keeps the long subtitles (the ones carrying the floor rule and the error-bar rule) inside
        # a narrow figure instead of running off both edges.
        subtitle = "\n".join(textwrap.wrap(subtitle, width=max(40, int(14 * figure.get_figwidth()))))
        subtitle_lines = subtitle.count("\n")
        figure.text(0.5, 1.0 - 0.62 / height, subtitle, ha="center", va="top", fontsize=8, color="0.3")
    wrapped = footer
    if len(wrapped) > 160:
        # Two lines at most; the rest of the call stays in run.log where it came from.
        split = wrapped.rfind(" ", 0, 160)
        wrapped = wrapped[:split] + "\n" + wrapped[split + 1 :]
    figure.text(0.006, 0.006, wrapped, ha="left", va="bottom", fontsize=5, color="0.45")
    top = 1.0 - ((0.78 + 0.16 * subtitle_lines) if subtitle else 0.46) / height
    figure.tight_layout(rect=(0.0, 0.30 / height, 1.0, top))
    return figure


def log_x_rates(axis, rates):
    """`p` on a log axis, ticked at the values that were actually run and labelled `1e-3` style."""
    axis.set_xscale("log")
    axis.set_xticks(list(rates))
    axis.set_xticklabels([fmt_p(value) for value in rates])
    axis.minorticks_off()


def plot_rate(axis, x, y, lo, hi, floor, colour, marker, label, linestyle="-"):
    """A rate series with §6.4's floor treatment.

    A zero rate is a measurement -- "no shot escalated in `n`" -- so it is drawn at the resolution
    floor `1/n` with a hollow marker and a downward arrow rather than dropped off a log axis. Its
    error bar keeps the Wilson upper end, which is the only end that says anything there.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    floor = np.asarray(floor, dtype=float)
    drawn = np.where(y > 0, y, floor)
    positive = y > 0

    axis.plot(x, drawn, linestyle=linestyle, color=colour, linewidth=1.2, zorder=2, label=label)
    if positive.any():
        lower = np.clip(y[positive] - np.asarray(lo, dtype=float)[positive], 0, None)
        upper = np.clip(np.asarray(hi, dtype=float)[positive] - y[positive], 0, None)
        axis.errorbar(
            x[positive],
            y[positive],
            yerr=np.vstack([lower, upper]),
            fmt=marker,
            markersize=4.5,
            color=colour,
            elinewidth=0.9,
            capsize=2,
            linestyle="none",
            zorder=3,
        )
    if (~positive).any():
        zero_x = x[~positive]
        zero_floor = floor[~positive]
        upper = np.clip(np.asarray(hi, dtype=float)[~positive] - zero_floor, 0, None)
        axis.errorbar(
            zero_x,
            zero_floor,
            yerr=np.vstack([np.zeros_like(upper), upper]),
            fmt=marker,
            markersize=5.0,
            markerfacecolor="none",
            markeredgecolor=colour,
            ecolor=colour,
            elinewidth=0.9,
            capsize=2,
            linestyle="none",
            zorder=3,
        )
        for xi, yi in zip(zero_x, zero_floor):
            axis.annotate(
                "",
                xy=(xi, yi * 0.45),
                xytext=(xi, yi),
                arrowprops=dict(arrowstyle="->", color=colour, linewidth=0.8, shrinkA=0, shrinkB=0),
                zorder=3,
            )
        # The arrow points below every drawn value, so the axis has to be told to make room for it;
        # left to autoscale it would clip the one mark that says "this is an upper bound, not a
        # measurement" and the zero would read as a datum at the floor.
        floors = getattr(axis, "_zero_floors", [])
        floors.extend(float(value) for value in zero_floor if np.isfinite(value) and value > 0)
        axis._zero_floors = floors


def apply_floor_headroom(axis):
    """Leaves room under the lowest floor mark for its arrow. Call once, after every series."""
    floors = getattr(axis, "_zero_floors", None)
    if floors:
        axis.set_ylim(bottom=min(floors) * 0.30)


def shots_subtitle(frame) -> str:
    counts = sorted(set(int(value) for value in frame["shots"]))
    if len(counts) == 1:
        return f"shots = {counts[0]:,} per (d, p, T)"
    return f"shots = {counts[0]:,}-{counts[-1]:,} per (d, p, T)"


# ---------------------------------------------------------------------------------------------
# §4.1 Escalation
# ---------------------------------------------------------------------------------------------


def figure_q_vs_p(artifact: Artifact, style: Style):
    """§4.1.1 -- the headline. `q` against `p`, one panel per `T`, one line per `d`."""
    horizons = artifact.horizons()
    figure, axes = panel_grid(len(horizons))
    for axis, horizon in zip(axes, horizons):
        panel = artifact.summary[artifact.summary["T"] == horizon]
        for d in sorted(panel["d"].unique()):
            series = panel[panel["d"] == d].sort_values("p")
            plot_rate(
                axis,
                series["p"],
                series["q"],
                series["q_ci_lo"],
                series["q_ci_hi"],
                1.0 / series["shots"].astype(float),
                style.colour(d),
                style.marker(horizon),
                f"d={fmt_d(d)}",
            )
        axis.set_yscale("log")
        log_x_rates(axis, sorted(panel["p"].unique()))
        axis.set_xlabel("p")
        axis.set_ylabel("q = shots_escalated / shots")
        axis.set_title(f"T = {fmt_t(horizon)}", fontsize=10)
        axis.grid(True, which="major", alpha=0.25, linewidth=0.5)
        axis.legend(fontsize=7, loc="best")
        apply_floor_headroom(axis)
    return finish(
        figure,
        "Escalation rate q vs p",
        shots_subtitle(artifact.summary) + "; error bars are 95% Wilson; " + FLOOR_SENTENCE,
        footer_text(artifact),
    )


def figure_q_vs_d(artifact: Artifact, style: Style):
    """§4.1.2 -- `q` against `d`, one panel per `p`, one line per `T`."""
    rates = artifact.error_rates()
    figure, axes = panel_grid(len(rates))
    for axis, rate in zip(axes, rates):
        panel = artifact.summary[artifact.summary["p"] == rate]
        for horizon in sorted(panel["T"].unique()):
            series = panel[panel["T"] == horizon].sort_values("d")
            plot_rate(
                axis,
                series["d"],
                series["q"],
                series["q_ci_lo"],
                series["q_ci_hi"],
                1.0 / series["shots"].astype(float),
                # `d` is the x axis here, so it cannot also be the colour; `T` takes both channels.
                style.horizon_colour(horizon),
                style.marker(horizon),
                f"T={fmt_t(horizon)}",
            )
        axis.set_yscale("log")
        axis.set_xticks(sorted(panel["d"].unique()))
        axis.set_xlabel("d")
        axis.set_ylabel("q = shots_escalated / shots")
        axis.set_title(f"p = {fmt_p(rate)}", fontsize=10)
        axis.grid(True, which="major", alpha=0.25, linewidth=0.5)
        axis.legend(fontsize=7, loc="best")
        apply_floor_headroom(axis)
    return finish(
        figure,
        "Escalation rate q vs d",
        shots_subtitle(artifact.summary) + "; error bars are 95% Wilson; " + FLOOR_SENTENCE,
        footer_text(artifact),
    )


def figure_q_vs_t(artifact: Artifact, style: Style):
    """§4.1.3 -- `q` against `T`, one panel per `d`, one line per `p`."""
    distances = artifact.distances()
    figure, axes = panel_grid(len(distances))
    for axis, d in zip(axes, distances):
        panel = artifact.summary[artifact.summary["d"] == d]
        for rate in sorted(panel["p"].unique()):
            series = panel[panel["p"] == rate].sort_values("T")
            plot_rate(
                axis,
                series["T"],
                series["q"],
                series["q_ci_lo"],
                series["q_ci_hi"],
                1.0 / series["shots"].astype(float),
                style.rate_colour(rate),
                "o",
                f"p={fmt_p(rate)}",
            )
        axis.set_yscale("log")
        axis.set_xticks(sorted(panel["T"].unique()))
        axis.set_xlabel("T (lattice edge weights)")
        axis.set_ylabel("q = shots_escalated / shots")
        axis.set_title(f"d = {fmt_d(d)}", fontsize=10)
        axis.grid(True, which="major", alpha=0.25, linewidth=0.5)
        axis.legend(fontsize=7, loc="best")
        apply_floor_headroom(axis)
    return finish(
        figure,
        "Escalation rate q vs T",
        shots_subtitle(artifact.summary) + "; error bars are 95% Wilson; " + FLOOR_SENTENCE,
        footer_text(artifact),
    )


def figure_escalation_composition(artifact: Artifact, style: Style, horizon: float):
    """§4.1.4 -- what an escalation is made of, at one `T`, over the `(d, p)` grid.

    Left axis: the per-shot rate `q` beside the per-component rate
    `components_truncated / components_total`. Right axis: how many components truncate on a shot
    that escalates at all, which is the number that says whether escalation is one bad component or
    many.
    """
    panel = artifact.summary[artifact.summary["T"] == horizon].sort_values(["d", "p"])
    figure, axes = panel_grid(1, width=max(6.0, 0.9 * len(panel) + 3.0), height=4.2)
    axis = axes[0]
    positions = np.arange(len(panel), dtype=float)
    width = 0.38

    q = panel["q"].to_numpy(dtype=float)
    q_floor = 1.0 / panel["shots"].to_numpy(dtype=float)
    components_total = panel["components_total"].to_numpy(dtype=float)
    truncated = panel["components_truncated"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        component_rate = np.where(components_total > 0, truncated / components_total, 0.0)
        component_floor = np.where(components_total > 0, 1.0 / np.where(components_total > 0, components_total, 1.0), np.nan)

    for offset, values, floor, colour, label in (
        (-width / 2, q, q_floor, "#1f77b4", "q (per shot)"),
        (width / 2, component_rate, component_floor, "#ff7f0e", "components_truncated / components_total"),
    ):
        axis.bar(
            positions + offset,
            np.where(values > 0, values, 0.0),
            width=width,
            color=colour,
            label=label,
            zorder=2,
        )
        # A zero gets a capped marker at the resolution floor, not a bar reaching it. The two floors
        # here are 1/shots and 1/components_total, and the second is much the coarser -- a bar drawn
        # up to it would be the tallest thing on the page and would read as the largest rate
        # measured, which is the opposite of what "no component truncated" says.
        zero = values <= 0
        if zero.any():
            axis.errorbar(
                positions[zero] + offset,
                floor[zero],
                fmt="_",
                markersize=11,
                markeredgewidth=1.4,
                color=colour,
                linestyle="none",
                zorder=4,
            )
            for xi, yi in zip(positions[zero] + offset, floor[zero]):
                if np.isfinite(yi):
                    axis.annotate(
                        "",
                        xy=(xi, yi * 0.30),
                        xytext=(xi, yi),
                        arrowprops=dict(arrowstyle="->", color=colour, linewidth=1.0, shrinkA=0, shrinkB=0),
                        zorder=4,
                    )
            existing = getattr(axis, "_zero_floors", [])
            existing.extend(float(value) for value in floor[zero] if np.isfinite(value) and value > 0)
            axis._zero_floors = existing
    axis.set_yscale("log")
    axis.set_ylabel("rate")
    apply_floor_headroom(axis)
    # Headroom for the legend, which would otherwise sit on top of the tallest bar.
    bottom, top = axis.get_ylim()
    axis.set_ylim(bottom, top * 3.0)
    axis.set_xticks(positions)
    axis.set_xticklabels(
        [f"d={fmt_d(d)}\np={fmt_p(p)}" for d, p in zip(panel["d"], panel["p"])], fontsize=7
    )
    axis.grid(True, axis="y", which="major", alpha=0.25, linewidth=0.5)

    right = axis.twinx()
    per_escalated = panel["truncated_components_per_escalated_shot"].to_numpy(dtype=float)
    escalated = panel["shots_escalated"].to_numpy(dtype=float)
    shown = escalated > 0
    right.plot(
        positions[shown],
        per_escalated[shown],
        "k^",
        markersize=5,
        linestyle="none",
        label="truncated components per escalated shot",
        zorder=4,
    )
    right.set_ylabel("truncated components per escalated shot")
    if shown.any():
        right.set_ylim(0, max(2.0, float(np.nanmax(per_escalated[shown])) * 1.25))
    else:
        right.set_ylim(0, 2.0)
        right.text(
            0.5,
            0.5,
            "no shot escalated: the right axis has nothing to plot",
            transform=right.transAxes,
            ha="center",
            va="center",
            fontsize=7,
            color="0.45",
        )
    handles, labels = axis.get_legend_handles_labels()
    right_handles, right_labels = right.get_legend_handles_labels()
    handles.append(Line2D([], [], color="0.4", marker="_", markersize=9, linestyle="none"))
    labels.append("0 events: floor mark at 1/shots or 1/components_total")
    axis.legend(handles + right_handles, labels + right_labels, fontsize=7, loc="upper left", framealpha=0.92)
    return finish(
        figure,
        f"Escalation composition, T = {fmt_t(horizon)}",
        shots_subtitle(panel) + "; " + FLOOR_SENTENCE,
        footer_text(artifact),
    )


# ---------------------------------------------------------------------------------------------
# §4.2 Structure scalars
# ---------------------------------------------------------------------------------------------


def figure_structure_scalar(artifact: Artifact, style: Style, name: str, mean_column: str, max_column: str):
    """§4.2 -- one scalar's mean (line) and max (faint marker) vs `p`, one panel per `T`."""
    horizons = artifact.horizons()
    figure, axes = panel_grid(len(horizons))
    for axis, horizon in zip(axes, horizons):
        panel = artifact.summary[artifact.summary["T"] == horizon]
        for d in sorted(panel["d"].unique()):
            series = panel[panel["d"] == d].sort_values("p")
            colour = style.colour(d)
            axis.plot(
                series["p"],
                series[mean_column],
                marker=style.marker(horizon),
                markersize=4.5,
                color=colour,
                linewidth=1.2,
                label=f"d={fmt_d(d)} mean",
            )
            if max_column in series.columns:
                axis.plot(
                    series["p"],
                    series[max_column],
                    marker=style.marker(horizon),
                    markersize=4.0,
                    markerfacecolor="none",
                    color=colour,
                    linewidth=0.0,
                    alpha=0.45,
                    label=f"d={fmt_d(d)} max",
                )
        log_x_rates(axis, sorted(panel["p"].unique()))
        axis.set_xlabel("p")
        axis.set_ylabel(name)
        axis.set_title(f"T = {fmt_t(horizon)}", fontsize=10)
        axis.grid(True, which="major", alpha=0.25, linewidth=0.5)
        axis.legend(fontsize=6, loc="best", ncol=2)
    return finish(
        figure,
        f"{name}: mean and max vs p",
        "means over shots with at least one defect; max over the same shots",
        footer_text(artifact),
    )


def figure_component_size_mix(artifact: Artifact, style: Style, horizon: float):
    """§4.2 -- the singleton / pair / size>=3 mix as a fraction of `num_components`.

    A stacked area needs a panel to itself, so `d` moves from a line to a panel here; `T` stays the
    figure, as it is everywhere else in §4.2.
    """
    panel = artifact.summary[artifact.summary["T"] == horizon]
    distances = sorted(panel["d"].unique())
    figure, axes = panel_grid(len(distances))
    for axis, d in zip(axes, distances):
        series = panel[panel["d"] == d].sort_values("p")
        total = series["mean_num_components"].to_numpy(dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            fractions = [
                np.where(total > 0, series[column].to_numpy(dtype=float) / np.where(total > 0, total, 1.0), 0.0)
                for column in (
                    "mean_singleton_components",
                    "mean_pair_components",
                    "mean_components_size_ge3",
                )
            ]
        axis.stackplot(
            series["p"].to_numpy(dtype=float),
            *fractions,
            labels=["singleton", "pair", "size >= 3"],
            colors=["#4c78a8", "#f58518", "#54a24b"],
            alpha=0.85,
        )
        log_x_rates(axis, sorted(series["p"].unique()))
        axis.set_ylim(0, 1)
        axis.set_xlabel("p")
        axis.set_ylabel("mean count / mean num_components")
        axis.set_title(f"d = {fmt_d(d)}", fontsize=10)
        axis.legend(fontsize=7, loc="lower left")
    return finish(
        figure,
        f"Component size mix, T = {fmt_t(horizon)}",
        shots_subtitle(panel),
        footer_text(artifact),
    )


def figure_odd_component_lower_bound(artifact: Artifact, style: Style):
    """§4.2 -- the odd-without-boundary rate on `q`'s own axes.

    An odd component with no boundary edge cannot match within itself and has nowhere else to go, so
    it truncates with certainty: the fraction of shots carrying one is a lower bound on `q`. Drawing
    the two together is the whole point -- the gap between them is the escalation that structure
    alone does not explain.
    """
    horizons = artifact.horizons()
    figure, axes = panel_grid(len(horizons))
    used_count_fallback = False
    for axis, horizon in zip(axes, horizons):
        panel = artifact.summary[artifact.summary["T"] == horizon]
        for d in sorted(panel["d"].unique()):
            series = panel[panel["d"] == d].sort_values("p")
            colour = style.colour(d)
            rates, denominators, is_rate = [], [], True
            for _, row in series.iterrows():
                value, denominator, kind = odd_component_measure(artifact, row)
                rates.append(value)
                denominators.append(denominator)
                is_rate = is_rate and kind == "rate"
            rates = np.asarray(rates, dtype=float)
            denominators = np.asarray(denominators, dtype=float)
            if is_rate:
                lo, hi = wilson(np.round(rates * denominators), denominators)
                plot_rate(
                    axis,
                    series["p"],
                    rates,
                    lo,
                    hi,
                    1.0 / np.where(denominators > 0, denominators, 1.0),
                    colour,
                    style.marker(horizon),
                    f"d={fmt_d(d)} odd-no-boundary rate",
                )
            else:
                used_count_fallback = True
                axis.plot(
                    series["p"],
                    rates,
                    marker=style.marker(horizon),
                    color=colour,
                    linewidth=1.2,
                    label=f"d={fmt_d(d)} mean count",
                )
            # `q` beside it, in the same colour and dashed, so the bound and the thing it bounds are
            # read off one axis.
            plot_rate(
                axis,
                series["p"],
                series["q"],
                series["q_ci_lo"],
                series["q_ci_hi"],
                1.0 / series["shots"].astype(float),
                colour,
                style.marker(horizon),
                f"d={fmt_d(d)} q",
                linestyle="--",
            )
        axis.set_yscale("log")
        log_x_rates(axis, sorted(panel["p"].unique()))
        axis.set_xlabel("p")
        axis.set_ylabel("rate" if not used_count_fallback else "rate / count")
        axis.set_title(f"T = {fmt_t(horizon)}", fontsize=10)
        axis.grid(True, which="major", alpha=0.25, linewidth=0.5)
        axis.legend(fontsize=6, loc="best", ncol=2)
        if used_count_fallback:
            axis.text(
                0.02,
                0.98,
                "no per-shot vector and no rate column:\nthe solid series is a mean COUNT, not a rate",
                transform=axis.transAxes,
                ha="left",
                va="top",
                fontsize=6,
                color="0.35",
            )
        apply_floor_headroom(axis)
    return finish(
        figure,
        "Shots with an odd component that has no boundary edge, and q",
        "solid: fraction of defect-carrying shots with such a component (dashed: q, same colour); "
        + FLOOR_SENTENCE,
        footer_text(artifact),
    )


def odd_component_measure(artifact: Artifact, row):
    """`(value, denominator, "rate"|"count")` for `num_odd_components_without_boundary`.

    §4.2 asks for `mean(count > 0)` when the per-shot vector is in `hists.json`, and says to fall
    back to the mean count otherwise. The binary also writes the same rate as a column, so that sits
    between the two: it is the identical quantity over the identical denominator, read rather than
    recomputed. Only when neither exists does the count get drawn, with the note §4.2 asks for.
    """
    cell = artifact.cell(row["d"], row["p"], row["T"])
    per_shot = (cell or {}).get("per_shot")
    if per_shot and "num_odd_components_without_boundary" in per_shot and "n_defects" in per_shot:
        odd = np.asarray(per_shot["num_odd_components_without_boundary"], dtype=float)
        defects = np.asarray(per_shot["n_defects"], dtype=float)
        if odd.size == defects.size and odd.size:
            # Restricted to shots with defects, which is the denominator the CSV's own means use.
            non_empty = defects > 0
            denominator = float(non_empty.sum())
            if denominator > 0:
                return float((odd[non_empty] > 0).sum()) / denominator, denominator, "rate"
    if "odd_component_without_boundary_rate" in row.index and "shots_with_defects" in row.index:
        return (
            float(row["odd_component_without_boundary_rate"]),
            float(row["shots_with_defects"]),
            "rate",
        )
    return float(row.get("mean_odd_components_without_boundary", 0.0)), float(row.get("shots", 0)), "count"


# ---------------------------------------------------------------------------------------------
# §4.3 Histograms
# ---------------------------------------------------------------------------------------------


def hist_bin_edges(kind: str, values, cell):
    """`(x centres, bar width, x label)` for one histogram.

    Every kind but `wint_diameter_hist` is binned by count: bin `k` counts the value `k`, and the
    last bin is the overflow. `wint_diameter_hist` is binned in sixteenths of `T_int`, so its axis is
    put back into `w_int` units, which is what makes the `T_int` and `2*T_int` marks meaningful.
    """
    n = len(values)
    if kind == "wint_diameter_hist":
        bins_per_t = 16
        t_int = float(cell.get("T_int", 0) or 0)
        if t_int > 0:
            step = t_int / bins_per_t
            centres = (np.arange(n, dtype=float) + 0.5) * step
            return centres, step, "w_int diameter"
        centres = np.arange(n, dtype=float)
        return centres, 0.9, f"w_int diameter (T_int/{bins_per_t} bins; T_int unavailable)"
    centres = np.arange(n, dtype=float)
    return centres, 0.9, HIST_KINDS[kind][0]


def draw_hist_panel(axis, kind: str, values, cell, cap_label: str, annotate_counts: bool = True):
    values = np.asarray(values, dtype=float)
    total = values.sum()
    fraction = values / total if total > 0 else np.zeros_like(values)
    centres, width, xlabel = hist_bin_edges(kind, values, cell)

    colours = ["#4c78a8"] * len(values)
    colours[-1] = "#b279a2"  # the overflow bin, always the last bar
    axis.bar(centres, fraction, width=width, color=colours, zorder=2)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(f"fraction of {HIST_KINDS[kind][2]}")
    axis.set_title(kind, fontsize=9)
    axis.grid(True, axis="y", alpha=0.2, linewidth=0.5)
    if annotate_counts and total > 0:
        top = fraction.max()
        for centre, value, height in zip(centres, values, fraction):
            if value > 0:
                axis.text(
                    centre,
                    height + top * 0.015,
                    f"{int(value)}",
                    ha="center",
                    va="bottom",
                    fontsize=4.2,
                    rotation=90,
                    color="0.3",
                )
        axis.set_ylim(0, top * 1.28 if top > 0 else 1.0)
    # Which bar the overflow is, said in the corner rather than under the bar: at the right-hand
    # edge of the axes the label runs off the figure, and above the bar it lands on that bar's own
    # count.
    axis.text(
        0.985,
        0.985,
        f"last bar (■) is the overflow bin: {cap_label}",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=5.5,
        color="#b279a2",
    )

    if kind == "wint_diameter_hist":
        t_int = float(cell.get("T_int", 0) or 0)
        if t_int > 0:
            for multiple, style in ((1, "-"), (2, "--")):
                axis.axvline(
                    multiple * t_int,
                    color="0.25",
                    linewidth=0.9,
                    linestyle=style,
                    zorder=3,
                )
                axis.text(
                    multiple * t_int,
                    axis.get_ylim()[1] * 0.97,
                    f" {'' if multiple == 1 else '2·'}T_int",
                    ha="left",
                    va="top",
                    fontsize=6,
                    color="0.25",
                )
    uncomputed = int(cell.get("diameter_uncomputed_components", 0) or 0)
    if uncomputed and kind in ("hop_diameter_hist", "wint_diameter_hist"):
        # §4.3: printed in the corner when non-zero. These components are above --diameter-cap and
        # are in no bin of this panel, so the histogram is over fewer components than the cell has.
        axis.text(
            0.985,
            0.925,
            f"diameter_uncomputed_components = {uncomputed}",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=6,
            color="0.35",
        )


def cap_label_for(kind: str, caps: dict, length: int) -> str:
    if kind == "degree_hist":
        cap = caps.get("degree_cap", length - 1)
    elif kind == "wint_diameter_hist":
        return ">= 4·T_int"
    else:
        cap = caps.get("size_cap", length - 1)
    return f">= {cap}"


def figure_cell_hists(artifact: Artifact, d, p, t):
    """§4.3 -- the four histograms of one `(d, p, T)`, normalised, raw counts above the bars."""
    cell = artifact.cell(d, p, t)
    kinds = [kind for kind in HIST_KINDS if kind in cell]
    figure, axes = panel_grid(len(kinds), width=5.2, height=3.4)
    for axis, kind in zip(axes, kinds):
        values = cell[kind]
        draw_hist_panel(axis, kind, values, cell, cap_label_for(kind, artifact.caps, len(values)))
    subtitle = (
        f"components = {int(cell.get('components_total', 0)):,}; "
        f"shots = {int(cell.get('shots', 0)):,}; "
        f"T_int = {int(cell.get('T_int', 0)):,}"
    )
    return finish(
        figure,
        f"Component histograms, d = {fmt_d(d)}, p = {fmt_p(p)}, T = {fmt_t(t)}",
        subtitle,
        footer_text(artifact, d=d, p=p),
    )


def figure_hist_overlay(artifact: Artifact, style: Style, kind: str, horizon: float):
    """§4.3's overlay -- every `(d, p)` at one `T`, normalised, as log-`y` step plots.

    `wint_diameter_hist` is left on its `T_int/16` bin axis here rather than converted to `w_int`:
    `T_int` moves with `(d, p)`, so a `w_int` axis would put bins that mean the same thing at
    different places and the shapes would no longer be comparable, which is the only reason this
    figure exists.
    """
    panel = artifact.summary[artifact.summary["T"] == horizon].sort_values(["d", "p"])
    figure, axes = panel_grid(1, width=8.0, height=4.6)
    axis = axes[0]
    uncomputed_total = 0
    drawn = 0
    for _, row in panel.iterrows():
        cell = artifact.cell(row["d"], row["p"], row["T"])
        if cell is None or kind not in cell:
            continue
        values = np.asarray(cell[kind], dtype=float)
        total = values.sum()
        if total <= 0:
            continue
        drawn += 1
        uncomputed_total += int(cell.get("diameter_uncomputed_components", 0) or 0)
        axis.step(
            np.arange(len(values), dtype=float),
            values / total,
            where="mid",
            color=style.colour(row["d"]),
            linestyle={0: "-", 1: "--", 2: ":", 3: "-."}[style.rates.index(row["p"]) % 4],
            linewidth=1.2,
            label=f"d={fmt_d(row['d'])} p={fmt_p(row['p'])}",
        )
    if drawn == 0:
        plt.close(figure)
        return None
    axis.set_yscale("log")
    if kind == "wint_diameter_hist":
        axis.set_xlabel("w_int diameter, in bins of T_int/16 (T_int at bin 16, 2·T_int at bin 32)")
        for position in (16, 32):
            axis.axvline(position, color="0.25", linewidth=0.9, linestyle="--", zorder=1)
    else:
        axis.set_xlabel(HIST_KINDS[kind][0] + "  (last bin is the overflow bin)")
    axis.set_ylabel(f"fraction of {HIST_KINDS[kind][2]}")
    axis.grid(True, which="major", alpha=0.25, linewidth=0.5)
    axis.legend(fontsize=6, loc="best", ncol=2)
    if uncomputed_total and kind in ("hop_diameter_hist", "wint_diameter_hist"):
        axis.text(
            0.98,
            0.02,
            f"diameter_uncomputed_components over the grid = {uncomputed_total}",
            transform=axis.transAxes,
            ha="right",
            va="bottom",
            fontsize=6,
            color="0.35",
        )
    return finish(
        figure,
        f"{kind} across the (d, p) grid, T = {fmt_t(horizon)}",
        shots_subtitle(panel),
        footer_text(artifact),
    )


# ---------------------------------------------------------------------------------------------
# §4.4 Joint tables
# ---------------------------------------------------------------------------------------------


def figure_status_table(artifact: Artifact, d, p, t, table_name: str, bin_label: str):
    """§4.4 -- one joint table: the two status counts per bin, and `P(TRUNCATED | bin)` under it.

    A bin with no components is left blank in both panels. Drawing it as a zero would say the bin
    was measured and came out zero, which is a different statement from not having been sampled.
    """
    cell = artifact.cell(d, p, t)
    table = np.asarray(cell[table_name], dtype=float)
    complete = table[:, 0]
    truncated = table[:, 1]
    total = complete + truncated
    bins = np.arange(len(total), dtype=float)
    populated = total > 0

    figure, axes = panel_grid(2, width=6.2, height=3.6)
    top, bottom = axes[0], axes[1]

    width = 0.4
    top.bar(bins[populated] - width / 2, complete[populated], width=width, color="#4c78a8", label="COMPLETE")
    top.bar(bins[populated] + width / 2, truncated[populated], width=width, color="#e45756", label="TRUNCATED")
    top.set_yscale("symlog", linthresh=1.0)
    top.set_xlabel(bin_label)
    top.set_ylabel("components")
    top.set_title("counts by status", fontsize=9)
    top.legend(fontsize=7, loc="upper left", framealpha=0.92)
    top.grid(True, axis="y", alpha=0.2, linewidth=0.5)

    with np.errstate(divide="ignore", invalid="ignore"):
        probability = np.where(populated, truncated / np.where(populated, total, 1.0), np.nan)
    lo, hi = wilson(truncated, total)
    lower = np.clip(probability - lo, 0, None)
    upper = np.clip(hi - probability, 0, None)
    bottom.errorbar(
        bins[populated],
        probability[populated],
        yerr=np.vstack([lower[populated], upper[populated]]),
        fmt="o",
        markersize=4,
        color="#e45756",
        elinewidth=0.9,
        capsize=2,
        linestyle="none",
    )
    bottom.set_ylim(-0.03, 1.03)
    bottom.set_xlabel(bin_label)
    bottom.set_ylabel("P(TRUNCATED | bin)")
    bottom.set_title("truncation probability, 95% Wilson on the two counts", fontsize=9)
    bottom.grid(True, alpha=0.2, linewidth=0.5)
    for axis in (top, bottom):
        axis.set_xlim(-0.8, len(total) - 0.2)
        axis.axvline(len(total) - 1.5, color="0.7", linewidth=0.8, linestyle=":", zorder=1)
        axis.text(
            0.985,
            0.985,
            f"right of the dotted line: the overflow bin ({len(total) - 1}+)",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=5.5,
            color="0.45",
        )
    uncomputed = int(cell.get("diameter_uncomputed_components", 0) or 0)
    if uncomputed and table_name == "hop_diameter_x_status":
        top.text(
            0.985,
            0.925,
            f"diameter_uncomputed_components = {uncomputed} (not in this table)",
            transform=top.transAxes,
            ha="right",
            va="top",
            fontsize=6,
            color="0.35",
        )
    return finish(
        figure,
        f"{table_name}, d = {fmt_d(d)}, p = {fmt_p(p)}, T = {fmt_t(t)}",
        f"components = {int(total.sum()):,}; empty bins are left blank, not drawn as 0",
        footer_text(artifact, d=d, p=p),
    )


def figure_truncation_grid(artifact: Artifact, style: Style, horizon: float):
    """§4.4's grid summary -- `P(TRUNCATED | size)` at a few sizes, vs `p`, one line per `d`."""
    panel = artifact.summary[artifact.summary["T"] == horizon]
    size_cap = artifact.caps.get("size_cap")
    figure, axes = panel_grid(len(GRID_SIZE_BINS), width=4.2, height=3.4)
    any_drawn = False
    for axis, size_bin in zip(axes, GRID_SIZE_BINS):
        for d in sorted(panel["d"].unique()):
            series = panel[panel["d"] == d].sort_values("p")
            xs, ys, los, his, floors = [], [], [], [], []
            for _, row in series.iterrows():
                cell = artifact.cell(row["d"], row["p"], row["T"])
                if cell is None or "size_x_status" not in cell:
                    continue
                table = np.asarray(cell["size_x_status"], dtype=float)
                if abs(size_bin) >= len(table):
                    continue
                complete, truncated = table[size_bin]
                total = complete + truncated
                if total <= 0:
                    continue
                lo, hi = wilson(np.array([truncated]), np.array([total]))
                xs.append(float(row["p"]))
                ys.append(truncated / total)
                los.append(float(lo[0]))
                his.append(float(hi[0]))
                floors.append(1.0 / total)
            if not xs:
                continue
            any_drawn = True
            plot_rate(
                axis,
                xs,
                ys,
                los,
                his,
                floors,
                style.colour(d),
                style.marker(horizon),
                f"d={fmt_d(d)}",
            )
        axis.set_yscale("log")
        rates = sorted(panel["p"].unique())
        log_x_rates(axis, rates)
        axis.set_xlabel("p")
        axis.set_ylabel("P(TRUNCATED | size)")
        name = f"size >= {size_cap}" if size_bin == -1 else f"size = {size_bin}"
        axis.set_title(name, fontsize=10)
        axis.grid(True, which="major", alpha=0.25, linewidth=0.5)
        axis.legend(fontsize=6, loc="best")
        apply_floor_headroom(axis)
    if not any_drawn:
        plt.close(figure)
        return None
    return finish(
        figure,
        f"P(TRUNCATED | size) vs p, T = {fmt_t(horizon)}",
        "error bars are 95% Wilson on the two counts; " + FLOOR_SENTENCE + " -- here 1 / (components in the bin)",
        footer_text(artifact),
    )


# ---------------------------------------------------------------------------------------------
# Driving
# ---------------------------------------------------------------------------------------------


class Recorder:
    """Collects what was written and what was not, for the two lists §5 asks the script to print."""

    def __init__(self, out_dir: Path, extension: str, show: bool):
        self.out_dir = out_dir
        self.extension = extension
        self.show = show
        self.written = []
        self.skipped = []

    def save(self, figure, name: str, group: str, quantity: str, fixed: dict):
        if figure is None:
            return
        path = self.out_dir / f"{name}.{self.extension}"
        metadata = {}
        if self.extension == "pdf":
            metadata = {"CreationDate": None}
        elif self.extension == "svg":
            metadata = {"Date": None}
        elif self.extension == "png":
            # Pinned rather than left to matplotlib's version string, so two machines agree.
            metadata = {"Software": "plot_sparse_graph_stats"}
        figure.savefig(path, dpi=200, metadata=metadata)
        if not self.show:
            plt.close(figure)
        self.written.append((path.name, group, quantity, fixed))

    def skip(self, name: str, reason: str):
        self.skipped.append((name, reason))


def required_columns(artifact: Artifact, columns):
    return [column for column in columns if not artifact.has(column)]


def draw_escalation(artifact: Artifact, style: Style, recorder: Recorder):
    missing = required_columns(artifact, ["q", "q_ci_lo", "q_ci_hi", "shots"])
    if missing:
        recorder.skip("escalation (all figures)", f"summary.csv is missing {', '.join(missing)}")
        return

    if len(artifact.error_rates()) < 2:
        recorder.skip("q_vs_p", "fewer than two p values on the intended x axis")
    else:
        recorder.save(figure_q_vs_p(artifact, style), "q_vs_p", "escalation", "q", {"x": "p", "panels": "T"})

    if len(artifact.distances()) < 2:
        recorder.skip("q_vs_d", "fewer than two d values on the intended x axis")
    else:
        recorder.save(figure_q_vs_d(artifact, style), "q_vs_d", "escalation", "q", {"x": "d", "panels": "p"})

    if len(artifact.horizons()) < 2:
        recorder.skip("q_vs_T", "fewer than two T values on the intended x axis")
    else:
        recorder.save(figure_q_vs_t(artifact, style), "q_vs_T", "escalation", "q", {"x": "T", "panels": "d"})

    missing = required_columns(
        artifact,
        ["components_total", "components_truncated", "truncated_components_per_escalated_shot"],
    )
    if missing:
        recorder.skip("escalation_composition", f"summary.csv is missing {', '.join(missing)}")
        return
    for horizon in artifact.horizons():
        recorder.save(
            figure_escalation_composition(artifact, style, horizon),
            f"escalation_composition__{label_of(T=horizon)}",
            "escalation",
            "q and components_truncated / components_total",
            {"T": fmt_t(horizon)},
        )


def draw_structure(artifact: Artifact, style: Style, recorder: Recorder):
    if len(artifact.error_rates()) < 2:
        recorder.skip("structure (all figures)", "fewer than two p values on the intended x axis")
        return
    for name, mean_column, max_column in STRUCTURE_SCALARS:
        if not artifact.has(mean_column):
            recorder.skip(f"structure_{name}", f"summary.csv is missing {mean_column}")
            continue
        recorder.save(
            figure_structure_scalar(artifact, style, name, mean_column, max_column),
            f"structure_{name}",
            "structure",
            f"{name} (mean, max)",
            {"x": "p", "panels": "T"},
        )

    mix_columns = [
        "mean_num_components",
        "mean_singleton_components",
        "mean_pair_components",
        "mean_components_size_ge3",
    ]
    missing = required_columns(artifact, mix_columns)
    if missing:
        recorder.skip("component_size_mix", f"summary.csv is missing {', '.join(missing)}")
    else:
        for horizon in artifact.horizons():
            recorder.save(
                figure_component_size_mix(artifact, style, horizon),
                f"component_size_mix__{label_of(T=horizon)}",
                "structure",
                "singleton / pair / size>=3 as a fraction of num_components",
                {"T": fmt_t(horizon)},
            )

    if not artifact.has("mean_odd_components_without_boundary"):
        recorder.skip(
            "odd_component_lower_bound_vs_p",
            "summary.csv is missing mean_odd_components_without_boundary",
        )
    else:
        recorder.save(
            figure_odd_component_lower_bound(artifact, style),
            "odd_component_lower_bound_vs_p",
            "structure",
            "P(a shot has an odd component with no boundary edge), with q",
            {"x": "p", "panels": "T"},
        )


def draw_hists(artifact: Artifact, style: Style, recorder: Recorder):
    if not artifact.cells:
        recorder.skip("hists (all figures)", "no hists.json found under any --in dir")
        return
    for _, row in artifact.summary.iterrows():
        d, p, t = row["d"], row["p"], row["T"]
        cell = artifact.cell(d, p, t)
        if cell is None:
            recorder.skip(
                f"component_hists__{label_of(d=d, p=p, T=t)}",
                "hists.json has no cell for this (d, p, T)",
            )
            continue
        if not any(kind in cell for kind in HIST_KINDS):
            recorder.skip(
                f"component_hists__{label_of(d=d, p=p, T=t)}",
                "the cell holds none of the four histograms",
            )
            continue
        recorder.save(
            figure_cell_hists(artifact, d, p, t),
            f"component_hists__{label_of(d=d, p=p, T=t)}",
            "hists",
            "component_size / degree / hop_diameter / wint_diameter",
            {"d": fmt_d(d), "p": fmt_p(p), "T": fmt_t(t)},
        )

    for kind in HIST_KINDS:
        for horizon in artifact.horizons():
            name = f"overlay_{kind}__{label_of(T=horizon)}"
            figure = figure_hist_overlay(artifact, style, kind, horizon)
            if figure is None:
                recorder.skip(name, f"no cell at this T holds a non-empty {kind}")
                continue
            recorder.save(figure, name, "hists", f"{kind}, normalised, across (d, p)", {"T": fmt_t(horizon)})


def draw_joint(artifact: Artifact, style: Style, recorder: Recorder):
    if not artifact.cells:
        recorder.skip("joint (all figures)", "no hists.json found under any --in dir")
        return
    tables = [
        ("size_x_status", "component size (last bin overflows)"),
        ("hop_diameter_x_status", "hop diameter (last bin overflows)"),
    ]
    for _, row in artifact.summary.iterrows():
        d, p, t = row["d"], row["p"], row["T"]
        cell = artifact.cell(d, p, t)
        for table_name, bin_label in tables:
            name = f"{table_name}__{label_of(d=d, p=p, T=t)}"
            if cell is None or table_name not in cell:
                recorder.skip(name, f"hists.json has no {table_name} for this (d, p, T)")
                continue
            recorder.save(
                figure_status_table(artifact, d, p, t, table_name, bin_label),
                name,
                "joint",
                f"{table_name} and P(TRUNCATED | bin)",
                {"d": fmt_d(d), "p": fmt_p(p), "T": fmt_t(t)},
            )

    if len(artifact.error_rates()) < 2:
        recorder.skip("p_truncated_given_size", "fewer than two p values on the intended x axis")
        return
    for horizon in artifact.horizons():
        name = f"p_truncated_given_size__{label_of(T=horizon)}"
        figure = figure_truncation_grid(artifact, style, horizon)
        if figure is None:
            recorder.skip(name, "no size bin at this T holds a component")
            continue
        recorder.save(figure, name, "joint", "P(TRUNCATED | size) vs p", {"T": fmt_t(horizon)})


def write_index(recorder: Recorder, artifact: Artifact):
    """`figures_index.md`, byte-identical for the same artifacts: sorted, and with no clock in it."""
    lines = ["# `sparse_graph_stats` figures", ""]
    lines.append("Source dirs:")
    lines.append("")
    for directory in sorted(artifact.source_dirs):
        lines.append(f"- `{directory}`")
    lines.append("")
    grid = (
        f"d = {', '.join(fmt_d(value) for value in artifact.distances())}; "
        f"p = {', '.join(fmt_p(value) for value in artifact.error_rates())}; "
        f"T = {', '.join(fmt_t(value) for value in artifact.horizons())}"
    )
    lines.append(f"Grid drawn: {grid}.")
    lines.append("")
    lines.append("| figure | group | quantity | fixed |")
    lines.append("| --- | --- | --- | --- |")
    for name, group, quantity, fixed in sorted(recorder.written):
        fixed_text = ", ".join(f"{key}={value}" for key, value in sorted(fixed.items())) or "-"
        # `P(TRUNCATED | size)` carries a pipe, which would end the cell it sits in.
        escaped = quantity.replace("|", "\\|")
        lines.append(f"| [`{name}`]({name}) | {group} | {escaped} | {fixed_text} |")
    if recorder.skipped:
        lines.append("")
        lines.append("## Skipped")
        lines.append("")
        lines.append("| figure | reason |")
        lines.append("| --- | --- |")
        for name, reason in sorted(recorder.skipped):
            lines.append(f"| `{name}` | {reason} |")
    lines.append("")
    (recorder.out_dir / "figures_index.md").write_text("\n".join(lines))


def parse_float_list(text):
    if text is None:
        return None
    return [float(token) for token in text.split(",") if token.strip()]


def parse_int_list(text):
    if text is None:
        return None
    return [int(float(token)) for token in text.split(",") if token.strip()]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Draw the summary.csv / hists.json that sparse_graph_stats writes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--in",
        dest="in_dirs",
        nargs="+",
        required=True,
        metavar="DIR",
        help="one or more sparse_graph_stats output dirs; concatenated, must not collide on (d,p,T)",
    )
    parser.add_argument("--out", default=None, metavar="DIR", help="default: first --in dir + /figures")
    parser.add_argument("--format", default="png", choices=["png", "pdf", "svg"])
    parser.add_argument("--d", default=None, help="comma list; default = every d present")
    parser.add_argument("--p", default=None, help="comma list; default = every p present")
    parser.add_argument("--T", dest="T", default=None, help="comma list; default = every T present")
    parser.add_argument(
        "--which",
        default="all",
        help="comma list of all|escalation|structure|hists|joint",
    )
    parser.add_argument("--show", action="store_true", help="also open interactive windows")
    arguments = parser.parse_args(argv)

    groups = {token.strip() for token in arguments.which.split(",") if token.strip()}
    known = {"all", "escalation", "structure", "hists", "joint"}
    unknown = groups - known
    if unknown:
        parser.error(f"unrecognised --which {','.join(sorted(unknown))} (want {'|'.join(sorted(known))})")
    if "all" in groups:
        groups = known - {"all"}

    artifact = load_artifacts(arguments.in_dirs)
    artifact = filter_artifact(
        artifact,
        parse_int_list(arguments.d),
        parse_float_list(arguments.p),
        parse_float_list(arguments.T),
    )
    if artifact.summary.empty:
        raise SystemExit("error: the --d/--p/--T filters left no rows")

    out_dir = Path(arguments.out) if arguments.out else Path(arguments.in_dirs[0]) / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    style = Style(artifact)
    recorder = Recorder(out_dir, arguments.format, arguments.show)

    if "escalation" in groups:
        draw_escalation(artifact, style, recorder)
    if "structure" in groups:
        draw_structure(artifact, style, recorder)
    if "hists" in groups:
        draw_hists(artifact, style, recorder)
    if "joint" in groups:
        draw_joint(artifact, style, recorder)

    write_index(recorder, artifact)

    print(f"wrote {len(recorder.written)} figures to {out_dir}")
    for name, _, _, _ in sorted(recorder.written):
        print(f"  {name}")
    print("  figures_index.md")
    if recorder.skipped:
        print(f"skipped {len(recorder.skipped)}")
        for name, reason in sorted(recorder.skipped):
            print(f"  {name}: {reason}")

    if arguments.show:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
