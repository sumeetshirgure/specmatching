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

"""Draws component size vs run time from what `component_profiler` writes, and fits `ns = A*size^t`.

    python benchmarks/spec_matching/plot_component_times.py --in ../component_stats
        [--out DIR] [--format png|pdf|svg] [--x-scale log|linear] [--y-scale log|linear]
        [--d 17,21] [--p 1e-3] [--T 1.5] [--status pooled|complete|truncated|split] [--split-status]
        [--stat mean|median|geomean] [--fit-on auto|rows|mean|geomean|median]
        [--split T,status] [--series-by d,p] [--fit-lines per-p,series] [--show]

An input is a directory holding `components.csv` and (optionally) `run.log`, or a `components.csv`
path directly. Several inputs are concatenated; a `(d, p, T)` cell present in more than one of them
is reported and kept, since two campaigns of the same cell on the same machine are more samples of
the same quantity and dropping one silently would be a choice made without the reader.

`design/branch_component_profiler.md` §5 defines the input: one row per connected component of `H`,
`d,p,T,shot,comp_size,status,ns`, where `ns` is the *raw uncorrected* timer delta. §2 puts the
overhead subtraction here rather than in the binary, so this script does it: `run.log`'s
`timer_overhead_ns` is parsed and subtracted from every row, and the caption says by how much.
`--no-overhead-correction` turns that off; `--timer-overhead-ns` overrides the parsed value for a
log that predates the line.

Two things about the record that the figures state rather than leave to be inferred:

  * `run.log`'s `min_logged_comp_size` (3 in the campaigns written so far) means the smallest
    components are *absent from the table*, not absent from the run -- their solve time is of the
    order of the timer overhead, so they were deliberately not logged. The x axis therefore starts
    where the record starts, and the caption says so;
  * `COMPLETE` and `TRUNCATED` are different quantities at the same size (§1): a truncated
    component ran to the horizon and stopped, a completed one ran until it was done. They are
    nonetheless pooled into one figure by default, since the cost of deciding a component is what
    the run pays whichever way that component ended; the caption then says how many of each went
    in. `--split-status` (equivalently `--status split`) draws them as separate figures instead,
    and `--status complete` / `--status truncated` keeps one of them alone.

Every figure carries one fitted curve per physical error rate and no other, each dotted and labelled
with its own `A` and `t` in the legend; a `p` is told apart by its dash pattern, since colour is
spent on the series. No curve is drawn over the whole panel: a single line through several error
rates at once is an average of things the figure is drawing separately, and it read as a thirteenth
series. The pooled fit is still computed, and still appears in the caption, the fit table and the
report, where it is labelled as the pooled scope rather than left to be inferred from a line.
`--fit-lines` chooses which of the two drawn families -- `per-p`, `series` -- appear.

The fit is ordinary least squares on `log10(ns)` against `log10(comp_size)`, which is what
`ns = A * comp_size^t` *is* -- there is no non-linear solver here and no scipy dependency. What is
fitted is by default the same per-size statistic that is drawn, so the line on the page is the line
the numbers describe. Reported with `t`: a heteroskedasticity-robust standard error (the residual
spread is strongly size-dependent, so the textbook SE would be wrong), the same interval clustered
by component size when the fit is over raw rows, R^2, the residual scatter in dex, a split-half
exponent comparison, and a runs test on the residual signs. The last two are the ones that matter:
a power law that is not one shows up as an exponent that drifts between the small and large halves
of the range and as residuals that change sign in long systematic sweeps, and neither of those
disturbs R^2 much.

This script computes nothing the binary could have computed, imports nothing from this repo, and
draws no conclusion in the figure text beyond the arithmetic it did.
"""

from __future__ import annotations

import argparse
import math
import re
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
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# The same artifacts must produce identical figures: `svg.hashsalt` fixes the element ids an SVG
# would otherwise seed from the process, and `FIG_METADATA` drops the creation date out of PDF and
# SVG. Nothing in this file reads a clock.
matplotlib.rcParams["svg.hashsalt"] = "component_times"
matplotlib.rcParams["figure.max_open_warning"] = 0

FIG_METADATA = {
    "png": {"Software": "plot_component_times.py"},
    "pdf": {"Creator": "plot_component_times.py", "Producer": "matplotlib", "CreationDate": None},
    "svg": {"Creator": "plot_component_times.py", "Date": None},
}

#: `d` and `p` are ordered magnitudes, not identities, so the series that carries colour gets a
#: single-hue ordinal ramp rather than a categorical wheel. Both ramps below are the blue ramp
#: sampled inside the band whose light end still clears the surface; the sampling in `ramp()` keeps
#: the lightness steps far enough apart to be told apart, and one hue means a colour-blind reader
#: loses nothing, since lightness alone carries the order.
THEMES = {
    "light": {
        "surface": "#fcfcfb",
        "text_primary": "#0b0b0b",
        "text_secondary": "#52514e",
        "text_muted": "#7a7975",
        "grid": "#dcdcd8",
        "axis": "#b4b3ae",
        "ramp": ("#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
                 "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"),
        "band": "#9ec5f4",
        "fit": "#0b0b0b",
        "accent": "#eb6834",
    },
    "dark": {
        "surface": "#1a1a19",
        "text_primary": "#ffffff",
        "text_secondary": "#c3c2b7",
        "text_muted": "#9a9a90",
        "grid": "#3a3a37",
        "axis": "#55554f",
        "ramp": ("#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec",
                 "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab"),
        "band": "#256abf",
        "fit": "#ffffff",
        "accent": "#d95926",
    },
}

#: The second series key gets marker and dash instead of a second hue, so that colour stays
#: one-dimensional and the legend can be read as two independent lists.
MARKERS = ("o", "s", "^", "D", "v", "P", "X", "*")

#: One dash pattern per physical error rate, for the fitted curve of each `p` (pooled over every
#: other key). They are all dotted and all drawn in one accent colour for the same reason as above:
#: `p` is told apart by the dash, not by a hue that would collide with the `d` ramp. Every pattern
#: is built out of the same round dot and differs only in its rhythm -- the gap between dots and how
#: many dots go in a group -- so no entry in the list stops reading as a dotted line.
DASHES = ((1, 2.2), (1, 5.5), (1, 2.2, 1, 5.5), (1, 2.2, 1, 2.2, 1, 5.5),
          (1, 9.0), (1, 2.2, 1, 9.0), (1, 1.3), (1, 1.3, 1, 1.3, 1, 6.5))

#: The caption block: where its last line sits in figure coordinates and how big it is set. Both are
#: read back when the bottom margin is worked out, so the two cannot drift apart.
CAPTION_BOTTOM = 0.012
CAPTION_FONT_SIZE = 7.6

#: Every column `components.csv` is required to have (§5 of the design). `shot` is read only to be
#: counted; nothing here is per-shot.
REQUIRED_COLUMNS = ("d", "p", "T", "comp_size", "status", "ns")

#: The keys a figure may be split or seriesed on, in the order they are printed in a title.
KEYS = ("d", "p", "T", "status")

#: 1.959964 is the two-sided 95% normal quantile; `student_q` widens it for a small point count.
Z95 = 1.959963984540054


# ---------------------------------------------------------------------------------------------
# run.log


@dataclass
class RunMeta:
    """The provenance lines of `run.log`. Every field is optional: a `components.csv` handed over
    without its log still draws, with the caption saying what it could not state."""

    path: Path | None = None
    timer_overhead_ns: float | None = None
    min_logged_comp_size: int | None = None
    timer_name: str | None = None
    timer_is_thread_scoped: str | None = None
    git_hash: str | None = None
    build_flags: str | None = None
    platform: str | None = None
    seed: str | None = None
    shots_per_cell: str | None = None
    warmup_shots_per_cell: str | None = None
    threads: str | None = None
    rows_total: str | None = None
    status: str | None = None

    def describe(self) -> str:
        bits = []
        if self.git_hash:
            bits.append(f"git {self.git_hash[:12]}")
        if self.build_flags:
            bits.append(self.build_flags)
        if self.platform:
            bits.append(self.platform)
        if self.timer_name:
            scope = " (thread-scoped)" if self.timer_is_thread_scoped == "1" else ""
            bits.append(f"timer {self.timer_name}{scope}")
        if self.threads:
            bits.append(f"threads={self.threads}")
        if self.seed:
            bits.append(f"seed={self.seed}")
        if self.shots_per_cell:
            warm = f"+{self.warmup_shots_per_cell} warm-up" if self.warmup_shots_per_cell else ""
            bits.append(f"{self.shots_per_cell} shots/cell{warm}")
        return "; ".join(bits)


#: `key=value` where the value may carry a trailing parenthesised gloss the binary writes for a
#: human, e.g. `timer_overhead_ns=9.25228 (10^6 back-to-back ...)`. The gloss is dropped.
_KV = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)=([^(\n]*)")


def read_run_log(path: Path) -> RunMeta:
    meta = RunMeta(path=path)
    if not path.is_file():
        return RunMeta()
    fields = {f: None for f in RunMeta.__dataclass_fields__ if f != "path"}
    for line in path.read_text(errors="replace").splitlines():
        match = _KV.match(line)
        if match is None:
            continue
        key, value = match.group(1), match.group(2).strip()
        if key in fields and fields[key] is None:
            fields[key] = value
    for key, value in fields.items():
        if value is None:
            continue
        if key == "timer_overhead_ns":
            try:
                meta.timer_overhead_ns = float(value)
            except ValueError:
                pass
        elif key == "min_logged_comp_size":
            try:
                meta.min_logged_comp_size = int(value)
            except ValueError:
                pass
        else:
            setattr(meta, key, value)
    return meta


# ---------------------------------------------------------------------------------------------
# loading


def resolve_inputs(paths: list[str]) -> list[tuple[Path, RunMeta]]:
    """Each argument is a results directory or a `components.csv`; returns `(csv, run.log meta)`."""
    resolved: list[tuple[Path, RunMeta]] = []
    for raw in paths:
        path = Path(raw).expanduser()
        if path.is_dir():
            csv = path / "components.csv"
            if not csv.is_file():
                raise FileNotFoundError(f"{path} holds no components.csv")
            resolved.append((csv, read_run_log(path / "run.log")))
        elif path.is_file():
            resolved.append((path, read_run_log(path.parent / "run.log")))
        else:
            raise FileNotFoundError(f"no such file or directory: {path}")
    return resolved


def load(inputs: list[tuple[Path, RunMeta]], notes: list[str]) -> tuple[pd.DataFrame, list[RunMeta]]:
    frames, metas = [], []
    seen_cells: dict[tuple, Path] = {}
    for csv, meta in inputs:
        frame = pd.read_csv(
            csv,
            usecols=lambda name: name in REQUIRED_COLUMNS,
            dtype={"d": "int32", "p": "float64", "T": "float64",
                   "comp_size": "int32", "status": "category", "ns": "float64"},
        )
        missing = [c for c in REQUIRED_COLUMNS if c not in frame.columns]
        if missing:
            raise ValueError(f"{csv}: missing column(s) {', '.join(missing)}")
        if frame.empty:
            notes.append(f"{csv} holds no rows and was skipped.")
            continue
        for cell in frame.groupby(["d", "p", "T"], observed=True).groups:
            if cell in seen_cells:
                notes.append(
                    f"cell d={cell[0]} p={cell[1]:g} T={cell[2]:g} appears in both "
                    f"{seen_cells[cell]} and {csv}; both were kept and pooled."
                )
            else:
                seen_cells[cell] = csv
        frame["source"] = str(csv)
        frames.append(frame)
        metas.append(meta)
    if not frames:
        raise ValueError("every input was empty")
    return pd.concat(frames, ignore_index=True), metas


def correct_overhead(frame: pd.DataFrame, overhead: float, notes: list[str]) -> pd.DataFrame:
    """Subtracts the calibrated timer overhead, per §2 of the design (rows are stored uncorrected).

    A row that goes non-positive is a component whose whole solve was inside the noise of one timer
    read; it cannot be placed on a log axis and is dropped, loudly."""
    if overhead == 0.0:
        return frame
    frame = frame.copy()
    frame["ns"] = frame["ns"] - overhead
    bad = int((frame["ns"] <= 0.0).sum())
    if bad:
        notes.append(
            f"{bad} row(s) ({100.0 * bad / len(frame):.3g}% of the selection) fell to <= 0 ns once "
            f"the {overhead:g} ns timer overhead was subtracted and were dropped: their whole solve "
            f"was inside one timer read."
        )
        frame = frame[frame["ns"] > 0.0]
    return frame


def parse_list(text: str | None, cast):
    if text is None:
        return None
    out = []
    for item in text.replace(";", ",").split(","):
        item = item.strip()
        if not item:
            continue
        try:
            out.append(cast(item))
        except ValueError as error:
            raise ValueError(f"cannot read {item!r} as {cast.__name__}") from error
    return out or None


def apply_filters(frame: pd.DataFrame, args, notes: list[str]) -> pd.DataFrame:
    for column, wanted in (("d", args.d), ("p", args.p), ("T", args.T)):
        if wanted is None:
            continue
        present = set(np.unique(frame[column]))
        if column == "d":
            missing = [w for w in wanted if w not in present]
        else:  # float compare against what is actually in the file
            missing = [w for w in wanted if not any(math.isclose(w, v, rel_tol=1e-9) for v in present)]
        if missing:
            notes.append(f"--{column} asked for {missing}, which the input does not hold.")
        keep = np.zeros(len(frame), dtype=bool)
        for want in wanted:
            keep |= np.isclose(frame[column].to_numpy(dtype=float), float(want), rtol=1e-9)
        frame = frame[keep]
    if args.status == "complete":
        frame = frame[frame["status"] == "COMPLETE"]
    elif args.status == "truncated":
        frame = frame[frame["status"] == "TRUNCATED"]
    if args.min_size is not None:
        frame = frame[frame["comp_size"] >= args.min_size]
    if args.max_size is not None:
        frame = frame[frame["comp_size"] <= args.max_size]
    return frame


# ---------------------------------------------------------------------------------------------
# aggregation


def aggregate(frame: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    """One row per `(*keys, comp_size)`: the count, the three centres, the quartiles, and the mean
    and variance of `log10(ns)`.

    The last two are what makes the row-level fit cheap: least squares on every raw row in log space
    is *identical* to least squares on these per-size means weighted by the counts, and the row-level
    sums of squares and the cluster-robust covariance both fall out of `(count, mean, var)`. A
    three-million-row regression is never formed."""
    logns = np.log10(frame["ns"].to_numpy())
    work = pd.DataFrame({
        "ns": frame["ns"].to_numpy(),
        "logns": logns,
        "comp_size": frame["comp_size"].to_numpy(),
    })
    for key in keys:
        work[key] = frame[key].to_numpy()
    grouped = work.groupby(keys + ["comp_size"], observed=True, sort=True)
    out = grouped.agg(
        count=("ns", "size"),
        mean=("ns", "mean"),
        median=("ns", "median"),
        q25=("ns", lambda s: s.quantile(0.25)),
        q75=("ns", lambda s: s.quantile(0.75)),
        logmean=("logns", "mean"),
        logvar=("logns", lambda s: s.var(ddof=0)),
    ).reset_index()
    out["geomean"] = np.power(10.0, out["logmean"].to_numpy())
    return out


def pool(agg: pd.DataFrame) -> pd.DataFrame:
    """Collapses an aggregate over its series keys onto `comp_size` alone.

    `count`, `mean` and the log moments combine exactly (the variance by the parallel-axis identity);
    the quartiles and the median do not, so they are dropped rather than approximated -- a caller
    that needs a pooled median re-aggregates the raw rows."""
    weights = agg["count"].to_numpy(dtype=float)
    work = agg.assign(
        _sum=agg["mean"].to_numpy() * weights,
        _logsum=agg["logmean"].to_numpy() * weights,
        _logsq=(agg["logvar"].to_numpy() + agg["logmean"].to_numpy() ** 2) * weights,
    )
    out = work.groupby("comp_size", observed=True, sort=True).agg(
        count=("count", "sum"), _sum=("_sum", "sum"), _logsum=("_logsum", "sum"),
        _logsq=("_logsq", "sum"),
    ).reset_index()
    n = out["count"].to_numpy(dtype=float)
    out["mean"] = out["_sum"].to_numpy() / n
    out["logmean"] = out["_logsum"].to_numpy() / n
    out["logvar"] = np.maximum(out["_logsq"].to_numpy() / n - out["logmean"].to_numpy() ** 2, 0.0)
    out["geomean"] = np.power(10.0, out["logmean"].to_numpy())
    out["median"] = np.nan
    out["q25"] = np.nan
    out["q75"] = np.nan
    return out.drop(columns=["_sum", "_logsum", "_logsq"])


# ---------------------------------------------------------------------------------------------
# the fit


@dataclass
class Fit:
    """`ns = A * comp_size^t`, least squares in log10 space."""

    t: float
    t_se: float
    t_se_kind: str
    log10_a: float
    log10_a_se: float
    r2: float
    r2_adj: float
    rmse_dex: float
    max_resid_dex: float
    max_resid_at: int
    points: int
    rows: int
    size_lo: int
    size_hi: int
    basis: str
    weighting: str
    split: tuple | None = None          # (t_lo, se_lo, t_hi, se_hi, z, cut)
    runs_z: float | None = None
    runs_p: float | None = None
    against: dict = field(default_factory=dict)   # null exponent -> (z, p)
    sensitivity: dict = field(default_factory=dict)  # basis -> t

    @property
    def ci(self) -> tuple[float, float]:
        half = student_q(self.points - 2) * self.t_se
        return self.t - half, self.t + half

    def curve(self, sizes: np.ndarray) -> np.ndarray:
        return np.power(10.0, self.log10_a + self.t * np.log10(sizes))

    def headline(self) -> str:
        lo, hi = self.ci
        return f"t = {self.t:.3f} +/- {self.t_se:.3f} (95% CI {lo:.3f}..{hi:.3f}), R^2 = {self.r2:.4f}"


def student_q(df: int) -> float:
    """Two-sided 95% quantile. scipy if it is installed, a Welch-style widening of `Z95` if not."""
    if df < 1:
        return float("nan")
    try:
        from scipy import stats  # noqa: PLC0415  (optional; the fallback below is exact enough)

        return float(stats.t.ppf(0.975, df))
    except Exception:
        # Within 0.3% of the true quantile for df >= 8 and conservative below it.
        return Z95 * (1.0 + 1.4 / df + 2.6 / df**2)


def normal_p(z: float) -> float:
    """Two-sided p for a standard normal statistic."""
    return math.erfc(abs(z) / math.sqrt(2.0))


def _wls(x: np.ndarray, y: np.ndarray, w: np.ndarray):
    design = np.column_stack([np.ones_like(x), x])
    dtw = design.T * w
    normal = dtw @ design
    beta = np.linalg.solve(normal, dtw @ y)
    resid = y - design @ beta
    return beta, resid, design, np.linalg.inv(normal)


def _sandwich(design, inv_normal, w, resid, clusters: np.ndarray | None, n_obs: int):
    """HC1 (or, given clusters, cluster-robust CR1) covariance of the WLS coefficients.

    The residual spread of a component-solve time grows with the component, so the constant-variance
    standard error would understate `t`; and every point of a per-size aggregate stands for many
    correlated rows, so when the fit claims to be over rows the cluster form is the honest one."""
    scores = design * (w * resid)[:, None]
    if clusters is None:
        meat = scores.T @ scores
        n, k = len(resid), design.shape[1]
        factor = n / max(n - k, 1)
    else:
        order = np.argsort(clusters, kind="stable")
        edges = np.flatnonzero(np.r_[True, np.diff(clusters[order]) != 0])
        summed = np.add.reduceat(scores[order], edges, axis=0)
        meat = summed.T @ summed
        groups, k = len(edges), design.shape[1]
        factor = (groups / max(groups - 1, 1)) * ((n_obs - 1) / max(n_obs - k, 1))
    return inv_normal @ meat @ inv_normal * factor


def _runs_z(resid: np.ndarray) -> tuple[float, float] | tuple[None, None]:
    """Wald-Wolfowitz runs test on the residual signs taken in size order.

    A power law that has curvature in it leaves residuals that sweep negative-positive-negative in
    long blocks rather than alternating, and that is invisible to R^2."""
    signs = np.sign(resid)
    signs = signs[signs != 0]
    n = len(signs)
    pos = int((signs > 0).sum())
    neg = n - pos
    if pos == 0 or neg == 0 or n < 8:
        return None, None
    runs = 1 + int((np.diff(signs) != 0).sum())
    mu = 2.0 * pos * neg / n + 1.0
    var = 2.0 * pos * neg * (2.0 * pos * neg - n) / (n * n * (n - 1.0))
    if var <= 0.0:
        return None, None
    z = (runs - mu) / math.sqrt(var)
    return z, normal_p(z)


def fit_power_law(agg: pd.DataFrame, basis: str, weighting: str, nulls: tuple[float, ...]) -> Fit | None:
    """`basis` picks what is regressed: `rows` (every raw row, via the per-size log moments), or the
    per-size `mean`, `geomean` or `median`."""
    frame = agg[agg["count"] > 0]
    if basis == "median":
        frame = frame[frame["median"].notna() & (frame["median"] > 0)]
    if len(frame) < 3:
        return None

    sizes = frame["comp_size"].to_numpy(dtype=float)
    counts = frame["count"].to_numpy(dtype=float)
    x = np.log10(sizes)

    if basis == "rows":
        y = frame["logmean"].to_numpy()
        w = counts
    else:
        y = np.log10(frame[{"mean": "mean", "geomean": "geomean", "median": "median"}[basis]].to_numpy())
        w = counts if weighting == "count" else np.ones_like(counts)

    beta, resid, design, inv_normal = _wls(x, y, w)
    n_rows = int(counts.sum())

    if basis == "rows":
        # Least squares over the raw rows: identical coefficients, but the sums of squares and the
        # covariance must be the row-level ones, rebuilt from the per-size (count, mean, var).
        var_within = frame["logvar"].to_numpy()
        sse = float((counts * (var_within + resid**2)).sum())
        grand = float((counts * y).sum() / counts.sum())
        sst = float((counts * (var_within + (y - grand) ** 2)).sum())
        cov = _sandwich(design, inv_normal, w, resid, frame["comp_size"].to_numpy(), n_rows)
        se_kind = "cluster-robust by component size, over every row"
        dof_n = n_rows
    else:
        sse = float((w * resid**2).sum())
        grand = float((w * y).sum() / w.sum())
        sst = float((w * (y - grand) ** 2).sum())
        cov = _sandwich(design, inv_normal, w, resid, None, len(resid))
        se_kind = "HC1 robust, over the per-size points"
        dof_n = len(resid)

    r2 = 1.0 - sse / sst if sst > 0 else float("nan")
    r2_adj = 1.0 - (1.0 - r2) * (dof_n - 1) / max(dof_n - 2, 1)
    worst = int(np.argmax(np.abs(resid)))

    fit = Fit(
        t=float(beta[1]),
        t_se=float(math.sqrt(max(cov[1, 1], 0.0))),
        t_se_kind=se_kind,
        log10_a=float(beta[0]),
        log10_a_se=float(math.sqrt(max(cov[0, 0], 0.0))),
        r2=r2,
        r2_adj=r2_adj,
        rmse_dex=float(np.sqrt((w * resid**2).sum() / w.sum())),
        max_resid_dex=float(resid[worst]),
        max_resid_at=int(sizes[worst]),
        points=len(resid),
        rows=n_rows,
        size_lo=int(sizes.min()),
        size_hi=int(sizes.max()),
        basis=basis,
        weighting=weighting if basis != "rows" else "count (rows)",
    )

    fit.runs_z, fit.runs_p = _runs_z(resid)
    for null in nulls:
        if fit.t_se > 0:
            z = (fit.t - null) / fit.t_se
            fit.against[null] = (z, normal_p(z))

    # Split-half: the same fit on the lower and upper halves of the *log* size range. A genuine
    # power law gives the same exponent on both; curvature gives two that disagree by more than
    # their errors, which is the failure R^2 hides.
    cut = float(np.median(x))
    lo_mask, hi_mask = x <= cut, x > cut
    if lo_mask.sum() >= 3 and hi_mask.sum() >= 3:
        halves = []
        for mask in (lo_mask, hi_mask):
            b, r, dsg, inv = _wls(x[mask], y[mask], w[mask])
            clusters = frame["comp_size"].to_numpy()[mask] if basis == "rows" else None
            c = _sandwich(dsg, inv, w[mask], r, clusters, int(counts[mask].sum()) if basis == "rows" else int(mask.sum()))
            halves.append((float(b[1]), float(math.sqrt(max(c[1, 1], 0.0)))))
        (t_lo, se_lo), (t_hi, se_hi) = halves
        spread = math.sqrt(se_lo**2 + se_hi**2)
        z = (t_hi - t_lo) / spread if spread > 0 else float("nan")
        fit.split = (t_lo, se_lo, t_hi, se_hi, z, 10.0**cut)

    # Estimator sensitivity: if the exponent moves when the centre moves, the spread at a given
    # size is what is being reported, not the trend.
    for other in ("rows", "mean", "geomean", "median"):
        if other == basis:
            fit.sensitivity[other] = fit.t
            continue
        try:
            side = fit_power_law_bare(frame, other, weighting)
        except Exception:
            side = None
        if side is not None:
            fit.sensitivity[other] = side
    return fit


def fit_power_law_bare(frame: pd.DataFrame, basis: str, weighting: str) -> float | None:
    """Just the exponent, for the sensitivity block. No covariance, no diagnostics."""
    work = frame
    if basis == "median":
        work = work[work["median"].notna() & (work["median"] > 0)]
    if len(work) < 3:
        return None
    x = np.log10(work["comp_size"].to_numpy(dtype=float))
    counts = work["count"].to_numpy(dtype=float)
    if basis == "rows":
        y, w = work["logmean"].to_numpy(), counts
    else:
        y = np.log10(work[{"mean": "mean", "geomean": "geomean", "median": "median"}[basis]].to_numpy())
        w = counts if weighting == "count" else np.ones_like(counts)
    beta, _, _, _ = _wls(x, y, w)
    return float(beta[1])


# ---------------------------------------------------------------------------------------------
# labels


def label_of(key: str, value) -> str:
    if key == "d":
        return f"d = {int(value)}"
    if key == "p":
        return f"p = {value:g}"
    if key == "T":
        return f"T = {value:g}"
    return str(value)


def slug_of(key: str, value) -> str:
    if key == "d":
        return f"d{int(value)}"
    if key == "p":
        return f"p{value:g}"
    if key == "T":
        return f"T{value:g}"
    return str(value).lower()


def ramp(theme: dict, n: int) -> list[str]:
    """`n` steps sampled evenly out of the single-hue ordinal ramp, so the ends are always the ends
    and a series keeps its step as long as the series count does not change."""
    steps = theme["ramp"]
    if n <= 1:
        return [steps[len(steps) // 2]]
    idx = np.linspace(0, len(steps) - 1, n)
    return [steps[int(round(i))] for i in idx]


# ---------------------------------------------------------------------------------------------
# drawing


def style_axes(ax, theme, args, xlabel: str, ylabel: str):
    ax.set_facecolor(theme["surface"])
    ax.grid(True, which="major", color=theme["grid"], linewidth=0.7, alpha=0.9, zorder=0)
    if args.x_scale == "log" or args.y_scale == "log":
        ax.grid(True, which="minor", color=theme["grid"], linewidth=0.4, alpha=0.5, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(theme["axis"])
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=theme["text_secondary"], labelsize=9, width=0.8)
    if xlabel:
        ax.set_xlabel(xlabel, color=theme["text_secondary"], fontsize=10)
    ax.set_ylabel(ylabel, color=theme["text_secondary"], fontsize=10)
    ax.set_axisbelow(True)


def draw_figure(panel, args, theme, meta_line: str, notes: list[str]):
    """One figure: the per-size statistic for every series and the fitted power laws."""
    caption = wrap(panel.caption(meta_line, args))

    # The caption is as long as the run needs it to be -- a selection with many error rates says
    # more than one with a single rate -- so the bottom margin is measured off it rather than fixed,
    # and the x label keeps its room whatever the caption grows to. 1.5 is `linespacing` below.
    height_in = 6.2
    line_frac = CAPTION_FONT_SIZE * 1.5 / (height_in * 72.0)
    bottom = min(0.55, CAPTION_BOTTOM + (caption.count("\n") + 1) * line_frac + 0.055)

    fig = plt.figure(figsize=(10.0, height_in), facecolor=theme["surface"])
    grid = GridSpec(1, 1, figure=fig, left=0.085, right=0.98, top=0.86, bottom=bottom)
    ax = fig.add_subplot(grid[0])

    stat_name = {"mean": "mean", "median": "median", "geomean": "geometric mean"}[args.stat]
    y_label = f"{stat_name} run time per component (ns)"
    style_axes(ax, theme, args, "component size (defects in the component of H)", y_label)
    ax.set_xscale(args.x_scale)
    ax.set_yscale(args.y_scale)

    colour_key = panel.series_keys[0] if panel.series_keys else None
    style_key = panel.series_keys[1] if len(panel.series_keys) > 1 else None
    colours = ramp(theme, max(len(panel.colour_values), 1))

    # The interquartile band is only honest when there is one series on the axes; overlapping bands
    # from twelve series are a wash of colour that hides the very points they are drawn around.
    band_drawn = False
    if args.band != "none" and len(panel.series) == 1 and args.stat != "geomean":
        series = panel.series[0]
        lo, hi = series.agg["q25"].to_numpy(), series.agg["q75"].to_numpy()
        if np.isfinite(lo).all() and np.isfinite(hi).all():
            ax.fill_between(series.agg["comp_size"].to_numpy(), lo, hi, color=theme["band"],
                            alpha=0.28, linewidth=0, zorder=1,
                            label="interquartile range of the raw rows")
            band_drawn = True

    for series in panel.series:
        colour = colours[panel.colour_values.index(series.colour_value)] if colour_key else colours[0]
        marker = MARKERS[panel.style_values.index(series.style_value) % len(MARKERS)] if style_key else MARKERS[0]
        ax.plot(
            series.agg["comp_size"].to_numpy(),
            series.agg[args.stat].to_numpy(),
            linestyle="none",
            marker=marker,
            markersize=args.marker_size,
            markerfacecolor=colour,
            markeredgecolor=theme["surface"],
            markeredgewidth=0.4,
            alpha=0.85,
            zorder=3,
        )

    # Fit lines. One dotted curve per physical error rate, all in one accent colour and told apart
    # by their dash, since `p` already spends the marker and colour is spent on `d`; per-series lines
    # are drawn only when there are few enough of them to tell apart. Nothing is drawn over the whole
    # panel: the pooled exponent is an average across the very rates these curves separate, and it is
    # said in the box, the caption and the table instead of being drawn as one more line.
    fit_lines = []
    if "per-p" in args.fit_lines:
        # A panel already fixed at one `p` by `--split p` grows no per-p fits of its own, since its
        # pooled fit *is* that rate's fit; it is drawn here under its rate's name so that such a
        # figure is not left with no curve at all.
        p_fits = panel.p_fits or (
            [(panel.fixed["p"], panel.pooled_fit)]
            if "p" in panel.fixed and panel.pooled_fit is not None else []
        )
        for i, (value, fit) in enumerate(p_fits):
            dash = DASHES[i % len(DASHES)]
            span = np.geomspace(fit.size_lo, fit.size_hi, 200)
            ax.plot(span, fit.curve(span), color=theme["accent"], linewidth=1.7, linestyle=":",
                    dashes=dash, zorder=5, dash_capstyle="round")
            # Both coefficients now that no other curve carries them; the legend has the room the
            # pooled entry gave back.
            fit_lines.append(Line2D([], [], color=theme["accent"], linewidth=1.7, linestyle=":",
                                    dashes=dash,
                                    label=f"{label_of('p', value)}:  ns = "
                                          f"{10 ** fit.log10_a:.3g} x size^{fit.t:.3f}"))
    if "series" in args.fit_lines and len(panel.series) <= 4:
        for series in panel.series:
            if series.fit is None:
                continue
            colour = colours[panel.colour_values.index(series.colour_value)] if colour_key else colours[0]
            span = np.geomspace(series.fit.size_lo, series.fit.size_hi, 200)
            ax.plot(span, series.fit.curve(span), color=colour, linewidth=1.6, alpha=0.9, zorder=4)
    elif "series" in args.fit_lines:
        notes.append(
            f"{panel.title}: {len(panel.series)} series is too many to draw a fit line for each "
            f"legibly, so no per-series curve was drawn; their exponents are in the table."
        )

    if panel.pooled_fit is not None:
        annotate_fit(ax, panel.pooled_fit, theme, args, len(panel.p_fits) > 1)

    build_legend(ax, panel, theme, colours, colour_key, style_key, fit_lines, band_drawn)

    fig.suptitle(panel.title, x=0.085, ha="left", y=0.965, fontsize=13.5,
                 color=theme["text_primary"], fontweight="semibold")
    fig.text(0.085, 0.917, panel.subtitle, ha="left", fontsize=10, color=theme["text_secondary"])
    fig.text(0.085, CAPTION_BOTTOM, caption, ha="left", va="bottom", fontsize=CAPTION_FONT_SIZE,
             color=theme["text_muted"], linespacing=1.5)
    return fig


def annotate_fit(ax, fit: Fit, theme, args, many_rates: bool):
    """The pooled fit as text. No line on the axes answers to it -- `many_rates` is what says so, by
    naming the scope, so the box is not read as the legend of a curve that is not there."""
    lo, hi = fit.ci
    basis_note = "" if fit.basis == args.stat else f"  (fitted on the {fit.basis}, not the drawn {args.stat})"
    scope_note = ", every p pooled (not drawn)" if many_rates else ""
    lines = [
        f"$ns = A\\,\\cdot\\,size^{{\\,t}}${scope_note}{basis_note}",
        f"t = {fit.t:.3f} $\\pm$ {fit.t_se:.3f}",
        f"95% CI  {lo:.3f} – {hi:.3f}",
        f"$R^2$ = {fit.r2:.4f}",
        f"residual scatter {fit.rmse_dex:.3f} dex ($\\times${10 ** fit.rmse_dex:.2f})",
    ]
    if fit.split is not None:
        t_lo, _, t_hi, _, z, cut = fit.split
        lines.append(f"halves: {t_lo:.3f} | {t_hi:.3f}  ({z:+.1f}$\\sigma$, cut at size {cut:.0f})")
    ax.text(
        0.985, 0.035, "\n".join(lines), transform=ax.transAxes, ha="right", va="bottom",
        fontsize=9, color=theme["text_primary"], linespacing=1.55,
        bbox=dict(boxstyle="round,pad=0.5", facecolor=theme["surface"], edgecolor=theme["grid"],
                  linewidth=0.8, alpha=0.94),
        zorder=6,
    )


def build_legend(ax, panel, theme, colours, colour_key, style_key, fit_lines, band_drawn):
    handles, labels = [], []
    if colour_key and len(panel.colour_values) > 1:
        handles.append(Patch(facecolor="none", edgecolor="none"))
        labels.append(f"{colour_key}:")
        for i, value in enumerate(panel.colour_values):
            handles.append(Line2D([], [], color="none", marker="s", markersize=8,
                                  markerfacecolor=colours[i], markeredgecolor=theme["surface"]))
            labels.append(label_of(colour_key, value))
    if style_key and len(panel.style_values) > 1:
        handles.append(Patch(facecolor="none", edgecolor="none"))
        labels.append(f"{style_key}:")
        for i, value in enumerate(panel.style_values):
            handles.append(Line2D([], [], color=theme["text_secondary"], linestyle="none",
                                  marker=MARKERS[i % len(MARKERS)], markersize=7))
            labels.append(label_of(style_key, value))
    if band_drawn:
        handles.append(Patch(facecolor=theme["band"], alpha=0.28, edgecolor="none"))
        labels.append("interquartile range")
    if len(fit_lines) > 1:
        handles.append(Patch(facecolor="none", edgecolor="none"))
        labels.append("fitted power laws:")
    for line in fit_lines:
        handles.append(line)
        labels.append(line.get_label())
    if not handles:
        return
    legend = ax.legend(
        handles, labels, loc="upper left", frameon=True, fontsize=8.6, ncol=1,
        labelspacing=0.35, handletextpad=0.6, borderpad=0.6,
        facecolor=theme["surface"], edgecolor=theme["grid"],
    )
    legend.get_frame().set_linewidth(0.8)
    for text in legend.get_texts():
        text.set_color(theme["text_secondary"])


def wrap(text: str, width: int = 148) -> str:
    return "\n".join(
        "\n".join(textwrap.wrap(part, width=width)) if part.strip() else ""
        for part in text.split("\n")
    )


# ---------------------------------------------------------------------------------------------
# panels


@dataclass
class Series:
    label: str
    colour_value: object
    style_value: object
    agg: pd.DataFrame
    fit: Fit | None
    rows: int
    key_values: dict


@dataclass
class Panel:
    title: str
    subtitle: str
    slug: str
    fixed: dict
    series_keys: list
    series: list
    colour_values: list
    style_values: list
    pooled_fit: Fit | None
    pooled_agg: pd.DataFrame
    counts: dict
    rows: int
    p_fits: list = field(default_factory=list)   # [(p, Fit)], each pooled over every other key

    def caption(self, meta_line: str, args) -> str:
        parts = []
        parts.append(
            f"{self.rows:,} component solves over {len(self.pooled_agg):,} distinct component "
            f"sizes, {self.pooled_agg['comp_size'].min()} to {self.pooled_agg['comp_size'].max()}. "
            f"One raw row is one connected component of H, timed around sub-H construction, the "
            f"truncated solve to the horizon, and the harvest; escalation is neither run nor timed."
        )
        if args.overhead_used:
            parts.append(
                f"{args.overhead_used:g} ns of calibrated timer overhead has been subtracted from "
                f"every row (the table stores raw deltas; the subtraction belongs here)."
            )
        else:
            parts.append("Times are the raw uncorrected timer deltas: no overhead was subtracted.")
        if args.min_logged is not None:
            parts.append(
                f"Components of size < {args.min_logged} are absent from the record by design, not "
                f"missing from the run: their solve is of the order of one timer read, so the "
                f"profiler times them and logs no row. The x axis starts where the record does."
            )
        status_note = ", ".join(f"{v:,} {k}" for k, v in sorted(self.counts.items()))
        parts.append(f"Status: {status_note}.")
        if self.pooled_fit is not None:
            fit = self.pooled_fit
            basis = {"rows": "every raw row (equivalently, the per-size mean of log10 ns weighted "
                             "by the row count)",
                     "mean": "the per-size arithmetic mean",
                     "geomean": "the per-size geometric mean",
                     "median": "the per-size median"}[fit.basis]
            scope = (" pooled over every error rate on the axes, which no drawn curve answers to,"
                     if len(self.p_fits) > 1 else "")
            parts.append(
                f"Fit:{scope} least squares of log10(ns) on log10(size) over {basis}, sizes "
                f"{fit.size_lo}-{fit.size_hi}; the standard error on t is {fit.t_se_kind}."
            )
            if fit.split is not None:
                t_lo, _, t_hi, _, z, cut = fit.split
                verdict = "consistent with one exponent" if abs(z) < 2 else "not one exponent across the range"
                parts.append(
                    f"Split at size {cut:.0f}, the two halves give t = {t_lo:.3f} and t = {t_hi:.3f}, "
                    f"differing by {z:+.1f} sigma -- {verdict}."
                )
        if self.p_fits and "per-p" in args.fit_lines:
            per_p = "; ".join(
                f"p = {value:g}: t = {fit.t:.3f} +/- {fit.t_se:.3f} over {fit.rows:,} rows"
                for value, fit in self.p_fits
            )
            told_apart = ("" if len(self.p_fits) == 1 else
                          " They are told apart by their dash pattern rather than by a colour, "
                          "which the series already spend.")
            parts.append(
                f"Each dotted curve is the same fit over one physical error rate alone, pooled over "
                f"every other key -- {per_p}.{told_apart}"
            )
        parts.append(meta_line)
        return "  ".join(p for p in parts if p)


def build_panels(frame: pd.DataFrame, args, notes: list[str]) -> list[Panel]:
    split_keys = args.split
    series_keys = args.series_by
    group_keys = list(dict.fromkeys(split_keys + series_keys))

    agg = aggregate(frame, group_keys) if group_keys else aggregate(frame, [])
    if not group_keys:
        agg = aggregate(frame.assign(_all=0), ["_all"])

    panels = []
    split_groups = ([((), frame, agg)] if not split_keys else
                    [(values if isinstance(values, tuple) else (values,), sub, None)
                     for values, sub in frame.groupby(split_keys, observed=True, sort=True)])

    for values, sub, pre in split_groups:
        fixed = dict(zip(split_keys, values))
        sub_agg = pre if pre is not None else aggregate(sub, series_keys) if series_keys else \
            aggregate(sub.assign(_all=0), ["_all"])

        series_list, colour_values, style_values = [], [], []
        if series_keys:
            for svalues, sagg in sub_agg.groupby(series_keys, observed=True, sort=True):
                svalues = svalues if isinstance(svalues, tuple) else (svalues,)
                kv = dict(zip(series_keys, svalues))
                sagg = sagg[sagg["count"] >= args.min_count]
                if sagg.empty:
                    continue
                colour_values.append(svalues[0])
                style_values.append(svalues[1] if len(svalues) > 1 else None)
                series_list.append(Series(
                    label=", ".join(label_of(k, v) for k, v in kv.items()),
                    colour_value=svalues[0],
                    style_value=svalues[1] if len(svalues) > 1 else None,
                    agg=sagg,
                    fit=fit_power_law(restrict(sagg, args), args.fit_basis, args.fit_weights, args.against),
                    rows=int(sagg["count"].sum()),
                    key_values=kv,
                ))
        else:
            sagg = sub_agg[sub_agg["count"] >= args.min_count]
            series_list.append(Series("all", None, None, sagg,
                                      fit_power_law(restrict(sagg, args), args.fit_basis,
                                                    args.fit_weights, args.against),
                                      int(sagg["count"].sum()), {}))
            colour_values, style_values = [None], [None]

        if not series_list:
            continue
        colour_values = sorted(dict.fromkeys(colour_values), key=lambda v: (v is None, v))
        style_values = sorted(dict.fromkeys(style_values), key=lambda v: (v is None, v))

        pooled_agg = pool(sub_agg) if series_keys else sub_agg
        pooled_agg = pooled_agg[pooled_agg["count"] >= args.min_count]
        if args.fit_basis == "median" and series_keys:
            # A pooled median cannot be combined out of per-series aggregates, so it is re-derived
            # from the rows rather than approximated.
            pooled_agg = aggregate(sub.assign(_all=0), ["_all"])
            pooled_agg = pooled_agg[pooled_agg["count"] >= args.min_count]
        pooled_fit = fit_power_law(restrict(pooled_agg, args), args.fit_basis, args.fit_weights, args.against)
        if pooled_fit is None:
            notes.append(f"{describe(fixed) or 'the selection'}: too few distinct sizes to fit.")

        # One fit per physical error rate, pooled over every other key (`d` in particular), so the
        # figure can say whether the exponent is a property of the decoder or of the noise. It is
        # aggregated from the rows rather than combined out of `sub_agg`, so it is exact for the
        # median basis too. A panel already fixed at one `p` by `--split p` has its `p` in the title
        # and its fit in the pooled line, so it grows no second identical curve here.
        p_fits = []
        if "p" not in fixed:
            for value, p_sub in sub.groupby("p", observed=True, sort=True):
                p_agg = aggregate(p_sub.assign(_all=0), ["_all"])
                p_agg = p_agg[p_agg["count"] >= args.min_count]
                fit = fit_power_law(restrict(p_agg, args), args.fit_basis, args.fit_weights, args.against)
                if fit is None:
                    notes.append(
                        f"{describe(fixed) or 'the selection'}, {label_of('p', value)}: too few "
                        f"distinct sizes to fit, so no line is drawn for it."
                    )
                    continue
                p_fits.append((float(value), fit))

        counts = sub["status"].value_counts().to_dict()
        title = "Component run time against component size"
        subtitle_bits = [describe(fixed)] if fixed else []
        subtitle_bits.append(f"{int(sub_agg['count'].sum()):,} component solves")
        if series_keys:
            subtitle_bits.append(f"{len(series_list)} series by {' x '.join(series_keys)}")
        panels.append(Panel(
            title=title,
            subtitle="   ".join(b for b in subtitle_bits if b),
            slug="__".join(slug_of(k, v) for k, v in fixed.items()) or "all",
            fixed=fixed,
            series_keys=series_keys,
            series=series_list,
            colour_values=colour_values,
            style_values=style_values,
            pooled_fit=pooled_fit,
            pooled_agg=pooled_agg,
            counts={str(k): int(v) for k, v in counts.items() if v},
            rows=int(sub_agg["count"].sum()),
            p_fits=p_fits,
        ))
    return panels


def restrict(agg: pd.DataFrame, args) -> pd.DataFrame:
    """The fit range, which is allowed to be narrower than the drawn range."""
    out = agg
    if args.fit_min_size is not None:
        out = out[out["comp_size"] >= args.fit_min_size]
    if args.fit_max_size is not None:
        out = out[out["comp_size"] <= args.fit_max_size]
    if args.fit_min_count > 1:
        out = out[out["count"] >= args.fit_min_count]
    return out


def describe(fixed: dict) -> str:
    return ", ".join(label_of(k, v) for k, v in fixed.items())


# ---------------------------------------------------------------------------------------------
# tables


def entries_of(panel: Panel) -> list[tuple[str, Fit | None, int, dict]]:
    """`(label, fit, rows, key values)` for the pooled fit, then the per-p fits, then the series --
    the same order and the same three scopes the figure draws, so a reader can move between the
    figure, the csv and the report without re-deriving which line is which.

    A per-p row count is the count inside the fitted size range rather than the whole `p` slice,
    since that fit is what the row reports."""
    entries = [("pooled", panel.pooled_fit, panel.rows, {})]
    entries += [(f"{label_of('p', value)}, pooled", fit, fit.rows, {"p": value})
                for value, fit in panel.p_fits]
    entries += [(s.label, s.fit, s.rows, s.key_values) for s in panel.series]
    return entries


def fit_rows(panels: list[Panel], args) -> list[dict]:
    rows = []
    for panel in panels:
        for label, fit, n_rows, kv in entries_of(panel):
            if fit is None:
                continue
            lo, hi = fit.ci
            row = {
                **{k: v for k, v in panel.fixed.items()},
                **kv,
                "scope": label,
                "rows": n_rows,
                "points": fit.points,
                "size_lo": fit.size_lo,
                "size_hi": fit.size_hi,
                "t": fit.t,
                "t_se": fit.t_se,
                "t_ci_lo": lo,
                "t_ci_hi": hi,
                "A_ns": 10.0 ** fit.log10_a,
                "log10_A": fit.log10_a,
                "log10_A_se": fit.log10_a_se,
                "r2": fit.r2,
                "r2_adj": fit.r2_adj,
                "rmse_dex": fit.rmse_dex,
                "rmse_factor": 10.0 ** fit.rmse_dex,
                "max_resid_dex": fit.max_resid_dex,
                "max_resid_at_size": fit.max_resid_at,
                "runs_z": fit.runs_z,
                "runs_p": fit.runs_p,
                "basis": fit.basis,
                "weighting": fit.weighting,
                "se_kind": fit.t_se_kind,
            }
            if fit.split is not None:
                t_lo, se_lo, t_hi, se_hi, z, cut = fit.split
                row.update({"t_lower_half": t_lo, "t_lower_half_se": se_lo,
                            "t_upper_half": t_hi, "t_upper_half_se": se_hi,
                            "split_z": z, "split_at_size": cut})
            for null, (z, p) in fit.against.items():
                row[f"z_vs_t{null:g}"] = z
                row[f"p_vs_t{null:g}"] = p
            for basis, value in fit.sensitivity.items():
                row[f"t_on_{basis}"] = value
            rows.append(row)
    return rows


def fit_report(panels: list[Panel], args, meta_line: str, notes: list[str]) -> str:
    out = ["# Power-law fits of component run time against component size", ""]
    out.append(f"Model `ns = A * comp_size^t`, fitted as least squares of `log10(ns)` on "
               f"`log10(comp_size)` over **{args.fit_basis}**"
               f"{'' if args.fit_basis == 'rows' else f' (weights: {args.fit_weights})'}.")
    out.append("")
    out.append(meta_line)
    out.append("")
    out.append("Three scopes appear in every table below: **pooled**, every row of the selection at "
               "once -- reported here and in the caption, but drawn as no curve, since it averages "
               "across the error rates the figure separates; `p = ..., pooled`, one physical error "
               "rate with every other key pooled into it -- these are the dotted curves and their "
               "legend entries; and one row per drawn series.")
    out.append("")
    for panel in panels:
        out.append(f"## {describe(panel.fixed) or 'all rows'}")
        out.append("")
        header = ("| scope | rows | sizes | size range | t | 95% CI | A (ns) | R^2 | adj R^2 | "
                  "scatter (dex) | halves (lo / hi) | split Z | runs Z |")
        out.append(header)
        out.append("|" + "---|" * 13)
        entries = [(f"**{label}**" if label == "pooled" else label, fit, n_rows)
                   for label, fit, n_rows, _ in entries_of(panel)]
        for label, fit, n_rows in entries:
            if fit is None:
                out.append(f"| {label} | {n_rows:,} | — | — | not fitted (too few distinct sizes) "
                           "| | | | | | | | |")
                continue
            lo, hi = fit.ci
            halves = split_z = "—"
            if fit.split is not None:
                t_lo, _, t_hi, _, z, _ = fit.split
                halves = f"{t_lo:.3f} / {t_hi:.3f}"
                split_z = f"{z:+.1f}"
            runs = "—" if fit.runs_z is None else f"{fit.runs_z:+.1f} (p={fit.runs_p:.2g})"
            out.append(
                f"| {label} | {n_rows:,} | {fit.points} | {fit.size_lo}–{fit.size_hi} | "
                f"**{fit.t:.4f}** ± {fit.t_se:.4f} | {lo:.4f}–{hi:.4f} | {10 ** fit.log10_a:.4g} | "
                f"{fit.r2:.4f} | {fit.r2_adj:.4f} | {fit.rmse_dex:.4f} (×{10 ** fit.rmse_dex:.2f}) | "
                f"{halves} | {split_z} | {runs} |"
            )
        out.append("")
        fit = panel.pooled_fit
        if fit is not None:
            if fit.against:
                out.append("Pooled exponent against round nulls: " + "; ".join(
                    f"`t = {null:g}` → Z = {z:+.2f}, p = {p:.3g}"
                    for null, (z, p) in sorted(fit.against.items())) + ".")
            if fit.sensitivity:
                # `pool()` drops the quantiles rather than approximating them, so a pooled row set
                # has no median to refit on unless the median *is* the basis; say so instead of
                # letting the entry go missing.
                shown = [f"{basis} → t = {value:.4f}" for basis, value in fit.sensitivity.items()
                         if value is not None]
                if "median" not in fit.sensitivity or fit.sensitivity["median"] is None:
                    shown.append("median → n/a (a median does not combine across series; "
                                 "`--fit-on median` re-derives it from the rows)")
                out.append("Same fit with the centre changed: " + "; ".join(shown) + ".")
            out.append(
                f"Largest single-point residual {fit.max_resid_dex:+.3f} dex at size "
                f"{fit.max_resid_at}."
            )
            out.append("")
    out.append("## How to read the fit-quality columns")
    out.append("")
    out.extend([
        "* **t ± se** — the exponent and a standard error that does not assume constant residual "
        "variance (the spread of a solve time grows with the component). Over raw rows the error "
        "is additionally clustered by component size, since the rows at one size are not "
        "independent draws about the trend line.",
        "* **R^2** — share of the variance in `log10(ns)` the straight line explains. It is the "
        "weakest of these columns: a curved relation over two decades still scores high.",
        "* **scatter (dex)** — root-mean-square residual in powers of ten; `×f` is the same number "
        "as a multiplicative factor, i.e. a typical point sits within a factor of `f` of the fit.",
        "* **halves** — the exponent refitted on the lower and upper half of the *log* size range, "
        "and **split Z**, their difference in units of its own error. This is the column that "
        "decides whether one power law describes the range: |Z| under about 2 says yes, a large "
        "|Z| says the exponent drifts and a single `t` is an average of two different regimes.",
        "* **runs Z** — Wald–Wolfowitz runs test on the signs of the residuals taken in size order. "
        "Strongly negative means the signs come in long blocks, i.e. systematic curvature rather "
        "than scatter. Strongly positive means they alternate more than chance.",
        "* **t on …** — the exponent when the per-size centre is the mean, the geometric mean or "
        "the median instead. If these disagree, the size-dependence of the *spread* is being read "
        "as trend by at least one of them.",
    ])
    if notes:
        out.append("")
        out.append("## Notes emitted while reading the input")
        out.append("")
        out.extend(f"* {note}" for note in notes)
    return "\n".join(out) + "\n"


def stdout_table(panels: list[Panel]):
    header = f"{'selection':<34}{'scope':<22}{'rows':>12}  {'t':>18}  {'R^2':>8}  {'dex':>7}  {'splitZ':>7}"
    print(header)
    print("-" * len(header))
    for panel in panels:
        selection = describe(panel.fixed) or "all rows"
        entries = [(label, fit, n_rows) for label, fit, n_rows, _ in entries_of(panel)]
        for i, (label, fit, n_rows) in enumerate(entries):
            left = selection if i == 0 else ""
            if fit is None:
                print(f"{left:<34}{label:<22}{n_rows:>12,}  {'not fitted':>18}")
                continue
            split = "—" if fit.split is None else f"{fit.split[4]:+.1f}"
            print(f"{left:<34}{label:<22}{n_rows:>12,}  "
                  f"{fit.t:>9.4f} ± {fit.t_se:<6.4f}  {fit.r2:>8.4f}  {fit.rmse_dex:>7.3f}  {split:>7}")
        print()


# ---------------------------------------------------------------------------------------------
# main


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("paths", nargs="*", metavar="DIR_OR_CSV",
                        help="results directories holding components.csv (and run.log), or csv paths")
    parser.add_argument("--in", dest="inputs", action="append", default=[], metavar="DIR_OR_CSV",
                        help="same as a positional argument; repeatable")
    parser.add_argument("--out", default=None, metavar="DIR",
                        help="where the figures go; default is the first input + /figures")
    parser.add_argument("--format", default="png", choices=["png", "pdf", "svg"])
    parser.add_argument("--dpi", type=int, default=200)
    parser.add_argument("--theme", default="light", choices=sorted(THEMES))
    parser.add_argument("--show", action="store_true", help="also open interactive windows")

    axes = parser.add_argument_group("axes")
    axes.add_argument("--x-scale", default="log", choices=["log", "linear"],
                      help="component size axis (default log: the sizes span three decades)")
    axes.add_argument("--y-scale", default="log", choices=["log", "linear"],
                      help="run time axis (default log)")
    axes.add_argument("--linear", action="store_true",
                      help="shorthand for --x-scale linear --y-scale linear")
    axes.add_argument("--loglog", action="store_true",
                      help="shorthand for --x-scale log --y-scale log (the default)")

    select = parser.add_argument_group("selection")
    select.add_argument("--d", default=None, help="comma list; default = every d present")
    select.add_argument("--p", default=None, help="comma list; default = every p present")
    select.add_argument("--T", dest="T", default=None, help="comma list; default = every T present")
    select.add_argument("--status", default="pooled",
                        choices=["pooled", "split", "complete", "truncated"],
                        help="pooled (default): COMPLETE and TRUNCATED go on the same axes and the "
                             "caption says how many of each; split: a figure each, since they are "
                             "different quantities at the same size")
    select.add_argument("--split-status", action="store_true",
                        help="shorthand for --status split")
    select.add_argument("--min-size", type=int, default=None, help="drop sizes below this")
    select.add_argument("--max-size", type=int, default=None, help="drop sizes above this")
    select.add_argument("--min-count", type=int, default=1, metavar="N",
                        help="drop a size with fewer than N rows from the drawing (default 1)")

    layout = parser.add_argument_group("layout")
    layout.add_argument("--split", default=None, metavar="KEYS",
                        help="comma list from d,p,T,status: one figure per combination "
                             "(default T,status, or T alone if --status is not 'split')")
    layout.add_argument("--series-by", default="d,p", metavar="KEYS",
                        help="comma list from d,p,T,status: the series inside a figure. The first "
                             "key carries colour (a single-hue ordinal ramp), the second the marker. "
                             "Empty or 'none' pools everything into one series")
    layout.add_argument("--stat", default="mean", choices=["mean", "median", "geomean"],
                        help="the per-size statistic that is drawn (default mean)")
    layout.add_argument("--band", default="iqr", choices=["iqr", "none"],
                        help="interquartile band around a single drawn series (default iqr)")
    layout.add_argument("--marker-size", type=float, default=4.0)

    fitting = parser.add_argument_group("power-law fit  (ns = A * comp_size^t)")
    fitting.add_argument("--fit-on", dest="fit_basis", default="auto",
                         choices=["auto", "rows", "mean", "geomean", "median"],
                         help="what is regressed. 'auto' (default) follows --stat, so the fitted "
                              "line is the line through the drawn points; 'rows' regresses every "
                              "raw row in log space")
    fitting.add_argument("--fit-weights", default="count", choices=["count", "none"],
                         help="per-size weight when the basis is a per-size statistic (default count)")
    fitting.add_argument("--fit-lines", default="per-p", metavar="LIST",
                         help="comma list from per-p,series (or 'none'): which fitted curves are "
                              "drawn. Default 'per-p' -- one dotted accent curve per physical error "
                              "rate, each with its own legend entry giving A and t. 'series' adds a "
                              "curve per drawn series, and only when there are at most four of them. "
                              "No curve is drawn over the whole panel; the pooled fit is reported in "
                              "the box, the caption and the table. 'pooled' is accepted and ignored, "
                              "and 'both' is kept as a spelling of per-p,series")
    fitting.add_argument("--fit-min-size", type=int, default=None,
                         help="fit only sizes at or above this (the drawn range is unaffected)")
    fitting.add_argument("--fit-max-size", type=int, default=None, help="fit only sizes at or below this")
    fitting.add_argument("--fit-min-count", type=int, default=1, metavar="N",
                         help="fit only sizes backed by at least N rows (default 1)")
    fitting.add_argument("--against", default="1,2", metavar="LIST",
                         help="comma list of round exponents to test t against (default 1,2; "
                              "empty to test against none)")

    timer = parser.add_argument_group("timer")
    timer.add_argument("--timer-overhead-ns", type=float, default=None,
                       help="override the value run.log calibrated")
    timer.add_argument("--no-overhead-correction", dest="correct", action="store_false",
                       help="plot the raw uncorrected deltas the table stores")

    parser.add_argument("--table", action="store_true", help="print the fit table to stdout as well")
    parser.add_argument("--table-name", default="component_time_fits", metavar="STEM",
                        help="stem of the .md and .csv written beside the figures")

    args = parser.parse_args(argv)
    notes: list[str] = []

    paths = args.paths + args.inputs
    if not paths:
        parser.error("give at least one results directory (or --in DIR)")
    if args.linear and args.loglog:
        parser.error("--linear and --loglog contradict each other")
    if args.linear:
        args.x_scale = args.y_scale = "linear"
    if args.loglog:
        args.x_scale = args.y_scale = "log"
    if args.split_status:
        if args.status not in ("pooled", "split"):
            parser.error(f"--split-status contradicts --status {args.status}")
        args.status = "split"

    try:
        args.d = parse_list(args.d, int)
        args.p = parse_list(args.p, float)
        args.T = parse_list(args.T, float)
        args.against = tuple(parse_list(args.against, float) or ())
    except ValueError as error:
        parser.error(str(error))

    fit_lines = [k for k in (args.fit_lines or "").replace(" ", "").split(",") if k and k != "none"]
    if "both" in fit_lines:  # what the flag meant when it was a single choice
        fit_lines = ["per-p", "series"]
    unknown = [k for k in fit_lines if k not in ("pooled", "per-p", "series")]
    if unknown:
        parser.error(f"unknown --fit-lines value(s) {', '.join(unknown)}; pick from per-p, series, none")
    # `pooled` no longer draws anything. An old command line that asks for it is honoured to the
    # extent that it still can be -- the pooled fit is in the box, the caption and the table -- and
    # told so, rather than failing on a word that used to work.
    if "pooled" in fit_lines:
        notes.append("--fit-lines pooled no longer draws a curve: the panel-wide fit is reported in "
                     "the box on the axes, in the caption and in the fit table instead.")
        fit_lines = [k for k in fit_lines if k != "pooled"]
    args.fit_lines = fit_lines

    args.series_by = [k for k in (args.series_by or "").replace(" ", "").split(",")
                      if k and k != "none"]
    if args.split is None:
        args.split = ["T", "status"] if args.status == "split" else ["T"]
    else:
        args.split = [k for k in args.split.replace(" ", "").split(",") if k and k != "none"]
    for key in args.split + args.series_by:
        if key not in KEYS:
            parser.error(f"unknown key {key!r}; pick from {', '.join(KEYS)}")
    overlap = set(args.split) & set(args.series_by)
    if overlap:
        parser.error(f"{', '.join(sorted(overlap))} cannot both split figures and index series")
    # A selection narrowed to one status has nothing left to split on; an explicit `--split ...,status`
    # under the pooled default is honoured, since asking for it is the same request `--split-status` makes.
    if args.status in ("complete", "truncated") and "status" in args.split:
        args.split = [k for k in args.split if k != "status"]

    try:
        inputs = resolve_inputs(paths)
        frame, metas = load(inputs, notes)
    except (FileNotFoundError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1

    meta = next((m for m in metas if m.timer_overhead_ns is not None), metas[0] if metas else RunMeta())
    overheads = {m.timer_overhead_ns for m in metas if m.timer_overhead_ns is not None}
    if len(overheads) > 1:
        notes.append(
            f"the inputs calibrated different timer overheads ({', '.join(f'{o:g}' for o in sorted(overheads))} ns); "
            f"{meta.timer_overhead_ns:g} ns was subtracted from all of them."
        )
    overhead = 0.0
    if args.correct:
        overhead = args.timer_overhead_ns if args.timer_overhead_ns is not None else (meta.timer_overhead_ns or 0.0)
        if overhead == 0.0 and args.timer_overhead_ns is None:
            notes.append("no run.log timer_overhead_ns was found, so nothing was subtracted.")
    args.overhead_used = overhead
    args.min_logged = meta.min_logged_comp_size

    frame = apply_filters(frame, args, notes)
    if frame.empty:
        print("error: the selection is empty", file=sys.stderr)
        for note in notes:
            print(f"  {note}", file=sys.stderr)
        return 1
    frame = correct_overhead(frame, overhead, notes)
    if frame.empty:
        print("error: every row fell to <= 0 ns after the overhead subtraction", file=sys.stderr)
        return 1

    if args.fit_basis == "auto":
        args.fit_basis = args.stat

    panels = build_panels(frame, args, notes)
    if not panels:
        print("error: nothing to draw after --min-count", file=sys.stderr)
        return 1

    out_dir = Path(args.out) if args.out else (
        (Path(paths[0]).expanduser() if Path(paths[0]).expanduser().is_dir()
         else Path(paths[0]).expanduser().parent) / "figures"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    meta_line = meta.describe()
    if meta_line:
        meta_line = f"Run: {meta_line}."

    theme = THEMES[args.theme]
    written = []
    for panel in panels:
        fig = draw_figure(panel, args, theme, meta_line, notes)
        name = f"component_time_vs_size__{panel.slug}__{args.x_scale}x_{args.y_scale}y.{args.format}"
        path = out_dir / name
        fig.savefig(path, dpi=args.dpi, facecolor=theme["surface"],
                    metadata=FIG_METADATA.get(args.format))
        written.append((path, panel))
        if not args.show:
            plt.close(fig)

    rows = fit_rows(panels, args)
    table = pd.DataFrame(rows)
    # The selection keys first, then the fit, then the diagnostics: a reader opening the csv should
    # meet the columns in the order the report discusses them.
    lead = [c for c in ("d", "p", "T", "status", "scope", "rows", "points", "size_lo", "size_hi")
            if c in table.columns]
    table = table[lead + [c for c in table.columns if c not in lead]]
    csv_path = out_dir / f"{args.table_name}.csv"
    table.to_csv(csv_path, index=False)
    md_path = out_dir / f"{args.table_name}.md"
    md_path.write_text(fit_report(panels, args, meta_line, notes))

    index = ["# component_time_vs_size", "", meta_line, ""]
    for path, panel in written:
        fit = panel.pooled_fit
        line = f"* `{path.name}` — {describe(panel.fixed) or 'all rows'}, {panel.rows:,} solves"
        if fit is not None:
            line += f"; pooled {fit.headline()}"
        index.append(line)
    index += ["", f"* `{csv_path.name}` — every fit, one row each",
              f"* `{md_path.name}` — the same, with what the columns mean"]
    (out_dir / "component_times_index.md").write_text("\n".join(index) + "\n")

    print(f"wrote {len(written)} figure(s) to {out_dir}")
    for path, _ in written:
        print(f"  {path.name}")
    print(f"  {csv_path.name}")
    print(f"  {md_path.name}")
    print("  component_times_index.md")
    if notes:
        print("\nnotes:")
        for note in notes:
            print(f"  - {note}")
    if args.table:
        print()
        stdout_table(panels)
    if args.show:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
