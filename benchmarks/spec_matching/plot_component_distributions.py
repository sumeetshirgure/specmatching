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

"""Draws the component size distribution of each requested `(d, p, T)` as a filled log-log histogram.

    python benchmarks/spec_matching/plot_component_distributions.py --in RESULTS_DIR [MORE_DIRS ...]
        --cells d=17,p=1e-3,T=1.5 [d=21,p=1e-3,T=2 ...]
        [--out-dir DIR] [--formats png,pdf,svg] [--dpi N] [--theme light|dark] [--show]

An input is a `sparse_graph_stats` output directory holding `hists.json`. Its
`component_size_hist` bin `k` counts the components of size exactly `k`, so bin `k` is drawn as the
bar spanning `[k - 0.5, k + 0.5]` at height `count / components in the cell`.

One or two cells share an axes, in the two hues of the corpus palette that hold their separation
under every colour-vision check. Three or more are drawn as a grid of panels, one hue each, so a
comparison never needs a third hue.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

import matplotlib

# The backend has to be chosen before `pyplot` is imported, and `--show` is what changes it.
if "--show" not in sys.argv:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402  (must follow the backend selection above)

# The corpus palette, shared with the other plotters here so a reader moving between their output
# sees one visual language. Slot 0 against slot 1 is ΔE 24.7 and passes every colour-vision check
# against both surfaces; slot 1 against slot 2 is ΔE 4.4 under protanopia, which is why a third
# category becomes a panel rather than a third hue.
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

MAX_COLUMNS = 3


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

    `hists.json` writes its key with `%g`, which is six significant digits. Rounding a value given
    on the command line to that same precision is what lets `1e-3` and `0.001` name one cell.
    """
    return (int(float(d)), float(f"{float(p):.6g}"), float(f"{float(t):.6g}"))


def label_of(key):
    d, p, t = key
    return f"d={d}_p={fmt_p(p)}_T={t:g}"


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


def load_cells(in_dirs):
    """Every cell of every `hists.json`, keyed `(d, p, T)`."""
    cells = {}
    for directory in in_dirs:
        path = Path(directory) / "hists.json"
        if not path.is_file():
            raise SystemExit(f"error: {path} does not exist")
        document = json.loads(path.read_text())
        for raw_key, cell in document.get("cells", document).items():
            if not isinstance(cell, dict) or not raw_key.startswith("d="):
                continue
            fields = dict(token.split("=", 1) for token in raw_key.split(","))
            cells[cell_key(fields["d"], fields["p"], fields["T"])] = cell
    return cells


def distribution_of(cell):
    """`(bin edges, fraction of components, component total)` for sizes 1 and up."""
    counts = np.asarray(cell["component_size_hist"], dtype=float)[1:]
    total = float(counts.sum())
    edges = np.arange(0.5, len(counts) + 1.0)
    fraction = counts / total if total > 0 else counts
    return edges, fraction, total


def style_axes(ax, theme):
    """The recessive frame every panel wears: log-log, grid behind the data, two spines."""
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(color=theme["grid"], linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(theme["grid"])
    ax.tick_params(colors=theme["text_secondary"], labelsize=9)
    ax.set_xlabel("component size", color=theme["text_secondary"], fontsize=10)
    ax.set_ylabel("fraction of components", color=theme["text_secondary"], fontsize=10)


def draw_distribution(ax, edges, fraction, colour, label):
    """Filled at low alpha under a 2px edge, so an overlap is legible either way round."""
    ax.stairs(fraction, edges, color=colour, alpha=0.45, fill=True, label=label, zorder=2)
    ax.stairs(fraction, edges, color=colour, linewidth=2.0, zorder=3)


def draw_overlay(drawn, theme):
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    fig.patch.set_facecolor(theme["surface"])
    ax.set_facecolor(theme["surface"])
    for index, (key, edges, fraction, total) in enumerate(drawn):
        draw_distribution(
            ax,
            edges,
            fraction,
            theme["series"][index],
            f"{title_of(key)}  ({int(total):,} components)",
        )
    style_axes(ax, theme)
    ax.set_xlim(0.5, max(edges[-1] for _, edges, _, _ in drawn))
    # Above the axes rather than inside it: the two fills between them reach most corners of the
    # panel, and a legend over either one hides the shape it is naming.
    legend = ax.legend(frameon=False, fontsize=9, loc="lower left", bbox_to_anchor=(0.0, 1.01))
    for text in legend.get_texts():
        text.set_color(theme["text_primary"])
    fig.tight_layout()
    return fig


def draw_panels(drawn, theme):
    columns = min(len(drawn), MAX_COLUMNS)
    rows = math.ceil(len(drawn) / columns)
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(4.6 * columns, 3.9 * rows),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    fig.patch.set_facecolor(theme["surface"])
    flat = [axes[r][c] for r in range(rows) for c in range(columns)]
    for ax in flat[len(drawn) :]:
        ax.set_visible(False)
    for ax, (key, edges, fraction, total) in zip(flat, drawn):
        ax.set_facecolor(theme["surface"])
        draw_distribution(ax, edges, fraction, theme["series"][0], None)
        style_axes(ax, theme)
        # The shared axes hide the tick labels of every panel that has another one under it, and the
        # last row is short, so the panels above the gap would be left with an unlabelled axis.
        ax.tick_params(labelbottom=True, labelleft=True)
        ax.set_title(
            f"{title_of(key)}\n{int(total):,} components",
            color=theme["text_primary"],
            fontsize=10,
        )
    flat[0].set_xlim(0.5, max(edges[-1] for _, edges, _, _ in drawn))
    fig.tight_layout()
    return fig


def stem_of(keys):
    if len(keys) <= 3:
        return "component_size_distribution__" + "__".join(label_of(key) for key in keys)
    return f"component_size_distribution__{len(keys)}_cells"


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Draw component size distributions as filled log-log histograms.",
    )
    parser.add_argument(
        "--in",
        dest="in_dirs",
        nargs="+",
        required=True,
        metavar="DIR",
        help="one or more sparse_graph_stats output dirs, each holding a hists.json",
    )
    parser.add_argument(
        "--cells",
        nargs="+",
        required=True,
        metavar="d=..,p=..,T=..",
        help="the (d, p, T) settings to draw",
    )
    parser.add_argument("--out-dir", default=None, metavar="DIR", help="default: first --in dir + /figures")
    parser.add_argument("--formats", default="png", help="comma separated: png,pdf,svg")
    parser.add_argument("--dpi", type=int, default=160)
    parser.add_argument("--theme", choices=sorted(THEMES), default="light")
    parser.add_argument("--show", action="store_true", help="also open an interactive window")
    args = parser.parse_args(argv)

    available = load_cells(args.in_dirs)
    keys = [parse_cell(text) for text in args.cells]

    drawn = []
    for key in keys:
        cell = available.get(key)
        if cell is None or "component_size_hist" not in cell:
            raise SystemExit(f"error: no component_size_hist for {title_of(key)} under {', '.join(args.in_dirs)}")
        edges, fraction, total = distribution_of(cell)
        if total <= 0:
            raise SystemExit(f"error: {title_of(key)} holds no components")
        drawn.append((key, edges, fraction, total))

    theme = THEMES[args.theme]
    fig = draw_overlay(drawn, theme) if len(drawn) <= 2 else draw_panels(drawn, theme)

    out_dir = Path(args.out_dir) if args.out_dir else Path(args.in_dirs[0]) / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = stem_of([key for key, _, _, _ in drawn])
    for extension in [token.strip() for token in args.formats.split(",") if token.strip()]:
        path = out_dir / f"{stem}.{extension}"
        fig.savefig(path, dpi=args.dpi, facecolor=theme["surface"])
        print(f"wrote {path}")

    if args.show:
        plt.show()
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main())
