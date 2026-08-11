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

"""Runs the syndrome buffer as a queue driven by the latencies `two_phase_latency_profiler` logged.

    python benchmarks/two_phase/simulate_syndrome_buffer.py \
        benchmarks/two_phase/results/latency

A mean latency does not size a buffer. Syndrome blocks arrive on the hardware's clock whether or not
the previous one has been decoded, so what a deployment has to know is not "how fast is the decoder
on average" but "how far behind does it fall, how often, and for how long" — and the answer is set by
the *shape* of the latency distribution, the escalation tail included, not by its mean. This script
reads the same per-shot logs `plot_latency_histograms.py` draws, and puts the distribution through
the queue it would actually feed:

  * **Arrivals** are deterministic. One block of `rounds` syndrome-extraction rounds arrives every
    `rounds x --round-ns`; the default 1 us round is the usual superconducting figure and is the one
    number here that comes from neither the log nor the code, so it is a flag.
  * **Service** is one shot's decode latency, drawn from the log. The two series are the same two the
    histogram plots — stock exact decode on `G`, and the §M7 front end on the sparsified `H` with its
    escalation tail charged in — simulated separately against the identical arrival stream, so the
    difference between the two buffers is the front end's own doing.
  * **The server** is `--servers` decoders (one by default), FIFO, work-conserving, non-preemptive.

which is a D/G/1 queue (D/G/c for `--servers`). Blocks are never dropped: the buffer is simulated as
unbounded and the *sizing* question is answered from its occupancy tail, because a simulation that
drops what overflows cannot tell you what size would not have overflowed.

What comes out, per series:

  * utilisation rho, and the **critical round period** — how fast the syndrome clock can run before
    rho reaches 1 and the backlog stops draining. A run with rho >= 1 is reported as unstable rather
    than summarised: its "mean occupancy" is a function of how long the simulation was left running
    and nothing else, and none of the steady-state numbers beside it mean anything.
  * the fraction of *time* the buffer is empty, taken as an exact time average over the event stream
    rather than off a sampling grid. For one server this must come out at 1 - rho, and the table
    prints both so the simulation is checkable against the identity it has to satisfy.
  * the occupancy distribution: mean, max, and the tail `P(blocks resident > n)`, which is the buffer
    sizing curve — read off it the depth at which overflow drops below a target.
  * the sojourn tail `P(wait + decode > D)`, which is the deadline-miss curve.

A simulation of `--cycles` blocks cannot see a probability below `1/cycles`, and the interesting
deadline-miss targets are far below anything simulatable. So both tails are also extrapolated, using
the decay rate they are known to have: for a GI/G/1 queue the waiting time's tail is asymptotically
`C exp(-theta x)`, where `theta > 0` solves `E[exp(theta (S - A))] = 1` — the Lundberg exponent,
computed here from the empirical service distribution directly, with no fitted shape and no
distributional assumption. Only the prefactor `C` is fitted, and only against the part of the
simulated tail that has enough samples to be worth fitting. Extrapolated numbers are labelled as
such everywhere they appear, and the figures draw the fit over the simulated curve so the reader can
see how well the two agree before trusting the part that runs off the end.

Two figures per log file, plus a text table beside them:

    buffer_<stem>.png        (a) latency density  (b) latency tail  (c,d) occupancy against time
    buffer_tail_<stem>.png   (e) buffer sizing  (f) deadline miss  (g,h) sweep over the round period
    buffer_simulation.txt    every number in the figures, as text

Panels (a)-(d) are the shape `ppt/buffer_tradeoff.py` draws from a fitted lognormal; here they are
drawn from the measured distribution instead. Panels (e)-(h) are what that figure does not answer:
how deep the buffer has to be, how often a deadline is missed, and how much clock headroom is left.
"""

import argparse
import heapq
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Imports the reader rather than reimplementing it: `series_of` is where the definition of the
# sparsified series lives (which stages are charged, and that escalating shots pay their Phase-2
# re-decode), and a second copy of that definition would be a second thing to keep in step with the
# profiler. Also sets the Agg backend, before pyplot is imported below.
from plot_latency_histograms import THEMES, collect, load, series_of  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.ticker import FuncFormatter, NullFormatter  # noqa: E402

# Fixed slots, the same two the histogram script uses: `G` is always series 1 and the sparsified
# graph always series 2, in every figure either script draws.
STOCK_LABEL = "stock decode on G"
SPARSE_LABEL = "stock decode on sparsified H"


class Series:
    """One latency series and the queue it drives, with everything the figures and table report.

    Built once per (log file, series) at the nominal round period. The sweep panels re-run the queue
    at other periods through `run_queue` and keep only the handful of scalars they plot, because
    holding a full event stream per sweep point is what turns a campaign-sized run into a swap storm.
    """

    def __init__(self, key, label, colour, latency_us, arrival_us, args, rng):
        self.key = key
        self.label = label
        self.colour = colour
        # The measured distribution, kept whole: it is the population every draw comes from, and the
        # density panel plots it rather than the resampled copy.
        self.latency_us = latency_us
        self.arrival_us = arrival_us
        self.mean_us = float(latency_us.mean())
        self.rho = self.mean_us / (arrival_us * args.servers)
        self.stable = self.rho < 1.0
        # The round period at which rho hits 1: the clock this decoder can just barely keep up with,
        # which is the number a hardware schedule is actually chosen against.
        rounds_per_block = arrival_us / (args.round_ns / 1000.0)
        self.critical_round_ns = 1000.0 * self.mean_us / (args.servers * rounds_per_block)

        service = draw_service(latency_us, args.cycles, rng, args.sampling, args.block_len)
        arrivals, _, sojourn, finishes = run_queue(service, arrival_us, args.servers)
        self.times, self.occupancy = occupancy_events(arrivals, finishes)

        # Steady state is quoted after the burn-in only; the queue starts empty, which is not a state
        # a loaded one spends any time in. The occupancy window and the kept sojourns are cut at the
        # same arrival so the two views of the run cover the same stretch of it.
        burn_in = min(args.burn_in, sojourn.size - 1)
        self.window = (arrivals[burn_in], arrivals[-1])
        self.occ = occupancy_stats(self.times, self.occupancy, self.window)
        # Block periods, not microseconds: that is the unit a reader of these figures thinks in, and
        # it is also what keeps the exponent fit's numbers near 1.
        self.sojourn_kept = sojourn[burn_in:] / arrival_us
        self.samples = int(self.sojourn_kept.size)
        # The smallest probability worth quoting off this run. Not `1/samples`: a tail bin holding one
        # or two observations is noise, and a depth the buffer reached twice in the whole run is not
        # evidence that it does not reach it a third time.
        self.resolution = 20.0 / max(self.samples, 1)

        # The tail's decay rate comes from the service distribution, not from the simulation; only
        # the prefactors below are fitted, and only over the range the simulation resolves.
        self.theta = lundberg_exponent(latency_us / arrival_us) if self.stable else float("nan")
        self.sojourn_c = fit_prefactor(*tail_curve(self.sojourn_kept), self.theta, self.resolution)
        occ_support = np.arange(self.occ["ccdf"].size, dtype=np.float64)
        self.occ_c = fit_prefactor(occ_support, self.occ["ccdf"], self.theta, self.resolution)
        if not np.isfinite(self.occ_c) and np.isfinite(self.sojourn_c) and self.theta > 0:
            # A lightly loaded buffer holds at most three or four blocks even in a long run, which
            # leaves the occupancy tail too few resolved points to fit a prefactor to. It does not
            # need one: with one arrival per block period, the time the buffer spends deeper than `n`
            # is the sojourn tail integrated from `n` up — a block is the (n+1)th resident only while
            # the one that arrived `n` periods earlier is still being decoded —
            #
            #     P(resident > n) = E[(R - n)+] = (C / theta) exp(-theta n)
            #
            # so the sojourn fit already made carries over, divided by the rate. Checked against the
            # measured curve wherever both exist; the two agree to tens of percent over the resolved
            # decades, which is the accuracy this kind of extrapolation has anyway.
            self.occ_c = self.sojourn_c / self.theta

    # -- reported scalars -------------------------------------------------------------------------

    def miss_probability(self, deadline_cycles):
        """P(a block is still undecoded `deadline_cycles` block periods after it arrived)."""
        if not self.samples:
            return float("nan")
        return float(np.count_nonzero(self.sojourn_kept > deadline_cycles)) / self.samples

    def miss_extrapolated(self, deadline_cycles):
        """The same probability off the fitted exponential tail, for deadlines the run cannot reach.

        Returns NaN rather than a number when the queue is unstable or the fit had nothing to sit on:
        an extrapolation with no fit behind it is worse than a blank.
        """
        if math.isinf(self.theta):
            return 0.0
        if not (np.isfinite(self.theta) and np.isfinite(self.sojourn_c)):
            return float("nan")
        return min(1.0, float(self.sojourn_c * math.exp(-self.theta * deadline_cycles)))

    def deadline_for(self, target):
        """The sojourn budget, in block periods, at which the miss probability falls to `target`.

        The second return says whether the number was read off the fitted tail rather than off the
        simulation. Nothing that came from the fit is ever printed without that mark.
        """
        if math.isinf(self.theta):
            # No shot outruns its block period, so nothing ever waits and the worst case is exactly
            # the slowest decode in the log. Measured, not extrapolated.
            return float(self.latency_us.max() / self.arrival_us), False
        if not (np.isfinite(self.theta) and np.isfinite(self.sojourn_c)) or self.sojourn_c <= 0:
            return float("nan"), True
        return max(0.0, math.log(self.sojourn_c / target) / self.theta), True

    def depth_for(self, target):
        """Blocks of buffer needed for `P(resident > depth) <= target`.

        Answered from the simulated occupancy tail while that tail still has samples in it, and from
        the fitted exponential below the run's resolution — the targets that matter for a hard commit
        are several decades below anything a simulation of this length can show directly.
        """
        if not self.stable:
            # The occupancy CCDF of an unstable queue does eventually cross any target, at whatever
            # depth the backlog happened to reach before the run ended. That crossing is a property
            # of the run length and not of the decoder, and reporting it as a buffer depth would be
            # the most expensive mistake this script could make.
            return float("nan"), True
        ccdf = self.occ["ccdf"]
        # `ccdf[n]` is P(resident > n), so the first index that clears the target is the depth.
        below = np.nonzero(ccdf <= target)[0]
        if below.size and (ccdf[below[0]] > self.resolution or math.isinf(self.theta)):
            return float(below[0]), False
        if not (np.isfinite(self.theta) and np.isfinite(self.occ_c)) or self.occ_c <= 0:
            return float("nan"), True
        return max(0.0, math.log(self.occ_c / target) / self.theta), True


def draw_service(latency_us, cycles, rng, sampling, block_len):
    """`cycles` service times drawn from the logged shots.

    Three ways to turn a finite campaign into an arbitrarily long arrival stream, and they answer
    slightly different questions:

      * `iid` (default) treats the log as a distribution and forgets the order. Shot-to-shot
        correlation — a burst of hard shots, a cache that stays warm — is discarded, which is the
        conservative choice for a *distribution* question and the wrong one if bursts are the point.
      * `replay` keeps the logged order and tiles it if the simulation is longer than the campaign.
        Every correlation the profiler captured survives, including the ones that are artefacts of
        the machine it ran on.
      * `block` resamples contiguous runs of `--block-len` shots, which keeps correlation out to the
        block length and randomises beyond it.
    """
    n = latency_us.size
    if sampling == "replay":
        index = np.arange(cycles) % n
    elif sampling == "block":
        starts = rng.integers(0, n, (cycles + block_len - 1) // block_len)
        index = (starts[:, None] + np.arange(block_len)[None, :]).ravel()[:cycles] % n
    else:
        index = rng.integers(0, n, cycles)
    return latency_us[index]


def run_queue(service_us, arrival_us, servers):
    """FIFO queue, deterministic arrivals, `servers` decoders. Returns arrivals, wait, sojourn, end.

    One server is the Lindley recursion `W[i+1] = max(0, W[i] + S[i] - A)`, which is a running
    minimum of the partial sums of `S - A` and so runs vectorised in one pass — the difference
    between a second and several minutes at campaign length:

        W[i] = C[i] - min(C[0..i]),   C[i] = sum of (S - A) over the arrivals before i

    More than one server has no such closed form and is simulated the explicit way, against a heap
    of per-decoder free times. FIFO with `c` servers sends each block to the decoder that frees up
    first, which is what the heap's root is.
    """
    n = service_us.size
    arrivals = np.arange(n, dtype=np.float64) * arrival_us
    if servers == 1:
        partial = np.concatenate(([0.0], np.cumsum(service_us - arrival_us)))[:n]
        wait = partial - np.minimum.accumulate(partial)
    else:
        wait = np.empty(n, dtype=np.float64)
        # All-equal, so it is already a heap; no heapify needed before the first pop.
        free = [0.0] * servers
        for i in range(n):
            arrival = arrivals[i]
            earliest = heapq.heappop(free)
            start = arrival if arrival > earliest else earliest
            heapq.heappush(free, start + service_us[i])
            wait[i] = start - arrival
    sojourn = wait + service_us
    return arrivals, wait, sojourn, arrivals + sojourn


def occupancy_events(arrivals, finishes):
    """The occupancy step function as an event stream: `(time, blocks resident after the event)`.

    Exact, and cheaper than the obvious alternative. Sampling occupancy on a fine time grid — the way
    a strip chart wants it — costs a grid fine enough to resolve the shortest decode and still gets
    the time averages slightly wrong; the event stream gets them exactly, at two events per block.

    Ties are broken with the departure first. Latencies are integer nanoseconds and arrivals land on
    integer multiples of the block period, so a departure and an arrival landing on the same
    timestamp is not the measure-zero event it would be in a continuous model. Ordering the arrival
    first would report an occupancy of n+1 over an interval of zero length: invisible in every time
    average, and wrong in exactly one place — the reported maximum, which is the number a buffer gets
    sized against.
    """
    times = np.concatenate((arrivals, finishes))
    delta = np.concatenate((np.ones(arrivals.size, np.int64), -np.ones(finishes.size, np.int64)))
    order = np.lexsort((delta, times))
    return times[order], np.cumsum(delta[order])


def occupancy_stats(times, occupancy, window):
    """Time-weighted occupancy over `window`: the empty fraction, the mean, the max and the tail.

    Weighted by how long the buffer *held* each depth, not by how many events left it there. The two
    are not the same distribution and it is the first one that says how much of the time a buffer of
    a given depth would have been enough.
    """
    start, end = window
    lo = np.clip(times[:-1], start, end)
    hi = np.clip(times[1:], start, end)
    duration = hi - lo
    live = duration > 0
    depth = occupancy[:-1][live]
    duration = duration[live]
    total = float(duration.sum())
    if total <= 0:
        return {"empty": float("nan"), "mean": float("nan"), "max": 0, "ccdf": np.zeros(1), "drain": float("nan")}
    weight = np.bincount(depth, weights=duration) / total
    support = np.arange(weight.size, dtype=np.float64)
    # P(resident > n) at n = 0, 1, ...: one minus the CDF, clipped at zero so float error in the last
    # place cannot produce a negative probability on a log axis.
    ccdf = np.clip(1.0 - np.cumsum(weight), 0.0, 1.0)
    # A busy period is a maximal stretch with the buffer non-empty. Occupancy moves by one at a time,
    # so the busy stretches are exactly what separates the idle ones, and counting idle spells counts
    # them — a zero-length idle interval, from a departure and an arrival at the same instant, is
    # already dropped above and so does not split a busy period in two.
    idle_spells = int(np.count_nonzero(depth == 0))
    busy_time = total * float(1.0 - weight[0])
    return {
        "empty": float(weight[0]),
        "mean": float((support * weight).sum()),
        "max": int(weight.size - 1),
        "ccdf": ccdf,
        "drain": busy_time / idle_spells if idle_spells else float("nan"),
    }


def lundberg_exponent(service_cycles):
    """`theta > 0` solving `E[exp(theta (S - A))] = 1`, with the block period as the unit of time.

    The decay rate of the waiting-time tail of a GI/G/1 queue, taken straight from the empirical
    service distribution: no fitted lognormal, no assumed shape, and in particular no assumption that
    the escalation tail looks like the bulk. It is what lets the deadline-miss curve be quoted at
    probabilities a simulation of any affordable length cannot reach.

    `g(theta) = log E[exp(theta (S - 1))]` is convex, zero at the origin, and has slope `E[S] - 1`
    there, so a stable queue gives it a unique positive root; bracket by doubling and bisect. The
    expectation is taken in log space because `exp(theta (S - 1))` overflows long before the root
    does when a single escalated shot is tens of block periods long.
    """
    excess = service_cycles - 1.0
    if excess.mean() >= 0:
        return float("nan")
    if excess.max() <= 0:
        # No shot ever outruns its block period, so the buffer never holds anything and there is no
        # tail to decay. Reported as infinite rather than as a very large finite root.
        return float("inf")

    offset = math.log(excess.size)

    def g(theta):
        scaled = theta * excess
        peak = scaled.max()
        return peak + math.log(np.exp(scaled - peak).sum()) - offset

    high = 1e-3
    for _ in range(200):
        if g(high) > 0:
            break
        high *= 2.0
    else:
        return float("nan")
    low = 0.0
    for _ in range(200):
        mid = 0.5 * (low + high)
        if g(mid) > 0:
            high = mid
        else:
            low = mid
    return 0.5 * (low + high)


def tail_curve(values, points=700):
    """`(x, P(X > x))` sampled log-spaced from the top, so every decade of the tail gets equal detail.

    A linear subsample of a sorted series spends all its points in the bulk and draws the last three
    decades — the only part a deadline is set from — with a handful of them.
    """
    if values.size == 0:
        return np.zeros(0), np.zeros(0)
    ordered = np.sort(values)
    n = ordered.size
    rank = np.unique(np.round(np.logspace(0, math.log10(n), points)).astype(np.int64))
    rank = rank[rank < n]
    return ordered[n - rank - 1], rank / n


def fit_prefactor(x, tail, theta, resolution):
    """`C` in `C exp(-theta x)`, fitted where the simulated tail still has samples behind it.

    Only the prefactor is fitted; the rate is `theta`, which came from the service distribution. The
    fit range stops at the run's resolution floor — the last few points of an empirical CCDF are one
    or two observations each and would drag the fit around — and starts below the bulk, where the
    asymptotic form does not hold yet.
    """
    if not np.isfinite(theta) or theta <= 0 or tail.size == 0:
        return float("nan")
    inside = (tail > resolution) & (tail < 0.1) & (x > 0)
    if np.count_nonzero(inside) < 3:
        return float("nan")
    # Median of the per-point implied prefactor rather than a least-squares line: the range spans
    # decades of probability and one heavy point should not tilt it.
    return float(np.exp(np.median(np.log(tail[inside]) + theta * x[inside])))


def sweep_round_periods(series, args, rng):
    """Re-runs the queue across syndrome clocks, for the two panels that ask how much headroom is left.

    Two scalars are kept per point and nothing else. Mean occupancy comes from Little's law — with one
    arrival per block period, the mean number resident is exactly the mean sojourn measured in block
    periods — and the empty fraction of a single-server queue is exactly `1 - rho`. So neither panel
    needs an event stream, and a sweep point costs one Lindley pass.
    """
    periods = np.geomspace(args.round_ns * args.sweep_low, args.round_ns * args.sweep_high, args.sweep_points)
    rows = []
    for round_ns in periods:
        arrival_us = series.arrival_us * (round_ns / args.round_ns)
        service = draw_service(series.latency_us, args.sweep_cycles, rng, args.sampling, args.block_len)
        _, _, sojourn, _ = run_queue(service, arrival_us, args.servers)
        sojourn_cycles = sojourn[min(args.burn_in, sojourn.size - 1):] / arrival_us
        rho = series.mean_us / (arrival_us * args.servers)
        rows.append(
            {
                "round_ns": float(round_ns),
                "occupancy": float(sojourn_cycles.mean()),
                # `1 - rho` is exact for one decoder. With more than one it is only a bound, so the
                # panel draws nothing rather than a number nothing measured.
                "empty": max(0.0, 1.0 - rho) if args.servers == 1 else float("nan"),
            }
        )
    return rows


# -- figures ---------------------------------------------------------------------------------------


def style(ax, theme):
    ax.set_facecolor(theme["surface"])
    ax.grid(True, color=theme["grid"], linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(theme["grid"])
    ax.tick_params(colors=theme["text_secondary"], labelsize=9)
    ax.xaxis.label.set_color(theme["text_secondary"])
    ax.yaxis.label.set_color(theme["text_secondary"])


def panel_title(ax, text, theme, colour=None):
    ax.set_title(text, loc="left", fontsize=10.5, fontweight="bold", color=colour or theme["text_primary"], pad=8)


# Below this the extrapolation is arithmetic rather than a statement about a decoder. The
# exponential tail keeps going as far as it is asked to, and a miss probability of 1e-90 is not a
# claim anyone should have to evaluate: nothing in a fault-tolerant architecture is specified below
# the low tens of decades, and the asymptotic form was never fitted anywhere near there.
QUOTE_FLOOR = 1e-15


def unstable_note(ax, missing, consequence, theme):
    """Says which series a tail panel left out, and why. An absent curve is never left unexplained."""
    if not missing:
        return
    ax.text(0.5, 0.5, "\n".join(f"{label}: ρ ≥ 1, {consequence}" for label in missing),
            transform=ax.transAxes, ha="center", va="center", fontsize=9, color=theme["text_primary"])


def marked(value, fitted, shape):
    """A cell for the summary card: `~` if the number came off the fitted tail, `--` if there is none."""
    if not np.isfinite(value):
        return "--"
    return ("~" if fitted else "") + shape.format(value)


def late_cell(series, deadline):
    """P(late) for the card, which is usually zero in the run and never zero in reality.

    A deadline a few block periods out is missed far too rarely for a simulation of any affordable
    length to see it happen, and printing the `0` it observed would say the opposite of what the
    figure is for. When the run saw no late block the fitted tail answers instead, marked `~`.
    """
    observed = series.miss_probability(deadline)
    if observed > 0:
        return f"{observed:.0e}"
    extrapolated = series.miss_extrapolated(deadline)
    if np.isfinite(extrapolated) and extrapolated < QUOTE_FLOOR:
        return f"~<{QUOTE_FLOOR:.0e}"
    return marked(extrapolated, True, "{:.0e}")


def mark_clock(ax, args, theme):
    """The swept round period, on a log axis matplotlib would otherwise label with overlapping minors.

    Every decade of a two-and-a-half decade sweep gets a 1-2-5 tick and a plain integer label; the
    minor ticks keep their gridlines and lose their labels, which is what makes the axis readable at
    this width.
    """
    ax.axvline(args.round_ns, color=theme["text_secondary"], linewidth=1.0)
    ax.text(args.round_ns, 0.965, f" {args.round_ns:g} ns clock", transform=ax.get_xaxis_transform(), fontsize=8,
            color=theme["text_secondary"], ha="left", va="top")
    low = args.round_ns * args.sweep_low
    high = args.round_ns * args.sweep_high
    ax.set_xscale("log")
    ax.set_xlim(low, high)
    ticks = [value * 10.0 ** power for power in range(-2, 8) for value in (1, 2, 5)]
    ax.set_xticks([tick for tick in ticks if low <= tick <= high])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("syndrome round period  (ns)")


def strip_window(series_list, width, search=3000.0):
    """A time window for the occupancy strips, chosen once and used for both series.

    Centred on the busiest moment in the first `search` block periods after the burn-in, taken from a
    series that is actually stable. Two choices, both deliberate:

      * on a spike, because a window picked at random from a lightly loaded queue is a flat line at
        zero and says nothing about a buffer whose whole story is what the escalation tail does to it;
      * on a spike in a bounded early slice rather than the largest spike in the whole run, which
        would be a once-in-a-campaign excursion drawn as if it were what the buffer looks like.

    The same window goes to both strips — two strips on different windows are not a comparison.
    """
    reference = next((s for s in series_list if s.stable), series_list[-1])
    start, end = reference.window
    end = min(end, start + search * reference.arrival_us)
    times, occupancy = reference.times, reference.occupancy
    margin = width * reference.arrival_us
    inside = np.nonzero((times >= start + margin) & (times <= end - margin))[0]
    if inside.size:
        centre = times[inside[np.argmax(occupancy[inside])]] / reference.arrival_us
    else:
        centre = 0.5 * (start + end) / reference.arrival_us
    return centre - width / 2.0, centre + width / 2.0


def draw_overview(log, series_list, args, theme, deadline):
    """Panels (a)-(d): the two latency distributions and the two buffers they produce."""
    fig = plt.figure(figsize=(12.4, 9.2))
    fig.patch.set_facecolor(theme["surface"])
    # Three columns on the top row so the summary card has a column of its own rather than a corner
    # of a panel to sit in; the strips below span all three.
    grid = GridSpec(3, 3, figure=fig, height_ratios=[1.3, 0.9, 0.9], width_ratios=[1.25, 1.25, 0.62],
                    left=0.065, right=0.98, top=0.86, bottom=0.085, hspace=0.55, wspace=0.26)

    # (a) the measured latency densities, in block periods, which is the unit that decides whether
    # the buffer drains: a shot to the right of one block period leaves the queue deeper than it
    # found it.
    ax = fig.add_subplot(grid[0, 0])
    style(ax, theme)
    pooled = np.concatenate([s.latency_us / s.arrival_us for s in series_list])
    x_max = max(float(np.quantile(pooled, args.x_max_percentile / 100.0)),
                1.2 * max(s.mean_us / s.arrival_us for s in series_list))
    edges = np.linspace(0.0, max(x_max, 1.05), args.bins + 1)
    for series in series_list:
        counts, _ = np.histogram(series.latency_us / series.arrival_us, bins=edges, density=True)
        ax.stairs(counts, edges, color=series.colour, alpha=0.42, fill=True, zorder=2,
                  label=f"{series.label}  (ρ={series.rho:.2f})")
        ax.stairs(counts, edges, color=series.colour, linewidth=2.0, zorder=3)
    for series in series_list:
        ax.axvline(series.mean_us / series.arrival_us, color=series.colour, linewidth=1.4, linestyle=(0, (4, 3)),
                   zorder=4)
    ax.axvline(1.0, color=theme["text_secondary"], linewidth=1.0, zorder=4)
    # Up the line rather than across it: the label has to sit at x = 1 wherever that falls, and a
    # horizontal one lands on the legend at one grid point and on the mode at the next.
    ax.text(1.0, 0.02, " 1 block period", transform=ax.get_xaxis_transform(), va="bottom", ha="left",
            rotation=90, fontsize=7.5, color=theme["text_secondary"], zorder=5)
    ax.set_xlim(0, edges[-1])
    # Headroom for the legend, which sits over the top right corner and would otherwise land on the
    # taller of the two distributions.
    ax.set_ylim(0, ax.get_ylim()[1] * 1.22)
    ax.set_xlabel("decode time per block  (block periods)")
    ax.set_ylabel("density")
    panel_title(ax, "(a)  Measured per-block decode latency", theme)
    legend = ax.legend(frameon=False, loc="upper right", fontsize=8.4)
    for text in legend.get_texts():
        text.set_color(theme["text_primary"])

    # (b) the same two series as tails, on their own x range rather than (a)'s: an escalation tail one
    # shot in a thousand long is invisible in a density and is the whole reason the buffer ever fills,
    # so the panel that exists to show it is not cropped at the density's 99.5th percentile.
    ax = fig.add_subplot(grid[0, 1])
    style(ax, theme)
    tail_max = 1.05
    for series in series_list:
        x, tail = tail_curve(series.latency_us / series.arrival_us)
        ax.plot(x, tail, color=series.colour, linewidth=1.9)
        tail_max = max(tail_max, float(np.quantile(series.latency_us / series.arrival_us, 0.9999)))
    ax.set_yscale("log")
    ax.set_xlim(0, tail_max * 1.05)
    ax.set_ylim(0.5 / max(s.latency_us.size for s in series_list), 1.5)
    ax.axvline(1.0, color=theme["text_secondary"], linewidth=1.0)
    ax.set_xlabel("decode time per block  (block periods)")
    ax.set_ylabel("P(decode > x)")
    panel_title(ax, "(b)  Latency tail", theme)

    # (c), (d) the buffers themselves, on one window shared by both so the two strips compare.
    low, high = strip_window(series_list, args.strip_cycles)
    # Scaled to the deepest *stable* series in the window. An unstable one is tens of thousands of
    # blocks deep by this point in the run and would flatten every other strip to a line at the axis;
    # it runs off the top of the panel instead, and is annotated with where it actually is.
    ceiling = 1
    for series in series_list:
        cycles = series.times / series.arrival_us
        inside = np.nonzero((cycles >= low) & (cycles <= high))[0]
        if inside.size and (series.stable or all(not other.stable for other in series_list)):
            ceiling = max(ceiling, int(series.occupancy[inside].max()))
    ceiling = max(3, ceiling + 2)
    for row, series in enumerate(series_list):
        ax = fig.add_subplot(grid[1 + row, :])
        style(ax, theme)
        cycles = series.times / series.arrival_us
        first = max(int(np.searchsorted(cycles, low)) - 1, 0)
        last = min(int(np.searchsorted(cycles, high)) + 1, cycles.size)
        window_t, window_occ = cycles[first:last], series.occupancy[first:last]
        ax.fill_between(window_t, window_occ, step="post", color=series.colour, alpha=0.22)
        ax.step(window_t, window_occ, where="post", color=series.colour, linewidth=1.4)
        if window_occ.size and window_occ.min() > ceiling:
            # Off the top of a scale set by the series that stays bounded. Says where it actually is,
            # so the panel is not just a solid block of colour.
            ax.text(0.5, 0.5, f"{int(window_occ.min()):,}–{int(window_occ.max()):,} blocks deep here"
                              " and climbing — off this scale",
                    transform=ax.transAxes, ha="center", va="center", fontsize=9.5,
                    color=theme["text_primary"])
        else:
            ax.axhline(series.occ["mean"], color=series.colour, linestyle="--", linewidth=1.0)
            ax.annotate(f"mean {series.occ['mean']:.2f}", xy=(0.998, series.occ["mean"]),
                        xycoords=("axes fraction", "data"), xytext=(0, 3), textcoords="offset points",
                        color=theme["text_primary"], fontsize=8.4, ha="right", va="bottom")
        ax.set_xlim(low, high)
        ax.set_ylim(0, ceiling)
        ax.set_ylabel("blocks\nresident")
        if row == len(series_list) - 1:
            ax.set_xlabel("block period")
        if series.stable:
            headline = (f"{series.label} (ρ={series.rho:.2f}): empty {series.occ['empty'] * 100:.0f}% of the"
                        f" time, busy for {series.occ['drain']:.1f} block periods at a stretch")
        else:
            headline = (f"{series.label} (ρ={series.rho:.2f}): UNSTABLE — the backlog grows by"
                        f" {series.rho - 1:.2f} blocks per block period and never drains")
        panel_title(ax, f"({'cd'[row]})  Buffer occupancy — {headline}", theme, series.colour)

    card = fig.add_subplot(grid[0, 2])
    card.axis("off")
    tightest = min(args.targets)
    # An unstable series has occupancy numbers, and every one of them is a statement about how long
    # this simulation ran rather than about the decoder. They are left out of the card entirely; the
    # strip's own title says what happened instead.
    rows = [("", [s.key for s in series_list]),
            ("rho", [f"{s.rho:.2f}" for s in series_list]),
            ("empty %", [f"{100.0 * s.occ['empty']:.1f}" if s.stable else "--" for s in series_list]),
            ("mean resident", [f"{s.occ['mean']:.2f}" if s.stable else "--" for s in series_list]),
            ("max resident", [f"{s.occ['max']:d}" if s.stable else "--" for s in series_list]),
            ("busy period", [f"{s.occ['drain']:.1f}" if s.stable else "never" for s in series_list]),
            (f"P(late>{deadline:g}cyc)", [late_cell(s, deadline) for s in series_list]),
            (f"depth @{tightest:g}", [marked(*s.depth_for(tightest), "{:.0f}") for s in series_list]),
            (f"budget @{tightest:g}", [marked(*s.deadline_for(tightest), "{:.1f}") for s in series_list])]
    card.text(0.0, 1.0, "\n".join(f"{label:<14}" + "".join(f"{value:>9}" for value in values)
                                  for label, values in rows),
              transform=card.transAxes, ha="left", va="top", fontsize=7.6, family="monospace",
              color=theme["text_primary"],
              bbox=dict(boxstyle="round,pad=0.6", fc=theme["surface"], ec=theme["grid"], lw=0.9))
    card.text(0.0, 0.31, f"depth and budget are what\nit takes to hold overflow\nand lateness under"
                         f" {tightest:g}.\n~ marks a number read off\nthe fitted tail rather than\nseen in the run:"
                         " panels (e)\nand (f) draw both.",
              transform=card.transAxes, ha="left", va="top", fontsize=7.4, color=theme["text_secondary"])

    fig.suptitle(f"Syndrome buffer — {log.title}", x=0.065, ha="left", fontsize=14, fontweight="bold",
                 color=theme["text_primary"])
    fig.text(0.065, 0.895, subtitle(series_list, args), ha="left", fontsize=9.2, color=theme["text_secondary"])
    return fig


def draw_tails(log, series_list, sweeps, args, theme, deadline):
    """Panels (e)-(h): how deep the buffer has to be, how often it is late, and the clock headroom."""
    fig = plt.figure(figsize=(12.0, 8.4))
    fig.patch.set_facecolor(theme["surface"])
    grid = GridSpec(2, 2, figure=fig, left=0.07, right=0.98, top=0.855, bottom=0.145, hspace=0.45, wspace=0.2)

    # (e) the sizing curve. Read a target overflow probability off the y axis and the depth that meets
    # it off the x axis; the dashed continuation is the fitted exponential, which is the only thing
    # that can answer a target below the run's own resolution.
    # An unstable series is left out of both tail panels rather than drawn and ignored. Its occupancy
    # reached six figures and its sojourn five, so on a shared axis it is not one curve among two —
    # it is the axis, with the series that has an answer compressed into the first pixel. There is no
    # depth and no budget that suffices for it, which is what the note in each panel says.
    drawn = [series for series in series_list if series.stable]
    missing = [series.label for series in series_list if not series.stable]

    ax = fig.add_subplot(grid[0, 0])
    style(ax, theme)
    # Both panels stop their y axis a decade below the tightest target asked for, and both run their
    # dashed continuation exactly to that floor: the depth (or budget) at which the fit reaches the
    # bottom of the panel is the last one anything on the page is claiming.
    floor = min(args.targets) / 10.0
    depth_limit = 4
    for series in drawn:
        ccdf = series.occ["ccdf"]
        # Stop where the run stops resolving rather than drawing the CCDF's collapse to zero: a step
        # to the floor of a log axis is a statement that the buffer never goes deeper, and what the
        # run actually saw is that it did not go deeper *often enough to count* in this many blocks.
        resolved = int(np.count_nonzero(ccdf > series.resolution))
        ax.step(np.arange(resolved + 1), ccdf[:resolved + 1], where="post", color=series.colour,
                linewidth=1.9, label=series.label)
        depth_limit = max(depth_limit, resolved + 2)
        if np.isfinite(series.theta) and np.isfinite(series.occ_c):
            # Same step geometry as the measured curve, so the two are compared at like points rather
            # than across half a step.
            extended = np.arange(0.0, math.log(series.occ_c / floor) / series.theta + 1.0)
            ax.plot(extended, series.occ_c * np.exp(-series.theta * extended), color=series.colour,
                    linewidth=1.2, linestyle=(0, (3, 3)), alpha=0.85, drawstyle="steps-post")
            depth_limit = max(depth_limit, int(extended[-1]) if extended.size else 0)
    for target in args.targets:
        ax.axhline(target, color=theme["grid"], linewidth=0.8)
        ax.text(0.995, target, f"{target:g} ", transform=ax.get_yaxis_transform(), ha="right", va="bottom",
                fontsize=7.5, color=theme["text_secondary"])
    ax.set_yscale("log")
    ax.set_xlim(0, depth_limit)
    ax.set_ylim(floor, 1.5)
    ax.set_xlabel("buffer depth  (blocks resident)")
    ax.set_ylabel("P(resident > depth)")
    panel_title(ax, "(e)  Buffer sizing", theme)
    if drawn:
        legend = ax.legend(frameon=False, loc="lower left", fontsize=8.4)
        for text in legend.get_texts():
            text.set_color(theme["text_primary"])
    unstable_note(ax, missing, "no depth is enough", theme)

    # (f) the deadline-miss curve: the same picture for time rather than depth.
    ax = fig.add_subplot(grid[0, 1])
    style(ax, theme)
    span = 4.0
    for series in drawn:
        x, tail = tail_curve(series.sojourn_kept)
        # Same truncation as (e): the last handful of points of an empirical tail are one observation
        # each, and drawing them is drawing the run's noise as if it were the distribution's shape.
        resolved = tail >= series.resolution
        x, tail = x[resolved], tail[resolved]
        ax.plot(x, tail, color=series.colour, linewidth=1.9, label=series.label)
        if x.size:
            span = max(span, float(x.max()))
        if np.isfinite(series.theta) and np.isfinite(series.sojourn_c):
            reach = math.log(series.sojourn_c / floor) / series.theta
            extended = np.linspace(0.0, reach, 200)
            ax.plot(extended, np.minimum(1.0, series.sojourn_c * np.exp(-series.theta * extended)),
                    color=series.colour, linewidth=1.2, linestyle=(0, (3, 3)), alpha=0.85)
            span = max(span, reach)
    ax.axvline(deadline, color=theme["text_secondary"], linewidth=1.0)
    ax.text(deadline, 0.98, f" deadline {deadline:g}", transform=ax.get_xaxis_transform(), fontsize=8,
            color=theme["text_secondary"], ha="left", va="top")
    ax.set_yscale("log")
    ax.set_xlim(0, max(span * 1.08, deadline * 1.3))
    ax.set_ylim(floor, 1.5)
    ax.set_xlabel("latency budget from arrival to decoded  (block periods)")
    ax.set_ylabel("P(block still undecoded)")
    panel_title(ax, "(f)  Deadline miss", theme)
    unstable_note(ax, missing, "misses every deadline eventually", theme)

    # (g), (h) the sweep. The syndrome clock is the hardware's knob, so it is the x axis; the coloured
    # vertical is where each series stops keeping up, and the grey one is the clock assumed above.
    ax = fig.add_subplot(grid[1, 0])
    style(ax, theme)
    for series in series_list:
        rows = sweeps[series.key]
        ax.plot([row["round_ns"] for row in rows], [max(row["occupancy"], 1e-3) for row in rows],
                color=series.colour, linewidth=1.9, label=series.label)
        ax.axvline(series.critical_round_ns, color=series.colour, linewidth=1.1, linestyle=(0, (2, 2)))
    mark_clock(ax, args, theme)
    ax.set_yscale("log")
    # Past rho = 1 the mean is set by how long the simulation ran and by nothing else, so the axis
    # stops where the numbers stop meaning something and the curve leaves the top of the panel.
    ax.set_ylim(0.05, 200)
    ax.set_ylabel("mean blocks resident")
    panel_title(ax, "(g)  Clock headroom — coloured vertical is ρ = 1", theme)

    ax = fig.add_subplot(grid[1, 1])
    style(ax, theme)
    for series in series_list:
        rows = sweeps[series.key]
        ax.plot([row["round_ns"] for row in rows], [100.0 * row["empty"] for row in rows], color=series.colour,
                linewidth=1.9, label=series.label)
        ax.axvline(series.critical_round_ns, color=series.colour, linewidth=1.1, linestyle=(0, (2, 2)))
    mark_clock(ax, args, theme)
    ax.set_ylim(0, 100)
    ax.set_ylabel("time the buffer is empty  (%)")
    panel_title(ax, "(h)  Idle fraction", theme)
    if args.servers > 1:
        ax.text(0.5, 0.5, f"not drawn for {args.servers} decoders:\n1 − ρ is only a bound on the idle time"
                          " once\nmore than one of them can be idle at once",
                transform=ax.transAxes, ha="center", va="center", fontsize=9, color=theme["text_primary"])
    legend = ax.legend(frameon=False, loc="upper left", fontsize=8.4)
    for text in legend.get_texts():
        text.set_color(theme["text_primary"])

    fig.suptitle(f"Buffer sizing and deadline miss — {log.title}", x=0.075, ha="left", fontsize=14,
                 fontweight="bold", color=theme["text_primary"])
    fig.text(0.07, 0.895, subtitle(series_list, args), ha="left", fontsize=9.2, color=theme["text_secondary"])
    fig.text(0.012, 0.012,
             "Solid in (e) and (f): simulated, drawn only as far as this run resolves — a probability under"
             f" {series_list[0].resolution:.0e} is not something {args.cycles:,} blocks can show.\nDashed"
             " continuation: C·exp(−θx), with θ the Lundberg exponent of the measured service distribution and only"
             " C fitted, to the resolved part of the curve it continues.\nEvery depth and budget quoted below"
             f" {series_list[0].resolution:.0e} comes off those dashed lines. In (g) and (h) the coloured verticals"
             " are where each series reaches ρ = 1.",
             color=theme["text_secondary"], fontsize=8.2, va="bottom", linespacing=1.5)
    return fig


def subtitle(series_list, args):
    reference = series_list[0]
    return (f"{args.cycles:,} blocks simulated · one block = {reference.arrival_us:.2f} us"
            f" ({args.round_ns:g} ns x {reference.arrival_us * 1000.0 / args.round_ns:.0f} rounds)"
            f" · {args.servers} decoder{'s' if args.servers > 1 else ''}"
            f" · service drawn from the logged shots ({args.sampling})")


# -- text ------------------------------------------------------------------------------------------


def table_lines(results, args):
    """Every number the figures carry, as text, so the plots are never the only read."""
    lines = [
        f"{'d':>4} {'p':>7} {'T':>4} {'series':>7} {'shots':>8} {'mean us':>9} {'block us':>9} {'rho':>8}"
        f" {'empty%':>8} {'1-rho%':>8} {'mean occ':>10} {'max occ':>8} {'busy':>8} {'p99 lat':>10}"
        f" {'p99.99':>10} {'':>8}"
    ]
    for result in results:
        meta = result["log"].meta
        for series in result["series"]:
            late = series.sojourn_kept
            drain = series.occ["drain"]
            lines.append(
                f"{meta.get('d', '?'):>4} {meta.get('p', '?'):>7} {meta.get('T', '?'):>4} {series.key:>7}"
                f" {series.latency_us.size:>8,} {series.mean_us:>9.3f} {series.arrival_us:>9.3f}"
                f" {series.rho:>8.3f} {series.occ['empty'] * 100:>8.2f}"
                f" {max(0.0, 1.0 - series.rho) * 100:>8.2f} {series.occ['mean']:>10.2f} {series.occ['max']:>8d}"
                f" {(f'{drain:.2f}' if np.isfinite(drain) else 'never'):>8} {np.quantile(late, 0.99):>10.2f}"
                f" {np.quantile(late, 0.9999):>10.2f} {('' if series.stable else 'UNSTABLE'):>8}"
            )
    lines += [
        "",
        "  `block us` is one arrival period: `rounds` syndrome rounds at the assumed round period. `rho` is mean",
        "  decode over that period (over that period times the decoder count, with --servers). `empty%` is the",
        "  fraction of *time* the buffer holds nothing, measured; `1-rho%` is what it must equal for one decoder,",
        "  printed beside it as a check on the simulation. `busy` is the mean busy period in block periods, and",
        "  `p99`/`p99.99 lat` are sojourn percentiles — arrival to decoded — also in block periods. On an UNSTABLE",
        "  row the occupancy and latency columns are what this run reached before it ended and would be larger in",
        "  a longer one; only `mean us`, `block us` and `rho` mean anything there.",
        "",
        f"{'d':>4} {'p':>7} {'T':>4} {'series':>7} {'crit ns':>9} {'theta':>7}"
        f" {f'late@{args.deadline_cycles:g}':>11}"
        + "".join(f"{f'depth {t:g}':>13}" for t in args.targets)
        + "".join(f"{f'budget {t:g}':>13}" for t in args.targets),
    ]
    for result in results:
        meta = result["log"].meta
        for series in result["series"]:
            depths = []
            for target in args.targets:
                depth, fitted = series.depth_for(target)
                depths.append("--" if not np.isfinite(depth) else f"{math.ceil(depth):d}{'*' if fitted else ''}")
            budgets = []
            for target in args.targets:
                budget, fitted = series.deadline_for(target)
                budgets.append("--" if not np.isfinite(budget) else f"{budget:.1f}{'*' if fitted else ''}")
            late = late_cell(series, args.deadline_cycles).replace("~", "")
            if series.miss_probability(args.deadline_cycles) <= 0:
                late += "*"
            lines.append(
                f"{meta.get('d', '?'):>4} {meta.get('p', '?'):>7} {meta.get('T', '?'):>4} {series.key:>7}"
                f" {series.critical_round_ns:>9.1f} {(f'{series.theta:.3f}' if series.stable else '--'):>7}"
                f" {late:>11}" + "".join(f"{value:>13}" for value in depths + budgets)
            )
    resolution = min((series.resolution for result in results for series in result["series"]), default=0.0)
    lines += [
        "",
        "  `crit ns` is the syndrome round period at which rho reaches 1: below it this decoder stops keeping up.",
        "  `theta` is the Lundberg exponent of the measured service distribution, the decay rate of the queue's",
        f"  tail. `late@{args.deadline_cycles:g}` is the chance a block is still undecoded {args.deadline_cycles:g}"
        " block periods after it arrived. `depth E` is",
        "  the buffer depth, in blocks, at which the overflow probability falls to E, and `budget E` the",
        "  arrival-to-decoded budget, in block periods, at which the miss probability falls to E.",
        "",
        "  A `*` marks a value read off the fitted exponential tail rather than off the simulation, which is every",
        f"  value under this run's resolution of {resolution:.0e} — including, usually, the miss probability itself:",
        "  a deadline missed once in a million blocks is not an event a simulation of this length observes, and the",
        "  `0` it would otherwise report is the one number on the page that would be actively misleading.",
        f"  `<{QUOTE_FLOOR:.0e}` is where quoting the fit stops: it keeps going, but an exponential extrapolated"
        " thirty decades",
        "  past anything it was fitted to is arithmetic and not a statement about a decoder.",
        "  `--` is an unstable queue: no depth and no budget is enough, because the backlog is not bounded.",
    ]
    for result in results:
        meta = result["log"].meta
        for series in result["series"]:
            if not series.stable:
                lines.append(
                    f"  UNSTABLE: d={meta.get('d', '?')} p={meta.get('p', '?')} {series.key} has rho ="
                    f" {series.rho:.2f}; the backlog grows by {series.rho - 1:.2f} blocks per block period."
                )
    return lines


def write_tables(results, args):
    by_dir = {}
    for result in results:
        out_dir = args.out_dir or os.path.dirname(result["log"].path) or "."
        by_dir.setdefault(out_dir, []).append(result)
    written = []
    for out_dir, group in sorted(by_dir.items()):
        out_path = os.path.join(out_dir, args.table_name)
        header = [
            "Syndrome buffer simulation: a D/G/1 queue driven by the measured per-shot decode latencies.",
            f"{len(group)} log file(s); {args.cycles:,} blocks per series; arrivals every `rounds` x"
            f" {args.round_ns:g} ns; {args.servers} decoder(s); service drawn {args.sampling}.",
            "",
        ]
        with open(out_path, "w") as handle:
            handle.write("\n".join(line.rstrip() for line in header + table_lines(group, args)) + "\n")
        written.append(out_path)
    return written


# -- driver ------------------------------------------------------------------------------------------


def process(path, args, theme):
    log = load(path)
    stock, sparse, escalated, dropped = series_of(log, args.include_contaminated, args.include_excluded)
    if not stock.size:
        print(f"  {log.stem}: every shot was contaminated; nothing to simulate")
        return None

    # One block is the window the profiler decoded: `rounds` rounds of syndrome extraction, arriving
    # at the hardware's rate. Taken from the log so a campaign at another depth needs no flag; the
    # round period itself is not in the log and is the flag.
    rounds = float(log.meta.get("rounds") or log.meta.get("d") or 1)
    arrival_us = args.arrival_ns / 1000.0 if args.arrival_ns else rounds * args.round_ns / 1000.0

    rng = np.random.default_rng(args.seed)
    series_list = [
        Series("stock", STOCK_LABEL, theme["series"][0], stock, arrival_us, args, rng),
        Series("sparse", SPARSE_LABEL, theme["series"][1], sparse, arrival_us, args, rng),
    ]
    sweeps = {series.key: sweep_round_periods(series, args, rng) for series in series_list}

    out_dir = args.out_dir or os.path.dirname(log.path) or "."
    written = []
    for name, figure in (
        (f"buffer_{log.stem}", draw_overview(log, series_list, args, theme, args.deadline_cycles)),
        (f"buffer_tail_{log.stem}", draw_tails(log, series_list, sweeps, args, theme, args.deadline_cycles)),
    ):
        for extension in args.formats:
            out_path = os.path.join(out_dir, f"{name}.{extension}")
            figure.savefig(out_path, dpi=args.dpi, facecolor=theme["surface"])
            written.append(out_path)
        plt.close(figure)
    return {"log": log, "series": series_list, "escalated": escalated, "dropped": dropped, "written": written}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("paths", nargs="+", metavar="PATH",
                        help="log files, or directories to take every latency_*.csv from")
    parser.add_argument("--out-dir", help="where the figures go; default is beside each log file")
    parser.add_argument("--formats", default="png", help="comma-separated: png,pdf,svg")
    parser.add_argument("--dpi", type=int, default=160)
    parser.add_argument("--round-ns", type=float, default=1000.0, metavar="NS",
                        help="one syndrome-extraction round; a block of `rounds` of them arrives every"
                             " `rounds` x this. The one number here that is neither measured nor derived"
                             " (default 1000, the usual superconducting figure)")
    parser.add_argument("--arrival-ns", type=float, default=0.0, metavar="NS",
                        help="the block arrival period outright, overriding --round-ns x rounds; for a decoder fed"
                             " by something other than one d-round window per d rounds")
    parser.add_argument("--servers", type=int, default=1, metavar="K",
                        help="parallel decoders sharing the one FIFO buffer (default 1)")
    parser.add_argument("--cycles", type=int, default=400_000, metavar="N",
                        help="blocks to simulate; also the resolution floor of every simulated tail, since a"
                             " probability under 1/N cannot be observed in N draws (default 400000)")
    parser.add_argument("--burn-in", type=int, default=2000, metavar="N",
                        help="blocks dropped from the statistics at the start; the queue begins empty, which is"
                             " not a state a loaded one spends time in (default 2000)")
    parser.add_argument("--sampling", choices=("iid", "replay", "block"), default="iid",
                        help="how the logged shots become an arbitrarily long arrival stream: `iid` forgets the"
                             " order, `replay` keeps it, `block` keeps it out to --block-len (default iid)")
    parser.add_argument("--block-len", type=int, default=64, metavar="N", help="run length for --sampling block")
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument("--deadline-cycles", type=float, default=10.0, metavar="D",
                        help="the arrival-to-decoded budget the miss probability is quoted against, in block"
                             " periods. A placeholder: set it from the reaction time your architecture actually"
                             " has (default 10)")
    parser.add_argument("--targets", default="1e-3,1e-6,1e-9,1e-12",
                        help="overflow / miss probabilities the table sizes the buffer and the budget for")
    parser.add_argument("--sweep-points", type=int, default=21, help="round periods in panels (g) and (h)")
    parser.add_argument("--sweep-low", type=float, default=0.15, metavar="F",
                        help="lowest swept round period, as a multiple of --round-ns")
    parser.add_argument("--sweep-high", type=float, default=8.0, metavar="F",
                        help="highest swept round period, as a multiple of --round-ns")
    parser.add_argument("--sweep-cycles", type=int, default=80_000, metavar="N",
                        help="blocks per sweep point; the sweep plots means, which converge long before a tail"
                             " does (default 80000)")
    parser.add_argument("--strip-cycles", type=float, default=190.0, metavar="N",
                        help="width of the occupancy strips in panels (c) and (d), in block periods")
    parser.add_argument("--bins", type=int, default=140, help="bins in the density panel")
    parser.add_argument("--x-max-percentile", type=float, default=99.5, metavar="Q",
                        help="where the latency axis stops, as a percentile of both series pooled")
    parser.add_argument("--include-excluded", action="store_true",
                        help="charge intersect + H build + Mwpm(H) build to the sparsified series, i.e. reject"
                             " §M2's reading that they pipeline out")
    parser.add_argument("--include-contaminated", action="store_true",
                        help="keep shots the scheduler interfered with; dropped from both series by default")
    parser.add_argument("--theme", choices=sorted(THEMES), default="light")
    parser.add_argument("--table", action="store_true", help="print the summary table to stdout as well")
    parser.add_argument("--table-name", default="buffer_simulation.txt",
                        help="filename of the text table written beside the figures")
    args = parser.parse_args()
    args.formats = [item.strip() for item in args.formats.split(",") if item.strip()]
    args.targets = sorted((float(item) for item in args.targets.split(",") if item.strip()), reverse=True)
    if args.cycles <= args.burn_in:
        print(f"  --burn-in {args.burn_in:,} leaves nothing of --cycles {args.cycles:,}")
        return 1

    files = collect(args.paths)
    if not files:
        print("  nothing to simulate.")
        return 1
    if args.out_dir:
        os.makedirs(args.out_dir, exist_ok=True)

    theme = THEMES[args.theme]
    results = []
    for path in files:
        try:
            result = process(path, args, theme)
        except (OSError, ValueError, KeyError) as error:
            print(f"  skipping {path}: {error}")
            continue
        if result is None:
            continue
        results.append(result)
        for out_path in result["written"]:
            print(f"  wrote {out_path}")
        for series in result["series"]:
            if series.stable:
                state = (f"empty {series.occ['empty'] * 100:.1f}% of the time, mean {series.occ['mean']:.2f} /"
                         f" max {series.occ['max']} blocks resident")
            else:
                state = f"UNSTABLE — backlog grows by {series.rho - 1:.2f} blocks per block period"
            print(f"    {result['log'].stem} {series.key:>6}: rho {series.rho:.3f}, {state}")

    if not results:
        print("  nothing simulated.")
        return 1
    for out_path in write_tables(results, args):
        print(f"  wrote {out_path}")
    if args.table:
        print()
        print("\n".join(table_lines(results, args)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
