# M1 exit report — truncated timeline + harvest

Covers `design/pyrematching_design.md` §M1 (horizon gating, `process_timeline_until_horizon`,
harvest, exposed root blossoms) plus §M5.2, which the M1 exit checkpoint pulls forward.

Everything new lives under `src/pyrematching/two_phase/`, `tests/two_phase/` and
`benchmarks/two_phase/`. The six vendored files touched are enumerated in
[`upstream_edits.md`](upstream_edits.md).

## What landed

| Design item | Where |
|---|---|
| M1.1 horizon gating at the event funnel | `sparse_blossom/tracker/queued_event_tracker.h`, `flooder/graph_flooder.{h,cc}` |
| M1.2 `process_timeline_until_horizon` + RAII sentinel guard | `two_phase/truncation/truncated_timeline.{h,cc}` |
| M1.2 shared timeline preamble (`pm::begin_timeline`) | `sparse_blossom/driver/mwpm_decoding.{h,cc}` |
| Unit conversion, one place only | `two_phase/truncation/horizon.h` |
| M1.3 harvest, both flavours | `two_phase/truncation/harvest.{h,cc}` |
| M1.4 exposed root blossom extraction | `two_phase/truncation/exposed_blossom.{h,cc}` |
| M5.2 per-shot profile + aggregate stats + `summarize` | `two_phase/perf/two_phase_profile.h` |
| M1 exit artifact harness | `benchmarks/two_phase/m1_exit_artifact.cc`, `summarize_m1_artifact.py` |

## Test coverage — design §M1.5

| # | Requirement | Test |
|---|---|---|
| 1 | `T = ∞` identity, bit-exact, both flavours | `TwoPhaseTruncation.InfiniteHorizonIsBitExact{ObsFlavor,MatchEdgesFlavor}` |
| 2 | Separation invariant `Y(u) == T`; `d(u,v) >= 2T`, `bcost(u) >= T` | `TwoPhaseHarvest.ResidualDefectsHaveDualSumExactlyEqualToTheHorizon`, `…AreMutuallySeparatedByTwiceTheHorizon` |
| 3 | Exposed root blossoms yield their base | `TwoPhaseExposedBlossom.*` (2 hand-built fixtures, 1 flavour cross-check, 1 corpus sweep) |
| 4 | `T = 0` commits nothing | `TwoPhaseTruncation.ZeroHorizonCommitsNothing` |
| 5 | Residual size monotonic in `T` | `TwoPhaseTruncation.ResidualSizeIsMonotonicInHorizon` |
| 6 | No throw, instance reusable | `TwoPhaseTruncation.TruncationNeverThrowsAndLeavesInstanceReusable`, `…HarvestLeavesTheQueueEmpty`, `TwoPhaseFuzz.InterleavedTruncatedAndExactDecodesAgreeWithFreshInstances` |
| 7 | Upper-bound sanity | `TwoPhaseHarvest.CommittedPlusExactResidualUpperBoundsTheExactOptimum` |
| 8 | Negative-weight DEMs | the negative-weight corpus is in tests 1 and 2 |
| 9 | Tight-edge invariant | `TwoPhaseHarvest.TreeCommittedPairsAreTight` |
| 10 | Stock-path perf gate | measured A/B, below — a benchmark, not a gtest |
| 11 | Truncation fuzz | `TwoPhaseFuzz.RandomHorizonsPreserveHarvestInvariants` |

Debug-build invariants 1–5 all have coverage. Invariant 2 ("nothing past the horizon is ever
enqueued") is discharged by `pm::horizon_gate_stats`, checked in
`TwoPhaseTruncation.NothingPastTheHorizonIsEverEnqueued` and `…SentinelHorizonRejectsNothing`.
Invariants 6–9 belong to M2–M4.

Two invariants beyond the design's list are asserted because they turned out to be cheap and sharp:

- **The dual sum splits exactly.** `committed.weight + Σ_u Y(u) == Σ_S y_S` over every region alive
  at truncation. Every region's dual either paid for a committed pair or belongs to a residual
  defect; nothing is double-counted or dropped. `TwoPhaseHarvest.DualSumSplitsExactly…`.
- **Feasibility.** `Σ_S y_S <= exact optimum` on every shot, which is what makes the per-shot
  certificate of M3.4 meaningful.

### Baselines

- Inherited C++ suite: **116 tests**, green. With the two-phase suites: **136 tests**, green.
- Inherited Python suite: **105 passed** (`pytest tests/ src/`). This is the M0 baseline referred
  to by every later milestone; `docs/upstream_provenance.md` does not exist in this tree, so it is
  recorded here instead.
- The gtest target is now built with `-UNDEBUG` so the debug-build invariants are actually live in
  it even when the tree is configured `Release`.

## Stock-path perf gate (test 10)

There is no `vendor-base` tag in this repository, so the comparison is against commit `93f4212`,
the tip of the vendored history before any M1 work. Both binaries were built `Release` from the
same sources otherwise, and run **interleaved** on the same machine — a first attempt at
sequential A/B showed a 10% "regression" that was entirely thermal drift, and disappeared under
interleaving.

`pyrematching_perf`, shots/s (higher is better), four interleaved rounds:

| Benchmark | base (93f4212) | with horizon gating |
|---|---|---|
| `Decode_surface_r11_d11_p1000` | 210, 200, 210, 210 k | 210, 200, 210, 220 k |
| `Decode_surface_r21_d21_p10000` | 220, 220, 220, 220 k | 230, 220, 230, 230 k |
| `Decode_surface_r21_d21_p1000` | 11–12 ms | 11 ms |
| `Decode_surface_r21_d21_p100` | 570–580 shots/s | 570–580 shots/s |
| `Decode_surface_r11_d11_p100` | 8.7–8.8 k | 8.6–8.8 k |

**Green.** The sentinel horizon costs nothing measurable, so the fallback of templating
`GraphFlooder` on a horizon policy is not needed.

## Deviations from the design, and why

1. **`run_until_next_mwpm_notification` gained a horizon stop.** The design says shrink events are
   exempt from the gate *and* that "the queue drains by itself since nothing beyond `T` was
   inserted". Those two cannot both hold: an exempt shrink event can be scheduled past `T`, and
   processing it would carry `queue.cur_time` — and therefore every surviving dual — past the
   horizon, breaking `Y(u) == T`. The exemption is kept as written and the timeline stops instead.
   Documented as a fourth row in `upstream_edits.md`.

2. **The residual defect is *a* member at the horizon, not provably a unique one.** The design says
   an exposed blossom's residual defect must be its base and "an arbitrary member will not do".
   That is right, and the implementation finds the base by descending on the maximal nested dual
   sum. What the corpus sweep showed is that the base is not always *unique*: on a degenerate graph
   like the surface code, a region can be matched and pulled back into a tree at the same instant,
   spending zero time frozen, so its dual also tracks the clock and reaches `T`. Every such member
   is a legitimate exposed defect — all the odd cycle's edges are tight, so removing any of them
   leaves a tight perfect matching on the rest, and both `Y(u) == T` and `d(u,v) >= 2T` follow from
   the dual value alone, not from which member was picked. Ties break to the lower detector index,
   for determinism. The tests assert the harvested defect is a member at the horizon, and the
   hand-built fixtures (where the base *is* unique) pin the exact index.

3. **`Mwpm::reset()` is not called after harvest.** The design says to call it. Doing so destroys
   both arenas, which frees the region pool and makes the next shot malloc every region again —
   measured at roughly 3× the cost of the truncated timeline itself at `d = 5`. Harvest already
   hands every region and tree node back to its arena, so `reset_for_next_shot` instead clears just
   what a truncated timeline leaves behind: the queue, and the queued-event bookkeeping of the
   nodes with events still in it. Debug builds verify the result is indistinguishable from a fresh
   instance, node by node.

## Exit artifact

`d ∈ {11,13,15,17,19,21}` rotated-memory-X, `d` rounds, `p ∈ {0.001, 0.003, 0.005}`, 500 shots per
configuration, `T` swept over `{0.25 … 8}` × the median discretised edge weight of the graph.

```
./build/two_phase_m1_artifact --distances 11,13,15,17,19,21 --error-rates 0.001,0.003,0.005 \
    --shots 500 --horizon-multiples 0.25,0.5,0.75,1.0,1.25,1.5,2.0,3.0,4.0,8.0 \
    --diameter-shots 10 --dijkstra-budget 1500 --csv benchmarks/two_phase/results/m1_exit_artifact.csv
python benchmarks/two_phase/summarize_m1_artifact.py --error-rate 0.003
```

Raw data: `benchmarks/two_phase/results/m1_exit_artifact.csv`. `T` is quoted in median edge
weights throughout; the surface code's discretised weights span a factor of ~1.3, so one unit is
essentially one lattice step. All times are microseconds per shot.

### Residual vs `T`, at `p = 0.003`

| `T`/edge | | d=11 | d=15 | d=21 | | | |
|---|---|---|---|---|---|---|---|
| | | `q` / residual | `q` / residual | `q` / residual | Phase-1 ÷ exact (d=21) | exposed-root-blossom shots (d=21) | |
| 0.25 | | 1.00 / 57.5 | 1.00 / 151.6 | 1.00 / 422.7 | 0.44 | 0% | |
| 0.50 | | 1.00 / 14.7 | 1.00 / 37.1 | 1.00 / 99.3 | 0.55 | 0% | |
| 0.75 | | 0.89 / 2.5 | 0.98 / 6.1 | 1.00 / 14.8 | 0.73 | 71% | |
| 1.00 | | 0.21 / 0.30 | 0.51 / 0.94 | 0.81 / 2.33 | 0.77 | 39% | |
| 1.50 | | 0.02 / 0.02 | 0.06 / 0.07 | 0.13 / 0.14 | 0.81 | 5% | |
| 2.00 | | 0.01 / 0.01 | 0.01 / 0.01 | 0.01 / 0.01 | 0.81 | 0.4% | |
| ≥3.00 | | 0 / 0 | 0 / 0 | 0 / 0 | 0.81 | 0% | |

### The three things the artifact says

**1. Locality holds.** At fixed `T`, the mean residual grows in proportion to the space-time
volume, so the residual *density* is flat in `d`: at `T = 1` edge and `p = 0.003` it is 0.0053,
0.0061, 0.0062, 0.0055, 0.0057, 0.0055 for `d = 11 … 21`. This is the property the whole design
needs, and it is the strongest positive result here. `max_S y_S` never exceeds `T` and the cluster
weight-diameter it lives in is 2.3–4.7 × `T`, so no dual outgrows its cluster.

**2. Sparse blossom has almost no timeline to truncate.** The work is front-loaded: at `d = 21`,
`p = 0.003`, cutting the timeline at `T = 0.25` edges saves 46% of Phase 1, at `T = 0.5` saves 32%,
at `T = 0.75` saves 10%, at `T = 1.0` saves 5%, and beyond `T = 1.5` nothing at all. A `d = 21`
shot's whole matching
resolves within about two lattice steps of *time*, because at these error rates almost every defect
pair is adjacent. There is no long tail of late events to cut.

**3. There is no "rare fallback" regime.** `q` goes from 1.0 to 0 between `T = 0.75` and `T = 1.5`
edges, and it moves the wrong way with `d`: at `T = 1` it is 0.21 at `d = 11` but 0.81 at `d = 21`,
so falling back becomes the common case exactly where the speedup is supposed to matter. The
fallback cost ratio is 1.0–1.3, i.e. truncated shots are not a heavy tail — they are the norm.
Empty-residual shots do exist and are cheap, but they are 19% of shots at `d = 21, T = 1`.

Tail numbers (d=21, p=0.003, µs): at `T = 0.5`, Phase 1 mean/p99 = 207 / 265 against exact
262 / 368, so 1.27× on the mean and 1.39× at p99 — with 99 residual defects still to solve. At
`T = 0.75`: 239 / 311 against 262 / 372, so 1.10× / 1.20× with 14.8 residual defects.

## Go / no-go read

**The `(T, c)` feasibility window looks empty as specified, and this needs an explicit decision
before M2 is written.**

The design's own M2 note says compression is only real when `rho` exceeds a few lattice
edge-weights, i.e. roughly `T >= 3c`. The M1 data says truncation only buys Phase-1 time at
`T <~ 1` edge weight. With `rho = T / c` and `c >= 2` (the smallest `c` that makes the `1/c` stub
error bound worth anything), `rho` would be well under one edge weight — every node is its own
portal, and the portal graph is the outer graph. Conversely at `T >= 3` edge weights, where a
portal net could compress, `q = 0`: the residual is always empty, Phase 2 never runs, and Phase 1
is a full decode.

Three ways forward, in the order I would rank them:

1. **Reframe the claim around problem-size reduction rather than fallback rarity.** The amortised
   mean is `C_phase1 + q * C_phase2` with `q ≈ 1`, not `q ≈ 0`. At `T = 0.5–0.75` edges, Phase 1
   costs 55–73% of an exact decode and hands Phase 2 a residual that is 3.5%–23% of the original
   defects with a constant density in `d`. That leaves a real budget — 24–55 µs/shot at `d = 21` —
   for Phase 2 to break even, and everything below that is speedup. This keeps the design intact
   but makes the M2 compression target harder: `C_phase2` has to be small *always*, not rarely.
2. **Check whether the residual needs the portal graph at all.** With 2–15 residual defects at
   `d = 21`, mutually separated by `>= 2T`, an exact blossom run on the residual alone may already
   be cheap — in which case M2 is an optimisation of an already-working M3 rather than a
   prerequisite, and the milestone order should change.
3. **Accept the worst-case framing.** The data does not support it well: the fallback cost ratio is
   ~1.0–1.3, so there is no heavy tail being tamed either.

M1 itself is done and green regardless of which is chosen: the truncated timeline, the harvest, and
the separation invariant are the foundation for all three.

