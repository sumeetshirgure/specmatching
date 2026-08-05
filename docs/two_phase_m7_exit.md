# M7 exit report — stock blossom on `H` with a max-dual certificate

Covers `design/stock_on_h_maxdual_design.md`. An alternative Phase-1 front end: run **stock
(untruncated) sparse blossom on the ball graph `H`** and decide exactness from a **terminal max-dual
certificate** rather than from a compiled-in horizon. The landed truncated-`H` path is untouched and
keeps its own M1 oracle; this is a front-end swap behind a flag, selected with
`TwoPhaseConfig::stock_on_h`.

Output is exact MWPM on every shot, measured rather than argued: **0 disagreements against stock
exact decode on `G`** over 1 500 000 shot-decodes.

## What this run measures, and what it no longer measures

`benchmarks/two_phase/results/m7_exit_artifact.csv`. Grid: `d ∈ {13, 17, 21, 25, 29}` with
`rounds = d`, `p ∈ {1e-4, 5e-4, 1e-3}`, `T = 2` edge weights, `T_max = 2`, `R = 4`, `SCAN` build,
**100 000 shots per point** in both sections — 1.5 M shot-decodes of identity checking and 1.5 M of
timing.

Two things changed in the harness since the milestone landed, and both narrow what the artifact
produces:

- **The A/B against the landed truncated-`H` path is gone.** The harness builds one front end. The
  read is now stock exact blossom on `G` against §M2's critical path of the *same stock solver* on
  `H` — both sides the same solver, one on a big graph and one on a small one, so the ratio is the
  ball graph's own win and nothing else.
- **The deep single-point escalation campaign is gone**, along with its `--q-*` flags. `q_this` now
  rides along on every `(d, p, T)` of the operating grid with `shots`; `q_current_on_same_corpus` is
  no longer produced.

Both of those had already been measured before the harness was cut down. Their results are kept
below under **Recorded results**, marked as such — the harness revision that produced them and the
CSVs it wrote are not in the tree, so this document is their only record. Nothing here re-states
them as if this run produced them.

The error-rate grid runs **low**, which is where the ball graph is meant to look best — fewer
defects means a smaller `H` and a shorter solve, against a stock decode whose cost is set by the
whole detector graph. It is also why §M7.6 level 3 is read as a disagreement count and not as an LER
fit: every LER in the campaign is zero (§M6.4's vacuous-fit warning). The distance grid runs **high**
— `d = 29` is 2.2× the largest distance any prior milestone measured — because the previous run's
headline was a scaling claim, and a scaling claim read off five points ending at `d = 13` is a claim
about the small-`d` regime. Extending it is the main new content here, and it turns the previous
conclusion over: see result 2.

## What landed

| Design item | Where |
|---|---|
| §M7.2 non-throwing completion entry point | `pm::CompletionStatus`, `pm::process_timeline_until_completion_or_report` in `sparse_blossom/driver/mwpm_decoding.{h,cc}` — the one vendored edit |
| §M7.0/§M7.2 the certificate | `two_phase/certificate/max_dual.{h,cc}` — `max_nested_dual`, `max_nested_dual_by_chain_walk`, `run_stock_and_certify` |
| §M7.1 `T` as a runtime scalar | `BallConfig::stock_on_h`, `TwoPhaseConfig::stock_on_h`; the flooder keeps `pm::NO_HORIZON` throughout |
| §M7.3 pipeline | the `stock_on_h` branch in `BallDecoder::decode_impl`; escalation and extraction are M3's and stock's, unchanged |
| §M7.5 invariants 1–6 | asserts in `run_stock_and_certify` and `decode_impl`; tests C1–C10 |
| §M7.6 tests | `tests/two_phase/stock_on_h_certificate.test.cc` (18 cases), `tests/two_phase/two_phase_bindings_test.py` (3 cases) |
| §M7.7 profile / `summarize` | `dual_scan_ns`, `certified`, `h_no_perfect_matching`, `max_dual_at_completion`; `q_this` |
| §M7.8 the read | `benchmarks/two_phase/m7_exit_artifact.cc`, `summarize_m7_artifact.py`, `results/m7_exit_artifact.csv` |

The vendored budget is unchanged: `mwpm_decoding.{h,cc}` was already on the edited list from M1.2
and M5, and `process_timeline_until_completion` is reimplemented on top of the reporting form rather
than copied, so there is still one timeline loop and one tree-survival test in the tree. The radix
heap and the timeline loop are untouched, and `process_timeline_until_completion`'s behaviour —
exception type and message included — is unchanged.

## Result 1 — the silent escalation trigger is 92% of escalations, and its share grows with defect count

§M7.0's central warning is that catching the no-perfect-matching throw is **not** sufficient, because
`H` can admit a *suboptimal* perfect matching whose terminal `max_u Y(u) > T`. Such a run completes,
throws nothing, and returning its matching is a silent wrong answer.

Over the 1 500 000-shot identity campaign, 270 shots escalated. 21 of them announced themselves:

| `d` | escalations | `H` had no perfect matching (announces itself) | `H` **completed** over `T` (silent) |
|---|---|---|---|
| 13 | 11 | 0 | 11 — 100% |
| 17 | 13 | 1 | 12 — 92.3% |
| 21 | 50 | 5 | 45 — 90.0% |
| 25 | 75 | 4 | 71 — 94.7% |
| 29 | 121 | 11 | 110 — 90.9% |
| **all** | **270** | **21 — 7.8%** | **249 — 92.2%** |

At `T = 2` on this grid, the throw catches **less than one escalation in twelve**. A decoder built on
"run stock on `H`, catch the throw, escalate on throw" would have returned a suboptimal matching on
249 of these 1.5 M shots and nothing anywhere in the system would have noticed.

The trend is the same one the previous campaign found against `p`, now found against `d`: the silent
trigger's share rises with defect density however the density is raised. The cheap, wrong scheme
looks least bad exactly where the decoder is least stressed. Debug invariant 4 and C1's hand-built
fixture exist for this case, and the fixture is worth reading: a star graph with one shortcut where
`H`'s only available perfect matching costs 16 against the true optimum's 12.

This is the result that justifies the milestone. The certificate is a certificate, not an exception
handler.

## Result 2 — the read: 1.44–2.69× on the critical path, and it **peaks at `d ≈ 25`**

The ratio §M2 is argued on, with both sides running stock blossom:

```
crit_ns = blossom_on_h_ns + dual_scan_ns + harvest_ns
speedup = stock_ns_on_G / crit_ns
```

`crit_ns` is what sits between a syndrome arriving and a correction leaving. Ball intersect, `H`
build and `Mwpm(H)` build are charged to nothing, on the assumption that they run ahead of the
syndrome they serve or overlap the previous shot's matching — the same assumption
`summarize_m2_artifact.py::critical_ns` makes. The certificate joins the critical path rather than
the pipelined-out set: the terminal scan gates whether the answer may be emitted at all, so it cannot
hide behind the next shot. Stock on `G` is timed on the same shot, in the same loop, immediately
after the two-phase decode.

| `d` | `p` | `stock_us` | `crit_us` | `e2e_us` | **`speedup`** | `e2e` |
|---|---|---|---|---|---|---|
| 13 | 1e-4 | 0.908 | 0.631 | 2.940 | **1.44** | 0.31 |
| 13 | 5e-4 | 4.829 | 3.108 | 14.357 | **1.55** | 0.34 |
| 13 | 1e-3 | 10.455 | 6.923 | 31.209 | **1.51** | 0.34 |
| 17 | 1e-4 | 2.672 | 1.326 | 6.541 | **2.02** | 0.41 |
| 17 | 5e-4 | 13.486 | 7.089 | 34.679 | **1.90** | 0.39 |
| 17 | 1e-3 | 29.238 | 15.641 | 76.364 | **1.87** | 0.38 |
| 21 | 1e-4 | 6.081 | 2.523 | 12.300 | **2.41** | 0.49 |
| 21 | 5e-4 | 32.452 | 13.874 | 67.825 | **2.34** | 0.48 |
| 21 | 1e-3 | 68.277 | 30.383 | 147.803 | **2.25** | 0.46 |
| 25 | 1e-4 | 12.814 | 4.767 | 21.330 | **2.69** | 0.60 |
| 25 | 5e-4 | 65.488 | 27.454 | 119.707 | **2.39** | 0.55 |
| 25 | 1e-3 | 140.687 | 63.822 | 265.101 | **2.20** | 0.53 |
| 29 | 1e-4 | 22.453 | 9.486 | 36.129 | **2.37** | 0.62 |
| 29 | 5e-4 | 119.903 | 57.514 | 206.126 | **2.08** | 0.58 |
| 29 | 1e-3 | 280.451 | 161.362 | 497.870 | **1.74** | 0.56 |

min 1.44, mean 2.05, max 2.69.

**The curve turns over, and that is the result.** The previous campaign, on `d ∈ {5..13}`, reported
this ratio rising monotonically with `d` and "still climbing at the largest distance measured", and
read that as the design's `d³`-volume argument showing up in the data. Extended to `d = 29` it does
not climb. It peaks at `d = 25` at `p = 1e-4` and at `d = 21` at `p = 1e-3`, then falls: from `d = 25`
to `d = 29` at `p = 1e-3`, `stock_us` grows 1.99× while `crit_us` grows 2.53×.

The scaling argument does not survive the extension, and it is worth saying why rather than filing it
as noise. The argument was that stock's cost is set by `G`, which grows as `d³` in spacetime volume,
while `H`'s solve is set by the defect count. But at **fixed** `p` the defect count is itself
proportional to the spacetime volume — both sides grow with `d³` together. There is no asymptotic
separation in this regime. What is left is a **constant-factor** advantage from matching on the
defect manifold instead of the detector graph, and this run measures that constant at **1.4–2.7×**,
best in the middle of the grid. The rise from `d = 5` to `d = 13` that the previous run saw was the
tail of the fixed-per-shot-cost regime, where building and solving a tiny `H` is dominated by costs
that do not scale, not the beginning of a growing win.

Two honest caveats, because this is still the most flattering number in the report:

- It is an **upper bound on the end-to-end win, not the win** — see result 3, which measures the
  same campaign at 0.31–0.62×.
- It charges `H` nothing for escalation. `escal_us = q_this · C_escalation` is reported beside it and
  is at most 0.336 µs per shot on this grid — 0.21% of `crit_us` at the worst point
  (`d = 29, p = 1e-3`). Small, not zero; see result 4.

## Result 3 — end to end it is 0.31–0.62×, and `intersect` is what decides that

The same campaign, undiscounted: `stock_us / e2e_us` over the whole two-phase shot runs **0.31 to
0.62**, mean 0.47 — the decoder costs 1.6 to 3.2 times a stock exact decode on this machine. This is
M6 result 1's measurement (0.28–0.42× on `d ∈ {5..13}`), re-taken on this front end and extended to
`d = 29`, and unlike the critical-path ratio it **does** rise monotonically with `d`: 0.31 at
`d = 13` to 0.62 at `d = 29`, with no sign of turning over. The two-phase shot grows more slowly than
stock's does even though the critical path inside it does not.

Where the shot goes, as a share of the whole two-phase decode:

| `d` | `p` | `isect%` | `hbuild%` | `mwpm%` | `blossom%` | `dscan%` | `hrvst%` | `charged%` |
|---|---|---|---|---|---|---|---|---|
| 13 | 1e-4 | 60.12 | 1.44 | 7.62 | 15.61 | 1.369 | 4.47 | 21.45 |
| 13 | 5e-4 | 66.15 | 1.86 | 7.47 | 17.86 | 0.612 | 3.18 | 21.65 |
| 13 | 1e-3 | 64.60 | 2.98 | 8.28 | 18.78 | 0.459 | 2.95 | 22.18 |
| 17 | 1e-4 | 67.28 | 1.02 | 6.53 | 15.77 | 0.881 | 3.62 | 20.27 |
| 17 | 5e-4 | 68.68 | 1.92 | 7.03 | 16.96 | 0.456 | 3.02 | 20.44 |
| 17 | 1e-3 | 66.36 | 3.07 | 8.66 | 17.38 | 0.377 | 2.72 | 20.48 |
| 21 | 1e-4 | 69.62 | 0.88 | 5.83 | 16.02 | 0.662 | 3.83 | 20.52 |
| 21 | 5e-4 | 68.99 | 1.99 | 6.97 | 16.98 | 0.421 | 3.06 | 20.46 |
| 21 | 1e-3 | 66.68 | 3.23 | 8.30 | 17.25 | 0.367 | 2.94 | 20.56 |
| 25 | 1e-4 | 69.03 | 0.88 | 5.33 | 17.10 | 0.537 | 4.71 | 22.35 |
| 25 | 5e-4 | 66.96 | 2.00 | 6.82 | 17.84 | 0.384 | 4.71 | 22.93 |
| 25 | 1e-3 | 63.08 | 3.18 | 8.58 | 19.25 | 0.346 | 4.48 | 24.07 |
| 29 | 1e-4 | 66.31 | 0.83 | 4.71 | 19.62 | 0.449 | 6.19 | 26.26 |
| 29 | 5e-4 | 62.75 | 1.87 | 6.40 | 21.15 | 0.356 | 6.40 | 27.90 |
| 29 | 1e-3 | 53.60 | 2.79 | 10.06 | 26.52 | 0.294 | 5.60 | 32.41 |

`charged%` — what the read in result 2 bills — is **20–32% of the shot**. The other 68–80% is
discounted as pipelined out, and **`intersect` alone is 54–70% of the whole decode at every single
point**. That is the entire distance between the two readings, and it is one stage.

This says where to push, and it says it more sharply than any previous milestone's split did.
`blossom_on_h` is a minority of the shot everywhere (15.6–26.5%), so making the solve faster cannot
move the end-to-end number much. Ball intersect can. It is also the stage whose pipelining assumption
is doing the most work: at `d = 21, p = 1e-4`, 69.6% of the shot is being assumed away, and if the
pipelining does not hold in a real system, essentially all of the result-2 gain comes back.

The `charged%` share climbing at `d = 29` (26–32%, against 20–22% in the middle of the grid) is the
same rollover as result 2 seen from the other side: `blossom_on_h` and `harvest` are taking a growing
share of the shot as the defect count grows, which is exactly what shrinks the critical-path ratio.

## Result 4 — the certificate's own cost, and the escalation charge

`dual_scan_ns` is the only cost the certificate itself adds. It is **0.29%–1.37% of the shot** and
*shrinking* as a fraction in both `d` and `p` — 1.37% at `d = 13, p = 1e-4`, 0.29% at
`d = 29, p = 1e-3`. It is one pass over the shot's defects against a solve that grows faster than
linearly in them, so it dilutes exactly where it would matter most.

Escalation, charged separately because it is a Phase-2 cost outside the Phase-1 critical path:

| `p` | `q_this` range over `d ∈ {13..29}` | `escal_us` | worst share of `crit_us` |
|---|---|---|---|
| 1e-4 | 0 at every `d` (below this campaign's `1/shots = 1e-5` floor) | 0.0000 | 0% |
| 5e-4 | 1.0e-5 – 1.6e-4 | 0.0005–0.0269 | 0.069% |
| 1e-3 | 1.0e-4 – 1.05e-3 | 0.0053–0.3359 | 0.208% |

`q_this` rises with `d` at fixed `p` — up to small-count wobble at `p = 5e-4`, where the `d = 25`
point saw 9 escalations against `d = 21`'s 12 — as it must: more defects means more chance that some
pair of them needs an edge `H` does not carry. It stays at or below 1.05e-3 across the grid, and the
amortised charge stays at or below 0.21% of the critical path and 0.07% of the end-to-end shot.
§M3.2 already capped what escalation can be worth at ~0.02% amortised at its own operating point;
nothing here disturbs that. Note how thin the counting is: five of the fifteen points saw **no**
escalation at all and sit on this campaign's `1/shots` floor, and seven more saw fewer than 30
events. The summarizer prints the event count and the `±√n/shots` band beside each `q_this` rather
than letting the figure pass as resolved.

## Result 5 — Level 1 pairing is a tie, not an identity

§M7.6 level 1 asks for the matched pair set to equal stock-on-`G`'s, "hard equality, no tolerance".
It does not, and the reason is the one §M2.6 level 1 already recorded: **the optimum is degenerate**.
When more than one minimum-weight perfect matching exists, `H` and `G` need not land on the same one
— in `G` a growing region's flood is blocked by its neighbours' territory, so `G` never observes some
tight collisions that `H`, where every pair within `2T` is a direct edge, does observe. Both are
correct blossom implementations of the same metric.

Measured over the level-1 corpora: 0–9 pairing ties per 40 shots, rising with defect density; the
toric corpus had none. What the certificate actually claims — that `H`'s matching is *a* global MWPM
— is checked and is **not** weakened: `H`'s weight is the sum of the ball tables' exact `d_G` values
over its own pairs, and it equals stock's optimum on every shot with no tolerance. A differing
pairing at equal weight is a second optimum; a differing pairing at a different weight would have
failed the weight assertion first. The test additionally asserts that the two pairings have the same
size and cover the same defects, which a tie cannot move.

**This is a deviation from the design as written**, and it is recorded here rather than silently
absorbed: level 1's weight half is a hard equality and stays one; level 1's pair half is a hard
equality *up to the choice among optima*, exactly as §M2.6 already treats it.

Level 2 came out cleaner than the design allowed for: **0 observable-mask divergences** across every
corpus and horizon. The licence to differ by a homologically trivial cycle was never exercised.

## Recorded results — measured, no longer produced by the artifact

Both of these were taken by the earlier revision of the harness — the one that built both front ends
and carried the `--q-*` flags — on `d ∈ {5..13}` at `T ∈ {1.5, 2.0}`. They stand as measurements;
they are not re-runnable from the current artifact, and no claim below is re-derived from the CSV in
the tree.

**The horizon gate is a net saving on the common path, not an overhead.** §M7.8 predicted the gain
over the landed truncated path would be "the gating overhead only", bounded above by M6 result 2 at
~20%. Measured with both front ends alternating **per shot** on the same corpus, the stock-on-`H`
front end came out **1–11% slower** (min 0.8933, mean 0.9536, max 0.9853) — inside the bound in
magnitude, opposite in sign. The stage split located it exactly: `intersect`, `h_build`, `mwpm_build`
and extraction were a noise floor at ratios of 0.97–1.03, while **`blossom_on_h` ran 1.04–1.56× the
truncated path's**. The M1 horizon gate lives in `QueuedEventTracker::set_desired_event` and rejects
an event *before* the queue insertion, so on every shot — not only on the ones that truncate — it
suppresses radix-heap pushes for collisions that would occur past `T`. Removing it gives that work
back, and the sparser the syndrome the larger that share is, because a lone defect floods further
before it meets anything.

This is the one surprise worth keeping in the design's memory. Any future proposal to remove gating
"because it is overhead" should be pointed at that measurement first.

**`q_this == q_current`, exactly, at 10⁶ shots.** §M7.0's corollary says the scheme escalates no more
often than the landed one. It escalated *exactly* as often, on a paired replay of the landed scheme's
decision on identical shots (`TwoPhaseConfig::measure_truncated_reference`):
`q_this = q_current = 7.343e-3` at `d = 13, p = 1e-3, T = 1.5` and `1.060e-4` at the same point with
`T = 2.0`, over 10⁶ shots each. Both points have plenty of events (7 343 and 106 escalations), so
these are measurements and not resolution floors. Not one shot in two million decodes differed in its
decision. The corollary's inequality is tight in practice, which is a stronger and less interesting
statement than the design allowed for — §M3.2 already capped the value of any reduction at ~0.02%
amortised, and the measured reduction is 0. `measure_truncated_reference` still exists and is still
exercised by `tests/two_phase/stock_on_h_certificate.test.cc`.

## Invariants and tests

All six §M7.5 invariants are asserted on every corpus shot in the debug (test) build.

| # | statement | where |
|---|---|---|
| 1 | no production horizon | `run_stock_and_certify` asserts the flooder sentinel on entry and exit; C5 additionally asserts `pm::horizon_gate_stats` is *never touched* across a 60-shot campaign, which says no gating code ran rather than that the field was restored |
| 2 | certificate soundness | `max_nested_dual` is cross-checked against an independent `blossom_parent` chain walk on every completing shot; output equality against stock-on-`G` is the level-1/level-3 identity count |
| 3 | escalation trigger exactness | asserted in `decode_impl` against the certificate's recomputed dual; swept over four horizons in C1 |
| 4 | no suboptimal `H`-matching escapes | asserted in `decode_impl`; C1's hand-built fixture is the positive case; result 1 is what it guards against, 249 times in this campaign alone |
| 5 | no surviving-tree artefacts on certified shots | `any_alternating_tree_survives` and the full sweep both asserted; `base_descents` asserted unchanged across a certified shot's extraction |
| 6 | output exactness | 0 disagreements over 1 500 000 shot-decodes here, plus the recorded 150 000 at `T ∈ {1.5, 2.0}` on `d ∈ {5..13}`, plus the C-series |

C1–C10 are all green, in the order the design's §M7.10 asks for (C3 first). The full suite — 207 C++
tests and 114 Python tests, the inherited PyMatching suite included — passes at its M0 baseline.

The stock-path perf gate is green, measured rather than assumed: `pyrematching_perf` run three times
on this tree and three times on `cf0cde0` with the change stashed gives identical throughput on
`Decode_surface_r11_d11_p1000` (230 kshots/s both), `Decode_surface_r21_d21_p1000` (20–21 kshots/s
both) and `Decode_surface_r11_d11_p100` (8.7–8.8 kshots/s both). Expected: the vendored edit is a
refactor of one function into two with no behavioural change, and the certificate is unreachable
from the stock path.

Two configurations are refused rather than silently reinterpreted:

- `stock_on_h` with `verify_against_g` / `full_harvest_for_verification`. The M1 oracle compares
  against a *truncated* harvest, and this path produces none. Its oracle is stock-on-`G`.
- `BallDecoder::decode_phase1` / `decode_phase1_to_match_edges` under `stock_on_h`. They return only
  a `HarvestResult`, and their callers read the escalation decision off `residual.empty()` — which is
  wrong in exactly the 92.2% case of result 1. Leaving a signature whose obvious use is a silent
  wrong answer was not worth the convenience.

## §M7.8's checkpoint, box by box

| box | status |
|---|---|
| Level 1 identity green against stock-on-`G`, full corpus, every `(d, p, T)`, no tolerance, certified and escalated | **green on weight**, over 1 500 000 shot-decodes at `T = 2` here and 150 000 at `T ∈ {1.5, 2.0}` on the recorded `d ∈ {5..13}` run. **Amended on pairs**: see result 5 — the optimum is degenerate, so the pairing is an identity up to the choice among optima, as §M2.6 already records |
| Level 2 divergences all homologically trivial; rate recorded | **green, vacuously**: `mask_divergence_rate = 0` on every corpus and horizon |
| Level 3 read as a per-shot disagreement count (expect 0) at `T ∈ {1.5, 2.0}`, `p ∈ {1e-3, 3e-3, 5e-3}` | **green: 0**, and `d_eff/d` is therefore exact MWPM's by construction rather than by fit. **Deviation**: this run covers `T = 2.0` and `p ≤ 1e-3` only; the `T = 1.5` and high-`p` halves are the recorded earlier campaigns, not this CSV |
| C1–C10 green; inherited suite at the M0 baseline; stock-path perf gate green | **green** — 207 C++ / 114 Python tests, perf gate measured against `cf0cde0` |
| Invariants 1–6 asserted on every corpus shot | **green** |
| `q_this` and `q_current_on_same_corpus` with `shots`, at `T ∈ {1.5, 2.0}`, 10⁶ shots | **green when taken: reduction exactly 0** — recorded above. **Deviation**: the paired 10⁶-shot campaign is no longer in the harness. `q_this` is reported with `shots` on every `(d, p, T)` of the operating grid instead — result 4 |
| **The read** — speedup vs the landed truncated-`H` path with the stage split; record whether the gain is inside the ~20% gating bound | **taken, and then retired.** Recorded above at 0.89–0.99 — a 1–11% slowdown, inside the bound in magnitude and opposite in sign. **Deviation**: the harness no longer races the two front ends. The read it takes now is §M2's `stock_ns / crit_ns` framing at **1.44–2.69×**, which turns over at `d ≈ 25` rather than climbing — result 2 — with the end-to-end number beside it at 0.31–0.62× — result 3 |

## How to reproduce

```
cmake --build build --target two_phase_m7_artifact

./build/two_phase_m7_artifact \
    --distances 13,17,21,25,29 --error-rates 0.0001,0.0005,0.001 --horizons 2.0 \
    --shots 100000 --identity-shots 100000 \
    --csv benchmarks/two_phase/results/m7_exit_artifact.csv

python benchmarks/two_phase/summarize_m7_artifact.py benchmarks/two_phase/results/m7_exit_artifact.csv
```

A bare `./build/two_phase_m7_artifact` runs the smaller default grid (`d ∈ {5..13}`,
`T ∈ {1.5, 2.0}`, 20 000 shots) instead; the flags above are what produced the CSV in the tree. The
binary exits non-zero if any identity disagreement is found, so it is usable as a gate.

The recorded results of the previous section are **not** reproducible from this tree. They need the
earlier revision of `m7_exit_artifact.cc` that built both front ends and carried the `--q-*` flags,
and neither it nor the CSVs it wrote survive in history.

This run was made on an ordinary desktop with no core pinning and no frequency-scaling lockout, like
the M6 run. Stock on `G` is timed on the same shot in the same loop, immediately after the two-phase
decode, in a fixed position rather than a rotation: a rotation is what an A/B between two candidate
front ends needs, and this is not one — the two sides are measured on different scales (a whole
decode against a subset of one), a fixed position is the same at every point in the grid, and an
earlier three-way rotation was measured to *move* the operating point rather than to centre it.

## Where this leaves the scheme

Keep it, as a flag, and do not make it the default.

- Its correctness story is genuinely better: the oracle is direct equality against stock-on-`G`
  rather than bit-exactness against a truncated harvest, and §M1.4 — the design's own highest-risk
  item — is off the path entirely rather than behind a bypass flag. Result 1 is the concrete payoff:
  92% of escalations on this grid are the trigger only a certificate can see.
- It costs 1–11% of Phase 1 against the landed path (recorded) and returns nothing measurable in
  escalation rate (recorded, reduction exactly 0). The certificate's own scan is 0.3–1.4% of the shot
  and shrinking with `d`.
- Its one surprise is worth keeping in the design's memory: **the horizon gate pays for itself on
  the common path**, by up to 1.56× of the solve at low `p`.

What this run says about the *design*, as against this front end, is the correction to carry forward.
The previous campaign read `stock_ns / crit_ns` rising monotonically to `d = 13` and took it as the
`d³`-volume argument appearing in the data. At `d = 29` the curve has turned over. At fixed `p` the
defect count scales with the spacetime volume too, so both sides grow together and there is no
asymptotic separation on offer — only a constant factor, measured at **1.4–2.7× on the critical path
and 0.31–0.62× end to end**. The end-to-end ratio, at least, is still improving with `d` and has not
turned over.

The place to push is not the solver. `intersect` is 54–70% of every shot in this campaign and
`blossom_on_h` never exceeds 27%, so the two-phase decoder's end-to-end cost is a ball-intersect cost
with a matching problem attached. Everything result 2 claims over result 3 rests on assuming
`intersect` pipelines away.
