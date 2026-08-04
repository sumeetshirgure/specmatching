# M4 and M5 exit report — the end-to-end driver, both flavours

Covers `design/pyrematching_design.md` §M4 and §M5. Both milestones are short: with M3's compression
removed, M4 is glue and M5 is a lift through storage the ball tables already have. They are reported
together because their tests share a file and neither has a measurement campaign of its own — M4's
exit checkpoint reads `q` and the latency distribution off the M3 and M6 artifacts.

## What landed

| Design item | Where |
|---|---|
| §M4.1 `TwoPhaseConfig` / `TwoPhaseDecoder` | `two_phase/driver/two_phase_decoding.{h,cc}` |
| §M4.1 batch entry | `TwoPhaseDecoder::decode_batch`, accumulating into `TwoPhaseAggregateStats` |
| §M4.1 `phase1_on_ball_graph = false` oracle path | same file; compiles **no** ball tables at all |
| §M4.2 tests | `tests/two_phase/two_phase_driver.test.cc` |
| M5 edges flavour | `TwoPhaseDecoder::decode_to_edges` |
| M5 `need_edge_lift` retired | renamed to `BallParams::store_paths`; `BALL_ARTIFACT_VERSION` bumped to 2 |
| M5 tests | same file |

`escalate_to_exact` is gone from the config, as the design says: escalation is not optional and has
no failure mode to guard.

## The driver

Per shot, literally §M4.1's block:

```
Phase 1                          # §M2.5 on H, or M1's decode on G in the oracle configuration
if the timeline completed:       # the common case, q ~ 5e-6 .. 1.3e-4 at T = 2 (§M3 exit)
    emit Phase 1's committed observables and weight
else:
    escalate_to_stock(...)       # §M3.1: re-decode the whole shot on G, discard Phase 1
```

Three things are worth writing down because they are not in the design.

**There is exactly one `pm::Mwpm` on `G`.** The design says `TwoPhaseDecoder` owns it, and
`BallDecoder` already owned one for the negative-weight preamble, the boundary-node mask and the
oracle path. Building a second for escalation would have doubled the largest allocation in the
decoder for no reason: the three uses are at disjoint times within a shot. `g_mwpm()` names it, and
in the oracle configuration — where there is no `BallDecoder` — it is a direct member instead. X6
asserts that in the common case that instance is never touched at all: its region and node arenas
have handed nothing out and its clock is still zero.

**`T = infinity` is rejected at construction on the ball graph**, rather than being left to fail
invariant 6, which is what §M4.1 asks for. It routes through `phase1_on_ball_graph = false`, where
no ball tables are compiled — with a finite `R` they could not supply the edges an unbounded run
needs, so compiling them would be waste as well as wrong.

**More than 64 observables works on both flavours.** Above 64 the `obs_int` mask is unusable, so the
obs flavour switches to the match-edges datapath and reads the committed observables out of the ball
tables' observable *id* lists — the same lists the masks are built from, by a route that does not
care how many observables there are (§0). Every corpus in `data/` is a single-observable memory
experiment, so this has a test of its own on a 71-observable chain, run at two horizons so that both
the committed readout and the escalation readout are exercised.

## The edges flavour

Non-escalating shots: Phase 1's committed match edges, lifted to `G`'s edges through the ball
tables' stored canonical paths. Escalating shots: stock `decode_detection_events_to_edges` on the
full syndrome, with Phase 1's edges discarded exactly as its observables are. Output format is the
inherited one — detector id pairs, `-1` for the boundary — so downstream tooling is untouched.

`need_edge_lift` is retired **by rename, not by removal**, and this is a place where the design and
the code as landed disagreed. §M5 says the flag "gated the storage of outer edge-index paths for the
portal lift" and that the edges flavour "needs only the ball tables' own path storage, which
`ball_path_offsets` / `ball_path_edges` already provide" — but as landed those were the same thing:
`need_edge_lift` is what gates `ball_path_nodes`. What §M5 actually retires is the *portal* lift, so
the honest resolution is to keep the flag and stop it naming a mechanism that no longer exists.
`BallParams::store_paths` says what it does. `BALL_ARTIFACT_VERSION` went from 1 to 2 with it, so an
artifact compiled before the rename is rejected at load rather than silently reinterpreted, which is
the hard-error contract of §M2.2.

Cancellation is the one piece of the lift that needed thought. A correction is a *subset* of `G`'s
edges, so an edge crossed by two committed pairs' paths is crossed zero times by the correction.
Stock does this by flipping a marker bit on the search graph and dropping whatever ends up unflipped;
the ball front end has the whole multiset in hand, so it sorts and drops the even multiplicities,
which needs no search graph and is deterministic by construction. `G`'s own negative-weight edges go
into the multiset first, exactly as stock adds them before cancelling.

One consequence, documented on the binding rather than left to be discovered: **`weight` is the
matching weight — the sum over matched pairs — not the weight of the emitted edge set.** When two
paths share an edge, the edge cancels out of the correction but the pairs still cost what they cost.
This is the same property stock has, and it is why stock's `to_edges` reports no weight at all; here
the number is reported because it is the same integer the obs flavour returns, which the tests
assert on every shot of both branches.

Getting that consistent took a correction. The escalating branch first read the weight back off the
emitted edge set, which is the *other* number, and the oracle branch read
`HarvestResult::committed.weight` — which the match-edges flavour never fills, so it was reading
zero. Both now go through `weight_of_match_edges`, which is `Mwpm::extract_paths_from_match_edges`
over the matched pairs. The tests that would have caught it — an equality against stock's obs-flavour
weight on the escalating and the unbounded-horizon paths — were added with the fix, because the ones
that existed compared edge *sets* and never the weight.

The **oracle** front end cannot lift through ball tables, because it has none. It reuses stock's own
expansion and cancellation pass instead, which was factored out of `decode_detection_events_to_edges`
into `pm::expand_match_edges_to_edges` — the M1.2 `begin_timeline` precedent applied a second time.
Copying it would have put the cancellation in two places, and the cancellation is the part of that
function with no independent check on it. See `docs/upstream_edits.md`.

## Test results — all green

| Test | What it establishes |
|---|---|
| `UnboundedHorizonIsIdenticalToStock` | `T = infinity` through the driver on the oracle front end is bit-exact against stock, on the surface-code and negative-weight corpora. M1 test 1 restated one layer up |
| `UnboundedHorizonOnTheBallGraphIsRejected` | The unsupported combination throws at construction |
| **`ExactnessEndToEnd`** | **The deliverable.** Obs bytes and weight bit-exact against stock on **every** shot of three corpora — surface code, negative weights, toric code — escalating or not. This replaces the old "LER vs exact, gate on no catastrophic regression": the output is exact MWPM, so the gate is equality |
| `PerShotCertificate` | `weight_out >= dual_sum_at_truncation` on every non-escalating shot. Measured ratio: **exactly 1.0000, mean and max, over 79 shots** — the truncated dual is tight, which is what §M4.2 expects on a completed shot. On escalating shots the dual is zero and the certificate is skipped rather than computed against a discarded number |
| `BatchDeterminism` | `decode_batch` is reproducible, agrees with the per-shot entry point, and its accumulator reconciles with `summarize()` |
| `EdgesFlavorOuterSyndromeMatchesTheInput` | `H @ edges` equals the input syndrome, computed independently in the test, over the full corpus and not a sample. Both branches exercised |
| `EdgesFlavorObservablesMatchTheObsFlavor` | Observables read off the edge set agree with the obs flavour, and so do the weights. The two go through different readouts — id lists against stored paths — so this is a real check |
| `EdgesFlavorUnboundedHorizonIsIdenticalToStock` | Edge-flavour identity at `T = infinity` through the oracle front end |
| `EdgesFlavorEscalationIsIdenticalToStock` | Edge-flavour escalation identity against stock (§M3.3 X1, edges flavour) |
| `EdgesFlavorRequiresItsConfiguration` | The flavour refuses to run without `store_paths` or without the search graph, rather than emitting a half-lifted edge set |
| `WideObservableCount` | 71 observables, both readouts, bit-exact against stock |

## M4 exit checkpoint

- [x] Identity, exactness, certificate and determinism tests green.
- [x] `q` measured and reported at the operating point, with shot counts — `docs/two_phase_m3_exit.md`,
      200 000 shots per point, `q` from `5e-6` (one event) to `1.3e-4` at `T = 2`.
- [x] Empty-residual shots measurably cost no more than Phase 1 alone. Asserted structurally rather
      than by timing (X6): `escalation_ns` is zero and `G`'s arenas have handed nothing out, so the
      cost is not "small", it is absent.
- [x] Headline latency distribution vs stock produced, per `(d, p)`, with the escalation spike
      visible at p99.99 rather than averaged away — `benchmarks/two_phase/results/m6_exit_artifact.csv`,
      section `latency`. **The escalation spike is not the largest thing in the tail**: at
      `q <= 1.3e-4` the max is 1.5–2.2x p99.99 and at four of eight M3 grid points the worst shot is a
      heavy *non-escalating* one. The tail is set by the defect-count distribution at least as much
      as by the fallback, which is worth knowing before anyone sizes a buffer against `q`.

## M5 exit checkpoint

- [x] The lifted edge set is a valid correction, over the full corpus.
- [x] Observables derived from the edge set match the obs-flavour output.
- [x] Edge-flavour identity at `T = infinity`.
- [x] Edge-flavour escalation identity against stock.
- [x] `need_edge_lift` retired; no other consumer existed. `BallParams` changed, so
      `BALL_ARTIFACT_VERSION` was bumped and the change is recorded in `docs/upstream_edits.md`.
