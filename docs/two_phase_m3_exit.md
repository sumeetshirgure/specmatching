# M3 exit report — stock escalation, and the harvest bypass

Covers `design/pyrematching_design.md` §M3 and §M3.4. When Phase 1 leaves a non-empty residual the
partial result is discarded and the **whole shot** is re-decoded with stock exact decode on `G`.
The decoder's output is therefore exact MWPM on every shot, by construction: there is no
approximation budget anywhere in the pipeline, no `eps_total`, and no accuracy sweep.

## What landed

| Design item | Where |
|---|---|
| §M3.1 `escalate_to_stock`, obs and edges flavours | `two_phase/escalation/escalate.{h,cc}` |
| §M3.4 extract-only production path | `Harvester::extract_only_to_obs` / `…_to_match_edges`, `truncation/harvest.cc` |
| §M3.4 abandon-on-truncate | `two_phase::abandon_shot`, called from `BallDecoder::decode_phase1_production` |
| §M3.4 O(1) `COMPLETE`/`TRUNCATED` branch | `two_phase::any_alternating_tree_survives`, one `empty()` on the node arena's live vector |
| Full harvest retained as the oracle | `TwoPhaseConfig::full_harvest_for_verification`, and `BallDecoder::decode_phase1` unchanged |
| §M3.3 X1–X9 | `tests/two_phase/escalation_escalate.test.cc` |
| Exit artifact | `benchmarks/two_phase/m3_exit_artifact.cc`, results in `benchmarks/two_phase/results/` |

Nothing was deleted. `BallDecoder::decode_phase1` and `decode_phase1_to_match_edges` still run
§M1.3's harvest in full on every shot; they are what §M2.6 level 1 compares against M1, what X8
compares the bypass against, and what debug invariants 3, 4, 18 and 19 read. Deleting the machinery
would delete the oracle, which the design says twice and the risk table says a third time.

## The bypass, and the one thing the design did not say

§M3.4's production control flow is implemented literally:

```
process_timeline_until_horizon(H, dets, T)
  TRUNCATED  -> abandon the shot, then escalate_to_stock(...)
  COMPLETE   -> extract_only_to_obs(...)          # shatter + reduce
```

The design's table says jobs 2 and 3 — tree commits and M1.4's base descent — are dead on the
production path because they run exactly when the shot escalates and their output is discarded.
That is true, and both are now counted rather than argued: `HarvestCounters::tree_nodes_visited` and
`base_descents` are asserted identically zero on every production shot (X9). **The real return is
that M1.4 — the design's highest-risk item, with a silent failure mode — leaves the production
path**, exactly as §M3.4 says.

What the design does not mention is that skipping the harvest still leaves the `Mwpm` on `H` holding
live regions and tree nodes, and the next shot's preamble requires a clean instance. That is
`abandon_shot`, which is `Mwpm::reset` — the same teardown the stock driver runs when it cannot find
a perfect matching. It is the expensive kind of reset: it sweeps `H`'s node vector and frees the
arena pools, so the next shot re-allocates them. Two bounds make that the right trade rather than a
regression:

- it is paid on `q <= 7.3e-3` of shots (measured below, and `<= 1.3e-4` at the operating horizon);
- it is proportional to `H`'s node count, which is the shot's defect count, not to `G`'s.

Writing a cheaper "hand everything back without extracting it" would duplicate the extraction's
ownership rules, which is the one part of this code with a silent failure mode, and §M3.2 caps what
any escalation-path optimisation can be worth at ~0.02% amortised. So it was not written.

`dual_sum_at_truncation` is the one thing extract-only keeps that a strict reading would drop. It is
a flat reduction over the live-region walk extraction already performs — one integer read and one
add per region — and it is the §M4.2 per-shot certificate, which is the cheap independent check on
the primal that the escalation design removed every other one of. X8 asserts the two paths agree on
it exactly.

### The O(1) branch

`any_alternating_tree_survives` is one `empty()` on `node_arena.live`, the vector §M2.9.1 landed
inside `Arena<T>`. The design asked for the counter to be "either already available or nearly free —
reuse it rather than adding a second mechanism", and it was already there. Debug builds cross-check
it on every shot of every test against two slower answers to the same question: the free-list
comparison the timeline used to read, and a full sweep over live regions looking for an
`alt_tree_node`.

## Test results — X1 to X9, all green

| Test | What it establishes |
|---|---|
| X5 | At `T = 0.25` edge weights every shot with a defect escalates, and the decoder is bit-identical to stock across the corpus. Written first, per the implementation order |
| X1 | Bit-exact against stock on **every** shot at the operating horizon, escalating or not. 74 of 120 shots escalated on the chosen corpus, so the branch is genuinely exercised |
| X2 | `escalated == truncated == (residual_size > 0)`, asserted against `HarvestResult` rather than a recomputed predicate. Invariant 12 |
| X3 | Stock is handed the shot's **raw** syndrome. Run on the negative-weight corpus, where the raw and post-preamble sets differ, and the test additionally checks that the two inputs give *different* answers before concluding that agreeing with the raw one is evidence |
| X4 | A shot decoded after an escalating shot matches a fresh decoder, in both directions. Both `Mwpm` instances come out clean |
| X6 | At an operating point where nothing escalates, `escalation_ns` is zero and `G`'s region and node arenas have never handed anything out. Asserted structurally, not by timing: "cheap" and "absent" are different claims |
| X7 | Bit-reproducible across runs and across a save/load of the ball artifact |
| X8 | Extract-only and the full harvest agree bit for bit on every non-escalating shot, in both flavours, including `dual_sum_at_truncation` and the commit counts. **This is what licenses the bypass** |
| X9 | The production path performs no per-tree work: `tree_nodes_visited == 0` and `base_descents == 0` on every shot, and an escalating shot runs no harvest at all |

The inherited PyMatching suite is at its M0 baseline; the whole C++ suite is green at 189 tests
(169 before this milestone), and the python suite at 111 (105 before).

## Measured results — do not re-derive

Two campaigns.

**The sweep** — `benchmarks/two_phase/results/m3_exit_artifact.{csv,log}`. `d ∈ {11, 13}` with
`rounds = d`, `p ∈ {5e-4, 1e-3}`, `T ∈ {1.5, 2.0}` edge weights, **200 000 shots per point**, `SCAN`
build, `T_max = T` and `R = 2 T_max` per point. Preemption detection on; contaminated shots are
excluded from percentiles and counted. Everything below is from this campaign unless stated.

**The operating point** — `benchmarks/two_phase/results/m3_operating_point.{csv,log}`. The same
harness at `d ∈ {11, 13}`, `p = 1e-3`, `T = 2`, **10⁶ shots per point**, which is what the M3 exit
checkpoint asks for and what the sweep's 200 000 does not reach. Its `q` is quoted in §1 beside the
sweep's.

Running 10⁶ shots needed a change to the harness, not just a larger `--shots`: `generate`
materialised every shot, and at `d = 13` a shot's detection events cost around half a kilobyte, so
10⁶ of them is most of a gigabyte before the decoder allocates anything. `ShotSampler` samples and
decodes in chunks (`--chunk`, default 50 000), so the shot count is a question of time rather than
memory, and the stream depends on the seed alone. §M6.4's escalation budget asks for 10⁷ at
`p = 1e-4`; that is now a matter of leaving it running.

### 1. `q`, with its denominator

| `d` | `p` | `T` | shots | escalated | `q` |
|---|---|---|---|---|---|
| 11 | 5e-4 | 1.5 | 200 000 | 257 | 1.29e-3 |
| 11 | 5e-4 | 2.0 | 200 000 | 1 | 5.0e-6 |
| 11 | 1e-3 | 1.5 | 200 000 | 953 | 4.77e-3 |
| 11 | 1e-3 | 2.0 | 200 000 | 14 | 7.0e-5 |
| 13 | 5e-4 | 1.5 | 200 000 | 346 | 1.73e-3 |
| 13 | 5e-4 | 2.0 | 200 000 | **0** | **at or below the 5e-6 floor** |
| 13 | 1e-3 | 1.5 | 200 000 | 1466 | 7.33e-3 |
| 13 | 1e-3 | 2.0 | 200 000 | 26 | 1.30e-4 |

And the operating point at 10⁶ shots, where the checkpoint actually asks for it:

| `d` | `p` | `T` | shots | escalated | `q` | the sweep's `q` at 200 000 |
|---|---|---|---|---|---|---|
| 11 | 1e-3 | 2.0 | **1 000 000** | 85 | **8.5e-5** | 7.0e-5 (14 events) |
| 13 | 1e-3 | 2.0 | **1 000 000** | 112 | **1.12e-4** | 1.30e-4 (26 events) |

The two campaigns agree inside their Poisson bars, which is the point of running both: 85 and 112
events give `±11%` and `±9%`, where 14 and 26 give `±27%` and `±20%`. Quote the 10⁶ numbers.

The zero in the sweep is reported as a floor, not as a rate, which is the discipline §M3.0's risk
row demands. Three sweep points rest on fewer than 30 events and their error bars are `sqrt(k)/N`;
they are quoted with the count so a reader can compute that rather than trust the ratio.

`q` here is larger at `T = 1.5` than §M3.0's `T = 2` table by one to two orders, which is M1
result 2 restated: `lambda` falls 2.8–4x per quarter edge weight, and half an edge weight of horizon
is worth about one and a half orders of `q`. **The lever is `T`.**

### 2. `mean_stock_ns_on_escalated` — §M3.2's estimate replaced, and it was low by 2.8–5.5x

§M3.2 estimated the `stock` term of an escalating shot as the mean stock decode over *all* shots,
and flagged it as a lower bound because escalating shots are harder than average. The escalation
path measures the true figure for free. It is a much larger correction than "a lower bound" suggests:

| `d` | `p` | `T` | stock on escalated (µs) | stock over all shots (µs) | ratio |
|---|---|---|---|---|---|
| 11 | 5e-4 | 2.0 | 23.5 | 2.89 | 8.1x |
| 11 | 1e-3 | 2.0 | 34.7 | 6.30 | 5.5x |
| 13 | 1e-3 | 2.0 | 36.9 | 13.1 | 2.8x |
| 11 | 1e-3 | 1.5 | 13.8 | 5.80 | 2.4x |
| 13 | 1e-3 | 1.5 | 19.6 | 11.2 | 1.7x |

The 10⁶-shot campaign confirms the two `T = 2` rows on 85 and 112 escalating shots rather than 14
and 26: **5.16x** at `d = 11` (31.9 µs against 6.19 µs) and **3.05x** at `d = 13` (41.5 µs against
13.6 µs).

The ratio grows as `T` rises and as `q` falls, which is what it should do: at a higher horizon only
the very hardest shots survive to `T`, so the escalating population is a further-out slice of the
same tail. **Any future cost model for escalation must use this number and not a mean over all
shots.** The two sweep points at `T = 2, p = 5e-4` rest on 1 and 0 escalating shots respectively and
are not measurements; they are in the CSV with their counts.

### 3. Amortised cost, and the reconciliation

The amortised penalty `q * (C_escalation / C_phase1)` is **0.0013%–0.74%** across the sweep, and
**0.0013%–0.014% at `T = 2`**; the 10⁶-shot campaign puts the operating point at **0.0147%**
(`d = 11`) and **0.0133%** (`d = 13`). `summarize()`'s amortised mean reconciles with the measured
mean to within **0.24%–1.4%**; that gap is the unattributed cost the profile does not name, and
noticing it is what the reconciliation is for.

### 4. Latency with escalation live, and contamination

| `d` | `p` | `T` | mean (µs) | p99 | p999 | **p99.99** | max | contaminated |
|---|---|---|---|---|---|---|---|---|
| 11 | 1e-3 | 2.0 | 19.5 | 43.7 | 59.0 | 79.2 | 108.2 | 0.156% |
| 13 | 1e-3 | 2.0 | 35.3 | 75.6 | 106.4 | 143.0 | 222.6 | 0.367% |
| 13 | 5e-4 | 2.0 | 15.9 | 37.2 | 51.2 | 71.0 | 133.2 | 0.208% |
| 13 | 1e-3 | 1.5 | 19.8 | 44.6 | 67.8 | 92.8 | 163.4 | 0.436% |

And at the operating point over 10⁶ shots, where the tail has five times as many chances to show
itself:

| `d` | `p` | `T` | mean (µs) | p99 | p999 | **p99.99** | max | contaminated |
|---|---|---|---|---|---|---|---|---|
| 11 | 1e-3 | 2.0 | 18.6 | 42.6 | 59.9 | 86.4 | **209.7** | 0.183% |
| 13 | 1e-3 | 2.0 | 35.2 | 76.9 | 113.3 | 161.5 | **320.5** | 0.338% |

The design predicted that at `q ~ 3e-4` the escalation spike sits near p99.97 and that p99 and p999
would not show it. At the horizons here `q` is `8.5e-5` to `1.1e-4`, so the spike sits nearer p99.99
still — and the max is 2.0–2.4x p99.99, which is where an escalating shot actually shows up.
**Quoting p99 alone would have hidden the entire fallback.** Note also that the max grew by ~1.5x
between 200 000 and 10⁶ shots at both points: this tail is not converged at either shot count, and a
budget taken from a short campaign is a budget taken from a truncated distribution.

`contaminated_shot_rate` runs 0.006%–0.44%. At 0.44% the scheduler touches roughly one shot in 230,
which is *four times* the p999 tail: without excluding those shots, p999 would have been a
measurement of the scheduler and not of the decoder. This is not a quiet machine — no core pinning,
no frequency-scaling lockout — and the rate is reported so that a reader can discount the
percentiles accordingly rather than being invited to trust them.

### 5. Streaming budget

At a 1 µs syndrome-extraction round and `rounds = d`, a shot's budget is `d` µs. The worst observed
shot overruns it by 39–210 µs, so the input buffer that absorbs it is **4–17 shots deep**:

| `d` | `p` | `T` | shot budget (µs) | worst shot (µs) | overrun (µs) | buffer (shots) |
|---|---|---|---|---|---|---|
| 11 | 5e-4 | 1.5 | 11 | 50.0 | 39.0 | 4 |
| 11 | 1e-3 | 2.0 | 11 | 108.2 | 97.2 | 9 |
| 13 | 5e-4 | 2.0 | 13 | 133.2 | 120.2 | 10 |
| 13 | 1e-3 | 2.0 | 13 | 222.6 | 209.6 | 17 |

At the operating point over 10⁶ shots the buffer is deeper, because the tail is: **19 shots** at
`d = 11` (worst shot 209.7 µs) and **24 shots** at `d = 13` (320.5 µs).

Three caveats belong next to this table. The worst shot is not always an escalating one — at four of
the eight sweep points the maximum is a heavy non-escalating shot, though at both 10⁶-shot points it
is the escalating one — so this budget is set by the defect-count tail as well as by the fallback,
and which of the two dominates depends on how long you look. The maximum grew by half again between
200 000 and 10⁶ shots, so neither number is a converged worst case. And `q` is extensive in
spacetime volume (M1 result 3), so a deployment must re-run this at its intended `(d, rounds)`;
scaling this table is not valid.

### 6. The decoder is slower than stock on this machine, and this milestone does not change that

`C_phase1` runs 5.4–35.2 µs against a stock exact decode of 2.5–13.1 µs on the same shots, i.e.
**`speedup_vs_stock` is 0.29–0.46**. That is the M2 exit read, unchanged: the `beta * V` floor that
`H` was built to remove is not present in M1 as landed, so re-hosting Phase 1 onto the defect
manifold is not cheaper as CPU throughput. M3 does not address it and was never going to — it makes
the decoder *exact*, and it costs `<= 0.014%` of the mean at the operating horizon to do so.

The governing metric from the M2 exit read onward is the critical path in dependent operations, not
wall time on one machine (see `docs/two_phase_m29_exit.md`). The wall-time numbers above are
reported because the M3 and M6 exit checkpoints ask for them by name, and they are reported as
measured.

## Decision

M3 lands. The decoder is exact MWPM on every shot, the fallback costs at most 0.014% of the mean at
`T = 2`, and M1.4 is off the production path. §M3.4's bypass is licensed by X8 and its branch by X9.

**Known-bad, recorded so they are not re-proposed.** Escalating on the residual rather than the full
syndrome — a correctness bug, caught by X3 and invariant 14, not an optimisation. Deleting the full
harvest once the bypass landed — it is the oracle for §M2.6 level 1 and X8. Optimising the
escalation path: §M3.2 caps it at ~0.02% amortised, and measurement 2 above makes the case stronger
rather than weaker, since the escalating population is harder than the estimate assumed and the only
lever that removes them is a larger `T`.
