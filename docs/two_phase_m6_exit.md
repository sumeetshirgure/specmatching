# M6 exit report — profiling, bindings, benchmarks

Covers `design/pyrematching_design.md` §M6. Three deliverables: a per-shot profile and a campaign
accumulator that reconcile against each other, python bindings that make the decoder usable end to
end, and §M6.4's benchmark list.

## What landed

| Design item | Where |
|---|---|
| §M6.1 per-shot profile, reworked for escalation | `two_phase/perf/two_phase_profile.h`, `TwoPhaseProfile` |
| §M6.1 timer backends | `HiResTimer`, `steady_clock` by default, `rdtsc` under `-DPYREMATCHING_USE_RDTSC=ON` |
| §M6.2 aggregate stats and `summarize()` | same file, `TwoPhaseAggregateStats` / `TwoPhaseSummary` |
| §M6.3 bindings | `two_phase/driver/two_phase.pybind.{h,cc}`, registered in `pyrematching.pybind.cc` |
| §M6.3 python entry point | `pyrematching.two_phase_decoder`, `src/pyrematching/_two_phase.py` |
| §M6.4 benchmarks | `benchmarks/two_phase/m6_exit_artifact.cc`, plus `m3_exit_artifact.cc` for the escalation and streaming budgets |
| §M6.4 timer discipline | `PreemptionProbe`, `TwoPhaseConfig::detect_preemption` |
| Python tests | `tests/two_phase/two_phase_bindings_test.py` |

### The profile, and what came out of it

The portal fields went with M3's compression: `inject_ns`, `inner_blossom_ns`, `lift_ns`,
`portal_collisions`, `inner_detection_events` and `sum_phase2_ns` are gone, and `escalation_ns`
replaces them. `c_phase2` and `fallback_cost_ratio` became `c_escalation` and
`escalation_cost_ratio`. The landed M1 exit artifact's CSV column names were **not** renamed with
them: both fields are identically zero on that harness — M1 has no second phase to run — so renaming
the columns would have desynchronised the generator from its own recorded artifact for no
information gained. There is a comment saying so at the emit site.

Two fields are new rather than renamed. `stock_ns` times the stock decode *inside* an escalation, so
that §M3.2's `mean_stock_ns_on_escalated` is a measurement rather than a mean over all shots — see
`docs/two_phase_m3_exit.md`, where it turns out to be a 1.7–8x correction. `contaminated` marks
shots the scheduler interfered with.

### Timer discipline — this was a real bug, not a style note

`HiResTimer` used `std::chrono::high_resolution_clock`, which on libstdc++ is a *typedef for*
`system_clock`: not monotonic, and free to step backwards under NTP. §M6.4 says to name
`steady_clock` explicitly and it now does, in both the default and the `rdtsc` calibration path.

`PreemptionProbe` reads `getrusage(RUSAGE_THREAD)` once before and once after each shot and marks
the shot contaminated if either context-switch counter moved. Contaminated shots are excluded from
the percentile vectors and counted, and `contaminated_shot_rate` is reported beside every
percentile. On the M6 grid that rate runs **0% to 3.4%**, rising with `d` and `p` — at `d = 13,
p = 5e-3` the scheduler touches one shot in thirty. Without the exclusion, p999 at that point would
have been a measurement of the scheduler and nothing else. The probe is Linux-specific; elsewhere it
reports "not contaminated" and the rate reads 0, which the artifact says explicitly rather than
implying the machine was quiet.

This run was made on an ordinary desktop with no core pinning and no frequency-scaling lockout. The
design asks for both; the contamination rate is what says how much to discount the numbers, and it
is why it is printed next to them rather than in a footnote.

### Bindings

`pyrematching.two_phase_decoder(dem, T=...)` builds a decoder from a `stim.DetectorErrorModel`.
The DEM crosses into C++ as text, not as a `stim` pybind object: the extension links `libstim` but
does not depend on `stim`'s own bindings, and coupling them would mean the module could only be
imported alongside a matching `stim` build. The inherited `Matching.from_detector_error_model` round-
trips through a string for the same reason.

`decode_batch(shots, profile=True)` returns the per-shot profile as a dict of arrays row-aligned
with the shots — columns rather than a list of objects, because every consumer of it is a plot or a
percentile. `aggregate=True` keeps accumulating across calls, which is what a long campaign that
should not materialise per-shot arrays wants: leave `profile=False` and read `get_aggregate_stats()`.

Weights cross the boundary in DEM float units, divided by the graph's `normalising_constant`, which
is what the inherited `Matching.decode(return_weight=True)` reports. The python helper deliberately
does **not** name a `num_distinct_weights` default of its own: naming a different one from
`Matching`'s would silently give the two decoders different integer metrics.

The helper defaults `R` to `2 * T_max`, which is invariant 6 with no slack. That exposed a real
papercut: `to_time_units` rounds `R` and `T_max` independently, and `round(2x)` can be
`2 * round(x) - 1`, so the derived time-unit form of the invariant could fail by one unit at the
exact-equality setting. `compile_ball_tables` now widens `r_int` to `2 * t_max_int` instead of
throwing. Widening is the safe direction — a ball one unit too big costs memory, one unit too small
silently omits a reachable pair — and `BallParams::validate` still rejects `R < 2 * T_max` in weight
units, which is the check that actually protects the theorem.

## Measured results — do not re-derive

`benchmarks/two_phase/results/m6_exit_artifact.{csv,log}`. Grid: `d ∈ {5, 7, 9, 11, 13}` with
`rounds = d`, `p ∈ {5e-4, 1e-3, 3e-3, 5e-3}`, `T = 2` edge weights, `T_max = 2`, `R = 4`, 20 000
shots per point, `SCAN` build. Escalation and streaming budgets are in `docs/two_phase_m3_exit.md`
and are not repeated.

### 1. Latency against stock exact — the decoder is 2.4–3.6x slower

| `d` | `p` | mean (µs) | p99 | p999 | p99.99 | stock (µs) | speedup | contaminated |
|---|---|---|---|---|---|---|---|---|
| 5 | 5e-4 | 0.82 | 4.20 | 7.66 | 11.5 | 0.23 | 0.280 | 0% |
| 7 | 1e-3 | 3.91 | 13.1 | 21.2 | 36.1 | 1.18 | 0.302 | 0.010% |
| 9 | 1e-3 | 9.30 | 25.6 | 34.8 | 44.1 | 2.82 | 0.303 | 0.035% |
| 11 | 1e-3 | 18.6 | 41.0 | 56.0 | 71.6 | 6.12 | 0.329 | 0.375% |
| 13 | 1e-3 | 33.7 | 71.2 | 105 | 149 | 12.8 | 0.378 | 0.735% |
| 13 | 5e-3 | 250 | 480 | 686 | 925 | 95.8 | 0.383 | 3.38% |

`speedup_vs_stock` runs **0.280–0.416** across all twenty points: the two-phase decoder costs 2.4 to
3.6 times a stock exact decode on this machine. The ratio is worst — 0.28 to 0.30 — at low `p` for
`d <= 9`, and settles near 0.38 at `d = 13` roughly independently of `p`. The gap is therefore
*widest* at the low-`p` end, which is the opposite of what §M2's exit checkpoint predicted before its
read fired.

This is the M2 exit read, unchanged and now measured end to end. The `beta * V` floor that `H` was
built to remove is not present in M1 as landed, so re-hosting Phase 1 onto the defect manifold is not
cheaper as CPU throughput. **M3–M6 do not address that and were never going to.** What they deliver
is exactness on every shot, a fallback costing `<= 0.014%` of the mean at `T = 2`, and the
correctness machinery off the production path.

The governing metric from the M2 exit read onward is the critical path in dependent operations, not
wall time on one machine (`docs/two_phase_m29_exit.md`). The wall-time numbers here are reported
because the M6 exit checkpoint asks for them by name, and they are reported as measured.

### 2. Front-end comparison — three curves

Stock exact on `G`, M1's Phase 1 on `G`, M2's Phase 1 on `H`, all on the same shots at `T = 2`:

| `d` | `p` | stock (ns) | M1 on `G` (ns) | M2 on `H` (ns) | `speedup_vs_m1` |
|---|---|---|---|---|---|
| 5 | 5e-4 | 231 | 279 | 725 | 0.385 |
| 5 | 5e-3 | 2 505 | 2 658 | 6 168 | 0.431 |
| 7 | 5e-4 | 590 | 687 | 2 165 | 0.317 |
| 9 | 1e-3 | 2 818 | 2 988 | 9 477 | 0.315 |
| 11 | 5e-3 | 47 901 | 50 973 | 135 631 | 0.376 |
| 13 | 1e-3 | 12 751 | 13 323 | 34 440 | 0.387 |

`speedup_vs_m1` is **0.315–0.431** across the grid, against `0.48–0.76` at the M2 exit. The direction
is the same and the conclusion is unchanged; the difference is partly that this campaign sizes the
ball tables at `R = 2T` per point rather than at a fixed `R`, and partly that this is a different
machine-day. Two things worth reading off the table rather than the ratio:

- **M1's Phase 1 on `G` costs almost exactly what a full stock decode costs** — 0.99x to 1.21x
  across the grid, and above 1 at seventeen of the twenty points. Truncation does not save Phase-1
  time (M1 result 5), and at `T = 2` it saves essentially none of it.
- **The gap does not close as `p` falls.** §M2's exit checkpoint expected the win to grow at low `p`,
  because the `beta * V` floor dominates at low defect density. It does not, because that floor is
  not there.

### 3. Effective distance — identical predictions, not merely an indistinguishable fit

| `p` | `d` points with signal | slope, two-phase | slope, exact MWPM | shots where the two disagreed |
|---|---|---|---|---|
| 5e-4 | 1 | — | — | **0** |
| 1e-3 | 2 | — | — | **0** |
| 3e-3 | 2 | — | — | **0** |
| 5e-3 | 5 | `-0.2421 ± 0.0254` | `-0.2421 ± 0.0254` | **0** |

The slope is of `ln(LER)` against `d` at fixed `p`; `c2 = d_eff/d` is proportional to it, with
constant `ln(p / p_th)`, so the two decoders are compared on the slope rather than on an absolute
`c2` that would need a threshold estimate.

**The number to read is the last column.** The two decoders made the same prediction on all 400 000
shots of the campaign, which is exactly what §M4.2's bit-exactness test asserts and what §M3
guarantees by construction. The LER fit can only ever be a weaker, statistical restatement of that,
and the identical slopes are not independent evidence — they are the same numbers twice.

The fit also reproduces the design's own warning verbatim: at `p <= 3e-3` fewer than three distances
have a non-zero LER at 20 000 shots and the fit has no signal at all. §M6.4 says to run at
`p ∈ {3e-3, 5e-3}` as well as the operating grid, "which has now produced a vacuous Level 3 result
three campaigns running". Only `p = 5e-3` produced a fit here. Anyone wanting a Level 3 result at the
operating grid needs shot counts several orders larger, and should read the disagreement count
instead.

### 4. Ball table cost and the structural counters

`T_max = 2`, `R = 4` edge weights, `SCAN` build, at `p = 5e-3`:

| `d` | nodes | entries | table bytes | mean \|B\| | max \|B\| | isect bytes/defect | `hbld_edges_written` | `mwpm_init_elements` |
|---|---|---|---|---|---|---|---|---|
| 5 | 120 | 11 238 | 288 704 | 93.7 | 116 | 849 | 33.9 | 69.9 |
| 7 | 336 | 41 848 | 1 066 524 | 124.5 | 182 | 1 127 | 124.1 | 257.6 |
| 9 | 720 | 101 040 | 2 582 004 | 140.3 | 225 | 1 260 | 305.9 | 637.3 |
| 11 | 1 320 | 200 440 | 5 172 660 | 151.8 | 241 | 1 363 | 617.3 | 1 289.3 |
| 13 | 2 184 | 351 136 | 9 196 176 | 160.8 | 242 | 1 438 | 1 087.0 | 2 274.7 |

Table bytes grow linearly in the node count with a slowly rising constant — 2.4 kB per node at
`d = 5` against 4.2 kB at `d = 13` — because the mean ball size is still growing towards its
asymptote at these distances. **`isect_bytes_per_defect` is the headline number for the local-memory
hardware argument**: 0.85–1.4 kB per defect, in `SCAN` mode, which is what a per-defect processing
element would have to be able to stream. `BITSET` reads roughly four times fewer bytes per defect;
the row says which mode it is, and the artifact carries the other mode's derived count beside it.

The mean ball size saturating near 160 while `max |B|` saturates near 242 is the number to watch if
`R` is ever raised: both are `R³`-ish in the continuum and the memory is the binding constraint on
`T`, which §M3.2 identifies as the only worthwhile lever on `q`.

### 5. Heavy-shot drill-down

The top 0.1% of shots by defect count, plus every escalating shot, with the full per-stage
breakdown, are in the `heavy` section of the CSV — 1 005 shots on this grid. The pattern in them is
that a heavy shot is heavy in Phase 1 and not in harvest: at `d = 5` the worst shot in the sample
spends 32.8 µs in Phase 1 against 1.5 µs in harvest, and its stock reference is 11.6 µs. That is the
same story as the mean, more sharply.

## M6 exit checkpoint

- [x] Bindings usable from python end to end, including artifact save/load — six tests in
      `tests/two_phase/two_phase_bindings_test.py`, covering decode, batch, profile alignment,
      edges, `summarize()`, and a ball artifact round trip.
- [x] `summarize()` reconciles the amortised mean against the measured mean. Gap **0.24%–1.4%** on
      the M3 grid; the reconciliation is reported as a number (`amortisation_gap`) rather than left
      as a judgement call.
- [x] `d_eff/d` equals exact MWPM's, measured where the LER is non-zero — and, more strongly, the
      two decoders disagreed on **0 of 400 000 shots**.
- [x] Escalation and streaming budgets reported — `docs/two_phase_m3_exit.md`.
- [x] `contaminated_shot_rate` reported beside every percentile. It is 0%–3.4% on this machine,
      which is large enough to matter at p999.
- [x] Stock-path perf gate green — **run against `vendor-base` for the first time**, not asserted.
      The tag did not exist (the M2.9 report records the same gap); it was created on `66d3746`,
      which is unambiguously the vendoring commit, and the gate was run from a worktree at it. See
      `docs/upstream_provenance.md` for what could and could not be reconstructed about it.
      The 189-test C++ suite and the 111-test python suite are green.

### The perf gate, as actually run

`pyrematching_perf` built from `vendor-base` and from `HEAD`, four interleaved rounds each,
`--target_seconds 0.4`, medians over rounds:

| benchmark | `vendor-base` | `HEAD` | ratio | base's own round-to-round spread |
|---|---|---|---|---|
| `Decode_surface_r11_d11_p1000` | 195 kshots/s | 190 | 0.974 | 15.4% |
| `Decode_surface_r21_d21_p1000` | 18 kshots/s | 17 | 0.944 | 11.1% |
| `Decode_surface_r21_d21_p10000` | 200 kshots/s | 190 | 0.950 | 10.0% |

`HEAD` is 2.6–5.6% below `vendor-base`, and `vendor-base` is 10–15% below *itself* between rounds on
this machine. The difference is inside the noise floor, so the gate is green — but "green" here
means "not resolvable", not "identical". Two caveats a reader should have: the harness quantises its
throughput figures to two significant figures, which is itself ~5% at these magnitudes; and this is
an unpinned desktop, which the 0–3.4% contamination rate above independently attests to. A gate that
needs to resolve a few percent needs a pinned core and more rounds than four.

Every stock-path edit across M1 through M6 is enumerated in `docs/upstream_edits.md`. M6's is one
factoring — `expand_match_edges_to_edges` out of `decode_detection_events_to_edges` — which adds a
call and changes no work.

### A build note that costs an afternoon if unwritten

`cmake --build build` with no target fails, and has nothing to do with this milestone: it tries to
build **stim's own** `stim_pybind` target, and vendored stim v1.16.0 does not compile against the
pinned pybind11 v3.0.2 (`inconsistent types ... deduced for lambda return type` in
`tableau_simulator.pybind.cc`). Nothing in pyrematching links that target. Build the targets by name
instead:

```
cmake --build build --target pyrematching pyrematching_tests pyrematching_perf \
    libpyrematching _cpp_pyrematching \
    two_phase_m1_artifact two_phase_m2_artifact two_phase_m29_artifact \
    two_phase_m3_artifact two_phase_m6_artifact
```

## Decision

M6 lands, and with it the last milestone the design gates on. The decoder is exact MWPM on every
shot, usable from python, and instrumented in the terms the design asks for.

**What it is not is faster.** `speedup_vs_stock` is 0.28–0.42 and `speedup_vs_m1` is 0.32–0.43. The
M2 exit read already established why, and nothing in M3–M6 was aimed at it. The honest summary of
the project as it stands: the two-phase construction is *correct*, its correctness is cheap to
maintain (`<= 0.014%` amortised for exactness at `T = 2`), and its cost model on a CPU is worse than
stock by a factor of three. The case for it rests on the critical-path and local-memory arguments —
`isect_bytes_per_defect`, `hbld_edges_written`, `mwpm_init_elements`, and the harvest/solve depth
model of §M2.9.6 — and those are the numbers to carry forward, not these.

§M2.8 (the shell ladder) remains the one unimplemented section of the design. It is optional, it is
last, and nothing downstream depends on it.
