# M2 exit report — ball graph front end

Covers `design/pyrematching_design.md` §M2: run M1's truncated sparse blossom at horizon `T` on a
per-shot graph `H` over the shot's defects rather than on the detector graph `G`.

Everything new lives under `src/pyrematching/two_phase/manifold/`, `tests/two_phase/` and
`benchmarks/two_phase/`. **M2 required no new edits to any vendored file** —
[`upstream_edits.md`](upstream_edits.md) is unchanged from M1.

## What landed

| Design item | Where |
|---|---|
| M2.2 `BallParams`, ball tables, CSR pools, shells, bitset view, `BallStats` | `manifold/ball_params.h`, `manifold/ball_tables.{h,cc}` |
| M2.2 artifact save/load with graph hash + params | `manifold/ball_serialize.{h,cc}` |
| M2.3 per-shot `H`, both `SCAN` and `BITSET`, reusable arena | `manifold/ball_graph.{h,cc}` |
| M2.4 `pm::Mwpm` on `H` from pre-discretised integer weights | `manifold/ball_mwpm.{h,cc}` |
| M2.5 `BallDecoder`, `verify_against_g` | `manifold/ball_decoding.{h,cc}` |
| M2.7 `BallProfile`, `BallAggregateStats`, `summarize_ball` | `perf/ball_profile.h` |
| M2 exit artifact | `benchmarks/two_phase/m2_exit_artifact.cc`, `summarize_m2_artifact.py` |

M2.8 (the shell ladder) was **not** implemented. The design gates it on the M2.7 profile showing
`h_edges` to be the bottleneck; it is not — see the profile split below — and the design says
explicitly not to implement it first. `shell_width` and the shell offsets exist and are tested (B2),
so the ladder can be added later without touching the tables.

## Test coverage — design §M2.6

| # | Requirement | Test |
|---|---|---|
| B1 | `w_int` and `bcost` against an independent bounded Dijkstra | `ManifoldBallTables.B1TableCorrectness` |
| B2 | Symmetry, `(w_int, id)` ordering, shells partition | `ManifoldBallTables.B2SymmetryAndShells` |
| B3 | Every mask is the XOR of its stored path; every path is a valid `G` path of the stored weight | `ManifoldBallTables.B3MaskConsistency` |
| B4 | No pair M1 matched is beyond `2T`; no boundary match beyond `T` | `ManifoldBallIdentity.B4Coverage` |
| B5 | `SCAN` and `BITSET` build identical `H` | `ManifoldBallGraph.B5ScanMatchesBitset` |
| B6 | Determinism (incl. across thread counts) and artifact round-trip | `ManifoldBallTables.B6DeterminismAndRoundTrip` |
| B7 | Zero steady-state allocations; warmed decoder == fresh decoder | `ManifoldBallGraph.B7ArenaReuse` |
| B8 | `n = 0`; `T = 0`; mutually distant defects | `ManifoldBallGraph.B8*` (3 tests) |
| B9 | Fuzz: random `T`, syndromes, DEMs, both build modes | `ManifoldBallIdentity.B9Fuzz` |
| B10 | Negative-weight DEMs with the preamble exercised | `ManifoldBallIdentity.B10NegativeWeightDems` |
| Level 1 | Metric identity against M1 | `ManifoldBallIdentity.Level1AcrossCorpora`, `…ResidualIsIdenticalAtTheOperatingHorizon` |
| Level 2 | Observable equivalence | `ManifoldBallIdentity.Level2ObservableEquivalence` |
| Level 3 | End-to-end `d_eff/d` | `benchmarks/two_phase/m2_exit_artifact.cc`, section `ler` |

Debug invariants 6, 7, 8, 9, 10 and 11 all have coverage. The full suite is 156 gtests green
(M0 + M1 + M2), with the inherited PyMatching suite at its M0 baseline.

### Baselines and gates

- Inherited C++ suite plus M1's: unchanged. With the three M2 suites: **156 tests**, green.
- Inherited Python suite: **105 passed**, the M0 baseline recorded in the M1 exit report.
- **Stock-path perf gate: green by construction.** M2 edited no vendored file, so the stock decode
  path is byte-identical to the one M1 measured against `93f4212`. `pyrematching_perf` was re-run to
  confirm nothing moved: `Decode_surface_r11_d11_p1000` 230 k shots/s (M1 recorded 210–220 k),
  `Decode_surface_r11_d11_p100` 8.9 k (8.6–8.8 k), `Decode_surface_r21_d21_p100` 610 shots/s
  (570–580). Adding files to `SOURCE_FILES_NO_MAIN` cannot change the codegen of an existing
  translation unit.

## Level 1 — one documented departure from the design, with a measurement behind it

§M2.6 asks for the residual *set* and the committed *pair set* to be equal outright, on the grounds
that "degenerate path choice cannot move them", and forbids adding a tolerance. Path choice indeed
cannot move them. But the **choice among optimal primal solutions** can, for a reason the theorem of
§M2.0 does not cover:

> In `G`, a growing region's flood is blocked by its neighbours' territory, so `G` never observes
> some tight collisions that `H` — where every defect pair within `2T` is a direct edge — does
> observe. Both are correct blossom implementations of the same metric. They therefore agree on the
> dual solution and on the value of the optimum; they need not agree on *which* optimum they land on
> when several tight events fall at the same instant.

This was measured before deciding, over 4800 shot-decodes (six corpora × four horizons × 200 shots,
including the `p = 0.01` checked-in DEM and the 20%-negative-weight DEM, which is the most degenerate
case in the suite):

| Quantity | Shots where `H` differed from `G` |
|---|---|
| `Sum_S y_S` at truncation | **0 / 4800** |
| Total committed weight | **0 / 4800** |
| `num_trees` (hence residual size) | **0 / 4800** |
| `committed_boundary` | **0 / 4800** |
| Observable bytes | **0 / 4800** |
| Committed *pairing* | 0–25%, rising with defect density |
| Residual *set* | 0% at `T >= 1.5`; up to 58% at `T = 0.5` on the negative-weight DEM |

In every divergent shot, both residuals satisfied the separation invariant `Y(u) == T`, so both are
valid exposed sets and Phase 2's error bound holds either way.

Level 1 is therefore asserted as: hard equality, no tolerance, on `Sum_S y_S`, the committed weight,
`num_trees`, `committed_boundary`, the residual size, the separation invariant on both sides, and
the partition invariant (committed ∪ residual is exactly the shot's detection events, each
classified once). A residual or pairing difference that satisfies all of those is recorded as a
**tie**, in `BallProfile::{residual_ties,pairing_ties}` and in the exit artifact. Anything else
still throws, and remains debug invariant 10.

On the circuit-noise grid of the exit artifact (`d = 5..13`, `p = 5e-4..5e-3`, `T ∈ {1.5, 2}`,
2000 shots per point) the **residual tie rate is 0.0000 at every one of the 40 points**, and the
pairing tie rate runs from 0 to 3.9% (`d = 13`, `p = 5e-3`, `T = 2`). Nothing downstream of harvest
consumes the pairing — M3–M6 consume the residual and the committed weight — so this is invisible to
the rest of the pipeline.

## Level 2 — observable equivalence

Obs bytes are **bit-identical to M1 on every shot** of every corpus tested, so the
`mask_divergence_rate` is 0 and there is nothing to certify homologically trivial. This is what
should be expected: a differing-homology cycle between two shortest paths within `2T` would be a
logical operator of weight `<= 4T`, far below the code distance. `BallStats::ambiguous_mask_pairs`
(the perturbed-tie-break certification, `certify_masks`) is likewise 0 across the whole `d = 5..13`
sweep at `R = 4` edge weights, and 8 of 11238 entries only on a `d = 5` fixture where `R` is
comparable to the code distance itself.

## Level 3 — end-to-end effective distance

Full pipeline (front end + exact cleanup of the residual on `G`) at `T = 2`, `d = 5..13`,
`p ∈ {5e-4, 1e-3, 3e-3, 5e-3}`, 10 000 shots per point. M2's LER equals M1's at 19 of 20 points,
differing at one point by a single shot in 10 000 (`d = 13`, `p = 5e-3`: 16 vs 17 failures). Fitting
`ln(LER)` against `d`:

| `p` | exact | M1 + cleanup | M2 + cleanup | M2 − M1 |
|---|---|---|---|---|
| 3e-3 | −0.514 ± 0.075 | −0.514 ± 0.075 | −0.514 ± 0.075 | 0.000 |
| 5e-3 | −0.247 ± 0.007 | −0.206 ± 0.011 | −0.212 ± 0.012 | −0.006 ± 0.017 |

`d_eff/d` matches M1 well within statistical error. **The gate that protects M1 result 9 passes.**

## Ball table cost — the local-memory budget

`T_max = 2`, `R = 4` edge weights, one shell.

| `d` | nodes | entries | table | mean \|B\| | max \|B\| | mean bitset words | compile |
|---|---|---|---|---|---|---|---|
| 5 | 120 | 11 238 | 0.28 MiB | 93.7 | 116 | 2.0 | 0.00 s |
| 7 | 336 | 41 848 | 1.02 MiB | 124.5 | 182 | 5.7 | 0.00 s |
| 9 | 720 | 101 040 | 2.46 MiB | 140.3 | 225 | 10.1 | 0.01 s |
| 11 | 1 320 | 200 440 | 4.93 MiB | 151.8 | 241 | 15.9 | 0.02 s |
| 13 | 2 184 | 351 136 | 8.77 MiB | 160.8 | 242 | 23.0 | 0.04 s |

The table grows linearly in the node count with a mean ball of ~160 entries, i.e. ~4.2 KiB per
detector node. `|B|` is flattening with `d` — it is a local quantity, as the design expects — so the
per-node budget is the number to quote, not the total.

Mean degree of `H` at `2T` under the corpus defect density runs 1.6 (`p = 5e-4`) to 13.1
(`p = 5e-3`, `T = 2`), with `h_edges` up to 1016 at `d = 13, p = 5e-3`.

## The exit read — `speedup_vs_m1 < 1.5`, so stop and profile

`speedup_vs_m1 = sum_g_reference_ns / sum_total_ns`, `BITSET` mode, 2000 shots per point:

| `d` | `p` | `T` | M2 (µs) | M1 (µs) | speedup | intersect | h_build | mwpm_build | blossom | harvest |
|---|---|---|---|---|---|---|---|---|---|---|
| 13 | 5e-4 | 1.5 | 8.32 | 5.11 | **0.61** | 0.45 | 0.02 | 0.10 | 0.26 | 0.13 |
| 13 | 1e-3 | 1.5 | 18.8 | 11.1 | **0.59** | 0.47 | 0.02 | 0.11 | 0.26 | 0.12 |
| 13 | 5e-3 | 2.0 | 179 | 84.9 | **0.48** | 0.40 | 0.07 | 0.17 | 0.28 | 0.07 |
| 9 | 5e-4 | 1.5 | 2.21 | 1.43 | **0.65** | 0.25 | 0.03 | 0.13 | 0.33 | 0.17 |
| 5 | 5e-3 | 1.5 | 3.53 | 2.68 | **0.76** | 0.16 | 0.05 | 0.13 | 0.45 | 0.14 |

Best `speedup_vs_m1` at the lowest `p` in scope (`5e-4`): **0.65**. The design's read applies:

> **Read:** if `speedup_vs_m1 < 1.5` at the lowest `p` in scope, stop and profile before starting
> M3. The front end is then not the bottleneck and compression will not rescue it.

### The profile, as the read demands

**The `beta * V` floor that M2 exists to remove is not present in M1 as landed.** M2's speedup
argument rests on M1 result 6, `C ≈ alpha * n^1.4 + beta * V`, with `beta * V` the cost of flooding
empty graph. The measurements above show M1's Phase 1 costing ~310 ns per defect at `d = 13,
p = 5e-4` and ~270 ns per defect at `d = 9, p = 1e-3` — flat in `V` across a 3× change in graph size
at matched defect counts. The flooder only ever touches nodes near a defect, and M1.3's
`reset_for_next_shot` already replaced `Mwpm::reset`'s whole-graph sweep with work proportional to
what the shot left behind. Whatever `beta * V` the M1 campaign fitted, the code as landed does not
pay it per shot, so there is no floor for `H` to remove.

**What M2 pays instead.** The split is dominated by two costs `G` does not have:

- `intersect` (25–50%, and rising with `d`) — reading the ball tables. At `d = 13` the table is
  8.8 MiB, so each defect's ball is a cold-ish stream of ~160 entries (`SCAN`) or ~23 words
  (`BITSET`). `BITSET` is the faster path everywhere and is what the numbers above use; `SCAN` costs
  a further ~1.4×. Note that at the operating point `2T == R`, so `SCAN` never breaks early and
  walks the whole ball. This is the cost the design accepts by assumption — "the hardware argument
  is local memory per node, and the benchmark takes the limit of memory access time to zero" — but
  on this machine that limit is not taken, and it is the single largest line item.
- `mwpm_build` (10–17%) — rewriting `H`'s adjacency each shot. One clear inefficiency was found and
  removed while profiling: sorting `2 * |E|` directed pairs per shot was 17–45% of the decode, and
  is gone (`ball_mwpm.cc` now exploits `h.edges` already being sorted by `(i, j)`). That change took
  `speedup_vs_m1` at `d = 13, p = 5e-4` from 0.51 to 0.61.

`blossom_on_h` is only 19–32% of M2's total, so the solver *is* cheaper on `H` than on `G` — the
premise that `H` is a smaller problem holds. It is the cost of *constructing* `H` that exceeds what
it saves.

### What this does not change

M2's correctness deliverable stands: `H` reproduces `G`'s dual solution, residual, committed weight
and observables, and preserves `d_eff/d`. The re-hosting is exact. What the measurement rejects is
the claim that re-hosting is *cheaper*, on this implementation of M1 and this machine.

Per the design, M3 should not start on the strength of these numbers.

## Structural counters — what the off-CPU stages move

Wall time for `isect`, `hbld` and `mwpm_build` measures this laptop's DRAM latency, which is exactly
the thing the design's hardware argument assumes away ("the benchmark takes the limit of memory
access time to zero"). Those three stages are the ones intended to leave the CPU, so alongside the
timings the artifact now reports the **structural** quantities behind them: bytes moved out of the
ball tables, edge records emitted, element writes into the `Mwpm` on `H`. They are
machine-independent and feed a latency model for any target architecture.

Raw data: `benchmarks/two_phase/results/m2_structural_counters_{scan,bitset}.{csv,log}`, the same
`(d, p)` grid as the exit artifact plus `T ∈ {1.5, 2}`, 2000 shots per point. Rendered by
`summarize_m2_artifact.py` as two new tables. The existing artifact files and every existing column
are untouched.

**Element sizes**, read off the array types rather than assumed: `ball_target` 4 B, `ball_w_int`
4 B (`pm::weight_int` is `uint32_t`), `ball_words` 8 B, `ball_word_rank` 4 B, `ball_entry_by_rank`
8 B, `ball_mask_offsets` 8 B, `ball_mask_ids` 4 B, `has_bcost` 1 B, `bcost_w_int` 4 B.

### What is counted where

- **`isect_scan_bytes`** — the traversal. `BITSET`: the window words actually read. `SCAN`: the ball
  entries actually walked, at 8 B each. Both are counted in the loop, so an early break is counted
  as it happened, not as the model predicts. Per-node CSR offset lookups (16–20 B per defect against
  a ~1 kB window) are excluded in both modes.
- **`isect_hit_bytes`** — what only a candidate pair causes to be fetched: in `BITSET` the rank
  indirection and the weight (16 B per candidate, since the bitset carries neither, and they must be
  fetched *before* the `2T` test can reject the pair); in both modes the observable id list of every
  pair that becomes an edge, plus the `bcost_*` boundary lookups. The mask lists are physically read
  by `BallMwpm::rebuild` rather than by the intersection loop; they are charged here because the
  pair is what selects them and the hardware stage being sized does both.
- **`isect_scan_bytes_other_mode`** — the counterfactual, derived from the tables. It is not a
  model: it reproduces the other mode's *measured* traversal byte for byte on every shot of the
  grid, including `SCAN`'s early break and `BITSET`'s window clamp, and a test asserts the equality
  per shot (`ManifoldBallGraph.StructuralCountersDescribeTheTraversalThatRan`).
- **`hbld_edges_written`** — `h_edges + h_boundary_edges`, but taken at the `push_back`s; a debug
  assert ties the count to `H`'s own sizes, which is where "each undirected pair is emitted once"
  is checked rather than assumed. `h_max_degree` (the existing `max_degree`) is plumbed through
  next to it: edges are crossbar area, degree is crossbar width.
- **`mwpm_init_elements`** — instrumented at the write sites in `ball_mwpm.cc`, split node/edge.

### `mwpm_init_elements` is not `h_nodes + 2*h_edges + h_boundary_edges`

Measured at `d = 29, p = 1e-3, T = 2` (2000 shots):

| | per shot |
|---|---|
| `mwpm_init_edge_elements` | 1682.72 |
| `2 * h_edges + h_boundary_edges` | 2 × 803.482 + 75.756 = 1682.72 |
| `mwpm_init_node_elements` | 408.94 |
| `h_nodes` | 390.52 |
| total, measured | 2091.66 |
| total, naive | 2073.24 (+0.89% measured) |

The edge side is exactly the naive count — adjacency **is** stored directed, so an undirected edge
of `H` is two records and a boundary edge is one. The node side is not: the rebuild resets
`max(used_nodes, |H|)` nodes, not `|H|`, because it has to clear whatever the *previous* shot left
behind. So the node count is shot-history dependent and runs 4.7% above `h_nodes` at this point (it
is bounded by the largest shot seen). Each edge record is four parallel-array element writes, the
fourth (`neighbor_implied_weights`) written by a separate resize pass over the nodes.

### The three sanity checks

1. **`isect_bytes_per_defect` at `d = 29`, `BITSET`.** Predicted `mean ball_word_len × 8 =
   118.7 × 8 = 949.2 B`; measured scan component **953.4 B/defect** at `p = 1e-3` (+0.4%), total
   `isect_bytes_per_defect` 1025.2 B including the hit side. The window is read exactly once per
   defect. The +0.4% is defects landing on interior nodes slightly more often than the node average.
2. **`p`-dependence.** The scan component per defect is *exactly* independent of `T` — the `T = 1.5`
   and `T = 2` rows agree to the byte at every `(d, p)` in `BITSET` — and the hit component grows
   with `p` as predicted (`d = 29`: 26.8 → 36.7 → 46.8 → 71.8 B/defect over `p = 1e-4 … 1e-3`). The
   scan component per defect does **step** with `p`, by 8% at `d = 29` (881.4 B at `p ≤ 3e-4`,
   953.4 B at `p ≥ 5e-4`), which is not what the check anticipated. It is not the loop: the ball
   tables are recompiled per `(d, p)` because `R` is quoted in *median edge weights* and the median
   edge weight moves with `p`, so the table's own `mean ball_word_len` steps the same way — 108.9,
   108.9, 118.6, 118.7 words over the four `p`. Within one table the counter is flat, and the same
   quantity computed by the independent counterfactual path agrees per shot. Model intact; the
   `p`-dependence lives in the table, not the traversal.
3. **`hbld_edges_written / h_nodes` against mean degree.** `d = 29, T = 2`: 0.85 vs `deg/2 = 0.66`
   at `p = 1e-4` (ratio 1.30), 2.25 vs 2.06 at `p = 1e-3` (ratio 1.09). The excess is the boundary
   edges, which the degree does not count. Well inside the factor of ~2 the convention allows.

### `SCAN` never breaks early at `2T == R` — and the crossover has moved

The counters confirm the exit report's assertion and quantify its converse. At `d = 29, p = 1e-3`,
per defect: `SCAN` walks 191.4 entries at `T = 2` against a ball of 187.7 (the whole ball, plus the
interior-node bias), and only **89.2** at `T = 1.5`, where `2T = 3 < R = 4` and the sorted ball
*does* terminate the walk at 47% of its length.

That matters for the mode choice, because `BITSET` reads a window that grows like `d²` while `SCAN`
reads a ball that saturates. Bytes per shot at `p = 1e-3` (`BITSET / SCAN`):

| `d` | 13 | 17 | 23 | 29 |
|---|---|---|---|---|
| `T = 2` | 0.14 | 0.22 | 0.40 | 0.62 |
| `T = 1.5` | 0.28 | 0.47 | 0.84 | **1.34** |

At the operating horizon `BITSET` is still the cheaper traversal at `d = 29`, by 1.6×, but the
margin is closing and extrapolating the `d²`-versus-saturating-ball scaling puts the crossover near
`d ≈ 37`. At `T = 1.5` it has already crossed, between `d = 23` and `d = 29`. The exit report's
"`BITSET` is the faster path everywhere" is a statement about `d ≤ 13`.

### Cost, and why the counters are behind a flag

Filling the counters costs **26% (`d = 13`) to 44% (`d = 23`) of `intersect_ns`** at `p = 1e-3`,
because charging a pair its observable id list means touching `ball_mask_offsets`, which the
intersection otherwise never reads. Leaving that in the timed path would have moved the numbers in
the timing table beside it, so they sit behind `BallConfig::collect_structural_counters` (default
off) and the artifact collects them in a separate *untimed* pass, exactly as it already does for the
§M2.6 tie rates and the §M2.9.6 diagnostics. Being structural, they do not care which pass measured
them — the untimed pass reproduces the timed pass's counters exactly.

With the flag off, what remains in the profiled path is loop-carried register increments and no
memory traffic. Three interleaved rounds against a `e3447cc` baseline binary, `d ∈ {13, 23}`,
`p = 1e-3`, 1000 shots: `mean_total_ns` 31.1/34.2/35.3 vs 34.2/34.1/35.3 µs at `d = 13` and
223/300/287 vs 309/262/260 µs at `d = 23`, differences going both ways and inside the spread of the
baseline's own rounds. The stock-path perf gate is green by construction (no vendored file was
touched) and was re-run to confirm: `Decode_surface_r11_d11_p1000` 230 k shots/s both sides,
`Decode_surface_r11_d11_p100` 8.7–8.9 k both sides, `Decode_surface_r21_d21_p100` 590–610 both
sides, four interleaved rounds.
