# M2.9 exit report — harvest restructuring

Covers `design/pyrematching_design.md` §M2.9: cut harvest's critical path without changing one byte
of its output. §M1.3's commit policy is unchanged and remains binding. What changed is *how the
committed set is enumerated*, and — for §M2.9.5 — the answer to whether blossoms have to be
shattered at all.

The governing metric from the M2 exit read onwards is the critical path in dependent operations,
not wall time on one machine, so the headline numbers below are counts and depths. Wall time
appears in exactly one place, the software A/B of §M2.9.6, which the design asks for by name.

## What landed

| Design item | Where |
|---|---|
| M2.9.1 direct enumeration of tree commits | `truncation/harvest.cc`, `Harvester::harvest_impl` phase B |
| M2.9.1/M2.9.2 arena liveness tracking | `sparse_blossom/arena.h`, plus one `uint32_t` on `AltTreeNode` and `GraphFillRegion` |
| M2.9.2 matched-region enumeration, dedup removed | `truncation/harvest.cc` phase A |
| M2.9.2 debug equivalence guard | `assert_live_top_regions_match_detection_event_sweep`, live in every debug harvest |
| M2.9.3 reductions, residual compaction unchanged | `truncation/harvest.cc` phases A and C |
| M2.9.4 base descent | **unchanged from M1.4** — see the decision below |
| M2.9.5 E1/E2 experiments | `truncation/split_experiment.{h,cc}`; **no deferred-shatter code exists** |
| M2.9.6 measurements 1–3 and the profile fields | `truncation/harvest.{h,cc}`, `perf/ball_profile.h`, `perf/two_phase_profile.h` |
| M2.9.6 measurement 4 (solve depth) | `truncation/truncated_timeline.{h,cc}`, `TimelineDepthModel` |
| M2.9.7 H1–H8 | `tests/two_phase/truncation_harvest_m29.test.cc`, `…/truncation_split_experiment.test.cc` |
| Exit artifact | `benchmarks/two_phase/m29_exit_artifact.cc`, results in `benchmarks/two_phase/results/` |

The M1.3 enumeration is still compiled in, behind `Harvester::use_legacy_enumeration`, as the H1
oracle — exactly as `verify_against_g` is the oracle for M2. It is also what the A/B measures
against, in-process.

## The enumeration, and why it is shaped this way

§M2.9.0's load-bearing observation holds: a root is `inner_region == nullptr`, one comparator, and
a non-root commits its own pair from its own `inner_to_outer_edge`. `AltTreeNode::find_root` and the
`children` descent are gone from harvest. So is the detection-event sweep, and with it the visited
stamp that deduplicated regions reached from several of their own defects.

Harvest is now three sweeps, none of which follows the structure it is reading:

- **A** — one scan of the region arena's live vector. Every region alive at truncation is visited
  once, so `dual_sum_at_truncation` and `max_region_dual` are flat reductions rather than a descent
  through blossom nesting, and "top-level" is one local `blossom_parent == nullptr` test.
- **B** — one scan of the node arena's live vector, each node classified from its own three fields,
  committing and freezing its own pair. Then one pass over the roots for the base descent.
- **C** — extraction, driven once per committed match from its representative endpoint.

### Enumerability: what was actually built, and why it is not the design's intrusive list

Upstream exposes no iterable container of live `AltTreeNode`s, so §M2.9.1's open question had to be
answered. The design's preference is an intrusive doubly-linked list, with a stated fallback: *"if
the intrusive list turns out to need edits at more than two call sites, prefer an arena sweep with a
liveness bit and keep the vendored tree untouched."*

Counting the sites: `AltTreeNode` is allocated at two places (`Mwpm::make_child`,
`Mwpm::create_detection_event`) and deleted at two (`Mwpm::shatter_descendants_into_matches_and_freeze`,
`AltTreeNode::prune_upward_path_stopping_before`) — four, across two files. §M2.9.2's list of
top-level matched regions is worse: four `add_match` calls, the direct boundary match, and the two
`Match::clear()`s in `handle_tree_hitting_match`, so seven, all in `mwpm.cc`. That is three times the
design's own threshold, and the design is right that a missed one is silent.

**Every one of those sites funnels through `Arena<T>`.** So liveness is maintained there instead:
one `std::vector<T*> live` per arena, membership by swap-remove, and one `uint32_t arena_live_index`
in the object so removal is O(1). There is no mutation site to miss, because allocation and
deallocation are the only two events and the arena *is* both of them. `mwpm.cc` still changes, but
only mechanically: two `alloc_unconstructed()` + placement-new pairs become `alloc_constructed(...)`,
because the link has to be installed after the object's lifetime begins or the placement-new would
clobber it.

Two deviations from the design's letter, both deliberate:

1. **A vector, not a doubly-linked list.** A `prev`/`next` pointer pair grew `GraphFillRegion` from
   128 bytes to 152 and `AltTreeNode` from 104 to 128, and that showed up as a **consistent 2–4%
   regression in the stock decode path** across four interleaved rounds — purely cache footprint. A
   single `uint32_t` slot instead costs nothing at all: it lands in padding that both structs
   already had, and both are byte-for-byte the size they were before. (On `GraphFillRegion` the hole
   is at offset 44, between `shrink_event_tracker` and `match`, which is why the field is declared
   there rather than at the end of the struct. There is a comment saying so; do not tidy it to the
   bottom.) A vector is also the better shape for what the list is *for*: enumeration is a linear
   scan of contiguous pointers, not a chase.
2. **Regions are tracked too, not just tree nodes.** This is what lets §M2.9.2's dedup disappear
   without seven maintenance sites, and it additionally turns the dual-sum reduction from a
   nesting descent into a flat loop.

Five vendored files gained rows in [`upstream_edits.md`](upstream_edits.md), against a design that
anticipated one or two. The trade taken is fewer *mutation sites*, not fewer files, and the design's
own guard rail is the site count — that is what the silent-failure risk scales with.

### The one genuinely new decision: which endpoint drives an extraction

A match is between regions, and both endpoints see the same pair, so exactly one has to drive the
extraction or the pair would be extracted twice. M1.3 made that choice implicitly — whichever
endpoint the detection-event sweep reached first. With the sweep gone it has to be explicit, and it
is made locally and canonically: `match.edge.loc_from` lies in this region and `match.edge.loc_to`
in the partner, so comparing the two node pointers compares their detector ids, and the endpoint
holding the lower one wins. A boundary match is always its own representative.

This is why H1's match-edges comparison is of the committed pair *set*, canonicalised so that
`from < to`: the obs flavour is bit-identical, and in the edges flavour a `CompressedEdge` can come
out reversed relative to M1.3 when M1.3 happened to reach the higher-id endpoint first. Every
consumer in the tree already normalises (`BallDecoder::sort_pairs`, the partition tests), and the
new rule is the more canonical of the two.

## H1 — the deliverable

**Bit-exact against M1.3's landed harvest, both flavours, no tolerance.** Compared per shot:
`residual`, `residual_dual_sum`, `committed` (obs mask and weight), `committed_pairs_frozen`,
`committed_pairs_tree`, `committed_pairs_blossom_cycle`, `committed_boundary`, `num_trees`,
`largest_tree_size`, `exposed_root_blossoms`, `dual_sum_at_truncation`, `max_region_dual`, and — in
the edges flavour — the committed pair set.

Corpus: five DEMs (`surface_code_d13_p0.01`, its negative-weight variant, `toric_code_d5`, and two
generated surface codes at `d = 5, 9`) × 40 shots × 6 horizons `T ∈ {0, 1, 2, 3, 5, 8}` edge
weights, both flavours. Green.

| # | Requirement | Test |
|---|---|---|
| H1 | Harvest identity against the M1.3 oracle | `TwoPhaseHarvestM29.H1ObsFlavourIsBitExactAgainstTheM1Oracle`, `…H1MatchEdgesFlavour…` |
| H2 | Maintained enumeration == detection-event sweep's region set | `…H2LiveRegionListEqualsTheDetectionEventSweep`, plus a debug assert inside **every** harvest |
| H3 | Committed ∪ residual partitions the shot's events | `…H3CommittedAndResidualStillPartitionTheDetectionEvents` |
| H4 | Separation invariant `Y(u) == T` | `…H4SeparationInvariantSurvivesTheRestructuring` |
| H5 | Base descent still yields the base | `…H5ExposedRootBlossomsStillYieldTheirBase` |
| H6 | E1/E2 reported | `TwoPhaseSplitExperiment.H6ReportE1AndE2OverTheCorpus` (+ evaluator self-check, + E2 structurally) |
| H7 | Reusability and determinism | `…H7ReusabilityAndDeterminism` |
| H8 | Fuzz | `…H8FuzzRandomHorizonsAndSyndromes` |

### Baselines and gates

- C++ suite: **168 tests green** — 156 at the M2 baseline plus 12 new. Inherited PyMatching tests
  unchanged.
- Inherited Python suite: **105 passed**, the M0 baseline.
- **Stock-path perf gate: green.** Four interleaved rounds of `pyrematching_perf` against a
  baseline binary built from `fb44df0` (the M2 tip). `Decode_surface_r11_d11_p1000` 230 k shots/s
  both sides; `Decode_surface_r21_d21_p10000` 210–240 k both sides; `Decode_surface_r21_d21_p1000`
  18–20 k base vs 19–20 k; `Decode_surface_r21_d21_p100` 500–630 both sides;
  `Decode_surface_r11_d11_p100` 8.8–8.9 k vs 8.6–8.9 k. Differences go both ways round to round.
  The earlier pointer-pair layout did **not** pass this gate, which is why it was replaced.

## §M2.9.6 — the measurements

Grid: `d ∈ {5,7,9,11,13}` × `p ∈ {1e-3, 3e-3, 5e-3}` × `T ∈ {0.75, 1.0, 1.5, 2.0}` edge weights ×
1000 shots, run on `H` through `BallDecoder`. Raw data in
`benchmarks/two_phase/results/m29_exit_artifact.{csv,log}`.

**1. `largest_tree_size`.** Maximum **8** over the whole grid; means below 1 everywhere. Any
tree-side readout that survives is bounded by a very small number. Note it is now a diagnostic
rather than a free by-product — direct enumeration never learns which tree a node is in, so
computing it costs a root id per node and it is filled only under
`Harvester::collect_diagnostics`.

**2. Blossom nesting depth and member count.** Over *all* blossoms alive at truncation, the maximum
depth reaches 9 and the maximum member count 25 — so the design's "expect depth 1–3" is not a bound.
But the distribution that matters for §M2.9.4 is the one restricted to **exposed root blossoms**,
which is the only nesting the base descent walks:

| max exposed-root-blossom depth in a shot | 1 | 2 | 3 | 4 | 5 | 6 | 8+ |
|---|---|---|---|---|---|---|---|
| shots | 3116 | 505 | 88 | 21 | 5 | 2 | 2 |
| share | 83.3% | 13.5% | 2.4% | 0.6% | 0.13% | 0.05% | 0.05% |

**99.2% at depth ≤ 3**, and depth 1 in five shots out of six. M1 result 7's prediction holds where
the design applies it; the long tail is in matched blossoms, which do not use this path.

**3. Blossom formations per shot.** 0.02 at `(d=5, p=1e-3)` rising to **3.06** at
`(d=13, p=5e-3, T=2)`. Base readouts — one per exposed root blossom — occur on 3739 of 60 000 shots,
i.e. **6.2%**. Formations therefore outnumber base readouts by one to two orders of magnitude.

**4. Dependent-event chain depth.** Modelled as a longest-chain count over detector nodes: an event
naming `u` and `v` has depth `1 + max(depth[u], depth[v])` and writes it back to both; a blossom
shatter takes the max over the nodes the blossom owns. Mean solve depth runs **0.6 → 5.2** across
the grid, against a mean of up to **100 events** per shot — the solve is very wide and very shallow.
Harvest's own modelled depth (`1 + max(ceil(log2(commits)), ceil(log2(residual)), nesting chase)`)
runs **1.3 → 7.9**.

### Harvest's share of the remaining critical path — the estimate this section was motivated by

§M2.9 opens with "harvest is 7–17% of the total and looks minor… once [the parallel stages] collapse,
harvest becomes roughly a third of [the critical path]". Measured, in depth:

> **harvest is 60–77% of the remaining critical path (median 66%)**, not a third.

The estimate understated it by a factor of two, and it understated it because the *solve* is far
more parallel than assumed: 100 events deep 5. The section's motivation is stronger than written,
and after the restructuring harvest is still the larger half. That is the number M3-onward planning
should carry.

Caveat, stated rather than buried: the two depths count different kinds of dependent operation —
the solve's are matcher events, harvest's are reduction and chase steps — so the ratio is
indicative, not a cycle count. Both are machine-independent, which is the property the metric was
chosen for.

### The software A/B

The design asks: "enumerate `AltTreeNode`s directly and measure how much of harvest's 7–17% was
enumeration rather than extraction." Run as two passes of one process over the same shots, differing
in one flag (`BallConfig::use_legacy_harvest_enumeration`), with a warmup shot each and the pass
order alternated across horizons:

> **Harvest is 2.2× faster (median over 60 points; 1.36–2.53).** Its share of the whole shot falls
> from 5.8–27.2% (median 15.6%) to 2.5–18.0% (median 8.3%).

Roughly half of M1.3's harvest was enumeration. So §M2.9.1–§M2.9.2 are a real CPU win as well as a
critical-path one, and the register-maintenance model is confirmed on hardware that already exists —
which is what the design wanted this A/B for.

One outlier: `d=5, p=1e-3, T=0.75` reports 0.81, i.e. the legacy path winning. It is the first point
of the run, harvest is ~200 ns there, and the three sibling horizons at the same `(d, p)` all report
1.36–1.37. Warm-up, not a result.

Stage split of the restructured harvest, as fractions of `harvest_ns`: enumerate 11–27%, base
descent 0.7–17%, shatter 12–74%, reduce and compact the rest. **Shattering is now the dominant
stage**, and it grows with `T` as more of the syndrome ends up matched. That is the thing to attack
next if harvest is attacked again — which makes §M2.9.5's result below the important one.

## §M2.9.4 — mechanism chosen: **M1.4's recursion, unchanged**

The design names this a valid outcome and asks for the measurement that justified it.

Measurement 3 is the test the design sets: *"eager is right iff formations are meaningfully rarer
than base readouts."* They are not — formations run up to 3.06 per shot while a base is read on 6.2%
of shots. The eager cache would pay maintenance on every formation to save work on a path taken two
orders of magnitude less often. Measurement 2 covers the second half of the test: depth is 1 in 83%
of the shots that have an exposed root blossom and ≤ 3 in 99.2%, so the recursion is one or two
dereferences in almost every case. And the stage is 0.7–17% of a harvest that is itself 2.5–18% of a
shot.

Path doubling was not re-proposed and is not implemented; the design's three reasons stand, and
reason 1 is now measured rather than inferred.

**One correction to the section, recorded because it changes the arithmetic.** §M2.9.4 describes the
base descent as a "nested pointer chase" of depth 1–3. M1.4 did not implement it as a chase:
`find_exposed_base_node` scans *every leaf* of the blossom nesting and takes the argmax of the
nested dual sum, because the base is identified by `Y(u) == T` rather than by a stored pointer. Its
cost is therefore O(members), not O(depth), and members reach 23. An eager cached base would save
more than the section's depth argument implies. It is still rejected, on measurement 3 — but if
harvest is revisited, this is a better lead than the depth distribution suggests, and it should be
weighed against the shatter stage, which is larger.

## §M2.9.5 — **closed. E1 fails; E2 passes.**

Both experiments were run before any deferred-shatter code was written, and no deferred-shatter code
exists. The machinery is a read-only evaluator (`evaluate_shatter`) that reproduces
`Mwpm::shatter_blossom_and_extract_matches` exactly as a pure function of the region structure, so
the same live state can be evaluated twice under two split rules. It is validated against the real
destructive extraction, pair by pair, over the whole corpus
(`TwoPhaseSplitExperiment.EvaluatorReproducesTheRealExtraction`) — without that, "no difference"
would only have meant the evaluator was not modelling extraction.

Over the full grid, 12 084 matches in which at least one side is a blossom:

| Question | Result |
|---|---|
| **E1** — committed obs mask invariant under a rotated split? | **No.** 141 / 12 084 moved (**1.17%**) |
| E1 in the form that decides the saving — is the *interior* pairing's obs contribution zero, i.e. would not shattering at all give the same answer? | **No.** 81 / 12 084 non-zero (**0.67%**) |
| **E2** — committed weight invariant under a rotated split? | **Yes.** 0 / 12 084 |
| The section's premise — are blossom cycles homologically trivial? | **Almost always.** 3 / 14 802 cycles have a non-zero total obs mask (**0.02%**) |

**E1 fails, so §M2.9.5 is closed and matched blossoms must continue to be shattered in both
flavours.** The obs flavour and the edges flavour therefore remain one datapath, and the design's
question of "which one the target hardware emits" does not arise.

Why it fails, since the design expected it to pass and the reasoning is worth keeping. §M2.9.5
argues that "two splits differ by the full blossom cycle" and that such a cycle is homologically
trivial far below the code distance. The second half is right — measured at 0.02% — but the first
half is not, with the external match edge held where it is. Splitting an odd cycle at `s` versus at
`s'` gives two matchings whose symmetric difference is the **alternating path between `s` and `s'`**,
an open arc, not a closed cycle. It closes into a cycle only if the external match edge moves from
`s` to `s'` as well — which is a different primal solution, not a different reading of the same one.
An open arc has no homological reason to carry a zero observable mask, and 1.17% of the time it does
not.

E2 passes, and structurally rather than statistically: the weight upstream's extraction accumulates
is `Sum_S y_S` over every region in the nesting, each counted exactly once, and no choice of split
point changes which regions exist. That is asserted directly as
`TwoPhaseSplitExperiment.E2CommittedWeightIsIndependentOfTheSplitPoint`, so a future change that
breaks it fails a test rather than silently reporting a different upper bound. The design's worry —
"if E2 fails… say so and stop" — does not arise.

The two hard exceptions are moot: M5's edges flavour was always going to shatter, and exposed root
blossoms were always outside E1/E2's scope.

## M2.9 exit checkpoint

- [x] H1 green on the full corpus, both flavours, no tolerance. Oracle path retained behind
      `Harvester::use_legacy_enumeration`.
- [x] H2–H5, H7, H8 green; inherited PyMatching suite at the M0 baseline (168 C++, 105 Python);
      stock-path perf gate green against a `fb44df0` baseline binary, four interleaved rounds.
- [x] Measurements 1–4 of §M2.9.6 recorded, in counts and depths.
- [x] §M2.9.4 mechanism chosen — **M1.4's recursion, unchanged** — with the depth distribution and
      the formation rate that decided it written down above.
- [x] E1 and E2 reported. **E1 fails; §M2.9.5 is closed**, with the reason recorded above.
- [x] `docs/upstream_edits.md` updated: five new rows, plus a note on the exceeded file budget.
- [x] Harvest's share of the remaining critical path re-stated against measurement 4: **60–77%,
      median 66%**, replacing the "roughly a third" estimate this section was motivated by.

## Deviations from the design, and why

1. **An arena-owned live vector, not an intrusive doubly-linked list, and regions tracked as well as
   tree nodes.** Site counts and a failed perf gate; both argued above.
2. **`largest_tree_size` is no longer free.** The design anticipates this and asks for it to be kept
   as a profiling-only diagnostic; it is, behind `Harvester::collect_diagnostics`, and every caller
   that consumed it (`TwoPhaseProfile::fill_from`, the exit artifacts) either sets the flag or reads
   zero deliberately.
3. **`Arena::~Arena` empties `live` by move rather than by `clear()`.** `Mwpm::reset()` invokes the
   arena destructor explicitly and then lets the object be destroyed again; the pre-existing members
   survive that only because they are moved-from. A `clear()` keeps capacity and double-frees on the
   second destruction — caught by ASan in
   `MwpmDecoding.NoValidSolutionForLineGraph`. Upstream's double-destruct is not fixed here; the new
   member is just made to survive it the same way the old ones do.
4. **E1 is reported in two forms.** The design words it as split-invariance; the form that actually
   decides whether the shatter can be dropped is whether the interior pairing contributes any
   observable flip at all. Both are measured, both fail, and the premise behind the design's
   expectation is measured separately and confirmed — the failure is in the step from the premise,
   not in the premise.
