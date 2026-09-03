// Copyright 2026 SpecMatching contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef SPECMATCHING_SPEC_MATCHING_MANIFOLD_BALL_DECODING_H
#define SPECMATCHING_SPEC_MATCHING_MANIFOLD_BALL_DECODING_H

#include <string>
#include <vector>

#include "specmatching/sparse_blossom/driver/mwpm_decoding.h"
#include "specmatching/sparse_blossom/driver/user_graph.h"
#include "specmatching/spec_matching/manifold/ball_graph.h"
#include "specmatching/spec_matching/manifold/ball_mwpm.h"
#include "specmatching/spec_matching/manifold/ball_tables.h"
#include "specmatching/spec_matching/perf/ball_profile.h"
#include "specmatching/spec_matching/truncation/harvest.h"
#include "specmatching/spec_matching/truncation/truncated_timeline.h"
#include "stim.h"

namespace pm {
namespace spec_matching {

struct BallConfig {
    /// The horizon, in DEM float weight units. Must be `<= ball.T_max`.
    ///
    /// On the truncated path this is the *truncation* horizon. On §M7's stock path it is a pure
    /// runtime scalar with exactly two uses: the `2 * T_int` / `T_int` filter that builds `H`, and
    /// the `max_u Y(u) >= T_int` escalation threshold. It is never handed to the flooder there.
    double T{2.0};
    BallParams ball;

    /// §M7 — run **stock** (untruncated) sparse blossom on `H` and decide exactness from a max-dual
    /// certificate, instead of truncating the timeline at `T` and harvesting.
    ///
    /// The `H` that is solved is the same graph either way: the `2 * T_int` / `T_int` ball filter is
    /// untouched, and so are the ball tables, `BallMwpm::rebuild`, and stock's own extraction. What
    /// changes is the front end. The flooder runs with its `pm::NO_HORIZON` sentinel — no funnel
    /// gating, no shrink exemption, no RAII horizon guard — and the shot is kept iff it completed on
    /// `H` with `max_u Y(u) <= T_int`, which certifies it a *global* MWPM (§M7.0). Everything else
    /// escalates to stock on `G`, exactly as the truncated path's residual does.
    ///
    /// The §M1.3 harvest and §M1.4's exposed-root-blossom base descent are not reachable on this
    /// path at all: a certified shot completed, so it has no surviving trees, no residual and no
    /// exposed root blossoms. The truncated path and its M1 oracle stay exactly as they are; this is
    /// an alternative front end, not a replacement (§M7 non-goals).
    bool stock_on_h{false};
    /// Run M1's path on `G` for every shot and assert §M2.6 level 1 inline. CI and benchmarks only;
    /// off by default because it costs more than the decode it is checking.
    bool verify_against_g{false};
    BallGraphBuildMode mode{BallGraphBuildMode::SCAN};
    /// Threads used to compile the ball tables. 0 means hardware concurrency.
    size_t compile_threads{0};
    /// Fill the §M2.9.6 measurements: harvest's stage split, the blossom counts and depths, and the
    /// dependent-event chain depth of the solve. Off by default, because collecting them costs more
    /// than some of the stages they are measuring — the M2.9 exit artifact takes them in a separate
    /// untimed pass, exactly as it already does for the §M2.6 tie rates.
    bool collect_harvest_diagnostics{false};
    /// Fill the §M2 structural counters — `isect_*_bytes`, `hbld_edges_written`,
    /// `mwpm_init_elements`. Off by default for the same reason as the flag above, and measured:
    /// the byte counting adds 26% (`d = 13`) to 44% (`d = 23`) to `intersect_ns` at `p = 1e-3`,
    /// because charging a pair its observable id list means touching `ball_mask_offsets`, which the
    /// intersection itself never reads. Leaving it on would move the numbers in the timing table
    /// next to it. The M2 exit artifact therefore collects the counters in a separate *untimed*
    /// pass, exactly as it already does for the §M2.6 tie rates and the §M2.9.6 diagnostics.
    ///
    /// The counters are structural, so the separate pass costs nothing in fidelity: they are a
    /// function of the shot and the tables, not of when they were measured.
    bool collect_structural_counters{false};
    /// §3.5.2 — the component and structural statistics of the shot's `H`: the connected
    /// components, their sizes, weighted and hop diameters, boundary structure and per-component
    /// status, plus the `H` edge-weight and boundary-cost distributions.
    ///
    /// **Profiling only, and untimed.** The analysis runs after the shot's timed window has closed
    /// and is charged to no latency number at all: the component work is assumed free
    /// (hardware-offloadable), so timing it would be measuring a stage that is not meant to be on
    /// this critical path. It reads the `H` the arena is still holding and changes nothing about
    /// the decode — with the flag on or off, every shot produces byte-identical output and the same
    /// set of shots escalates.
    ///
    /// On by default on this branch: the component structure is a first-class output now that the
    /// decode is per component, and `sparse_graph_stats` exists to collect it.
    bool collect_component_stats{true};

    /// Components above this size are counted rather than given a diameter (§3.3's
    /// `--diameter-cap`). Both diameters are `O(s^2)` walks over the component, and at `p = 1e-3` a
    /// component this large is already far outside the distribution the statistic is for. Read only
    /// under `collect_component_stats`.
    uint32_t diameter_cap{MAX_DIAMETER_COMPONENT_SIZE};

    /// §2.5 — debug/bench: additionally run the **monolithic** solve on the whole of `H`, on a
    /// separate instance, and assert that it agrees with the per-component one.
    ///
    /// What is asserted: the escalation predicate (`any component TRUNCATED` against the monolithic
    /// `TimelineStatus`), and on a completing shot the committed weight and the observable bytes.
    /// With the full harvest on, also that the sorted residual set is the union of the
    /// per-component residual sets and that `num_trees` is the sum.
    ///
    /// This is what §1's independence claim is checked as. Off by default and costing nothing on
    /// the production path, exactly like `verify_against_g`.
    bool verify_component_decomposition{false};

    /// §B — skip the §M2.1 negative-weight preamble on `G` when the DEM has no negative-weight
    /// edge, which is the overwhelmingly common case (`log((1-p)/p) > 0` for `p < 0.5`).
    ///
    /// The skip only ever fires when the decoder's own scan of `G` also says all-positive, in which
    /// case the preamble is provably a no-op: there are no negative-weight detection events to
    /// symmetric-difference in and the weight/observable offsets are zero. With a negative-weight
    /// edge anywhere in the DEM the preamble runs exactly as it does today. The flag exists to be
    /// turned off for the A/B, not because the skip is conditional on anything else.
    bool skip_negative_weight_preamble_when_positive{true};

    /// Harvest with M1.3's enumeration instead of §M2.9.1's. The output is identical either way
    /// (that is H1); this exists so that the A/B of §M2.9.6 can be run as two passes of one process
    /// rather than as two runs of two binaries, which is the only way the difference — a few
    /// percent of a shot — is measurable above run-to-run noise.
    bool use_legacy_harvest_enumeration{false};
};

/// One committed pair, in `G`'s detector ids. `to == -1` means matched to the boundary. This is the
/// form §M2.6 level 1 compares and the form M5 lifts to `G`'s edges.
struct CommittedPair {
    int64_t from;
    int64_t to;
    /// Index into the ball pools of the canonical shortest path behind this pair, or
    /// `NO_BALL_ENTRY` when the pair is a boundary match (in which case `bcost_path_*` of `from`
    /// holds the path instead). M5 reads this; the obs flavour does not.
    uint64_t ball_entry;

    static constexpr uint64_t NO_BALL_ENTRY = UINT64_MAX;
};

/// What §M3.4's production Phase 1 yielded.
///
/// `status` is the branch, and after §2 it is a **reduction over the components**: `TRUNCATED` iff
/// some component's own truncated solve did not finish by `T`, in which case the shot escalates and
/// whatever the other components produced is discarded. Reading the status rather than
/// `harvest.residual.empty()` is what keeps the two from drifting apart when the bypass is on.
struct Phase1Outcome {
    TimelineStatus status{TimelineStatus::COMPLETE};
    HarvestResult harvest;
    /// §2.6's tally, so the one escalation predicate has a count behind it.
    /// `status == TRUNCATED  <=>  components_truncated > 0`.
    int components_total{0};
    int components_truncated{0};
};

/// Phase 1 executed on the defect manifold (§M2.5), **one connected component at a time** (§2).
///
/// Per shot:
///
/// ```
/// negative-weight preamble on G -> intersect balls -> build H
///   -> union-find over H's defect-defect edges -> components, roots ascending
///   -> for each component C, in ascending-root order:
///          build sub-H(C) on the component instance          (§2.3, §2.4)
///          status_C = process_timeline_until_horizon(C, T)
///          if COMPLETE: extract obs_C, w_C                    (§M3.4)
///          reset the instance
///   -> escalate iff some status_C is TRUNCATED
///   -> otherwise obs = XOR_C obs_C, weight = SUM_C w_C
/// ```
///
/// There is **one** question asked of a component, whatever its size, and it is asked of truncated
/// sparse blossom: no lookup table, no closed-form rule and no static feasibility test decides a
/// commit or a residual (§5.1). What licenses the split is §1's independence property — regions in
/// different components cannot interact before `T`, because interaction needs `d_G <= 2T` and that
/// is exactly the condition for an `H` edge — so the per-component event sequences are the disjoint
/// union of the monolithic one, and the escalating set and the output are identical either way.
///
/// Everything it emits is in `G`'s detector ids, so M3–M6 consume it unchanged: nothing downstream
/// of harvest knows `H` exists, let alone that it was split.
struct BallDecoder {
    BallConfig config;
    /// The detector graph. Kept because the negative-weight preamble, the boundary-node mask and
    /// the oracle path all live on it.
    pm::Mwpm g_mwpm;
    BallTables tables;
    /// §2.3's component instance: one `pm::Mwpm`, built once with capacity for the largest
    /// component seen and reused across components and shots — write the component's edges, solve,
    /// extract, reset, next component. Zero per-shot allocation after a warm-up shot (invariant 11).
    BallMwpm h_mwpm;
    BallGraphArena arena;
    Harvester harvester;

    /// §2.5. A **separate** instance for the monolithic cross-check, so the two solves share no
    /// state, and its own harvester so the counters of one are not read as the other's. Both stay
    /// empty unless `verify_component_decomposition` is set.
    BallMwpm verify_mwpm;
    Harvester verify_harvester;

    /// `T` in the flooder's time units, converted once from `config.T` (§0 unit rule).
    horizon_int horizon{0};

    /// §B. Does `G` carry a negative-weight edge at all? Scanned once at construction off the
    /// matching graph's own negative-weight record. When this is false the §M2.1 preamble has
    /// nothing to do, which is what `skip_negative_weight_preamble_when_positive` skips.
    bool dem_has_negative_weights{false};

    static BallDecoder from_detector_error_model(
        const stim::DetectorErrorModel& dem,
        BallConfig config,
        pm::weight_int num_distinct_weights = pm::NUM_DISTINCT_WEIGHTS,
        const char* ball_artifact_path = nullptr);

    /// Builds a decoder around an already-constructed `pm::Mwpm` for `G`. The `Mwpm` is moved in.
    static BallDecoder from_mwpm(pm::Mwpm g_mwpm, BallConfig config, const char* ball_artifact_path = nullptr);

    /// Same contract and same struct as M1's harvest, in `G`'s detector ids.
    ///
    /// **Rejected under `stock_on_h`.** The two verification entry points return only a
    /// `HarvestResult`, and on the truncated path the caller reads the escalation decision off
    /// `residual.empty()`. Under the certificate that reading is wrong in precisely the case §M7.0
    /// calls the most important correctness point in the design: a shot that completes on `H` with a
    /// suboptimal matching and `max_u Y(u) > T` has an empty residual and must still escalate. So
    /// rather than leave a signature whose obvious use is a silent wrong answer, these throw.
    HarvestResult decode_phase1(const std::vector<uint64_t>& dets, BallProfile* prof = nullptr);

    /// As `decode_phase1`, and additionally the committed pairs in `G`'s detector ids. Works for
    /// any number of observables, which the obs flavour does not. Rejected under `stock_on_h` for
    /// the reason above.
    HarvestResult decode_phase1_to_match_edges(
        const std::vector<uint64_t>& dets, std::vector<CommittedPair>& committed_pairs, BallProfile* prof = nullptr);

    /// §M3.4's production Phase 1: extraction only when the timeline completes, and no harvest at
    /// all when it truncates — the shot escalates and Phase 1's partial result is discarded.
    ///
    /// The two entry points above stay exactly as they were and remain the verification path: they
    /// are what §M2.6 level 1 compares against M1, what §M3.3 X8 compares this against, and what
    /// debug invariants 3, 4, 18 and 19 read. **Neither is dead code and neither may be deleted**
    /// (§M3.4, "keep the full harvest compiled in").
    Phase1Outcome decode_phase1_production(const std::vector<uint64_t>& dets, BallProfile* prof = nullptr);
    Phase1Outcome decode_phase1_production_to_match_edges(
        const std::vector<uint64_t>& dets, std::vector<CommittedPair>& committed_pairs, BallProfile* prof = nullptr);

    /// The observable mask and weight that `G`'s negative-weight edges contribute, to be combined
    /// with `HarvestResult::committed` exactly as the stock decode path does.
    inline pm::obs_int negative_weight_obs_mask() const {
        return g_mwpm.flooder.negative_weight_obs_mask;
    }
    inline pm::total_weight_int negative_weight_sum() const {
        return g_mwpm.flooder.negative_weight_sum;
    }

    /// The post-preamble detection events of a shot: the shot's events symmetric-differenced with
    /// `G`'s negative-weight detection events, with user-graph boundary nodes removed, sorted
    /// ascending. This is the one place `H` cannot be derived from the raw syndrome (§M2.1).
    void compute_seeded_detection_events(const std::vector<uint64_t>& dets, std::vector<uint64_t>& out) const;

    /// M1's Phase 1 on `G` for the same shot — the oracle §M2.6 level 1 is compared against.
    HarvestResult reference_phase1_on_g(
        const std::vector<uint64_t>& dets, std::vector<CommittedPair>* committed_pairs = nullptr);

    /// §M7.7's `q_current_on_same_corpus`: would the **landed truncated scheme** have escalated this
    /// shot? Rebuilds the same `H` at the same `T` and runs `process_timeline_until_horizon` on it,
    /// then tears the instance down without harvesting.
    ///
    /// Benchmark mode only, and it costs a second full Phase 1 on `H`. Call it *after* the shot's
    /// real decode has finished and released `h_mwpm`, never around it. It exists so that the
    /// design's `q_this <= q_current` claim is replayed on identical shots instead of compared
    /// across campaigns.
    bool truncated_scheme_escalates(const std::vector<uint64_t>& dets);

    /// §3.5.2. Decomposes the `H` the arena is still holding — the one the shot just decoded on —
    /// into connected components, and fills `prof.components`, `histograms` and the two joint
    /// tables from it.
    ///
    /// **Call it after the shot's timed window has closed, never inside one.** It starts no timer
    /// of its own and must not be wrapped in one: this experiment's latency account is the solver
    /// and the harvest, and the component work is assumed free.
    ///
    /// The per-component `COMPLETE`/`TRUNCATED` statuses are **read off the decode**, not
    /// recomputed: `arena.split` still holds them, and both decompositions enumerate in ascending
    /// root order, so component `c` is the same component in both. That correspondence is asserted
    /// rather than assumed.
    ///
    /// `histograms` and the tables are cleared and refilled with *this shot's* distributions; the
    /// campaign accumulator adds them up. Valid until the next `build_ball_graph`, which is why the
    /// driver calls this before anything that rebuilds `H` — `truncated_scheme_escalates`, in
    /// particular.
    void analyze_last_shot_components(
        BallProfile& prof,
        ComponentHistograms& histograms,
        ComponentStatusTable& size_x_status,
        ComponentStatusTable& hop_diameter_x_status);

    /// Retargets an already-constructed decoder at a different horizon, **without** recompiling the
    /// ball tables (§3.4 step 2: compile once at `T_max = max(T list)`, then sweep `T`).
    ///
    /// Legal exactly when the new `T` is within the compiled `T_max`, which is the same check
    /// construction makes; `H`'s `2 * T_int` / `T_int` filter and the truncation horizon both read
    /// `horizon`, so moving it is all a `T` sweep is. Rejected rather than clamped, so a run cannot
    /// report a `T` it did not decode at.
    void set_horizon(double T);

    void save_ball_artifact(const std::string& path) const;

    /// Scratch, so a steady-state shot allocates nothing.
    std::vector<uint64_t> seeded_scratch;
    std::vector<uint64_t> h_dets_scratch;
    mutable std::vector<uint64_t> sort_scratch;
    /// The component instance's match edges for the component being solved right now. Drained into
    /// `CommittedPair`s the moment that component's extraction returns, because the next component
    /// rebuilds the instance underneath the `DetectorNode*`s these hold.
    std::vector<pm::CompressedEdge> match_edge_scratch;
    /// The shot's committed pairs, in `G`'s detector ids, accumulated across the components.
    /// Cleared every shot, and cleared again on an escalating one, which discards Phase 1 in full.
    std::vector<CommittedPair> committed_pair_scratch;
    /// §2.5's scratch: the whole of `H`'s detection events, for the monolithic solve. Untouched
    /// unless `verify_component_decomposition` is on.
    std::vector<uint64_t> verify_dets;
    /// §M2.9.6 measurement 4, over the component being solved. Only touched when
    /// `collect_harvest_diagnostics` is set.
    TimelineDepthModel depth_model;

    void finish_construction(const char* ball_artifact_path);

    /// `harvest_component(mwpm, dets, status)` decides what to do with one solved component: the
    /// verification entry points always harvest it in full, the production ones extract a
    /// `COMPLETE` component and abandon a `TRUNCATED` one.
    ///
    /// `want_pairs` asks for the committed pairs in `G`'s ids; the caller's harvest lambda is the
    /// one that fills `match_edge_scratch`, and this maps and drains it per component.
    template <typename HarvestComponent>
    Phase1Outcome decode_impl(
        const std::vector<uint64_t>& dets,
        BallProfile* prof,
        bool want_pairs,
        const HarvestComponent& harvest_component);

    /// Solves the component sitting in `arena.split.sub` on `h_mwpm` and returns its status —
    /// so the caller's `build_component_subgraph` is what selects which component this is. Fills
    /// `harvest` through the caller's lambda; leaves the instance reset and ready for the next
    /// component either way.
    template <typename HarvestComponent>
    TimelineStatus solve_component(
        BallProfile* prof, const HarvestComponent& harvest_component, HarvestResult& harvest);

    /// Turns the component instance's `CompressedEdge`s into `CommittedPair`s over `G` — ball entry
    /// and all — and appends them to `committed_pair_scratch`. Called once per component, while
    /// `arena.split.sub` is still that component's sub-`H`.
    void drain_component_match_edges();

    /// §2.5. The monolithic solve on the whole of `H`, on `verify_mwpm`, compared against what the
    /// per-component path produced. Called only under `verify_component_decomposition`, and after
    /// the per-component path has finished with its own instance.
    ///
    /// `full_harvest` says whether the per-component side harvested every component in full, which
    /// is what makes the residual set and `num_trees` comparable. `committed_comparable` says
    /// whether it took the **obs** flavour: the match-edges flavour leaves
    /// `HarvestResult::committed` untouched by design — the pairs go to the caller's vector instead
    /// — so comparing it there would compare a real weight against an unfilled zero. On that
    /// flavour the weight check falls back to `dual_sum_at_truncation`, which both flavours do
    /// fill, and which is the same reduction over the same regions.
    void verify_decomposition(const Phase1Outcome& outcome, bool full_harvest, bool committed_comparable);
    /// Asserts §M2.6 level 1 against M1 on `G`. `actual_pairs` may be null when the caller took the
    /// obs flavour, in which case the committed *pair set* is not part of the comparison.
    void verify_level1(
        const std::vector<uint64_t>& dets,
        const HarvestResult& actual,
        const std::vector<CommittedPair>* actual_pairs,
        BallProfile* prof);
};

}  // namespace spec_matching
}  // namespace pm

#endif  // SPECMATCHING_SPEC_MATCHING_MANIFOLD_BALL_DECODING_H
