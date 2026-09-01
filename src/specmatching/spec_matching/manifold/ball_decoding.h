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
    /// §C — the component and structural statistics of the shot's `H`: the connected components,
    /// their sizes, weighted diameters and boundary structure, the `H` edge-weight and boundary-cost
    /// distributions, and §A.3's classification of each component into "the trivial resolver would
    /// commit it", "it would leave a residual" and "it goes to the solver".
    ///
    /// **Profiling only, and untimed.** The analysis runs after the shot's timed window has closed
    /// and is charged to no latency number at all: the component work is assumed free
    /// (hardware-offloadable), so timing it would be measuring a stage that is not meant to be on
    /// this critical path. It reads the `H` the arena is still holding and changes nothing about
    /// the decode — with the flag on or off, every shot produces byte-identical output and the same
    /// set of shots escalates.
    ///
    /// Off by default, for the same reason `collect_structural_counters` is: the exit artifact
    /// collects it in its own pass rather than beside the timings it would otherwise perturb
    /// through the cache.
    bool collect_component_stats{false};

    /// §A — `k`: the largest connected component of `H` the small-component resolver takes. The
    /// solver sees only components of size `> k`; components of size `<= k` are resolved exactly,
    /// off the solver, before it runs.
    ///
    /// Components are disconnected by construction — no `H` edge crosses one, and the boundary is
    /// not a node, so it joins nothing — which means the timeline on `H` factorises over them.
    /// Removing a component the resolver has already settled therefore cannot change the solve on
    /// what is left. The resolver is exact for weight and observable at every `k <= 4`: it
    /// enumerates defect-pairings-with-boundary-fill under `H`'s own `<=` cutoffs and takes the
    /// minimum, breaking ties deterministically. Nothing routes to the solver except by being size
    /// `> k`; a component with no feasible matching becomes a residual and the shot escalates.
    ///
    /// `k = 0` resolves nothing and restores the un-pruned production path exactly, which is the
    /// oracle the A/B is run against. `k = 2` is the previous branch's production behaviour.
    /// **Capped at 4** — no resolver is defined above that — and rejected at construction above it.
    ///
    /// What moves with `k` is `blossom_on_h_ns` and `harvest_ns`, which run on a smaller node set —
    /// and, when every component resolves, do not run at all (§A.4).
    ///
    /// **Production path only.** The two verification entry points (`decode_phase1`,
    /// `decode_phase1_to_match_edges`) ignore this and always solve the whole of `H`: they are what
    /// §M2.6 level 1 and §M3.3 X8 compare against, and pruning what the oracle sees would defeat
    /// them.
    int prune_component_max_size{2};

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

/// What the small-component resolver did with the shot, in counts. Not a latency measurement and
/// not derived from one: these are the branch tallies the combine step of §A.5 needs, plus the two
/// harvest counters the off-solver commits belong in.
struct SmallCommitCounts {
    /// Components of size `<= k` with no feasible pairing-with-boundary-fill at all — an odd
    /// component with no legal boundary, or a "star" whose far members cannot pair. Non-zero forces
    /// escalation (§A.5).
    int residual{0};
    /// Committed pairs and boundary matches, mirrored into `HarvestResult`'s own counters so that
    /// the profile's commit tallies still add up to the shot's defect count.
    int pairs{0};
    int boundary{0};
    /// Defects the resolver settled off the solver this shot: the members of every committed
    /// component, over **all** sizes `1..k` and not merely the sizes above the previous branch's 2.
    int defects_resolved{0};
};

/// What §M3.4's production Phase 1 yielded.
///
/// `status` is the branch: `COMPLETE` means the harvest below is the answer, `TRUNCATED` means the
/// shot escalates and `harvest` is empty because it was never run. Reading the status rather than
/// `harvest.residual.empty()` is what keeps the two from drifting apart when the bypass is on.
struct Phase1Outcome {
    TimelineStatus status{TimelineStatus::COMPLETE};
    HarvestResult harvest;
};

/// Phase 1 executed on the defect manifold (§M2.5).
///
/// Per shot: negative-weight preamble on `G` -> intersect balls -> build `H` -> rebuild `Mwpm(H)`
/// -> `process_timeline_until_horizon(H, ..., T)` -> `harvest` (M1 code, unmodified) -> map `H`
/// indices back to detector ids.
///
/// Everything it emits is in `G`'s detector ids, so M3–M6 consume it unchanged: nothing downstream
/// of harvest knows `H` exists.
struct BallDecoder {
    BallConfig config;
    /// The detector graph. Kept because the negative-weight preamble, the boundary-node mask and
    /// the oracle path all live on it.
    pm::Mwpm g_mwpm;
    BallTables tables;
    BallMwpm h_mwpm;
    BallGraphArena arena;
    Harvester harvester;

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

    /// §C. Decomposes the `H` the arena is still holding — the one the shot just decoded on — into
    /// connected components, and fills `prof.components` and `histograms` from it.
    ///
    /// **Call it after the shot's timed window has closed, never inside one** (hard constraint 1).
    /// It starts no timer of its own and must not be wrapped in one: this experiment's latency
    /// account is the solver and the harvest, and the component work is assumed free.
    ///
    /// `histograms` is cleared and refilled with *this shot's* distributions; the campaign
    /// accumulator adds them up. Valid until the next `build_ball_graph`, which is why the driver
    /// calls this before anything that rebuilds `H` — `truncated_scheme_escalates`, in particular.
    void analyze_last_shot_components(BallProfile& prof, ComponentHistograms& histograms);

    void save_ball_artifact(const std::string& path) const;

    /// Scratch, so a steady-state shot allocates nothing.
    std::vector<uint64_t> seeded_scratch;
    std::vector<uint64_t> h_dets_scratch;
    mutable std::vector<uint64_t> sort_scratch;
    std::vector<pm::CompressedEdge> match_edge_scratch;
    /// §A's committed pairs, in `G`'s detector ids — the small-component resolver's contribution to
    /// the match-edge flavours. Cleared every shot, and cleared again on an escalating one, where
    /// §A.5 discards Phase 1 in full.
    std::vector<CommittedPair> resolved_pairs;
    /// The graph the last solve actually ran on: the whole of `H`, or §A.4's sub-`H`. Everything
    /// that maps a solver index back to a detector id reads it, so that the two cases go through
    /// one path.
    const BallGraph* solved_graph{nullptr};
    /// §M2.9.6 measurement 4, over `H`. Only touched when `collect_harvest_diagnostics` is set.
    TimelineDepthModel depth_model;

    void finish_construction(const char* ball_artifact_path);
    /// `harvest_on_h(mwpm, h_dets, status)` decides what to do with the solved timeline: the
    /// verification entry points always harvest in full, the production ones branch on the status.
    ///
    /// `allow_prune` is what keeps §A off the verification path: only the production entry points
    /// pass true, and even they defer to `BallConfig::prune_component_max_size`.
    template <typename HarvestOnH>
    Phase1Outcome decode_impl(
        const std::vector<uint64_t>& dets, BallProfile* prof, bool allow_prune, const HarvestOnH& harvest_on_h);
    /// §A's small-component resolver, run over every component of `h` of size `<= k`, and §A.4's
    /// induced sub-graph over what is left.
    ///
    /// Returns the graph the solver should be handed: the sub-`H` over the SOLVER set, or `h`
    /// itself when nothing was resolved away (in which case no copy is made). Accumulates the
    /// committed observables and weight into `resolved`, the committed pairs into `resolved_pairs`,
    /// and counts the components with no feasible matching.
    ///
    /// **Untimed by construction** — no timer is started here and none may be added. It is a serial
    /// pre-pass on the critical path, outside this branch's reported latency by scope.
    const BallGraph& resolve_small_components(
        const BallGraph& h, pm::MatchingResult& resolved, SmallCommitCounts& counts);
    /// Turns harvest's `CompressedEdge`s over the solved graph into `CommittedPair`s over `G`, ball
    /// entry and all, and merges in §A's trivially committed pairs. Shared by the verification and
    /// production match-edge entry points.
    void map_match_edges_to_committed_pairs(std::vector<CommittedPair>& committed_pairs) const;
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
