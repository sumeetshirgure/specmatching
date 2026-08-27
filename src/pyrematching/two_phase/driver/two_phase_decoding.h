// Copyright 2026 PyReMatching contributors
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

#ifndef PYREMATCHING_TWO_PHASE_DRIVER_TWO_PHASE_DECODING_H
#define PYREMATCHING_TWO_PHASE_DRIVER_TWO_PHASE_DECODING_H

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "pyrematching/two_phase/escalation/escalate.h"
#include "pyrematching/two_phase/manifold/ball_decoding.h"
#include "pyrematching/two_phase/perf/ball_profile.h"
#include "pyrematching/two_phase/perf/two_phase_profile.h"
#include "stim.h"

namespace pm {
namespace two_phase {

struct TwoPhaseConfig {
    /// The horizon, in DEM float weight units. A truncation horizon on the landed path; a pure
    /// runtime scalar — ball filter plus escalation threshold, never a flooder field — under
    /// `stock_on_h` (§M7.1).
    double T{2.0};
    BallParams ball;

    /// §M7 — swap Phase 1's front end for **stock** sparse blossom on `H` plus a max-dual
    /// certificate. Requires `phase1_on_ball_graph`. See `BallConfig::stock_on_h`.
    ///
    /// The output is exact MWPM on every shot either way. What this buys is a cleaner exact
    /// certificate — direct equality against stock-on-`G` rather than bit-exactness against a
    /// truncated harvest — plus the removed gating and harvest overhead, and an escalation rate
    /// that is no higher than the truncated scheme's (§M7.0's corollary). It is **not** expected to
    /// move the headline speedup: the common path does the same blossom work in both schemes.
    bool stock_on_h{false};

    /// §M7.7 benchmark mode: replay the landed truncated scheme's escalation decision on every shot
    /// so that `q_current_on_same_corpus` — and therefore `q_current - q_this` — is measured on
    /// identical shots rather than compared across campaigns. Only meaningful with `stock_on_h`.
    ///
    /// **Read `q` from such a campaign, never a time.** The replay costs a second full Phase 1 on
    /// `H`, and it tears its instance down with `Mwpm::reset`, which frees the arena pools — so the
    /// *next* shot re-allocates them and its `total_ns` is not a steady-state number. The replay
    /// itself runs outside the timed window, but that does not undo what it leaves behind. §M7.8's
    /// speedup read therefore runs two decoders alternately rather than turning this on.
    bool measure_truncated_reference{false};

    /// Run Phase 1 on the defect manifold `H` (§M2). `false` runs M1's decode on `G` instead and
    /// compiles no ball tables at all — the oracle path, and the only one that supports
    /// `T = infinity`, since with a finite `R` the tables cannot supply the edges an unbounded run
    /// would need (§M4.1).
    bool phase1_on_ball_graph{true};

    /// `T = infinity`: run the timeline to completion. Only legal with
    /// `phase1_on_ball_graph = false`; the constructor rejects the combination rather than letting
    /// invariant 6 fail first.
    bool unbounded_horizon{false};

    /// §M3.4. `false` — the production path — is extract-only on a completed timeline and runs no
    /// harvest at all on a truncated one, because that shot escalates and Phase 1's partial result
    /// is discarded. `true` runs §M1.3's full harvest on every shot; that is the verification path
    /// §M3.3 X8 compares against, and it is what §M2.6 level 1 and debug invariants 3, 4, 18 and 19
    /// read. Forced on by `verify_against_g`, which has nothing to compare without it.
    bool full_harvest_for_verification{false};

    /// Build the search graph on `G`. Required by `decode_to_edges` in both branches: stock's
    /// `decode_detection_events_to_edges` needs it on the escalating path, and the
    /// negative-weight edge list it reads lives there.
    bool edges_flavor{false};

    BallGraphBuildMode mode{BallGraphBuildMode::SCAN};
    size_t compile_threads{0};

    /// §M2.5's inline level-1 check against M1 on `G`, for CI and benchmarks. Off by default
    /// because it costs more than the decode it is checking.
    bool verify_against_g{false};
    bool collect_harvest_diagnostics{false};
    bool collect_structural_counters{false};
    /// §C — the component structure of every shot's `H`, and §A.3's classification of it. Profiling
    /// only, computed after the shot's timed window has closed, and charged to no latency number:
    /// see `BallConfig::collect_component_stats`. Ignored without `phase1_on_ball_graph`, where
    /// there is no `H` to decompose.
    bool collect_component_stats{false};

    /// §A — `k`: resolve `H`'s connected components of size `<= k` exactly, directly off the ball
    /// tables, and run the solver on the size-`> k` remainder only. Capped at 4 and rejected above
    /// it; `k = 0` is the un-pruned path. See `BallConfig::prune_component_max_size`.
    ///
    /// Ignored without `phase1_on_ball_graph` (there is no `H` to decompose) and on the
    /// verification path, which always solves the whole of `H`.
    int prune_component_max_size{2};
    /// §B — skip the §M2.1 negative-weight preamble on an all-positive DEM, where it is provably a
    /// no-op. See `BallConfig::skip_negative_weight_preamble_when_positive`.
    bool skip_negative_weight_preamble_when_positive{true};

    /// Benchmark mode: additionally run stock exact decode on every shot and record it as
    /// `exact_reference_ns`, so `speedup_vs_stock` is a paired measurement rather than two
    /// campaigns divided.
    bool measure_exact_reference{false};
    /// §M6.4 timer discipline: probe the thread once per shot and mark shots the scheduler
    /// interfered with. One preemption defines p999 outright at 2000 shots. The instrument is
    /// `getrusage(RUSAGE_THREAD)`'s switch counters where they exist and the wall-minus-thread-CPU
    /// gap where they do not (macOS); `PreemptionProbe::mechanism_name()` says which, and a campaign
    /// should record it beside the rate, since "no shot was flagged" and "no shot could be flagged"
    /// are different statements.
    bool detect_preemption{false};
};

/// M4/M5 — the end-to-end two-phase decoder.
///
/// Per shot (§M4.1):
///
/// ```
/// Phase 1                          # §M2.5 on H, or M1's decode on G in the oracle configuration
/// if the residual is empty:        # the common case, q ~ 1e-6 .. 3e-4 (§M3.0)
///     emit Phase 1's committed observables and weight
/// else:
///     escalate_to_stock(...)       # §M3.1: re-decode the whole shot on G, discard Phase 1
/// ```
///
/// The output is **exact MWPM on every shot**, by construction: Phase 1 is exact within the
/// horizon, and anything it could not finish is handed to stock in full. There is no approximation
/// budget, no `eps_total`, and nothing to sweep.
///
/// The empty-residual path costs nothing beyond Phase 1: no `Mwpm(G)` touch, no allocation, no
/// reset (§M3.3 X6).
struct TwoPhaseDecoder {
    TwoPhaseConfig config;

    /// Phase 1 on `H`, and the owner of the one `pm::Mwpm` on `G`. Null in the oracle
    /// configuration, where `oracle_g_mwpm` holds that instance instead — either way there is
    /// exactly one, and `g_mwpm()` names it.
    std::unique_ptr<BallDecoder> ball;
    pm::Mwpm oracle_g_mwpm;
    Harvester oracle_harvester;

    /// `T` in the flooder's time units, converted once (§0 unit rule). `pm::NO_HORIZON` when
    /// `config.unbounded_horizon`.
    horizon_int horizon{0};
    size_t num_observables{0};

    /// Accumulated over a campaign whenever a profile is passed to `decode_batch` (§M6.2).
    TwoPhaseAggregateStats stats;

    static TwoPhaseDecoder from_detector_error_model(
        const stim::DetectorErrorModel& dem,
        TwoPhaseConfig config,
        pm::weight_int num_distinct_weights = pm::NUM_DISTINCT_WEIGHTS,
        const char* ball_artifact_path = nullptr);

    /// The single `pm::Mwpm` on `G`: the negative-weight preamble, the escalation fallback and the
    /// oracle front end all run on it, at disjoint times.
    pm::Mwpm& g_mwpm();
    const pm::Mwpm& g_mwpm() const;

    /// The one constant that converts weights to time units (§0's unit rule). Every weight this
    /// decoder reports is an integer in the discretised metric; divide by this to get the DEM float
    /// units the inherited `Matching.decode` reports.
    inline double normalising_constant() const {
        return g_mwpm().flooder.graph.normalising_constant;
    }

    /// Obs flavour. `obs` must point at `num_observables` bytes; it and `weight` are **overwritten**.
    /// Works for any number of observables: above 64 the `obs_int` mask is unusable, so the
    /// committed pairs are read out through the ball tables' observable id lists instead (§0).
    void decode_to_obs(
        const std::vector<uint64_t>& dets, uint8_t* obs, pm::total_weight_int& weight, TwoPhaseProfile* prof = nullptr);

    /// Edges flavour (M5). `edges` is cleared and refilled with detector id pairs, `-1` for the
    /// boundary — the inherited `decode_detection_events_to_edges` contract, so downstream tooling
    /// is untouched. Treat it as a set: the ball front end emits it sorted, escalation emits it in
    /// stock's own flip order.
    ///
    /// `weight` receives Phase 1's committed weight, or stock's on an escalating shot.
    void decode_to_edges(
        const std::vector<uint64_t>& dets,
        std::vector<int64_t>& edges,
        pm::total_weight_int& weight,
        TwoPhaseProfile* prof = nullptr);

    /// Decodes a batch, accumulating into `stats` when `profile` is set. `obs_out`, when non-null,
    /// receives `shots.size() * num_observables` bytes; `weights_out`, when non-null, one weight
    /// per shot; `profiles_out`, when non-null, one profile per shot.
    void decode_batch(
        const std::vector<std::vector<uint64_t>>& shots,
        uint8_t* obs_out,
        pm::total_weight_int* weights_out,
        bool profile,
        std::vector<TwoPhaseProfile>* profiles_out = nullptr);

    void save_ball_artifact(const std::string& path) const;

    /// Scratch, so a steady-state shot allocates nothing.
    std::vector<CommittedPair> committed_pairs;
    std::vector<pm::CompressedEdge> oracle_match_edges;
    /// Somewhere for `extract_paths_from_match_edges` to XOR observables it crosses while computing
    /// a weight. Discarded; the caller already has the observables by another route.
    std::vector<uint8_t> obs_sink;
    BallProfile ball_profile;
    /// §C.2's distributions for the shot just decoded. Refilled per shot and folded into `stats` by
    /// `decode_batch`; a driver that decodes shot by shot instead reads it here and calls
    /// `TwoPhaseAggregateStats::accumulate_component_histograms` itself.
    ComponentHistograms component_histograms;

    /// The matching weight of a set of matched detection-event pairs: the sum of the shortest-path
    /// weights between them plus the graph's negative-weight offset.
    ///
    /// **Not** the weight of the correction's edge set, which is smaller whenever two paths share a
    /// `G` edge and cancel it out. The obs flavour reports the matching weight, so the edges flavour
    /// reports it too; stock sidesteps the question by reporting no weight from `to_edges` at all.
    pm::total_weight_int weight_of_match_edges(pm::Mwpm& mwpm, const std::vector<pm::CompressedEdge>& match_edges);

   private:
    /// Runs Phase 1 and returns the branch. On `COMPLETE`, `committed_pairs` is filled iff the
    /// caller asked for pairs (the edges flavour, or the wide-observable obs flavour).
    Phase1Outcome run_phase1(const std::vector<uint64_t>& dets, bool need_pairs, TwoPhaseProfile* prof);
    /// The committed observables and weight of a completed shot, in the wide-observable case, read
    /// off the ball tables' observable id lists. `obs` and `weight` are overwritten.
    void obs_from_committed_pairs(uint8_t* obs, pm::total_weight_int& weight) const;
    void copy_phase1_stats(const Phase1Outcome& outcome, TwoPhaseProfile* prof) const;
    /// §C. Decomposes the shot's `H` and mirrors the result onto `prof`. Runs **outside** the timed
    /// window, after `total_ns` has been read, and before anything that rebuilds `H` — the analysis
    /// reads the graph the arena is still holding.
    void record_component_stats(TwoPhaseProfile* prof);
    /// §M7.7's `q_current_on_same_corpus`. Runs **outside** the shot's timed window, after
    /// `total_ns` has been read, because it is a second full Phase 1 on `H`.
    void record_truncated_reference(const std::vector<uint64_t>& dets, TwoPhaseProfile* prof);
};

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_DRIVER_TWO_PHASE_DECODING_H
