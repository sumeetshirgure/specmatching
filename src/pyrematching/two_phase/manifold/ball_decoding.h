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

#ifndef PYREMATCHING_TWO_PHASE_MANIFOLD_BALL_DECODING_H
#define PYREMATCHING_TWO_PHASE_MANIFOLD_BALL_DECODING_H

#include <string>
#include <vector>

#include "pyrematching/sparse_blossom/driver/mwpm_decoding.h"
#include "pyrematching/sparse_blossom/driver/user_graph.h"
#include "pyrematching/two_phase/manifold/ball_graph.h"
#include "pyrematching/two_phase/manifold/ball_mwpm.h"
#include "pyrematching/two_phase/manifold/ball_tables.h"
#include "pyrematching/two_phase/perf/ball_profile.h"
#include "pyrematching/two_phase/truncation/harvest.h"
#include "pyrematching/two_phase/truncation/truncated_timeline.h"
#include "stim.h"

namespace pm {
namespace two_phase {

struct BallConfig {
    /// The truncation horizon, in DEM float weight units. Must be `<= ball.T_max`.
    double T{2.0};
    BallParams ball;
    /// Run M1's path on `G` for every shot and assert §M2.6 level 1 inline. CI and benchmarks only;
    /// off by default because it costs more than the decode it is checking.
    bool verify_against_g{false};
    BallGraphBuildMode mode{BallGraphBuildMode::SCAN};
    /// Threads used to compile the ball tables. 0 means hardware concurrency.
    size_t compile_threads{0};
};

/// One committed pair, in `G`'s detector ids. `to == -1` means matched to the boundary. This is the
/// form §M2.6 level 1 compares and the form M5 will lift.
struct CommittedPair {
    int64_t from;
    int64_t to;
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

    static BallDecoder from_detector_error_model(
        const stim::DetectorErrorModel& dem,
        BallConfig config,
        pm::weight_int num_distinct_weights = pm::NUM_DISTINCT_WEIGHTS,
        const char* ball_artifact_path = nullptr);

    /// Builds a decoder around an already-constructed `pm::Mwpm` for `G`. The `Mwpm` is moved in.
    static BallDecoder from_mwpm(pm::Mwpm g_mwpm, BallConfig config, const char* ball_artifact_path = nullptr);

    /// Same contract and same struct as M1's harvest, in `G`'s detector ids.
    HarvestResult decode_phase1(const std::vector<uint64_t>& dets, BallProfile* prof = nullptr);

    /// As `decode_phase1`, and additionally the committed pairs in `G`'s detector ids. Works for
    /// any number of observables, which the obs flavour does not.
    HarvestResult decode_phase1_to_match_edges(
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

    void save_ball_artifact(const std::string& path) const;

    /// Scratch, so a steady-state shot allocates nothing.
    std::vector<uint64_t> seeded_scratch;
    std::vector<uint64_t> h_dets_scratch;
    mutable std::vector<uint64_t> sort_scratch;
    std::vector<pm::CompressedEdge> match_edge_scratch;

    void finish_construction(const char* ball_artifact_path);
    template <typename HarvestOnH>
    HarvestResult decode_impl(const std::vector<uint64_t>& dets, BallProfile* prof, const HarvestOnH& harvest_on_h);
    /// Asserts §M2.6 level 1 against M1 on `G`. `actual_pairs` may be null when the caller took the
    /// obs flavour, in which case the committed *pair set* is not part of the comparison.
    void verify_level1(
        const std::vector<uint64_t>& dets,
        const HarvestResult& actual,
        const std::vector<CommittedPair>* actual_pairs,
        BallProfile* prof);
};

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_MANIFOLD_BALL_DECODING_H
