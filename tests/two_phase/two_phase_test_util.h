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

#ifndef PYREMATCHING_TESTS_TWO_PHASE_TEST_UTIL_H
#define PYREMATCHING_TESTS_TWO_PHASE_TEST_UTIL_H

#include <algorithm>
#include <cstdio>
#include <limits>
#include <queue>
#include <random>
#include <string>
#include <vector>

#include "pyrematching/sparse_blossom/driver/mwpm_decoding.h"
#include "pyrematching/sparse_blossom/driver/user_graph.h"
#include "pyrematching/two_phase/truncation/harvest.h"
#include "pyrematching/two_phase/truncation/truncated_timeline.h"
#include "stim.h"

namespace pm {
namespace two_phase {
namespace test {

const pm::weight_int NUM_DISTINCT_WEIGHTS = 10001;

inline std::string find_data_file(const char* name) {
    for (const std::string& directory : {"data/", "../data/", "../../data/"}) {
        std::string path = directory + name;
        FILE* f = fopen(path.c_str(), "r");
        if (f != nullptr) {
            fclose(f);
            return path;
        }
    }
    throw std::invalid_argument("Failed to find test data file " + std::string(name));
}

/// A detector error model plus a corpus of sampled syndromes to decode against it.
struct Corpus {
    std::string name;
    stim::DetectorErrorModel dem;
    std::vector<std::vector<uint64_t>> shots;

    pm::Mwpm to_mwpm(bool with_search_flooder = false) const {
        return pm::detector_error_model_to_mwpm(dem, NUM_DISTINCT_WEIGHTS, with_search_flooder);
    }
};

inline Corpus load_corpus(const char* name, const char* dem_file, const char* b8_file, size_t max_shots) {
    Corpus corpus;
    corpus.name = name;

    FILE* dem_handle = std::fopen(find_data_file(dem_file).c_str(), "r");
    if (dem_handle == nullptr)
        throw std::invalid_argument("Failed to open dem file");
    corpus.dem = stim::DetectorErrorModel::from_file(dem_handle);
    fclose(dem_handle);

    FILE* shots_handle = std::fopen(find_data_file(b8_file).c_str(), "r");
    if (shots_handle == nullptr)
        throw std::invalid_argument("Failed to open shots file");
    auto reader = stim::MeasureRecordReader<stim::MAX_BITWORD_WIDTH>::make(
        shots_handle,
        stim::SampleFormat::SAMPLE_FORMAT_B8,
        0,
        corpus.dem.count_detectors(),
        corpus.dem.count_observables());
    stim::SparseShot sparse_shot;
    while (corpus.shots.size() < max_shots && reader->start_and_read_entire_record(sparse_shot)) {
        corpus.shots.push_back(sparse_shot.hits);
        sparse_shot.clear();
    }
    fclose(shots_handle);
    return corpus;
}

inline Corpus load_surface_code_d13(size_t max_shots = 100) {
    return load_corpus(
        "surface_code_d13_p0.01",
        "surface_code_rotated_memory_x_13_0.01.dem",
        "surface_code_rotated_memory_x_13_0.01_1000_shots.b8",
        max_shots);
}

inline Corpus load_surface_code_d13_negative_weights(size_t max_shots = 100) {
    return load_corpus(
        "surface_code_d13_p0.01_negative_weights",
        "surface_code_rotated_memory_x_13_0.01_prob_0.2_negative.dem",
        "surface_code_rotated_memory_x_13_0.01_prob_0.2_negative_1000_shots.b8",
        max_shots);
}

inline Corpus load_toric_code_d5(size_t max_shots = 100) {
    return load_corpus(
        "toric_code_unrotated_memory_x_5_0.005",
        "toric_code_unrotated_memory_x_5_0.005.dem",
        "toric_code_unrotated_memory_x_5_0.005_1000.b8",
        max_shots);
}

/// Generates a surface code memory experiment and samples syndromes from it, so that tests are not
/// limited to the handful of DEMs checked into `data/`.
inline Corpus generate_surface_code_corpus(
    size_t distance, size_t rounds, double noise, size_t num_shots, uint64_t seed) {
    stim::CircuitGenParameters gen(rounds, distance, "rotated_memory_x");
    gen.after_clifford_depolarization = noise;
    gen.after_reset_flip_probability = noise;
    gen.before_measure_flip_probability = noise;
    stim::Circuit circuit = stim::generate_surface_code_circuit(gen).circuit;

    std::mt19937_64 rng(seed);
    size_t num_detectors = circuit.count_detectors();
    auto dets_obs = stim::sample_batch_detection_events<stim::MAX_BITWORD_WIDTH>(circuit, num_shots, rng);
    auto& dets = dets_obs.first;

    Corpus corpus;
    corpus.name = "generated_surface_code_d" + std::to_string(distance) + "_p" + std::to_string(noise);
    corpus.dem = stim::ErrorAnalyzer::circuit_to_detector_error_model(circuit, false, true, false, 0, false, false);
    for (size_t shot = 0; shot < num_shots; shot++) {
        corpus.shots.emplace_back();
        for (size_t d = 0; d < num_detectors; d++) {
            if (dets[d][shot])
                corpus.shots.back().push_back(d);
        }
    }
    return corpus;
}

/// The horizon, in flooder time units, corresponding to `num_edges` typical edges of the graph.
/// Uses the median discretised edge weight, which for a uniform-noise memory experiment is the
/// natural unit to express `T` in.
inline pm::cumulative_time_int median_edge_weight(const pm::MatchingGraph& graph) {
    std::vector<pm::weight_int> weights;
    for (const auto& node : graph.nodes) {
        for (size_t i = 0; i < node.neighbors.size(); i++)
            weights.push_back(node.neighbor_weights[i]);
    }
    if (weights.empty())
        return 2;
    std::sort(weights.begin(), weights.end());
    return (pm::cumulative_time_int)weights[weights.size() / 2];
}

/// One "edge weight" of the graph, expressed in the DEM float weight units that `BallParams` and
/// `BallConfig` take. Horizons are quoted in multiples of a lattice edge weight throughout the M1
/// results, and that is the normalisation under which every structure in the data is
/// `d`-independent — so it is the unit the ball tests size `T` and `R` in too.
inline double edge_weight_units(const pm::MatchingGraph& graph) {
    return (double)median_edge_weight(graph) / graph.normalising_constant;
}

/// Phase-1-only decode: truncated timeline, harvest, and the same negative-weight post-processing
/// that `pm::decode_detection_events` applies. At `pm::NO_HORIZON` this must be bit-identical to
/// the stock decoder.
struct Phase1Decode {
    TimelineStatus status;
    HarvestResult harvest;
    pm::obs_int obs_mask;
    pm::total_weight_int weight;
};

inline Phase1Decode phase1_decode_to_obs(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    pm::cumulative_time_int horizon,
    Harvester* harvester = nullptr) {
    Phase1Decode out;
    out.status = process_timeline_until_horizon(mwpm, detection_events, horizon);
    Harvester local;
    out.harvest = (harvester != nullptr ? *harvester : local).harvest_to_obs(mwpm, detection_events);
    out.obs_mask = out.harvest.committed.obs_mask ^ mwpm.flooder.negative_weight_obs_mask;
    out.weight = out.harvest.committed.weight + mwpm.flooder.negative_weight_sum;
    return out;
}

/// Exact single-phase reference for a syndrome, on a scratch `Mwpm` of the same graph.
inline pm::MatchingResult exact_decode(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    return pm::decode_detection_events_for_up_to_64_observables(mwpm, detection_events, false);
}

/// Shortest path distances from `src` over the matching graph, in flooder time units. Only for use
/// on graphs small enough that an all-pairs sweep is affordable.
inline std::vector<pm::total_weight_int> dijkstra(const pm::MatchingGraph& graph, size_t src) {
    const pm::total_weight_int infinity = std::numeric_limits<pm::total_weight_int>::max();
    std::vector<pm::total_weight_int> dist(graph.nodes.size(), infinity);
    using Entry = std::pair<pm::total_weight_int, size_t>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> frontier;
    dist[src] = 0;
    frontier.push({0, src});
    while (!frontier.empty()) {
        auto [d, u] = frontier.top();
        frontier.pop();
        if (d != dist[u])
            continue;
        const auto& node = graph.nodes[u];
        for (size_t i = 0; i < node.neighbors.size(); i++) {
            if (node.neighbors[i] == nullptr)
                continue;  // boundary half-edge
            size_t v = (size_t)(node.neighbors[i] - graph.nodes.data());
            pm::total_weight_int candidate = d + (pm::total_weight_int)node.neighbor_weights[i];
            if (candidate < dist[v]) {
                dist[v] = candidate;
                frontier.push({candidate, v});
            }
        }
    }
    return dist;
}

/// Cheapest path from `src` to the boundary, or -1 if the component has no boundary.
inline pm::total_weight_int boundary_cost(const pm::MatchingGraph& graph, size_t src) {
    auto dist = dijkstra(graph, src);
    const pm::total_weight_int infinity = std::numeric_limits<pm::total_weight_int>::max();
    pm::total_weight_int best = infinity;
    for (size_t v = 0; v < graph.nodes.size(); v++) {
        if (dist[v] == infinity)
            continue;
        const auto& node = graph.nodes[v];
        if (!node.neighbors.empty() && node.neighbors[0] == nullptr)
            best = std::min(best, dist[v] + (pm::total_weight_int)node.neighbor_weights[0]);
    }
    return best == infinity ? -1 : best;
}

/// How deeply blossoms are nested inside a region: 0 for a plain region, 1 for a blossom of plain
/// regions, 2 once a blossom contains a blossom, and so on.
inline int blossom_nesting_depth(const pm::GraphFillRegion& region) {
    if (region.blossom_children.empty())
        return 0;
    int deepest = 0;
    for (const auto& child : region.blossom_children)
        deepest = std::max(deepest, blossom_nesting_depth(*child.region));
    return deepest + 1;
}

/// The detection events a region owns, i.e. the sources of every leaf of its blossom nesting.
inline void collect_defects_in_region(
    const pm::GraphFillRegion& region, const pm::MatchingGraph& graph, std::vector<uint64_t>& out) {
    if (region.blossom_children.empty()) {
        if (!region.shell_area.empty())
            out.push_back((uint64_t)(region.shell_area[0] - graph.nodes.data()));
        return;
    }
    for (const auto& child : region.blossom_children)
        collect_defects_in_region(*child.region, graph, out);
}

/// The detection events that `begin_timeline` will actually create regions for: the shot's events
/// symmetric-differenced with the negative weight events, minus any boundary nodes.
inline std::vector<uint64_t> expected_seeded_detection_events(
    const pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    std::vector<uint64_t> negatives = mwpm.flooder.negative_weight_detection_events;
    std::vector<uint64_t> events = detection_events;
    std::sort(events.begin(), events.end());
    std::sort(negatives.begin(), negatives.end());
    std::vector<uint64_t> seeded;
    std::set_symmetric_difference(
        events.begin(), events.end(), negatives.begin(), negatives.end(), std::back_inserter(seeded));
    const auto& graph = mwpm.flooder.graph;
    seeded.erase(
        std::remove_if(
            seeded.begin(),
            seeded.end(),
            [&](uint64_t det) {
                return det < graph.is_user_graph_boundary_node.size() && graph.is_user_graph_boundary_node[det];
            }),
        seeded.end());
    return seeded;
}

}  // namespace test
}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TESTS_TWO_PHASE_TEST_UTIL_H
