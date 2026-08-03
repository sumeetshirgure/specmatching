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

/// The M1 exit artifact.
///
/// Sweeps the truncation horizon `T` across a set of code distances and error rates, and reports —
/// per configuration — the residual size distribution, the fallback rate `q`, the commit
/// composition, the largest surviving tree, the exposed-root-blossom rate, Phase-1 wall time
/// against a full exact decode, and `max_S y_S` against the cluster weight-diameter.
///
/// The go/no-go read at the end of M1 is made from `q` and the fallback cost ratio: if
/// `q * C_phase2` is not a meaningful share of the mean, the honest claim is a worst-case bound
/// rather than a speedup.
///
/// Usage:
///   two_phase_m1_artifact [--distances 11,13,15] [--error-rates 0.001,0.003]
///                         [--shots 500] [--horizon-multiples 1,2,4,8]
///                         [--diameter-shots 20] [--csv path]

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "pyrematching/sparse_blossom/driver/mwpm_decoding.h"
#include "pyrematching/sparse_blossom/driver/user_graph.h"
#include "pyrematching/sparse_blossom/matcher/alternating_tree.h"
#include "pyrematching/two_phase/perf/two_phase_profile.h"
#include "pyrematching/two_phase/truncation/harvest.h"
#include "pyrematching/two_phase/truncation/truncated_timeline.h"
#include "stim.h"

using namespace pm::two_phase;

namespace {

const pm::weight_int NUM_DISTINCT_WEIGHTS = 10001;

struct Options {
    std::vector<size_t> distances = {11, 13, 15, 17, 19, 21};
    std::vector<double> error_rates = {0.001, 0.003, 0.005};
    std::vector<double> horizon_multiples = {0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 8.0, 16.0};
    size_t shots = 500;
    size_t diameter_shots = 20;
    /// Cluster weight-diameters cost one Dijkstra per cluster member; this caps the total per
    /// configuration. Whatever the cap skips is reported, never silently dropped.
    size_t dijkstra_budget = 3000;
    uint64_t seed = 20260803;
    std::string csv_path;
};

template <typename T>
std::vector<T> parse_list(const std::string& text) {
    std::vector<T> out;
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        if (item.empty())
            continue;
        std::stringstream item_stream(item);
        T value;
        item_stream >> value;
        out.push_back(value);
    }
    return out;
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; i++) {
        std::string flag = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc)
                throw std::invalid_argument("missing value for " + flag);
            return argv[++i];
        };
        if (flag == "--distances") {
            options.distances = parse_list<size_t>(next());
        } else if (flag == "--error-rates") {
            options.error_rates = parse_list<double>(next());
        } else if (flag == "--horizon-multiples") {
            options.horizon_multiples = parse_list<double>(next());
        } else if (flag == "--shots") {
            options.shots = std::stoul(next());
        } else if (flag == "--diameter-shots") {
            options.diameter_shots = std::stoul(next());
        } else if (flag == "--dijkstra-budget") {
            options.dijkstra_budget = std::stoul(next());
        } else if (flag == "--seed") {
            options.seed = std::stoull(next());
        } else if (flag == "--csv") {
            options.csv_path = next();
        } else {
            throw std::invalid_argument("unrecognised flag " + flag);
        }
    }
    return options;
}

struct Experiment {
    stim::DetectorErrorModel dem;
    std::vector<std::vector<uint64_t>> shots;
};

Experiment generate(size_t distance, double noise, size_t num_shots, uint64_t seed) {
    stim::CircuitGenParameters gen(distance, distance, "rotated_memory_x");
    gen.after_clifford_depolarization = noise;
    gen.after_reset_flip_probability = noise;
    gen.before_measure_flip_probability = noise;
    stim::Circuit circuit = stim::generate_surface_code_circuit(gen).circuit;

    std::mt19937_64 rng(seed);
    size_t num_detectors = circuit.count_detectors();
    auto dets_obs = stim::sample_batch_detection_events<stim::MAX_BITWORD_WIDTH>(circuit, num_shots, rng);
    auto& dets = dets_obs.first;

    Experiment experiment;
    experiment.dem = stim::ErrorAnalyzer::circuit_to_detector_error_model(circuit, false, true, false, 0, false, false);
    for (size_t shot = 0; shot < num_shots; shot++) {
        experiment.shots.emplace_back();
        for (size_t d = 0; d < num_detectors; d++) {
            if (dets[d][shot])
                experiment.shots.back().push_back(d);
        }
    }
    return experiment;
}

pm::cumulative_time_int median_edge_weight(const pm::MatchingGraph& graph) {
    std::vector<pm::weight_int> weights;
    for (const auto& node : graph.nodes) {
        for (auto weight : node.neighbor_weights)
            weights.push_back(weight);
    }
    if (weights.empty())
        return 2;
    std::sort(weights.begin(), weights.end());
    return (pm::cumulative_time_int)weights[weights.size() / 2];
}

std::vector<pm::total_weight_int> dijkstra(const pm::MatchingGraph& graph, size_t src) {
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
                continue;
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

void collect_defects(const pm::GraphFillRegion& region, const pm::MatchingGraph& graph, std::vector<uint64_t>& out) {
    if (region.blossom_children.empty()) {
        if (!region.shell_area.empty())
            out.push_back((uint64_t)(region.shell_area[0] - graph.nodes.data()));
        return;
    }
    for (const auto& child : region.blossom_children)
        collect_defects(*child.region, graph, out);
}

/// The largest pairwise distance inside any region alive at truncation. This is the quantity
/// `max_S y_S` is meant to be read against: a dual that grows past a constant fraction of the
/// cluster it lives in is the signal that truncation is cutting into real structure.
///
/// Each measured cluster costs one Dijkstra per member, so the caller passes a budget. What the
/// budget skipped is counted and reported rather than silently dropped.
pm::total_weight_int cluster_weight_diameter(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    size_t& dijkstra_budget,
    uint64_t& clusters_measured,
    uint64_t& clusters_skipped) {
    // A cluster is a group of defects that actually interacted: a matched pair (each side possibly
    // a blossom), or a whole surviving alternating tree. Keyed by the matched partner with the
    // lower address, or by the tree root, so both sides land in the same group.
    std::map<const void*, std::vector<uint64_t>> clusters;
    std::vector<pm::GraphFillRegion*> tops;
    for (uint64_t det : detection_events) {
        if (det >= mwpm.flooder.graph.nodes.size())
            continue;
        auto& node = mwpm.flooder.graph.nodes[det];
        if (node.region_that_arrived == nullptr)
            continue;
        pm::GraphFillRegion* top = node.region_that_arrived_top;
        if (std::find(tops.begin(), tops.end(), top) != tops.end())
            continue;
        tops.push_back(top);

        const void* key = top;
        if (top->alt_tree_node != nullptr) {
            key = top->alt_tree_node->find_root();
        } else if (top->match.region != nullptr) {
            key = std::min<const void*>(top, top->match.region);
        }
        collect_defects(*top, mwpm.flooder.graph, clusters[key]);
    }

    pm::total_weight_int diameter = 0;
    for (auto& [key, defects] : clusters) {
        (void)key;
        if (defects.size() < 2)
            continue;
        if (defects.size() > dijkstra_budget) {
            clusters_skipped++;
            continue;
        }
        dijkstra_budget -= defects.size();
        clusters_measured++;
        for (size_t i = 0; i < defects.size(); i++) {
            auto distances = dijkstra(mwpm.flooder.graph, defects[i]);
            for (size_t j = i + 1; j < defects.size(); j++)
                diameter = std::max(diameter, distances[defects[j]]);
        }
    }
    return diameter;
}

struct ConfigResult {
    size_t distance;
    double error_rate;
    double horizon_multiple;
    pm::cumulative_time_int horizon;
    TwoPhaseAggregateStats stats;
    TwoPhaseSummary summary;
    pm::total_weight_int max_region_dual{0};
    pm::total_weight_int max_cluster_diameter{0};
    double mean_max_region_dual{0};
    uint64_t clusters_measured{0};
    uint64_t clusters_skipped{0};
    /// The same pipeline at `T = infinity` on the same shots. Isolates the effect of truncation
    /// from the fixed cost of harvesting, which the stock decoder does not pay.
    long long sum_untruncated_pipeline_ns{0};
};

double percentile(std::vector<long long> values, double fraction) {
    if (values.empty())
        return 0;
    std::sort(values.begin(), values.end());
    size_t index = (size_t)(fraction * (double)(values.size() - 1));
    return (double)values[std::min(index, values.size() - 1)];
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    try {
        options = parse_options(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return 1;
    }

    std::vector<ConfigResult> results;

    for (size_t distance : options.distances) {
        for (double error_rate : options.error_rates) {
            std::cerr << "generating d=" << distance << " p=" << error_rate << " (" << options.shots << " shots)\n";
            Experiment experiment = generate(distance, error_rate, options.shots, options.seed + distance);
            auto mwpm = pm::detector_error_model_to_mwpm(experiment.dem, NUM_DISTINCT_WEIGHTS);
            auto reference_mwpm = pm::detector_error_model_to_mwpm(experiment.dem, NUM_DISTINCT_WEIGHTS);
            auto unit = median_edge_weight(mwpm.flooder.graph);
            Harvester harvester;

            for (double multiple : options.horizon_multiples) {
                ConfigResult result;
                result.distance = distance;
                result.error_rate = error_rate;
                result.horizon_multiple = multiple;
                // Keep the horizon even, like every discretised edge weight, so that it lands on
                // the lattice collisions happen on.
                result.horizon = 2 * (pm::cumulative_time_int)((double)unit * multiple / 2);
                result.stats.keep_per_shot = true;

                // Untimed pre-pass: the cluster weight-diameters need the regions to still exist,
                // so they have to run between the timeline and the harvest, and each one costs a
                // pile of Dijkstras. Doing them in their own pass keeps them off the clock.
                size_t dijkstra_budget = options.dijkstra_budget;
                for (size_t shot_index = 0; shot_index < options.diameter_shots && shot_index < experiment.shots.size();
                     shot_index++) {
                    const auto& shot = experiment.shots[shot_index];
                    process_timeline_until_horizon(mwpm, shot, result.horizon);
                    result.max_cluster_diameter = std::max(
                        result.max_cluster_diameter,
                        cluster_weight_diameter(
                            mwpm, shot, dijkstra_budget, result.clusters_measured, result.clusters_skipped));
                    harvester.harvest_to_obs(mwpm, shot);
                }

                pm::total_weight_int sum_max_region_dual = 0;
                for (size_t shot_index = 0; shot_index < experiment.shots.size(); shot_index++) {
                    const auto& shot = experiment.shots[shot_index];
                    TwoPhaseProfile profile;
                    profile.num_defects = (int)shot.size();

                    HiResTimer total_timer;
                    HiResTimer phase_timer;
                    total_timer.start();

                    phase_timer.start();
                    process_timeline_until_horizon(mwpm, shot, result.horizon);
                    profile.phase1_ns = phase_timer.elapsed_ns();

                    phase_timer.start();
                    auto harvest = harvester.harvest_to_obs(mwpm, shot);
                    profile.harvest_ns = phase_timer.elapsed_ns();
                    profile.total_ns = total_timer.elapsed_ns();
                    profile.fill_from(harvest);
                    profile.weight_out = harvest.committed.weight;

                    // Stock exact decode, and the same pipeline without truncation. The first says
                    // what the two-phase decoder has to beat; the second separates "truncation did
                    // not help" from "harvesting costs something the stock path never pays".
                    HiResTimer exact_timer;
                    exact_timer.start();
                    pm::decode_detection_events_for_up_to_64_observables(reference_mwpm, shot, false);
                    profile.exact_reference_ns = exact_timer.elapsed_ns();

                    HiResTimer untruncated_timer;
                    untruncated_timer.start();
                    process_timeline_until_horizon(mwpm, shot, pm::NO_HORIZON);
                    harvester.harvest_to_obs(mwpm, shot);
                    result.sum_untruncated_pipeline_ns += untruncated_timer.elapsed_ns();

                    result.stats.accumulate(profile);
                    result.max_region_dual = std::max(result.max_region_dual, harvest.max_region_dual);
                    sum_max_region_dual += harvest.max_region_dual;
                }

                result.mean_max_region_dual =
                    (double)sum_max_region_dual / (double)std::max<size_t>(1, experiment.shots.size());
                result.summary = summarize(result.stats);
                results.push_back(result);

                std::cerr << "  T=" << result.horizon << " (" << multiple << " edges)"
                          << "  q=" << result.summary.q << "  mean_residual=" << result.summary.mean_residual_size
                          << "  phase1_ns=" << result.summary.c_phase1
                          << "  exact_ns=" << (double)result.stats.sum_exact_reference_ns / (double)result.stats.shots
                          << "\n";
            }
        }
    }

    std::ostringstream csv;
    csv << "distance,error_rate,horizon_multiple,horizon,shots,mean_defects,"
           "q_fallback_rate,mean_residual_size,max_residual_size,mean_residual_density,"
           "committed_pairs_frozen_per_shot,committed_pairs_tree_per_shot,committed_boundary_per_shot,"
           "max_largest_tree_size,exposed_root_blossom_rate,exposed_root_blossoms_per_shot,"
           "phase1_only_ns,harvest_ns,c_phase1_ns,c_phase2_ns,amortised_mean_ns,measured_mean_ns,"
           "fallback_cost_ratio,untruncated_pipeline_ns,"
           "exact_mean_ns,p50_ns,p99_ns,p999_ns,exact_p99_ns,speedup_vs_stock,"
           "max_region_dual,mean_max_region_dual,max_cluster_weight_diameter,"
           "clusters_measured,clusters_skipped_for_budget,residual_hist\n";
    for (const auto& result : results) {
        const auto& stats = result.stats;
        const auto& summary = result.summary;
        double shots = (double)stats.shots;
        csv << result.distance << "," << result.error_rate << "," << result.horizon_multiple << "," << result.horizon
            << "," << stats.shots << "," << (double)stats.sum_num_defects / shots << "," << summary.q << ","
            << summary.mean_residual_size << "," << stats.max_residual_size << "," << summary.mean_residual_density
            << "," << (double)stats.sum_committed_pairs_frozen / shots << ","
            << (double)stats.sum_committed_pairs_tree / shots << "," << (double)stats.sum_committed_boundary / shots
            << "," << stats.max_largest_tree_size << "," << summary.exposed_root_blossom_rate << ","
            << (double)stats.sum_exposed_root_blossoms / shots << "," << (double)stats.sum_phase1_ns / shots << ","
            << (double)stats.sum_harvest_ns / shots << "," << summary.c_phase1 << "," << summary.c_phase2 << ","
            << summary.amortised_mean_ns << "," << summary.measured_mean_ns << "," << summary.fallback_cost_ratio << ","
            << (double)result.sum_untruncated_pipeline_ns / shots << "," << (double)stats.sum_exact_reference_ns / shots
            << "," << summary.p50_total_ns << "," << summary.p99_total_ns << "," << summary.p999_total_ns << ","
            << percentile(stats.per_shot_exact_ns, 0.99) << "," << summary.speedup_vs_stock << ","
            << result.max_region_dual << "," << result.mean_max_region_dual << "," << result.max_cluster_diameter << ","
            << result.clusters_measured << "," << result.clusters_skipped << ",\"";
        for (size_t bin = 0; bin < stats.residual_size_hist.size(); bin++)
            csv << (bin ? " " : "") << stats.residual_size_hist[bin];
        csv << "\"\n";
    }

    if (options.csv_path.empty()) {
        std::cout << csv.str();
    } else {
        std::ofstream out(options.csv_path);
        out << csv.str();
        std::cerr << "wrote " << options.csv_path << "\n";
    }
    return 0;
}
