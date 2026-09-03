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

/// §3 — the `H`-structure profiler.
///
/// Measures the structure of the sparsified graph `H` and its matchability at horizon `T` for the
/// `d x d x d` rotated surface-code memory experiment under circuit-level noise at rate `p`. It
/// produces two things:
///
///  - the **escalation rate** `q` — the fraction of shots on which some connected component of `H`
///    fails to resolve within `T`, which is a property of `(DEM, shot, T)` alone; and
///  - the **connected-component statistics** of `H`: sizes, degrees, diameters, boundary structure,
///    and each component's own `COMPLETE`/`TRUNCATED` status.
///
/// The second comes for free from the first. Every component is decided by truncated sparse blossom
/// run on it in isolation (§2), so a per-component status is what the decode already computes; this
/// binary reads it rather than deriving anything of its own.
///
/// **No latency, no timers, no decoder output kept.** This file includes no clock header beyond
/// `<ctime>`, declares no `HiResTimer`, and starts nothing. It reads counters. The one clock it
/// touches is `std::time`, for the wall-clock start/end line in `run.log` — a run log entry, not a
/// measurement, and it appears in no other output, which is what keeps `summary.csv` and
/// `hists.json` byte-identical across runs at the same seed.
///
/// It also does **not** escalate. `escalate_to_stock` is never called: the escalation *trigger* is
/// the whole measurement, and the answer the fallback would produce is thrown away here. Observables
/// and weights are discarded on every shot.
///
/// ## Units
///
/// `--T` is in **multiples of one lattice edge weight**, the normalisation under which the `M1`
/// structural results are `d`-independent, and the one every horizon in `benchmarks/spec_matching/`
/// is quoted in. That is the unit that makes `T = 2` here the same `T = 2` as in the design's
/// tables; a horizon quoted in raw DEM float weight units would be a fraction of an edge at
/// `p = 1e-3` and would escalate every shot. Both readings are written out — `T` (edge weights),
/// `T_weight_units` (DEM float) and `T_int` (the flooder's integer time units, read off the decoder
/// rather than reconverted) — so no consumer has to infer which is which.
///
/// ## Usage
///
///   sparse_graph_stats [--d 5,7,9] [--p 1e-3,5e-4] [--T 1.5,2,2.5] [--shots N] [--seed S]
///                      [--build SCAN|BITSET] [--size-cap 32] [--degree-cap 32]
///                      [--diameter-cap 64] [--verify] [--out DIR]
///                      [--dem FILE --dets FILE]
///
/// Writes `summary.csv` (one row per `(d, p, T)`), `hists.json` (keyed `"d=..,p=..,T=.."`) and
/// `run.log` to `--out`. The binary prints no interpretation.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "benchmarks/spec_matching/profiler_util.h"
#include "specmatching/spec_matching/manifold/ball_decoding.h"

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

#ifndef SPECMATCHING_GIT_HASH
#define SPECMATCHING_GIT_HASH "unknown"
#endif
#ifndef SPECMATCHING_BUILD_FLAGS
#define SPECMATCHING_BUILD_FLAGS "unknown"
#endif
#ifndef SPECMATCHING_STIM_VERSION
#define SPECMATCHING_STIM_VERSION "unknown"
#endif

using namespace pm::spec_matching;
using namespace pm::spec_matching::profiler;

namespace {

/// Shots sampled and decoded per batch. Bounds the memory a 10^6-shot cell needs — a `d = 17`
/// shot's detection events cost most of a kilobyte — and is what lets every `T` see the *same*
/// shots without materialising the whole corpus: one batch is sampled, then decoded once per `T`.
constexpr size_t SHOT_BATCH = 4096;

/// Above this many shots the per-shot vectors are dropped from `hists.json` and only the aggregates
/// are written (§3.5.2). At 10^6 shots per cell the raw vectors would dominate the file.
constexpr size_t PER_SHOT_VECTOR_LIMIT = 100000;

struct Options {
    std::vector<size_t> distances = {5};
    std::vector<double> error_rates = {1e-3};
    std::vector<double> horizons = {2.0};
    size_t shots = 10000;
    uint64_t seed = 20260903;
    BallGraphBuildMode mode{BallGraphBuildMode::SCAN};
    size_t size_cap = 32;
    size_t degree_cap = 32;
    /// Components larger than this get no diameter and are counted instead. Note this is a
    /// *component size*, not a diameter value: both diameter walks are `O(s^2)` over the component.
    uint32_t diameter_cap = 64;
    bool verify = false;
    std::string out_dir = "benchmarks/spec_matching/results/sparse_graph_stats";
    /// The stim bypass of §3.3. Both or neither.
    std::string dem_path;
    std::string dets_path;
};

/// Whether the progress line may use carriage returns and the erase-to-end-of-line escape. A
/// redirected run gets one plain line per batch instead of a megabyte of control codes.
bool stderr_is_terminal() {
#if defined(__unix__) || defined(__APPLE__)
    if (isatty(fileno(stderr)) != 1)
        return false;
    const char* term = std::getenv("TERM");
    return term != nullptr && std::string(term) != "dumb";
#else
    return false;
#endif
}

const char* mode_name(BallGraphBuildMode mode) {
    return mode == BallGraphBuildMode::BITSET ? "BITSET" : "SCAN";
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
        if (flag == "--d" || flag == "--distances") {
            options.distances = parse_list<size_t>(next());
        } else if (flag == "--p" || flag == "--error-rates") {
            options.error_rates = parse_list<double>(next());
        } else if (flag == "--T" || flag == "--horizons") {
            options.horizons = parse_list<double>(next());
        } else if (flag == "--shots") {
            options.shots = std::stoul(next());
        } else if (flag == "--seed") {
            options.seed = std::stoull(next());
        } else if (flag == "--build" || flag == "--mode") {
            std::string value = next();
            if (value == "BITSET" || value == "bitset") {
                options.mode = BallGraphBuildMode::BITSET;
            } else if (value == "SCAN" || value == "scan") {
                options.mode = BallGraphBuildMode::SCAN;
            } else {
                throw std::invalid_argument("unrecognised --build " + value + " (want SCAN or BITSET)");
            }
        } else if (flag == "--size-cap") {
            options.size_cap = std::stoul(next());
        } else if (flag == "--degree-cap") {
            options.degree_cap = std::stoul(next());
        } else if (flag == "--diameter-cap") {
            options.diameter_cap = (uint32_t)std::stoul(next());
        } else if (flag == "--verify") {
            options.verify = true;
        } else if (flag == "--out" || flag == "--out-dir") {
            options.out_dir = next();
        } else if (flag == "--dem") {
            options.dem_path = next();
        } else if (flag == "--dets") {
            options.dets_path = next();
        } else {
            throw std::invalid_argument("unrecognised flag " + flag);
        }
    }
    if (options.distances.empty() || options.error_rates.empty() || options.horizons.empty())
        throw std::invalid_argument("--d, --p and --T each need at least one value");
    if (options.shots == 0)
        throw std::invalid_argument("--shots must be positive");
    if (options.size_cap == 0 || options.degree_cap == 0 || options.diameter_cap == 0)
        throw std::invalid_argument("the three caps must be positive");
    if (options.dem_path.empty() != options.dets_path.empty())
        throw std::invalid_argument("--dem and --dets go together: one without the other has nothing to read");
    for (double horizon : options.horizons) {
        if (!(horizon > 0))
            throw std::invalid_argument("--T values must be positive");
    }
    return options;
}

/// One `(d, p, T)` cell: §3.5.1's escalation counters and §3.5.2's structure, accumulated over the
/// cell's shots.
///
/// Everything here is a count or a sum of counts. Nothing is timed, and the struct holds no clock.
struct CellStats {
    // ---- §3.5.1, over every shot.
    uint64_t shots{0};
    uint64_t shots_zero_defects{0};
    uint64_t shots_escalated{0};
    uint64_t components_total{0};
    uint64_t components_truncated{0};
    uint64_t truncated_components_on_escalated{0};

    // ---- §3.5.2, over the non-empty shots. `shots_with_defects` is the divisor of every mean
    // below: a zero-defect shot stays in `q`'s denominator and is excluded from these.
    uint64_t shots_with_defects{0};
    uint64_t sum_n_defects{0};
    uint64_t sum_h_nodes{0};
    uint64_t sum_h_edges{0};
    uint64_t sum_h_boundary_edges{0};
    uint64_t sum_nodes_with_boundary_edge{0};
    uint64_t sum_num_components{0};
    uint64_t sum_largest_component_size{0};
    uint64_t sum_singleton_components{0};
    uint64_t sum_pair_components{0};
    uint64_t sum_components_size_ge3{0};
    uint64_t sum_odd_components_without_boundary{0};
    uint64_t sum_diameter_uncomputed_components{0};
    uint64_t shots_with_odd_component_without_boundary{0};

    uint64_t max_n_defects{0};
    uint64_t max_h_nodes{0};
    uint64_t max_h_edges{0};
    uint64_t max_h_boundary_edges{0};
    uint64_t max_nodes_with_boundary_edge{0};
    uint64_t max_num_components{0};
    uint64_t max_largest_component_size{0};
    uint64_t max_singleton_components{0};
    uint64_t max_pair_components{0};
    uint64_t max_components_size_ge3{0};
    uint64_t max_odd_components_without_boundary{0};
    uint64_t max_component_diameter_wint{0};
    uint64_t max_component_hop_diameter{0};

    ComponentHistograms hist;
    ComponentStatusTable size_x_status;
    ComponentStatusTable hop_diameter_x_status;

    /// §3.5.2's raw per-shot vectors, kept only when `shots <= PER_SHOT_VECTOR_LIMIT`.
    bool keep_per_shot{false};
    std::vector<int> per_shot_n_defects;
    std::vector<int> per_shot_h_edges;
    std::vector<int> per_shot_h_boundary_edges;
    std::vector<int> per_shot_num_components;
    std::vector<int> per_shot_largest_component;
    std::vector<int> per_shot_singleton_components;
    std::vector<int> per_shot_pair_components;
    std::vector<int> per_shot_components_size_ge3;
    std::vector<int> per_shot_odd_components_without_boundary;
    std::vector<int> per_shot_components_truncated;

    /// The cell's labels, filled once at the end of the run.
    double t_weight_units{0};
    int64_t t_int{0};

    void configure(const Options& options) {
        hist.configure(options.size_cap, options.degree_cap);
        size_x_status.configure(options.size_cap + 1);
        hop_diameter_x_status.configure(options.size_cap + 1);
        keep_per_shot = options.shots <= PER_SHOT_VECTOR_LIMIT;
    }

    void accumulate(
        const BallProfile& profile,
        const ComponentHistograms& shot_hist,
        const ComponentStatusTable& shot_size_x_status,
        const ComponentStatusTable& shot_hop_x_status) {
        shots++;
        if (profile.n_defects == 0)
            shots_zero_defects++;
        // §2.6's one predicate. Read off the decode, never re-derived: this binary never harvests a
        // residual and never runs the fallback, so there is nothing else it *could* be read from.
        if (profile.any_component_truncated) {
            shots_escalated++;
            truncated_components_on_escalated += (uint64_t)profile.components_truncated;
        }
        components_total += (uint64_t)profile.components_total;
        components_truncated += (uint64_t)profile.components_truncated;

        if (keep_per_shot) {
            // One entry per shot, empty ones included, so the vectors are index-aligned with the
            // shot number and a reader can join them to each other without a mask.
            per_shot_n_defects.push_back(profile.n_defects);
            per_shot_h_edges.push_back(profile.h_edges);
            per_shot_h_boundary_edges.push_back(profile.h_boundary_edges);
            per_shot_num_components.push_back(profile.components_total);
            per_shot_largest_component.push_back(profile.components.largest_component_size);
            per_shot_singleton_components.push_back(profile.components.num_singleton_components);
            per_shot_pair_components.push_back(profile.components.num_pair_components);
            per_shot_components_size_ge3.push_back(profile.components.num_components_size_ge3);
            per_shot_odd_components_without_boundary.push_back(
                profile.components.num_odd_components_without_boundary);
            per_shot_components_truncated.push_back(profile.components_truncated);
        }

        const ComponentStats& components = profile.components;
        // The distributions and the joint tables are over components and edges, so they are folded
        // in whatever the shot's defect count — an empty shot contributes nothing to them anyway.
        hist.add(shot_hist);
        size_x_status.add(shot_size_x_status);
        hop_diameter_x_status.add(shot_hop_x_status);

        // §3.4: zero-defect shots count in `shots_zero_defects`, stay in `q`'s denominator, and are
        // excluded from the per-shot structural means.
        if (profile.h_nodes == 0)
            return;
        shots_with_defects++;
        sum_n_defects += (uint64_t)profile.n_defects;
        sum_h_nodes += (uint64_t)profile.h_nodes;
        sum_h_edges += (uint64_t)profile.h_edges;
        sum_h_boundary_edges += (uint64_t)profile.h_boundary_edges;
        sum_nodes_with_boundary_edge += (uint64_t)components.nodes_with_boundary_edge;
        sum_num_components += (uint64_t)components.num_components;
        sum_largest_component_size += (uint64_t)components.largest_component_size;
        sum_singleton_components += (uint64_t)components.num_singleton_components;
        sum_pair_components += (uint64_t)components.num_pair_components;
        sum_components_size_ge3 += (uint64_t)components.num_components_size_ge3;
        sum_odd_components_without_boundary += (uint64_t)components.num_odd_components_without_boundary;
        sum_diameter_uncomputed_components += (uint64_t)components.diameter_uncomputed_components;
        if (components.num_odd_components_without_boundary > 0)
            shots_with_odd_component_without_boundary++;

        max_n_defects = std::max(max_n_defects, (uint64_t)profile.n_defects);
        max_h_nodes = std::max(max_h_nodes, (uint64_t)profile.h_nodes);
        max_h_edges = std::max(max_h_edges, (uint64_t)profile.h_edges);
        max_h_boundary_edges = std::max(max_h_boundary_edges, (uint64_t)profile.h_boundary_edges);
        max_nodes_with_boundary_edge =
            std::max(max_nodes_with_boundary_edge, (uint64_t)components.nodes_with_boundary_edge);
        max_num_components = std::max(max_num_components, (uint64_t)components.num_components);
        max_largest_component_size = std::max(max_largest_component_size, (uint64_t)components.largest_component_size);
        max_singleton_components = std::max(max_singleton_components, (uint64_t)components.num_singleton_components);
        max_pair_components = std::max(max_pair_components, (uint64_t)components.num_pair_components);
        max_components_size_ge3 = std::max(max_components_size_ge3, (uint64_t)components.num_components_size_ge3);
        max_odd_components_without_boundary =
            std::max(max_odd_components_without_boundary, (uint64_t)components.num_odd_components_without_boundary);
        max_component_diameter_wint =
            std::max(max_component_diameter_wint, (uint64_t)components.max_component_diameter_wint);
        max_component_hop_diameter =
            std::max(max_component_hop_diameter, (uint64_t)components.max_component_hop_diameter);
    }

    double q() const {
        return shots ? (double)shots_escalated / (double)shots : 0.0;
    }
    double mean_over_non_empty(uint64_t sum) const {
        return shots_with_defects ? (double)sum / (double)shots_with_defects : 0.0;
    }
};

/// The 95% Wilson score interval for `k` successes in `n` trials.
///
/// Wilson rather than the normal approximation because `q` is small — at `q ~ 2e-4` and `n = 10^6`
/// the normal interval is asymmetric enough to matter and can run below zero, which is not a
/// statement anyone should have to read past.
std::pair<double, double> wilson_interval(uint64_t successes, uint64_t trials) {
    if (trials == 0)
        return {0.0, 0.0};
    const double z = 1.959963984540054;  // two-sided 95%
    double n = (double)trials;
    double phat = (double)successes / n;
    double denominator = 1.0 + z * z / n;
    double centre = (phat + z * z / (2 * n)) / denominator;
    double half = (z / denominator) * std::sqrt(phat * (1 - phat) / n + z * z / (4 * n * n));
    double lo = centre - half;
    double hi = centre + half;
    // The two ends are exact at the extremes — `centre == half` when no shot escalated — and
    // `centre - half` is where the cancellation lands, so it comes back as `1e-18` rather than as
    // the 0 it is. Set them from the count instead of trying to round the subtraction.
    if (successes == 0)
        lo = 0.0;
    if (successes == trials)
        hi = 1.0;
    return {lo < 0 ? 0.0 : lo, hi > 1 ? 1.0 : hi};
}

/// The generator call, as one string, so `run.log` records what `(d, p)` actually meant.
///
/// It is `benchmarks/spec_matching/profiler_util.h`'s `ShotSampler::make`, which is the same call
/// the rest of the `benchmarks/spec_matching/` corpus is generated by — that is the requirement:
/// `(d, p)` has to mean here what it means in every table beside it.
///
/// Two deviations from §3.3's prose, recorded rather than silently taken. The basis is
/// `rotated_memory_x`, not `_z`; and three noise parameters are set, not four
/// (`before_round_data_depolarization` is left at 0). Both are what the existing corpus does, and
/// matching the existing corpus is the stated reason the shared call is required at all.
std::string generator_call(size_t distance, size_t rounds, double noise) {
    std::ostringstream out;
    out << "stim::CircuitGenParameters(rounds=" << rounds << ", distance=" << distance
        << ", task=\"rotated_memory_x\"); after_clifford_depolarization=" << noise
        << "; after_reset_flip_probability=" << noise << "; before_measure_flip_probability=" << noise
        << "; before_round_data_depolarization=0 (unset)"
        << "; stim::generate_surface_code_circuit(...).circuit"
        << "; stim::ErrorAnalyzer::circuit_to_detector_error_model(circuit, false, true, false, 0, false, false)";
    return out.str();
}

/// A corpus read off disk instead of sampled: a DEM text file and a b8 detection-event file.
///
/// b8 packs each shot into `ceil(num_detectors / 8)` bytes, bit `i` of byte `i / 8` at position
/// `i % 8` — stim's own convention, so a file written by `stim sample_dem`/`stim detect` reads back
/// directly.
struct FileCorpus {
    stim::DetectorErrorModel dem;
    std::vector<std::vector<uint64_t>> shots;
};

FileCorpus read_file_corpus(const std::string& dem_path, const std::string& dets_path, size_t max_shots) {
    FileCorpus corpus;
    FILE* dem_file = std::fopen(dem_path.c_str(), "r");
    if (dem_file == nullptr)
        throw std::invalid_argument("could not open --dem " + dem_path);
    corpus.dem = stim::DetectorErrorModel::from_file(dem_file);
    std::fclose(dem_file);

    size_t num_detectors = corpus.dem.count_detectors();
    size_t bytes_per_shot = (num_detectors + 7) / 8;
    if (bytes_per_shot == 0)
        throw std::invalid_argument("the DEM in " + dem_path + " declares no detectors");

    std::ifstream dets(dets_path, std::ios::binary);
    if (!dets.is_open())
        throw std::invalid_argument("could not open --dets " + dets_path);
    std::vector<char> buffer(bytes_per_shot);
    while (corpus.shots.size() < max_shots && dets.read(buffer.data(), (std::streamsize)bytes_per_shot)) {
        corpus.shots.emplace_back();
        for (size_t d = 0; d < num_detectors; d++) {
            if ((buffer[d / 8] >> (d % 8)) & 1)
                corpus.shots.back().push_back(d);
        }
    }
    if (corpus.shots.empty())
        throw std::invalid_argument("--dets " + dets_path + " held no complete shot of " +
                                    std::to_string(bytes_per_shot) + " bytes");
    return corpus;
}

/// Decodes one batch of shots at one horizon and folds the results into `cell`.
///
/// This is §3.4 step 4 in full: build `H`, split it, solve every component in isolation, read the
/// statuses, and analyse the structure. `escalate_to_stock` is not called and the observables and
/// weight are discarded — the trigger is all this binary needs.
void decode_batch_into(
    BallDecoder& decoder,
    const std::vector<std::vector<uint64_t>>& shots,
    ComponentHistograms& shot_hist,
    ComponentStatusTable& shot_size_x_status,
    ComponentStatusTable& shot_hop_x_status,
    CellStats& cell) {
    BallProfile profile;
    for (const std::vector<uint64_t>& shot : shots) {
        // The production Phase 1. Its `Phase1Outcome` — the committed observables and weight — is
        // deliberately dropped on the floor; the escalation predicate and the per-component
        // statuses ride out on the profile.
        decoder.decode_phase1_production(shot, &profile);
        // Untimed by construction, and after the decode rather than inside it. It reads the `H` the
        // arena is still holding and the statuses the decode just wrote.
        decoder.analyze_last_shot_components(profile, shot_hist, shot_size_x_status, shot_hop_x_status);
        cell.accumulate(profile, shot_hist, shot_size_x_status, shot_hop_x_status);
    }
}

void write_summary_header(std::ofstream& out) {
    // `corpus` says whether `d` and `p` describe the circuit that produced these shots or are
    // merely the labels the run was given. On `--dem`/`--dets` they are labels: the DEM came off
    // disk and this binary has no idea what generated it, so reading `d = 5` off a `corpus=file`
    // row would be inventing a fact.
    out << "corpus,d,p,T,T_weight_units,T_int,R,build_mode,seed,dem_has_negative_weights,"
           "shots,shots_zero_defects,shots_escalated,q,q_ci_lo,q_ci_hi,"
           "components_total,components_truncated,truncated_component_rate,"
           "truncated_components_per_escalated_shot,"
           "shots_with_defects,mean_n_defects,mean_h_nodes,mean_h_edges,mean_h_boundary_edges,"
           "mean_nodes_with_boundary_edge,mean_num_components,mean_largest_component_size,"
           "mean_singleton_components,mean_pair_components,mean_components_size_ge3,"
           "mean_odd_components_without_boundary,odd_component_without_boundary_rate,"
           "max_n_defects,max_h_nodes,max_h_edges,max_h_boundary_edges,max_nodes_with_boundary_edge,"
           "max_num_components,max_largest_component_size,max_singleton_components,max_pair_components,"
           "max_components_size_ge3,max_odd_components_without_boundary,"
           "max_component_diameter_wint,max_component_hop_diameter,"
           "diameter_uncomputed_components\n";
}

void write_summary_row(
    std::ofstream& out,
    const Options& options,
    size_t distance,
    double noise,
    double horizon,
    double ball_radius,
    bool dem_has_negative_weights,
    bool from_file,
    const CellStats& cell) {
    std::pair<double, double> ci = wilson_interval(cell.shots_escalated, cell.shots);
    double truncated_component_rate =
        cell.components_total ? (double)cell.components_truncated / (double)cell.components_total : 0.0;
    double truncated_per_escalated =
        cell.shots_escalated ? (double)cell.truncated_components_on_escalated / (double)cell.shots_escalated : 0.0;

    out << (from_file ? "file" : "stim") << "," << distance << "," << noise << "," << horizon << "," << cell.t_weight_units << "," << cell.t_int << ","
        << ball_radius << "," << mode_name(options.mode) << "," << options.seed << ","
        << (dem_has_negative_weights ? 1 : 0) << "," << cell.shots << "," << cell.shots_zero_defects << ","
        << cell.shots_escalated << "," << cell.q() << "," << ci.first << "," << ci.second << ","
        << cell.components_total << "," << cell.components_truncated << "," << truncated_component_rate << ","
        << truncated_per_escalated << "," << cell.shots_with_defects << ","
        << cell.mean_over_non_empty(cell.sum_n_defects) << "," << cell.mean_over_non_empty(cell.sum_h_nodes) << ","
        << cell.mean_over_non_empty(cell.sum_h_edges) << "," << cell.mean_over_non_empty(cell.sum_h_boundary_edges)
        << "," << cell.mean_over_non_empty(cell.sum_nodes_with_boundary_edge) << ","
        << cell.mean_over_non_empty(cell.sum_num_components) << ","
        << cell.mean_over_non_empty(cell.sum_largest_component_size) << ","
        << cell.mean_over_non_empty(cell.sum_singleton_components) << ","
        << cell.mean_over_non_empty(cell.sum_pair_components) << ","
        << cell.mean_over_non_empty(cell.sum_components_size_ge3) << ","
        << cell.mean_over_non_empty(cell.sum_odd_components_without_boundary) << ","
        << (cell.shots_with_defects
                ? (double)cell.shots_with_odd_component_without_boundary / (double)cell.shots_with_defects
                : 0.0)
        << "," << cell.max_n_defects << "," << cell.max_h_nodes << "," << cell.max_h_edges << ","
        << cell.max_h_boundary_edges << "," << cell.max_nodes_with_boundary_edge << "," << cell.max_num_components
        << "," << cell.max_largest_component_size << "," << cell.max_singleton_components << ","
        << cell.max_pair_components << "," << cell.max_components_size_ge3 << ","
        << cell.max_odd_components_without_boundary << "," << cell.max_component_diameter_wint << ","
        << cell.max_component_hop_diameter << "," << cell.sum_diameter_uncomputed_components << "\n";
}

/// A JSON array of integers, with the last element understood as the overflow bin. Written by hand
/// rather than through a library, because the whole file is arrays of integers and the shape is
/// what a reader needs stated, not the encoder.
void write_json_array(std::ofstream& out, const std::vector<uint64_t>& values) {
    out << "[";
    for (size_t i = 0; i < values.size(); i++)
        out << (i ? "," : "") << values[i];
    out << "]";
}

void write_json_array(std::ofstream& out, const std::vector<int>& values) {
    out << "[";
    for (size_t i = 0; i < values.size(); i++)
        out << (i ? "," : "") << values[i];
    out << "]";
}

/// The joint tables, as `[[complete, truncated], ...]` — one row per bin, last row overflowing.
void write_json_status_table(std::ofstream& out, const ComponentStatusTable& table) {
    out << "[";
    for (size_t bin = 0; bin < table.bins; bin++) {
        out << (bin ? "," : "") << "[" << table.at(bin, ComponentStatusTable::COMPLETE) << ","
            << table.at(bin, ComponentStatusTable::TRUNCATED) << "]";
    }
    out << "]";
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

    std::error_code dir_error;
    std::filesystem::create_directories(options.out_dir, dir_error);
    if (dir_error) {
        std::cerr << "error: could not create " << options.out_dir << ": " << dir_error.message() << "\n";
        return 1;
    }

    std::ofstream summary(options.out_dir + "/summary.csv");
    std::ofstream hists(options.out_dir + "/hists.json");
    std::ofstream log(options.out_dir + "/run.log");
    if (!summary.is_open() || !hists.is_open() || !log.is_open()) {
        std::cerr << "error: could not open the output files under " << options.out_dir << "\n";
        return 1;
    }

    // The one clock this binary touches. It is a run-log entry, not a measurement, and it appears in
    // neither `summary.csv` nor `hists.json` — which is what keeps those two byte-identical across
    // runs at the same seed (§3.7).
    auto stamp = []() {
        std::time_t now = std::time(nullptr);
        char buffer[64];
        std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S%z", std::localtime(&now));
        return std::string(buffer);
    };

    log << "sparse_graph_stats — the H-structure profiler (§3)\n";
    log << "git_hash=" << SPECMATCHING_GIT_HASH << "\n";
    log << "build_flags=" << SPECMATCHING_BUILD_FLAGS << "\n";
    log << "stim_version=" << SPECMATCHING_STIM_VERSION << "\n";
    log << "started=" << stamp() << "\n";
    log << "seed=" << options.seed << "\n";
    log << "shots_per_cell=" << options.shots << "\n";
    log << "build_mode=" << mode_name(options.mode) << "\n";
    log << "size_cap=" << options.size_cap << "\n";
    log << "degree_cap=" << options.degree_cap << "\n";
    log << "diameter_cap=" << options.diameter_cap << " (a component SIZE above which no diameter is computed)\n";
    log << "verify_component_decomposition=" << (options.verify ? 1 : 0) << "\n";
    log << "T_unit=one lattice edge weight (the median discretised edge of G, in DEM float units);"
           " summary.csv carries T, T_weight_units and T_int\n";
    log << "no_timing=1 (this binary contains no timing code and calls escalate_to_stock never)\n";
    log.flush();

    write_summary_header(summary);
    hists << "{\n";
    hists << "  \"caps\": {\"size_cap\": " << options.size_cap << ", \"degree_cap\": " << options.degree_cap
          << ", \"diameter_cap\": " << options.diameter_cap
          << ", \"weight_hist_bins_per_T\": " << ComponentHistograms::BINS_PER_T << "},\n";
    hists << "  \"note\": \"every histogram's last bin is the overflow bin; the joint tables are"
             " [[complete, truncated], ...] with one row per bin\",\n";
    // What a bin *is*, stated in the file rather than left to the reader. `component_size_hist`,
    // `degree_hist` and `hop_diameter_hist` are binned by count — bin `k` counts the value `k`.
    // `wint_diameter_hist` is binned in sixteenths of `T_int`, the §C.2 convention every weight
    // histogram in this project shares, so a consumer that wants coarser bins sums groups of
    // `weight_hist_bins_per_T` and gets width-`T` bins exactly.
    hists << "  \"bin_units\": {\"component_size_hist\": \"count\", \"degree_hist\": \"count\","
             " \"hop_diameter_hist\": \"count (H edges on the longest shortest hop path)\","
             " \"wint_diameter_hist\": \"T/16 in the flooder's integer time units; see T_int per cell\"},\n";
    hists << "  \"diameter_note\": \"both diameters are confined to the component's own subgraph of"
             " H. A pair whose G geodesic leaves the component is measured by the route that stays"
             " inside it, so these are H diameters and must not be read as distances in G.\",\n";
    hists << "  \"cells\": {\n";
    bool first_cell = true;

    bool interactive = stderr_is_terminal();
    bool ok = true;
    // The file bypass reads one corpus, not one per `(d, p)`; sweeping the labels over it would
    // write the same numbers under different names, so the first cell ends the run.
    bool stop = false;
    for (size_t distance : options.distances) {
        if (stop)
            break;
        for (double noise : options.error_rates) {
            // ---- §3.4 step 1. Circuit -> DEM, or the file bypass.
            stim::DetectorErrorModel dem;
            ShotSampler sampler;
            FileCorpus file_corpus;
            bool from_file = !options.dem_path.empty();
            if (from_file) {
                try {
                    file_corpus = read_file_corpus(options.dem_path, options.dets_path, options.shots);
                } catch (const std::exception& error) {
                    std::cerr << "error: " << error.what() << "\n";
                    return 1;
                }
                dem = file_corpus.dem;
                log << "corpus=file dem=" << options.dem_path << " dets=" << options.dets_path
                    << " shots_read=" << file_corpus.shots.size()
                    << " (--d and --p are labels only on this path)\n";
            } else {
                sampler = ShotSampler::make(distance, distance, noise, options.seed);
                dem = sampler.dem;
                log << "corpus=stim d=" << distance << " rounds=" << distance << " p=" << noise << "\n";
                log << "  generator_call=" << generator_call(distance, distance, noise) << "\n";
            }
            size_t cell_shots = from_file ? file_corpus.shots.size() : options.shots;
            char label[160];
            if (from_file) {
                std::snprintf(label, sizeof(label), "corpus=file");
            } else {
                std::snprintf(label, sizeof(label), "d=%zu p=%g", distance, noise);
            }

            // ---- §3.4 step 2. Ball tables once, at `T_max = max(T list)` and `R = 2 * T_max`.
            pm::Mwpm probe = pm::detector_error_model_to_mwpm(dem, NUM_DISTINCT_WEIGHTS, false);
            double unit = edge_weight_units(probe.flooder.graph);
            double max_horizon = *std::max_element(options.horizons.begin(), options.horizons.end());

            BallConfig ball_config;
            ball_config.T = max_horizon * unit;
            ball_config.ball.T_max = max_horizon * unit;
            ball_config.ball.R = 2.0 * ball_config.ball.T_max;
            ball_config.mode = options.mode;
            ball_config.collect_component_stats = true;
            ball_config.diameter_cap = options.diameter_cap;
            ball_config.verify_component_decomposition = options.verify;

            BallDecoder decoder = BallDecoder::from_detector_error_model(dem, ball_config, NUM_DISTINCT_WEIGHTS);
            log << "  dem_has_negative_weights=" << (decoder.dem_has_negative_weights ? 1 : 0)
                << " edge_weight_unit=" << unit << " T_max=" << ball_config.ball.T_max
                << " R=" << ball_config.ball.R << "\n";
            log.flush();

            std::vector<CellStats> cells(options.horizons.size());
            for (CellStats& cell : cells)
                cell.configure(options);

            ComponentHistograms shot_hist;
            shot_hist.configure(options.size_cap, options.degree_cap);
            ComponentStatusTable shot_size_x_status = ComponentStatusTable::with_bins(options.size_cap + 1);
            ComponentStatusTable shot_hop_x_status = ComponentStatusTable::with_bins(options.size_cap + 1);

            // ---- §3.4 steps 3-4. One batch of shots at a time, decoded once per `T`, so shot `i`
            // of `(d, p)` is the identical shot at every horizon.
            size_t decoded = 0;
            try {
                while (decoded < cell_shots) {
                    size_t batch = std::min(SHOT_BATCH, cell_shots - decoded);
                    const std::vector<std::vector<uint64_t>>* shots = nullptr;
                    std::vector<std::vector<uint64_t>> slice;
                    if (from_file) {
                        slice.assign(
                            file_corpus.shots.begin() + (ptrdiff_t)decoded,
                            file_corpus.shots.begin() + (ptrdiff_t)(decoded + batch));
                        shots = &slice;
                    } else {
                        sampler.sample(batch);
                        shots = &sampler.shots;
                    }
                    for (size_t t = 0; t < options.horizons.size(); t++) {
                        decoder.set_horizon(options.horizons[t] * unit);
                        cells[t].t_weight_units = options.horizons[t] * unit;
                        cells[t].t_int = (int64_t)decoder.horizon;
                        decode_batch_into(
                            decoder, *shots, shot_hist, shot_size_x_status, shot_hop_x_status, cells[t]);
                    }
                    decoded += batch;
                    // Between batches, and to stderr: it lands in no measurement — there are none —
                    // and it keeps stdout clean for the one line this binary prints.
                    std::fprintf(
                        stderr,
                        interactive ? "\r%s  %zu/%zu shots\x1b[K" : "%s  %zu/%zu shots\n",
                        label,
                        decoded,
                        cell_shots);
                    std::fflush(stderr);
                }
            } catch (const std::exception& error) {
                std::fprintf(stderr, "\n");
                std::cerr << "error while decoding d=" << distance << " p=" << noise << ": " << error.what() << "\n";
                log << "  ERROR=" << error.what() << "\n";
                ok = false;
                break;
            }
            if (interactive)
                std::fprintf(stderr, "\r\x1b[K");

            // ---- §3.4 step 5, and §3.6's two files.
            for (size_t t = 0; t < options.horizons.size(); t++) {
                const CellStats& cell = cells[t];
                write_summary_row(
                    summary,
                    options,
                    distance,
                    noise,
                    options.horizons[t],
                    ball_config.ball.R,
                    decoder.dem_has_negative_weights,
                    from_file,
                    cell);

                char key[128];
                std::snprintf(key, sizeof(key), "d=%zu,p=%g,T=%g", distance, noise, options.horizons[t]);
                hists << (first_cell ? "    \"" : ",\n    \"") << key << "\": {\n";
                first_cell = false;
                hists << "      \"shots\": " << cell.shots << ",\n";
                hists << "      \"shots_escalated\": " << cell.shots_escalated << ",\n";
                hists << "      \"components_total\": " << cell.components_total << ",\n";
                hists << "      \"components_truncated\": " << cell.components_truncated << ",\n";
                hists << "      \"T_int\": " << cell.t_int << ",\n";
                hists << "      \"component_size_hist\": ";
                write_json_array(hists, cell.hist.size_hist);
                hists << ",\n      \"degree_hist\": ";
                write_json_array(hists, cell.hist.degree_hist);
                hists << ",\n      \"hop_diameter_hist\": ";
                write_json_array(hists, cell.hist.hop_diameter_hist);
                hists << ",\n      \"wint_diameter_hist\": ";
                write_json_array(hists, cell.hist.diameter_hist);
                hists << ",\n      \"diameter_uncomputed_components\": " << cell.sum_diameter_uncomputed_components;
                hists << ",\n      \"size_x_status\": ";
                write_json_status_table(hists, cell.size_x_status);
                hists << ",\n      \"hop_diameter_x_status\": ";
                write_json_status_table(hists, cell.hop_diameter_x_status);
                if (cell.keep_per_shot) {
                    hists << ",\n      \"per_shot\": {\"n_defects\": ";
                    write_json_array(hists, cell.per_shot_n_defects);
                    hists << ", \"h_edges\": ";
                    write_json_array(hists, cell.per_shot_h_edges);
                    hists << ", \"h_boundary_edges\": ";
                    write_json_array(hists, cell.per_shot_h_boundary_edges);
                    hists << ", \"num_components\": ";
                    write_json_array(hists, cell.per_shot_num_components);
                    hists << ", \"largest_component_size\": ";
                    write_json_array(hists, cell.per_shot_largest_component);
                    hists << ", \"num_singleton_components\": ";
                    write_json_array(hists, cell.per_shot_singleton_components);
                    hists << ", \"num_pair_components\": ";
                    write_json_array(hists, cell.per_shot_pair_components);
                    hists << ", \"num_components_size_ge3\": ";
                    write_json_array(hists, cell.per_shot_components_size_ge3);
                    hists << ", \"num_odd_components_without_boundary\": ";
                    write_json_array(hists, cell.per_shot_odd_components_without_boundary);
                    hists << ", \"components_truncated\": ";
                    write_json_array(hists, cell.per_shot_components_truncated);
                    hists << "}";
                }
                hists << "\n    }";
            }
            summary.flush();
            hists.flush();
            if (from_file) {
                log << "  (--dem/--dets given: the d and p lists are not swept)\n";
                stop = true;
                break;
            }
        }
    }

    hists << "\n  }\n}\n";
    log << "finished=" << stamp() << "\n";
    if (!ok)
        log << "status=incomplete (see ERROR above)\n";
    else
        log << "status=ok\n";

    std::printf("wrote %s/{summary.csv,hists.json,run.log}\n", options.out_dir.c_str());
    return ok ? 0 : 1;
}
