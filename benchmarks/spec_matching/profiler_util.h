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

/// Benchmark scaffolding for the profilers under this directory: the surface-code corpus generator,
/// the argument parsing helpers and the CSV writer, kept in their own header so that each binary
/// reads as the measurement it is.

#ifndef SPECMATCHING_BENCHMARKS_SPEC_MATCHING_PROFILER_UTIL_H
#define SPECMATCHING_BENCHMARKS_SPEC_MATCHING_PROFILER_UTIL_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "specmatching/sparse_blossom/driver/mwpm_decoding.h"
#include "specmatching/sparse_blossom/driver/user_graph.h"
#include "stim.h"

namespace pm {
namespace spec_matching {
namespace profiler {

const pm::weight_int NUM_DISTINCT_WEIGHTS = 10001;

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

/// A surface-code memory experiment plus a corpus of sampled syndromes, and the true logical flips
/// so that a logical error rate can be measured against them.
struct Experiment {
    stim::DetectorErrorModel dem;
    std::vector<std::vector<uint64_t>> shots;
    std::vector<uint8_t> observable_flips;
};

inline Experiment generate(size_t distance, size_t rounds, double noise, size_t num_shots, uint64_t seed) {
    stim::CircuitGenParameters gen(rounds, distance, "rotated_memory_x");
    gen.after_clifford_depolarization = noise;
    gen.after_reset_flip_probability = noise;
    gen.before_measure_flip_probability = noise;
    stim::Circuit circuit = stim::generate_surface_code_circuit(gen).circuit;

    std::mt19937_64 rng(seed);
    size_t num_detectors = circuit.count_detectors();
    auto dets_obs = stim::sample_batch_detection_events<stim::MAX_BITWORD_WIDTH>(circuit, num_shots, rng);
    auto& dets = dets_obs.first;
    auto& obs = dets_obs.second;

    Experiment experiment;
    experiment.dem = stim::ErrorAnalyzer::circuit_to_detector_error_model(circuit, false, true, false, 0, false, false);
    for (size_t shot = 0; shot < num_shots; shot++) {
        experiment.shots.emplace_back();
        for (size_t d = 0; d < num_detectors; d++) {
            if (dets[d][shot])
                experiment.shots.back().push_back(d);
        }
        experiment.observable_flips.push_back(obs[0][shot] ? 1 : 0);
    }
    return experiment;
}

/// The same experiment, sampled in chunks instead of all at once.
///
/// `generate` materialises every shot, which is fine for the tens of thousands the timing sweeps
/// need and impossible for the campaigns the escalation budget needs: §M3's exit checkpoint asks
/// for 10⁶ shots at the operating point and §M6.4 asks for 10⁷, and at `d = 13` a shot's detection
/// events cost around half a kilobyte, so 10⁶ of them is most of a gigabyte before the decoder has
/// allocated anything. Chunking makes the shot count a question of time rather than of memory.
///
/// The stream is a function of the seed alone, so two runs at the same seed see the same shots
/// whatever the chunk size.
struct ShotSampler {
    stim::Circuit circuit;
    stim::DetectorErrorModel dem;
    std::mt19937_64 rng;
    size_t num_detectors{0};

    /// Refilled by every `sample` call.
    std::vector<std::vector<uint64_t>> shots;
    std::vector<uint8_t> observable_flips;

    static ShotSampler make(size_t distance, size_t rounds, double noise, uint64_t seed) {
        stim::CircuitGenParameters gen(rounds, distance, "rotated_memory_x");
        gen.after_clifford_depolarization = noise;
        gen.after_reset_flip_probability = noise;
        gen.before_measure_flip_probability = noise;

        ShotSampler sampler;
        sampler.circuit = stim::generate_surface_code_circuit(gen).circuit;
        sampler.dem =
            stim::ErrorAnalyzer::circuit_to_detector_error_model(sampler.circuit, false, true, false, 0, false, false);
        sampler.rng = std::mt19937_64(seed);
        sampler.num_detectors = sampler.circuit.count_detectors();
        return sampler;
    }

    void sample(size_t count) {
        auto dets_obs = stim::sample_batch_detection_events<stim::MAX_BITWORD_WIDTH>(circuit, count, rng);
        auto& dets = dets_obs.first;
        auto& obs = dets_obs.second;
        shots.assign(count, {});
        observable_flips.assign(count, 0);
        for (size_t shot = 0; shot < count; shot++) {
            for (size_t d = 0; d < num_detectors; d++) {
                if (dets[d][shot])
                    shots[shot].push_back(d);
            }
            observable_flips[shot] = obs[0][shot] ? 1 : 0;
        }
    }
};

/// One lattice edge weight in DEM float units — the normalisation under which every structure in
/// the M1 results is `d`-independent, and therefore the unit every horizon here is quoted in.
inline double edge_weight_units(const pm::MatchingGraph& graph) {
    std::vector<pm::weight_int> weights;
    for (const auto& node : graph.nodes) {
        for (auto weight : node.neighbor_weights)
            weights.push_back(weight);
    }
    if (weights.empty())
        return 1.0;
    std::sort(weights.begin(), weights.end());
    return (double)weights[weights.size() / 2] / graph.normalising_constant;
}

/// Rows go to stdout as a table and, optionally, to a CSV. Each row carries its section, so one
/// file holds the whole artifact.
struct CsvWriter {
    std::ofstream stream;
    bool enabled{false};

    void open(const std::string& path) {
        if (path.empty())
            return;
        stream.open(path);
        enabled = stream.is_open();
        if (!enabled)
            std::cerr << "warning: could not open " << path << " for writing\n";
    }
    void row(const std::string& text) {
        if (enabled)
            stream << text << "\n";
    }
};

inline double percentile(std::vector<long long> values, double fraction) {
    if (values.empty())
        return 0;
    std::sort(values.begin(), values.end());
    size_t index = (size_t)(fraction * (double)(values.size() - 1));
    return (double)values[std::min(index, values.size() - 1)];
}

/// A one-standard-error linear fit of `ln(LER)` against `d` at fixed `p`. `c2 = d_eff/d` is
/// proportional to the slope, and what the two front ends have to agree on is the slope.
struct Fit {
    double slope{0};
    double stderr_slope{0};
    size_t points{0};
};

inline Fit fit_log_ler(const std::vector<double>& distances, const std::vector<double>& lers) {
    Fit fit;
    std::vector<double> xs;
    std::vector<double> ys;
    for (size_t i = 0; i < lers.size(); i++) {
        if (lers[i] > 0) {
            xs.push_back(distances[i]);
            ys.push_back(std::log(lers[i]));
        }
    }
    fit.points = xs.size();
    if (xs.size() < 3)
        return fit;
    double n = (double)xs.size();
    double mean_x = 0;
    double mean_y = 0;
    for (size_t i = 0; i < xs.size(); i++) {
        mean_x += xs[i] / n;
        mean_y += ys[i] / n;
    }
    double sxx = 0;
    double sxy = 0;
    for (size_t i = 0; i < xs.size(); i++) {
        sxx += (xs[i] - mean_x) * (xs[i] - mean_x);
        sxy += (xs[i] - mean_x) * (ys[i] - mean_y);
    }
    if (sxx == 0)
        return fit;
    fit.slope = sxy / sxx;
    double residual = 0;
    for (size_t i = 0; i < xs.size(); i++) {
        double predicted = mean_y + fit.slope * (xs[i] - mean_x);
        residual += (ys[i] - predicted) * (ys[i] - predicted);
    }
    if (xs.size() > 2)
        fit.stderr_slope = std::sqrt(residual / (n - 2) / sxx);
    return fit;
}

}  // namespace profiler
}  // namespace spec_matching
}  // namespace pm

#endif  // SPECMATCHING_BENCHMARKS_SPEC_MATCHING_PROFILER_UTIL_H
