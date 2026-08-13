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

#include "pyrematching/two_phase/driver/two_phase.pybind.h"

#include <string>
#include <vector>

#include "pyrematching/sparse_blossom/driver/user_graph.h"
#include "pyrematching/sparse_blossom/driver/user_graph.pybind.h"
#include "pyrematching/two_phase/driver/two_phase_decoding.h"

using namespace pybind11::literals;
using namespace pm::two_phase;
using pm_pybind::vec_to_array;

namespace {

/// The DEM crosses the boundary as text rather than as a `stim` pybind object.
///
/// The extension links `libstim` but does not depend on `stim`'s own Python bindings, and coupling
/// the two would mean this module could only be imported alongside a matching `stim` build. The
/// inherited `Matching.from_detector_error_model` already round-trips through a string for the same
/// reason. Callers pass `str(dem)`; the python wrapper does it for them.
stim::DetectorErrorModel dem_from_text(const std::string& text) {
    return stim::DetectorErrorModel(text.c_str());
}

py::dict summary_to_dict(const TwoPhaseSummary& summary) {
    py::dict out;
    // Always beside `q`, never without it: a zero below `1 / shots` is a resolution floor rather
    // than a measurement (§M3.0).
    out["q"] = summary.q;
    // §M7.7 asks for the escalation rate under the name `q_this`, next to the landed scheme's
    // `q_current` on the same shots. It is the same number as `q`; both names are exposed so that
    // an M7 report and an M3–M6 report can each use their own vocabulary without a second field.
    out["q_this"] = summary.q;
    out["shots"] = summary.shots;
    out["shots_escalated"] = summary.shots_escalated;
    // `None` rather than 0 when the replay was not run: "not measured" and "measured, zero" are
    // different claims, and at these rates the second one is usually a resolution floor anyway.
    out["q_current_on_same_corpus"] = summary.q_current_on_same_corpus < 0
                                          ? py::cast<py::object>(py::none())
                                          : py::cast<py::object>(py::float_(summary.q_current_on_same_corpus));
    out["shots_with_truncated_reference"] = summary.shots_with_truncated_reference;
    out["certified_rate"] = summary.certified_rate;
    out["h_no_perfect_matching_rate"] = summary.h_no_perfect_matching_rate;
    out["mean_dual_scan_ns"] = summary.mean_dual_scan_ns;
    out["c_phase1_ns"] = summary.c_phase1;
    out["c_escalation_ns"] = summary.c_escalation;
    out["amortised_mean_ns"] = summary.amortised_mean_ns;
    out["measured_mean_ns"] = summary.measured_mean_ns;
    out["amortisation_gap"] = summary.amortisation_gap;
    out["amortised_penalty"] = summary.amortised_penalty;
    out["escalation_cost_ratio"] = summary.escalation_cost_ratio;
    out["mean_stock_ns_on_escalated"] = summary.mean_stock_ns_on_escalated;
    out["speedup_vs_stock"] = summary.speedup_vs_stock;
    out["p50_total_ns"] = summary.p50_total_ns;
    out["p99_total_ns"] = summary.p99_total_ns;
    out["p999_total_ns"] = summary.p999_total_ns;
    out["p9999_total_ns"] = summary.p9999_total_ns;
    out["max_total_ns"] = summary.max_total_ns;
    out["max_escalated_total_ns"] = summary.max_escalated_total_ns;
    out["contaminated_shot_rate"] = summary.contaminated_shot_rate;
    out["mean_residual_size"] = summary.mean_residual_size;
    out["mean_residual_density"] = summary.mean_residual_density;
    out["exposed_root_blossom_rate"] = summary.exposed_root_blossom_rate;

    // §C/§D's component structure. Always beside `shots_with_component_stats`, for the same reason
    // `q` is always beside `shots`: zero analysed shots means the campaign did not measure this,
    // which is not the same claim as "it measured zero".
    out["shots_with_component_stats"] = summary.shots_with_component_stats;
    out["shots_with_component_defects"] = summary.shots_with_component_defects;
    out["mean_components"] = summary.mean_components;
    out["mean_singleton_components"] = summary.mean_singleton_components;
    out["mean_pair_components"] = summary.mean_pair_components;
    out["mean_nontrivial_components"] = summary.mean_nontrivial_components;
    out["mean_largest_component_size"] = summary.mean_largest_component_size;
    out["max_component_size"] = summary.max_component_size;
    out["mean_max_component_diameter_wint"] = summary.mean_max_component_diameter_wint;
    out["max_component_diameter_wint"] = summary.max_component_diameter_wint;
    out["boundary_touching_component_fraction"] = summary.boundary_touching_component_fraction;
    out["diameter_uncomputed_component_fraction"] = summary.diameter_uncomputed_component_fraction;
    out["frac_defects_in_trivial_components"] = summary.frac_defects_in_trivial_components;
    out["frac_defects_to_solver"] = summary.frac_defects_to_solver;
    out["frac_defects_committed_trivially"] = summary.frac_defects_committed_trivially;
    out["frac_defects_residual_trivially"] = summary.frac_defects_residual_trivially;
    out["solver_set_empty_rate"] = summary.solver_set_empty_rate;
    out["trivial_residual_rate"] = summary.trivial_residual_rate;
    // The three weight histograms share one axis: bin `k` is `[k * T / bins_per_T, ...)`, last bin
    // overflowing. Carried so a plot can label it without restating the convention.
    out["weight_hist_bins_per_T"] = summary.weight_hist_bins_per_T;
    out["component_size_hist"] = summary.component_size_hist;
    out["component_diameter_hist"] = summary.component_diameter_hist;
    out["h_degree_hist"] = summary.h_degree_hist;
    out["h_edge_weight_hist"] = summary.h_edge_weight_hist;
    out["boundary_cost_hist"] = summary.boundary_cost_hist;
    return out;
}

py::dict stats_to_dict(const TwoPhaseAggregateStats& stats) {
    py::dict out;
    out["shots"] = stats.shots;
    out["shots_truncated"] = stats.shots_truncated;
    out["shots_escalated"] = stats.shots_escalated;
    out["shots_zero_defects"] = stats.shots_zero_defects;
    out["shots_contaminated"] = stats.shots_contaminated;
    out["shots_certified"] = stats.shots_certified;
    out["shots_h_no_perfect_matching"] = stats.shots_h_no_perfect_matching;
    out["shots_escalated_truncated_reference"] = stats.shots_escalated_truncated_reference;
    out["shots_with_truncated_reference"] = stats.shots_with_truncated_reference;
    out["sum_phase1_ns"] = stats.sum_phase1_ns;
    out["sum_harvest_ns"] = stats.sum_harvest_ns;
    out["sum_dual_scan_ns"] = stats.sum_dual_scan_ns;
    out["sum_escalation_ns"] = stats.sum_escalation_ns;
    out["sum_total_ns"] = stats.sum_total_ns;
    out["sum_total_ns_truncated"] = stats.sum_total_ns_truncated;
    out["sum_exact_reference_ns"] = stats.sum_exact_reference_ns;
    out["sum_stock_ns_on_escalated"] = stats.sum_stock_ns_on_escalated;
    out["max_total_ns"] = stats.max_total_ns;
    out["max_escalated_total_ns"] = stats.max_escalated_total_ns;
    out["max_residual_size"] = stats.max_residual_size;
    out["max_largest_tree_size"] = stats.max_largest_tree_size;
    out["sum_residual_size"] = stats.sum_residual_size;
    out["sum_num_defects"] = stats.sum_num_defects;
    out["residual_size_hist"] = stats.residual_size_hist;

    // §C.3's raw accumulators, so a consumer can re-derive `summarize`'s fractions or pool two
    // campaigns without going back to the per-shot rows.
    out["shots_with_component_stats"] = stats.shots_with_component_stats;
    out["shots_with_component_defects"] = stats.shots_with_component_defects;
    out["shots_solver_set_empty"] = stats.shots_solver_set_empty;
    out["shots_with_trivial_residual"] = stats.shots_with_trivial_residual;
    out["sum_num_components"] = stats.sum_num_components;
    out["sum_trivial_components"] = stats.sum_trivial_components;
    out["sum_singleton_components"] = stats.sum_singleton_components;
    out["sum_pair_components"] = stats.sum_pair_components;
    out["sum_nontrivial_components"] = stats.sum_nontrivial_components;
    out["sum_boundary_touching_components"] = stats.sum_boundary_touching_components;
    out["sum_diameter_uncomputed_components"] = stats.sum_diameter_uncomputed_components;
    out["sum_component_defects"] = stats.sum_component_defects;
    out["sum_defects_in_trivial_components"] = stats.sum_defects_in_trivial_components;
    out["sum_defects_committed_trivially"] = stats.sum_defects_committed_trivially;
    out["sum_defects_residual_trivially"] = stats.sum_defects_residual_trivially;
    out["sum_defects_to_solver"] = stats.sum_defects_to_solver;
    out["sum_largest_component_size"] = stats.sum_largest_component_size;
    out["max_component_size"] = stats.max_component_size;
    out["sum_max_component_diameter_wint"] = stats.sum_max_component_diameter_wint;
    out["max_component_diameter_wint"] = stats.max_component_diameter_wint;
    out["weight_hist_bins_per_T"] = (uint64_t)ComponentHistograms::BINS_PER_T;
    out["component_size_hist"] = stats.component_hist.size_hist;
    out["component_diameter_hist"] = stats.component_hist.diameter_hist;
    out["h_degree_hist"] = stats.component_hist.degree_hist;
    out["h_edge_weight_hist"] = stats.component_hist.edge_weight_hist;
    out["boundary_cost_hist"] = stats.component_hist.bcost_hist;
    return out;
}

py::dict ball_stats_to_dict(const BallStats& stats) {
    py::dict out;
    out["num_nodes"] = stats.num_nodes;
    out["total_entries"] = stats.total_entries;
    out["mean_ball_size"] = stats.mean_ball_size;
    out["max_ball_size"] = stats.max_ball_size;
    out["shell_entry_counts"] = stats.shell_entry_counts;
    out["nodes_with_boundary"] = stats.nodes_with_boundary;
    out["mean_ball_word_len"] = stats.mean_ball_word_len;
    out["max_ball_word_len"] = stats.max_ball_word_len;
    out["ambiguous_mask_pairs"] = stats.ambiguous_mask_pairs;
    out["bytes_targets"] = stats.bytes_targets;
    out["bytes_weights"] = stats.bytes_weights;
    out["bytes_masks"] = stats.bytes_masks;
    out["bytes_paths"] = stats.bytes_paths;
    out["bytes_words"] = stats.bytes_words;
    out["bytes_boundary"] = stats.bytes_boundary;
    out["bytes_offsets"] = stats.bytes_offsets;
    out["bytes_total"] = stats.bytes_total;
    out["compile_wall_seconds"] = stats.compile_wall_seconds;
    return out;
}

/// The per-shot profile as a dict of arrays, row-aligned with the shots (§M6.3). Columns rather
/// than a list of objects, because every consumer of this is a plot or a percentile.
py::dict profiles_to_dict(const std::vector<TwoPhaseProfile>& profiles) {
    size_t n = profiles.size();
    auto column_i64 = [&](long long TwoPhaseProfile::* field) {
        py::array_t<long long> array((py::ssize_t)n);
        auto view = array.mutable_unchecked<1>();
        for (size_t i = 0; i < n; i++)
            view((py::ssize_t)i) = profiles[i].*field;
        return array;
    };
    auto column_int = [&](int TwoPhaseProfile::* field) {
        py::array_t<int> array((py::ssize_t)n);
        auto view = array.mutable_unchecked<1>();
        for (size_t i = 0; i < n; i++)
            view((py::ssize_t)i) = profiles[i].*field;
        return array;
    };
    /// §C.1's fields live in a nested struct so that `BallProfile` and `TwoPhaseProfile` carry one
    /// definition between them; the columns are flattened here, where the dict keys are chosen.
    auto column_component = [&](int ComponentStats::* field) {
        py::array_t<int> array((py::ssize_t)n);
        auto view = array.mutable_unchecked<1>();
        for (size_t i = 0; i < n; i++)
            view((py::ssize_t)i) = profiles[i].components.*field;
        return array;
    };
    auto column_bool = [&](bool TwoPhaseProfile::* field) {
        py::array_t<bool> array((py::ssize_t)n);
        auto view = array.mutable_unchecked<1>();
        for (size_t i = 0; i < n; i++)
            view((py::ssize_t)i) = profiles[i].*field;
        return array;
    };

    py::dict out;
    out["phase1_ns"] = column_i64(&TwoPhaseProfile::phase1_ns);
    out["harvest_ns"] = column_i64(&TwoPhaseProfile::harvest_ns);
    out["dual_scan_ns"] = column_i64(&TwoPhaseProfile::dual_scan_ns);
    out["escalation_ns"] = column_i64(&TwoPhaseProfile::escalation_ns);
    out["stock_ns"] = column_i64(&TwoPhaseProfile::stock_ns);
    out["total_ns"] = column_i64(&TwoPhaseProfile::total_ns);
    out["exact_reference_ns"] = column_i64(&TwoPhaseProfile::exact_reference_ns);
    out["num_defects"] = column_int(&TwoPhaseProfile::num_defects);
    out["residual_size"] = column_int(&TwoPhaseProfile::residual_size);
    out["num_trees"] = column_int(&TwoPhaseProfile::num_trees);
    out["committed_pairs_frozen"] = column_int(&TwoPhaseProfile::committed_pairs_frozen);
    out["committed_pairs_tree"] = column_int(&TwoPhaseProfile::committed_pairs_tree);
    out["committed_boundary"] = column_int(&TwoPhaseProfile::committed_boundary);
    out["exposed_root_blossoms"] = column_int(&TwoPhaseProfile::exposed_root_blossoms);
    out["certified"] = column_int(&TwoPhaseProfile::certified);
    out["h_no_perfect_matching"] = column_int(&TwoPhaseProfile::h_no_perfect_matching);
    out["truncated"] = column_bool(&TwoPhaseProfile::truncated);
    out["escalated"] = column_bool(&TwoPhaseProfile::escalated);
    out["contaminated"] = column_bool(&TwoPhaseProfile::contaminated);
    out["truncated_reference_escalates"] = column_bool(&TwoPhaseProfile::truncated_reference_escalates);
    out["truncated_reference_measured"] = column_bool(&TwoPhaseProfile::truncated_reference_measured);

    // §C.1, row-aligned with the shots like everything else here. `components_measured` is 0 on
    // every row unless `collect_component_stats` was on, which is how a reader tells "not
    // collected" from "collected, and this shot had no defects".
    out["components_measured"] = column_component(&ComponentStats::measured);
    out["num_components"] = column_component(&ComponentStats::num_components);
    out["num_trivial_components"] = column_component(&ComponentStats::num_trivial_components);
    out["num_singleton_components"] = column_component(&ComponentStats::num_singleton_components);
    out["num_pair_components"] = column_component(&ComponentStats::num_pair_components);
    out["num_nontrivial_components"] = column_component(&ComponentStats::num_nontrivial_components);
    out["defects_in_trivial_components"] = column_component(&ComponentStats::defects_in_trivial_components);
    out["largest_component_size"] = column_component(&ComponentStats::largest_component_size);
    out["max_component_diameter_wint"] = column_component(&ComponentStats::max_component_diameter_wint);
    out["diameter_uncomputed_components"] = column_component(&ComponentStats::diameter_uncomputed_components);
    out["num_boundary_touching_components"] = column_component(&ComponentStats::num_boundary_touching_components);
    out["defects_committed_trivially"] = column_component(&ComponentStats::defects_committed_trivially);
    out["defects_residual_trivially"] = column_component(&ComponentStats::defects_residual_trivially);
    out["defects_to_solver"] = column_component(&ComponentStats::defects_to_solver);

    py::array_t<double> dual_sum((py::ssize_t)n);
    py::array_t<double> weight_out((py::ssize_t)n);
    py::array_t<double> max_dual((py::ssize_t)n);
    auto dual_view = dual_sum.mutable_unchecked<1>();
    auto weight_view = weight_out.mutable_unchecked<1>();
    auto max_dual_view = max_dual.mutable_unchecked<1>();
    for (size_t i = 0; i < n; i++) {
        dual_view((py::ssize_t)i) = (double)profiles[i].dual_sum_at_truncation;
        weight_view((py::ssize_t)i) = (double)profiles[i].weight_out;
        max_dual_view((py::ssize_t)i) = (double)profiles[i].max_dual_at_completion;
    }
    out["dual_sum_at_truncation"] = dual_sum;
    out["weight_out"] = weight_out;
    out["max_dual_at_completion"] = max_dual;
    return out;
}

}  // namespace

void pm_pybind::pybind_two_phase(py::module& m) {
    auto config = py::class_<TwoPhaseConfig>(m, "TwoPhaseConfig", R"pbdoc(
Configuration of the two-phase decoder.

`T` is the truncation horizon in DEM float weight units; `ball_T_max` and `ball_R` size the
compiled ball tables and must satisfy `ball_R >= 2 * ball_T_max` (§M2.0). `T` must not exceed
`ball_T_max`.
)pbdoc");
    config.def(py::init<>());
    config.def_readwrite("T", &TwoPhaseConfig::T);
    config.def_property(
        "ball_T_max",
        [](const TwoPhaseConfig& self) {
            return self.ball.T_max;
        },
        [](TwoPhaseConfig& self, double value) {
            self.ball.T_max = value;
        });
    config.def_property(
        "ball_R",
        [](const TwoPhaseConfig& self) {
            return self.ball.R;
        },
        [](TwoPhaseConfig& self, double value) {
            self.ball.R = value;
        });
    config.def_property(
        "shell_width",
        [](const TwoPhaseConfig& self) {
            return self.ball.shell_width;
        },
        [](TwoPhaseConfig& self, double value) {
            self.ball.shell_width = value;
        });
    config.def_property(
        "store_paths",
        [](const TwoPhaseConfig& self) {
            return self.ball.store_paths;
        },
        [](TwoPhaseConfig& self, bool value) {
            self.ball.store_paths = value;
        });
    config.def_property(
        "certify_masks",
        [](const TwoPhaseConfig& self) {
            return self.ball.certify_masks;
        },
        [](TwoPhaseConfig& self, bool value) {
            self.ball.certify_masks = value;
        });
    config.def_readwrite("phase1_on_ball_graph", &TwoPhaseConfig::phase1_on_ball_graph);
    config.def_readwrite("stock_on_h", &TwoPhaseConfig::stock_on_h);
    config.def_readwrite("measure_truncated_reference", &TwoPhaseConfig::measure_truncated_reference);
    config.def_readwrite("unbounded_horizon", &TwoPhaseConfig::unbounded_horizon);
    config.def_readwrite("full_harvest_for_verification", &TwoPhaseConfig::full_harvest_for_verification);
    config.def_readwrite("edges_flavor", &TwoPhaseConfig::edges_flavor);
    config.def_readwrite("compile_threads", &TwoPhaseConfig::compile_threads);
    config.def_readwrite("verify_against_g", &TwoPhaseConfig::verify_against_g);
    config.def_readwrite("collect_harvest_diagnostics", &TwoPhaseConfig::collect_harvest_diagnostics);
    config.def_readwrite("collect_structural_counters", &TwoPhaseConfig::collect_structural_counters);
    config.def_readwrite("collect_component_stats", &TwoPhaseConfig::collect_component_stats);
    config.def_readwrite("prune_trivial_components", &TwoPhaseConfig::prune_trivial_components);
    config.def_readwrite("trivial_component_max_size", &TwoPhaseConfig::trivial_component_max_size);
    config.def_readwrite(
        "skip_negative_weight_preamble_when_positive", &TwoPhaseConfig::skip_negative_weight_preamble_when_positive);
    config.def_readwrite("measure_exact_reference", &TwoPhaseConfig::measure_exact_reference);
    config.def_readwrite("detect_preemption", &TwoPhaseConfig::detect_preemption);
    config.def_property(
        "build_mode",
        [](const TwoPhaseConfig& self) {
            return self.mode == BallGraphBuildMode::BITSET ? "bitset" : "scan";
        },
        [](TwoPhaseConfig& self, const std::string& value) {
            if (value == "bitset") {
                self.mode = BallGraphBuildMode::BITSET;
            } else if (value == "scan") {
                self.mode = BallGraphBuildMode::SCAN;
            } else {
                throw std::invalid_argument("build_mode must be 'scan' or 'bitset'");
            }
        });

    auto stats = py::class_<TwoPhaseAggregateStats>(m, "TwoPhaseAggregateStats");
    stats.def("as_dict", &stats_to_dict);
    stats.def("reset", &TwoPhaseAggregateStats::reset);

    m.def(
        "summarize",
        [](const TwoPhaseAggregateStats& stats) {
            return summary_to_dict(summarize(stats));
        },
        "stats"_a,
        R"pbdoc(
The derived quantities of §M6.2, computed once in C++ rather than re-derived per benchmark script.

`q` is always returned beside `shots`: a zero below `1 / shots` is the resolution floor of the
campaign, not a measurement. `amortisation_gap` is the reconciliation the M6 exit checkpoint asks
for — `|amortised_mean - measured_mean| / measured_mean`, which a gap in unattributed cost shows up
in directly.

With `collect_component_stats`, the result also carries the component structure of `H`: the size,
weighted-diameter, degree, `H` edge-weight and boundary-cost histograms, and the fractions of the
defect set that sit in a trivial (size <= 2) component and that the solver would still have to
take. The three weight histograms share one axis — bin `k` counts `[k * T / weight_hist_bins_per_T,
(k + 1) * T / weight_hist_bins_per_T)`, last bin overflowing — and `shots_with_component_stats` is
0 when the campaign did not collect any of it, which is a different statement from a zero mean.
)pbdoc");

    auto decoder = py::class_<TwoPhaseDecoder>(m, "TwoPhaseDecoder", R"pbdoc(
Two-phase decoder: truncated sparse blossom on the defect manifold, with stock exact decode as the
residual fallback. The output is exact MWPM on every shot.
)pbdoc");
    decoder.def_static(
        "from_detector_error_model",
        [](const std::string& dem_text,
           const TwoPhaseConfig& config,
           pm::weight_int num_distinct_weights,
           const std::string& ball_artifact_path) {
            return TwoPhaseDecoder::from_detector_error_model(
                dem_from_text(dem_text),
                config,
                num_distinct_weights,
                ball_artifact_path.empty() ? nullptr : ball_artifact_path.c_str());
        },
        "detector_error_model"_a,
        "config"_a,
        "num_distinct_weights"_a = pm::NUM_DISTINCT_WEIGHTS,
        "ball_artifact_path"_a = std::string(),
        "The detector error model is passed as text; use `str(dem)`.");
    decoder.def_readonly("num_observables", &TwoPhaseDecoder::num_observables);
    decoder.def("save_ball_artifact", &TwoPhaseDecoder::save_ball_artifact, "path"_a);
    decoder.def(
        "ball_stats",
        [](const TwoPhaseDecoder& self) {
            if (self.ball == nullptr)
                throw std::invalid_argument("this decoder has no ball tables (phase1_on_ball_graph = False)");
            return ball_stats_to_dict(self.ball->tables.stats);
        },
        "Ball table cost: bytes per array, ball sizes, and compile wall time (§M2.2).");
    decoder.def(
        "get_aggregate_stats",
        [](TwoPhaseDecoder& self) {
            return self.stats;
        },
        "A copy of the campaign accumulator. Pass it to `summarize`.");
    decoder.def("reset_aggregate_stats", [](TwoPhaseDecoder& self) {
        self.stats.reset();
    });

    decoder.def(
        "decode_to_obs",
        [](TwoPhaseDecoder& self, const std::vector<uint64_t>& dets) {
            auto* obs = new std::vector<uint8_t>(self.num_observables, 0);
            pm::total_weight_int weight = 0;
            self.decode_to_obs(dets, obs->data(), weight, nullptr);
            return py::make_tuple(vec_to_array(obs), (double)weight / self.normalising_constant());
        },
        "detection_events"_a,
        "Returns `(observables, weight)`, with the weight in DEM float units as the inherited "
        "`Matching.decode` reports it.");
    decoder.def(
        "decode_to_edges",
        [](TwoPhaseDecoder& self, const std::vector<uint64_t>& dets) {
            auto* edges = new std::vector<int64_t>();
            pm::total_weight_int weight = 0;
            self.decode_to_edges(dets, *edges, weight, nullptr);
            py::array_t<int64_t> array = vec_to_array(edges);
            array.resize({(py::ssize_t)(array.size() / 2), (py::ssize_t)2});
            return py::make_tuple(array, (double)weight / self.normalising_constant());
        },
        "detection_events"_a,
        R"pbdoc(
Returns `(edges, weight)`, where `edges` is an `(n, 2)` array of detector ids with `-1` for the
boundary — the inherited `decode_to_edges` contract. Treat it as a set.

`weight` is the *matching* weight, i.e. the sum over committed pairs, which is the same number the
obs flavour reports. It is not the sum of the emitted edges' weights: two committed pairs whose
paths share an edge cancel it out of the correction, exactly as stock's own edge decode does.
)pbdoc");
    decoder.def(
        "decode_batch",
        [](TwoPhaseDecoder& self, const std::vector<std::vector<uint64_t>>& shots, bool profile, bool aggregate) {
            auto* obs = new std::vector<uint8_t>(shots.size() * self.num_observables, 0);
            auto* weights = new std::vector<double>(shots.size(), 0);
            std::vector<pm::total_weight_int> raw_weights(shots.size(), 0);
            std::vector<TwoPhaseProfile> profiles;
            if (profile)
                profiles.reserve(shots.size());

            if (!aggregate)
                self.stats.reset();
            self.decode_batch(shots, obs->data(), raw_weights.data(), profile, profile ? &profiles : nullptr);
            for (size_t i = 0; i < shots.size(); i++)
                (*weights)[i] = (double)raw_weights[i] / self.normalising_constant();

            py::array_t<uint8_t> obs_array = vec_to_array(obs);
            if (self.num_observables > 0)
                obs_array.resize({(py::ssize_t)shots.size(), (py::ssize_t)self.num_observables});
            py::object profile_object = py::none();
            if (profile)
                profile_object = profiles_to_dict(profiles);
            return py::make_tuple(obs_array, vec_to_array(weights), profile_object);
        },
        "shots"_a,
        py::kw_only(),
        "profile"_a = false,
        "aggregate"_a = false,
        R"pbdoc(
Decodes a batch and returns `(observables, weights, profile)`.

`observables` is `(shots, num_observables)`. `profile` is `None` unless requested, in which case it
is a dict of arrays row-aligned with the shots. `aggregate=True` keeps accumulating into the
decoder's campaign stats across calls instead of resetting them, which is what a long campaign that
should not materialise per-shot arrays wants: leave `profile=False` and read `get_aggregate_stats`.
)pbdoc");
}
