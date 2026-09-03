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

#include "specmatching/spec_matching/driver/spec_matching_decoding.h"

#include <algorithm>
#include <cassert>
#include <stdexcept>

#include "specmatching/sparse_blossom/driver/mwpm_decoding.h"

namespace pm {
namespace spec_matching {

namespace {

/// Cancels a correction's edge multiset over GF(2) and leaves it sorted.
///
/// A correction is a *subset* of `G`'s edges, so an edge crossed by two committed pairs' paths is
/// crossed zero times by the correction. Stock does this by flipping a marker bit on the search
/// graph and dropping whatever ends up unflipped; the ball front end has the whole multiset in
/// hand and can do it with a sort, which needs no search graph and is deterministic by
/// construction.
void cancel_edge_multiset(std::vector<int64_t>& edges) {
    size_t num_edges = edges.size() / 2;
    std::vector<std::pair<int64_t, int64_t>> pairs;
    pairs.reserve(num_edges);
    for (size_t i = 0; i < num_edges; i++) {
        int64_t u = edges[2 * i];
        int64_t v = edges[2 * i + 1];
        // The boundary end always sits second, so a boundary half-edge has one canonical form.
        if (v >= 0 && v < u)
            std::swap(u, v);
        pairs.emplace_back(u, v);
    }
    std::sort(pairs.begin(), pairs.end());

    edges.clear();
    size_t i = 0;
    while (i < pairs.size()) {
        size_t j = i;
        while (j < pairs.size() && pairs[j] == pairs[i])
            j++;
        if (((j - i) & 1) != 0) {
            edges.push_back(pairs[i].first);
            edges.push_back(pairs[i].second);
        }
        i = j;
    }
}

}  // namespace

pm::total_weight_int SpecMatchingDecoder::weight_of_match_edges(
    pm::Mwpm& mwpm, const std::vector<pm::CompressedEdge>& match_edges) {
    // The *matching* weight: the sum of the shortest-path weights between the matched detection
    // events, plus the graph's negative-weight offset. Not the weight of the correction's edge set,
    // which is smaller whenever two paths share an edge and cancel it.
    //
    // `extract_paths_from_match_edges` computes exactly that, and needs somewhere to XOR the
    // observables it crosses on the way; they are discarded here, because the caller has already
    // obtained them by whichever route its flavour uses.
    obs_sink.assign(std::max<size_t>(num_observables, 1), 0);
    pm::total_weight_int weight = 0;
    mwpm.extract_paths_from_match_edges(match_edges, obs_sink.data(), weight);
    return weight + mwpm.flooder.negative_weight_sum;
}

SpecMatchingDecoder SpecMatchingDecoder::from_detector_error_model(
    const stim::DetectorErrorModel& dem,
    SpecMatchingConfig config,
    pm::weight_int num_distinct_weights,
    const char* ball_artifact_path) {
    if (config.unbounded_horizon && config.phase1_on_ball_graph)
        throw std::invalid_argument(
            "T = infinity is not supported on the ball graph: with a finite R the ball tables cannot supply the "
            "edges an unbounded run needs (§M4.1). Route it through phase1_on_ball_graph = false.");
    if (config.verify_against_g && !config.phase1_on_ball_graph)
        throw std::invalid_argument("verify_against_g compares H against G; it needs phase1_on_ball_graph = true.");
    if (config.stock_on_h && !config.phase1_on_ball_graph)
        throw std::invalid_argument(
            "stock_on_h is a front end for the ball graph H (§M7); it needs phase1_on_ball_graph = true. Stock "
            "blossom on G with no horizon is simply the inherited decoder.");
    if (config.stock_on_h && (config.verify_against_g || config.full_harvest_for_verification))
        throw std::invalid_argument(
            "The §M7 stock-on-H path produces no truncated intermediate state, so the M1 harvest oracle has "
            "nothing to reproduce. Its oracle is stock exact decode on G (§M7.6 level 1).");
    if (config.measure_truncated_reference && !config.phase1_on_ball_graph)
        throw std::invalid_argument(
            "measure_truncated_reference replays the truncated scheme's decision on H; it needs "
            "phase1_on_ball_graph = true.");
    if (config.edges_flavor && config.phase1_on_ball_graph && !config.ball.store_paths)
        throw std::invalid_argument(
            "The edges flavour lifts a committed pair through the ball tables' stored path, so it needs "
            "BallParams::store_paths (§M5).");
    // The oracle is what level 1 is compared against; without the full harvest there is no residual
    // and no `num_trees` to compare (§M3.4, "keep the full harvest compiled in").
    if (config.verify_against_g)
        config.full_harvest_for_verification = true;

    SpecMatchingDecoder decoder;
    decoder.config = config;

    // The search graph is needed by the edges flavour, and by the ball tables whenever the
    // `obs_int` masks are unusable — above 64 observables they have to read observable ids off the
    // search graph instead.
    bool need_search_graph = config.edges_flavor || dem.count_observables() > sizeof(pm::obs_int) * 8;
    pm::Mwpm g_mwpm = pm::detector_error_model_to_mwpm(dem, num_distinct_weights, need_search_graph);
    decoder.num_observables = g_mwpm.flooder.graph.num_observables;
    double normalising_constant = g_mwpm.flooder.graph.normalising_constant;

    if (config.phase1_on_ball_graph) {
        BallConfig ball_config;
        ball_config.T = config.T;
        ball_config.ball = config.ball;
        ball_config.stock_on_h = config.stock_on_h;
        ball_config.verify_against_g = config.verify_against_g;
        ball_config.mode = config.mode;
        ball_config.compile_threads = config.compile_threads;
        ball_config.collect_harvest_diagnostics = config.collect_harvest_diagnostics;
        ball_config.collect_structural_counters = config.collect_structural_counters;
        ball_config.collect_component_stats = config.collect_component_stats;
        ball_config.diameter_cap = config.diameter_cap;
        ball_config.verify_component_decomposition = config.verify_component_decomposition;
        ball_config.skip_negative_weight_preamble_when_positive = config.skip_negative_weight_preamble_when_positive;
        decoder.ball =
            std::make_unique<BallDecoder>(BallDecoder::from_mwpm(std::move(g_mwpm), ball_config, ball_artifact_path));
        decoder.horizon = decoder.ball->horizon;
    } else {
        decoder.oracle_g_mwpm = std::move(g_mwpm);
        decoder.horizon = config.unbounded_horizon ? pm::NO_HORIZON : to_time_units(config.T, normalising_constant);
        decoder.oracle_harvester.collect_diagnostics = config.collect_harvest_diagnostics;
    }
    return decoder;
}

pm::Mwpm& SpecMatchingDecoder::g_mwpm() {
    return ball != nullptr ? ball->g_mwpm : oracle_g_mwpm;
}

const pm::Mwpm& SpecMatchingDecoder::g_mwpm() const {
    return ball != nullptr ? ball->g_mwpm : oracle_g_mwpm;
}

void SpecMatchingDecoder::save_ball_artifact(const std::string& path) const {
    if (ball == nullptr)
        throw std::invalid_argument("This decoder has no ball tables to save (phase1_on_ball_graph = false).");
    ball->save_ball_artifact(path);
}

Phase1Outcome SpecMatchingDecoder::run_phase1(const std::vector<uint64_t>& dets, bool need_pairs, SpecMatchingProfile* prof) {
    Phase1Outcome outcome;
    if (ball != nullptr) {
        BallProfile* bp = prof != nullptr ? &ball_profile : nullptr;
        if (config.full_harvest_for_verification) {
            // The verification path: §M1.3's harvest on every shot, residual and all, which is what
            // §M2.6 level 1 and §M3.3 X8 read. The status is then a *fact about the harvest* rather
            // than the O(1) branch, and the two are asserted equal below.
            HarvestResult harvest = need_pairs ? ball->decode_phase1_to_match_edges(dets, committed_pairs, bp)
                                               : ball->decode_phase1(dets, bp);
            outcome.status = harvest.residual.empty() ? TimelineStatus::COMPLETE : TimelineStatus::TRUNCATED;
            outcome.harvest = std::move(harvest);
            // §2.6's tally still comes from the decode, not from the residual: the residual is what
            // this branch checks the status *against*, so deriving both from it would make the
            // check vacuous.
            if (bp != nullptr) {
                outcome.components_total = ball_profile.components_total;
                outcome.components_truncated = ball_profile.components_truncated;
                assert(
                    (outcome.status == TimelineStatus::TRUNCATED) == (ball_profile.any_component_truncated != 0) &&
                    "the full harvest's residual and §2.6's escalation predicate disagree");
            }
        } else {
            outcome = need_pairs ? ball->decode_phase1_production_to_match_edges(dets, committed_pairs, bp)
                                 : ball->decode_phase1_production(dets, bp);
        }
        if (prof != nullptr) {
            // Phase 1 proper is everything the front end did except the extraction, which the
            // profile reports separately; `total_ns` is measured end to end and is not their sum.
            prof->phase1_ns = ball_profile.total_ns - ball_profile.harvest_ns;
            prof->harvest_ns = ball_profile.harvest_ns;
        }
        return outcome;
    }

    // The oracle front end: M1's decode on `G`. Always the full harvest — this path exists to be
    // the thing the ball path is checked against, so bypassing anything here would defeat it.
    HiResTimer step;
    if (prof != nullptr)
        step.start();
    outcome.status = process_timeline_until_horizon(oracle_g_mwpm, dets, horizon);
    if (prof != nullptr)
        prof->phase1_ns = step.elapsed_ns();

    if (prof != nullptr)
        step.start();
    if (need_pairs) {
        oracle_match_edges.clear();
        outcome.harvest = oracle_harvester.harvest_to_match_edges(oracle_g_mwpm, dets, oracle_match_edges);
    } else {
        outcome.harvest = oracle_harvester.harvest_to_obs(oracle_g_mwpm, dets);
    }
    if (prof != nullptr)
        prof->harvest_ns = step.elapsed_ns();

    assert((outcome.status == TimelineStatus::TRUNCATED) == !outcome.harvest.residual.empty());
    // The oracle front end solves `G` as one problem, so §2.6's tally reads it as the one component
    // it is. There is no `H` here and nothing to decompose; reporting zero components would make
    // the profile's `components_truncated > 0 <=> escalated` invariant false on this path alone.
    outcome.components_total = 1;
    outcome.components_truncated = outcome.status == TimelineStatus::TRUNCATED ? 1 : 0;
    return outcome;
}

void SpecMatchingDecoder::copy_phase1_stats(const Phase1Outcome& outcome, SpecMatchingProfile* prof) const {
    if (prof == nullptr)
        return;
    prof->fill_from(outcome.harvest);
    // §2.6's one predicate, carried out of Phase 1 rather than re-derived from a residual: on the
    // production path an escalating shot never builds one, so reading it here would report every
    // shot as complete. This is the thing invariant 12 is about.
    prof->any_component_truncated = outcome.status == TimelineStatus::TRUNCATED;
    prof->components_total = outcome.components_total;
    prof->components_truncated = outcome.components_truncated;
    assert(
        prof->any_component_truncated == (outcome.components_truncated > 0) &&
        "§2.6: Phase 1's status and its component tally disagree");
    if (ball != nullptr) {
        prof->blossom_formations = ball_profile.blossom_formations;
        prof->solve_dependent_depth = ball_profile.solve_dependent_depth;
        prof->solve_events = ball_profile.solve_events;
        // §M7.7. Zero on the truncated path, where there is no certificate to report.
        prof->dual_scan_ns = ball_profile.dual_scan_ns;
        prof->certified = ball_profile.certified;
        prof->h_no_perfect_matching = ball_profile.h_no_perfect_matching;
        prof->max_dual_at_completion = ball_profile.max_dual_at_completion;
    }
}

void SpecMatchingDecoder::configure_component_tables(size_t size_cap, size_t degree_cap) {
    component_histograms.configure(size_cap, degree_cap);
    component_size_x_status.configure(size_cap + 1);
    component_hop_diameter_x_status.configure(size_cap + 1);
    stats.configure_component_tables(size_cap, degree_cap);
}

void SpecMatchingDecoder::record_component_stats(SpecMatchingProfile* prof) {
    if (prof == nullptr || ball == nullptr || !config.collect_component_stats)
        return;
    // Deliberately after `total_ns` has been read, and deliberately untimed: the pruning experiment
    // this feeds assumes the component work is free, so charging it to a latency number here would
    // be measuring a stage that is not meant to run on this critical path. It also has to come
    // before `record_truncated_reference`, which rebuilds `H` in the same arena.
    ball->analyze_last_shot_components(
        ball_profile, component_histograms, component_size_x_status, component_hop_diameter_x_status);
    prof->components = ball_profile.components;
    assert(
        prof->components.components_truncated == prof->components_truncated &&
        "§2.6: the post-shot analysis and the decode disagree on how many components truncated");
}

void SpecMatchingDecoder::record_truncated_reference(const std::vector<uint64_t>& dets, SpecMatchingProfile* prof) {
    if (!config.measure_truncated_reference || ball == nullptr || prof == nullptr)
        return;
    // Deliberately after `total_ns` has been read: this is a second full Phase 1 on `H`, and folding
    // it into the timed window would make every benchmark number the sum of two schemes.
    prof->truncated_reference_measured = true;
    prof->truncated_reference_escalates = ball->truncated_scheme_escalates(dets);
}

void SpecMatchingDecoder::obs_from_committed_pairs(uint8_t* obs, pm::total_weight_int& weight) const {
    // Above 64 observables the `obs_int` mask is unusable, so the readout goes through the ball
    // tables' observable *id* lists — the same lists the masks were built from, so this is the same
    // answer by a route that does not care how many observables there are (§0).
    const BallTables& tables = ball->tables;
    std::fill(obs, obs + num_observables, (uint8_t)0);
    weight = 0;
    for (const CommittedPair& pair : committed_pairs) {
        if (pair.ball_entry == CommittedPair::NO_BALL_ENTRY) {
            size_t det = (size_t)pair.from;
            for (uint64_t k = tables.bcost_mask_offsets[det]; k < tables.bcost_mask_offsets[det + 1]; k++)
                obs[tables.bcost_mask_ids[k]] ^= 1;
            weight += (pm::total_weight_int)tables.bcost_w_int[det];
        } else {
            uint64_t e = pair.ball_entry;
            for (uint64_t k = tables.ball_mask_offsets[e]; k < tables.ball_mask_offsets[e + 1]; k++)
                obs[tables.ball_mask_ids[k]] ^= 1;
            weight += (pm::total_weight_int)tables.ball_w_int[e];
        }
    }
    for (uint64_t o : g_mwpm().flooder.negative_weight_observables)
        obs[o] ^= 1;
    weight += g_mwpm().flooder.negative_weight_sum;
}

void SpecMatchingDecoder::decode_to_obs(
    const std::vector<uint64_t>& dets, uint8_t* obs, pm::total_weight_int& weight, SpecMatchingProfile* prof) {
    PreemptionProbe before;
    HiResTimer total_timer;
    if (prof != nullptr) {
        prof->clear();
        prof->num_defects = (int)dets.size();
        if (config.detect_preemption)
            before.sample();
        total_timer.start();
    }

    bool wide = num_observables > sizeof(pm::obs_int) * 8;
    Phase1Outcome outcome = run_phase1(dets, wide, prof);
    copy_phase1_stats(outcome, prof);

    if (outcome.status == TimelineStatus::TRUNCATED) {
        // §M3.1. The whole syndrome, on `G`, with Phase 1's partial result discarded rather than
        // combined. Invariant 12: this fires iff the residual is non-empty, and nothing else.
        escalate_to_stock(g_mwpm(), dets, obs, num_observables, weight, prof);
    } else if (wide) {
        obs_from_committed_pairs(obs, weight);
    } else if (ball != nullptr) {
        pm::obs_int mask = outcome.harvest.committed.obs_mask ^ ball->negative_weight_obs_mask();
        std::fill(obs, obs + num_observables, (uint8_t)0);
        pm::fill_bit_vector_from_obs_mask(mask, obs, num_observables);
        weight = outcome.harvest.committed.weight + ball->negative_weight_sum();
    } else {
        pm::obs_int mask = outcome.harvest.committed.obs_mask ^ oracle_g_mwpm.flooder.negative_weight_obs_mask;
        std::fill(obs, obs + num_observables, (uint8_t)0);
        pm::fill_bit_vector_from_obs_mask(mask, obs, num_observables);
        weight = outcome.harvest.committed.weight + oracle_g_mwpm.flooder.negative_weight_sum;
    }

    if (prof != nullptr) {
        prof->weight_out = weight;
        prof->total_ns = total_timer.elapsed_ns();
        if (config.detect_preemption) {
            PreemptionProbe after;
            after.sample();
            prof->contaminated = after.switched_since(before);
        }
        // Invariant 12, and §2.6's "one escalation predicate, in one place": the branch that was
        // taken and the predicate it was taken on are the same fact under two names.
        assert(
            prof->escalated == prof->any_component_truncated &&
            "invariant 12: escalation fires iff some component of H did not resolve within T, and nothing else");
        record_component_stats(prof);
        record_truncated_reference(dets, prof);
    }
}

void SpecMatchingDecoder::decode_to_edges(
    const std::vector<uint64_t>& dets,
    std::vector<int64_t>& edges,
    pm::total_weight_int& weight,
    SpecMatchingProfile* prof) {
    if (!config.edges_flavor)
        throw std::invalid_argument(
            "decode_to_edges needs SpecMatchingConfig::edges_flavor, which builds the search graph.");

    PreemptionProbe before;
    HiResTimer total_timer;
    if (prof != nullptr) {
        prof->clear();
        prof->num_defects = (int)dets.size();
        if (config.detect_preemption)
            before.sample();
        total_timer.start();
    }

    Phase1Outcome outcome = run_phase1(dets, true, prof);
    copy_phase1_stats(outcome, prof);
    edges.clear();
    weight = 0;

    if (outcome.status == TimelineStatus::TRUNCATED) {
        escalate_to_stock_edges(g_mwpm(), dets, edges, prof);
        // Stock's `to_edges` reports no weight, so it is read back off the *matched pairs* the
        // expansion consumed — which `decode_detection_events_to_edges` leaves in place — and not
        // off the emitted edge set. Those are different numbers whenever two paths share an edge:
        // the edge cancels out of the correction but the pairs still cost what they cost. The
        // matching weight is what the obs flavour reports, so it is what this reports too.
        weight = weight_of_match_edges(g_mwpm(), g_mwpm().flooder.match_edges);
    } else if (ball != nullptr) {
        // Non-escalating shots: Phase 1's committed match edges, lifted to `G`'s edges through the
        // ball tables' own stored canonical paths (§M5). There is no portal lift and no stub-path
        // append — both went with M3's compression.
        const BallTables& tables = ball->tables;
        for (const CommittedPair& pair : committed_pairs) {
            if (pair.ball_entry == CommittedPair::NO_BALL_ENTRY) {
                size_t det = (size_t)pair.from;
                uint64_t begin = tables.bcost_path_offsets[det];
                uint64_t end = tables.bcost_path_offsets[det + 1];
                assert(end > begin && "a boundary-matched defect with no stored path; needs store_paths");
                for (uint64_t k = begin; k + 1 < end; k++) {
                    edges.push_back((int64_t)tables.bcost_path_nodes[k]);
                    edges.push_back((int64_t)tables.bcost_path_nodes[k + 1]);
                }
                // The last node of the stored path is the one carrying the boundary half-edge.
                edges.push_back((int64_t)tables.bcost_path_nodes[end - 1]);
                edges.push_back(-1);
                weight += (pm::total_weight_int)tables.bcost_w_int[det];
            } else {
                uint64_t e = pair.ball_entry;
                uint64_t begin = tables.ball_path_offsets[e];
                uint64_t end = tables.ball_path_offsets[e + 1];
                assert(end > begin && "a committed pair with no stored path; needs store_paths");
                for (uint64_t k = begin; k + 1 < end; k++) {
                    edges.push_back((int64_t)tables.ball_path_nodes[k]);
                    edges.push_back((int64_t)tables.ball_path_nodes[k + 1]);
                }
                weight += (pm::total_weight_int)tables.ball_w_int[e];
            }
        }
        // `G`'s own negative-weight edges are part of every correction, exactly as stock's
        // `decode_detection_events_to_edges` adds them before cancelling.
        for (const auto& neg : g_mwpm().search_flooder.graph.negative_weight_edges) {
            edges.push_back((int64_t)neg.first);
            edges.push_back(neg.second != SIZE_MAX ? (int64_t)neg.second : -1);
        }
        cancel_edge_multiset(edges);
        weight += g_mwpm().flooder.negative_weight_sum;
    } else {
        // Oracle front end. The match edges came off a *truncated* timeline, so stock's own
        // completed-decode entry point cannot be used — but its expansion and cancellation pass
        // can, and is (§M5, `expand_match_edges_to_edges`).
        oracle_g_mwpm.flooder.match_edges = oracle_match_edges;
        pm::expand_match_edges_to_edges(oracle_g_mwpm, edges);
        // `harvest_to_match_edges` does not fill `HarvestResult::committed` — that is the obs
        // flavour's field — so the weight comes off the match edges, the same way it does on the
        // escalating branch.
        weight = weight_of_match_edges(oracle_g_mwpm, oracle_match_edges);
    }

    if (prof != nullptr) {
        prof->weight_out = weight;
        prof->total_ns = total_timer.elapsed_ns();
        if (config.detect_preemption) {
            PreemptionProbe after;
            after.sample();
            prof->contaminated = after.switched_since(before);
        }
        // Invariant 12, and §2.6's "one escalation predicate, in one place": the branch that was
        // taken and the predicate it was taken on are the same fact under two names.
        assert(
            prof->escalated == prof->any_component_truncated &&
            "invariant 12: escalation fires iff some component of H did not resolve within T, and nothing else");
        record_component_stats(prof);
        record_truncated_reference(dets, prof);
    }
}

void SpecMatchingDecoder::decode_batch(
    const std::vector<std::vector<uint64_t>>& shots,
    uint8_t* obs_out,
    pm::total_weight_int* weights_out,
    bool profile,
    std::vector<SpecMatchingProfile>* profiles_out) {
    std::vector<uint8_t> scratch_obs;
    if (obs_out == nullptr)
        scratch_obs.resize(num_observables);
    SpecMatchingProfile local;
    pm::total_weight_int local_weight = 0;

    for (size_t i = 0; i < shots.size(); i++) {
        uint8_t* obs = obs_out != nullptr ? obs_out + i * num_observables : scratch_obs.data();
        pm::total_weight_int& weight = weights_out != nullptr ? weights_out[i] : local_weight;
        SpecMatchingProfile* prof = profile ? &local : nullptr;

        decode_to_obs(shots[i], obs, weight, prof);

        if (profile) {
            if (config.measure_exact_reference) {
                // A paired measurement on the same shot, so `speedup_vs_stock` is a ratio of two
                // numbers taken under the same conditions rather than of two campaigns.
                HiResTimer exact_timer;
                std::vector<uint8_t> reference(num_observables, 0);
                pm::total_weight_int reference_weight = 0;
                exact_timer.start();
                pm::decode_detection_events(g_mwpm(), shots[i], reference.data(), reference_weight, false);
                local.exact_reference_ns = exact_timer.elapsed_ns();
            }
            stats.accumulate(local);
            // §3.5.2. The scalars ride in on the profile; the distributions and the joint tables are
            // per component and per edge, so they are folded in from the decoder's own per-shot
            // buffers instead.
            if (config.collect_component_stats && local.components.measured) {
                stats.accumulate_component_histograms(component_histograms);
                stats.accumulate_component_status_tables(component_size_x_status, component_hop_diameter_x_status);
            }
            if (profiles_out != nullptr)
                profiles_out->push_back(local);
        }
    }
}

}  // namespace spec_matching
}  // namespace pm
