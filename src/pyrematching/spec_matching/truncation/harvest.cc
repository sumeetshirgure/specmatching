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

#include "pyrematching/spec_matching/truncation/harvest.h"

#include <algorithm>
#include <cassert>
#include <sstream>

#include "pyrematching/perf/thread_timer.h"
#include "pyrematching/sparse_blossom/flooder/graph.h"
#include "pyrematching/sparse_blossom/matcher/alternating_tree.h"
#include "pyrematching/sparse_blossom/tracker/flood_check_event.h"

namespace pm {
namespace spec_matching {

namespace {

/// The same thread-scoped stopwatch the rest of the profile uses, aliased locally so that
/// `harvest.h` stays free of the profile headers (`spec_matching_profile.h` includes it, not the other
/// way around). Harvest's stage split is summed against `harvest_ns`, which is timed by
/// `HiResTimer`, so the two must be the same clock or the reconciliation is meaningless.
using Stopwatch = pm::perf::ThreadTimer;

/// Nesting depth and member count of a region's blossom structure. Depth 0 is a plain region.
void measure_blossom(const pm::GraphFillRegion& region, int& depth_out, int& members_out) {
    if (region.blossom_children.empty()) {
        depth_out = 0;
        members_out = 1;
        return;
    }
    int deepest = 0;
    int members = 0;
    for (const pm::RegionEdge& child : region.blossom_children) {
        int child_depth = 0;
        int child_members = 0;
        measure_blossom(*child.region, child_depth, child_members);
        deepest = std::max(deepest, child_depth);
        members += child_members;
    }
    depth_out = deepest + 1;
    members_out = members;
}

/// `ceil(log2(n))`, the depth of a balanced reduction tree over `n` inputs. Zero for `n <= 1`.
int reduction_depth(size_t n) {
    int depth = 0;
    while ((size_t)1 << depth < n)
        depth++;
    return depth;
}

/// Which endpoint of a match drives its extraction.
///
/// A match is between *regions*, and both endpoints see the same pair, so exactly one of them has
/// to be chosen or the pair would be extracted twice (and the second attempt would touch regions
/// the first one already returned to the arena). M1.3 made that choice implicitly, by whichever of
/// the two endpoints the detection-event sweep reached first. With the sweep gone the choice has to
/// be explicit, and it is made locally and canonically: `match.edge.loc_from` lies in this region
/// and `match.edge.loc_to` in the partner, so comparing the two node pointers is comparing their
/// detector ids, and the endpoint holding the lower one wins. A boundary match has no partner and
/// is always its own representative.
inline bool is_extraction_representative(const pm::GraphFillRegion& region) {
    if (region.match.region == nullptr)
        return true;
    return region.match.edge.loc_from < region.match.edge.loc_to;
}

#ifndef NDEBUG
/// §M2.9.2's guard, debug builds only: the set of top-level regions the arena's live list yields
/// must be exactly the set M1.3's detection-event sweep would have produced.
///
/// This is the cheapest possible check against the failure mode the restructuring introduces, and
/// that failure mode is silent — a region missing from the enumeration drops committed pairs and
/// inflates the residual without any crash, and both flavours would still return a well-formed
/// answer. Called before anything is mutated, so both sides see the same primal state.
void assert_live_top_regions_match_detection_event_sweep(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    auto& flooder = mwpm.flooder;
    auto& nodes = flooder.graph.nodes;

    std::unordered_set<pm::GraphFillRegion*> from_live_list;
    for (pm::GraphFillRegion* region : flooder.region_arena.live) {
        if (region->blossom_parent == nullptr)
            from_live_list.insert(region);
    }

    std::unordered_set<pm::GraphFillRegion*> from_sweep;
    const std::vector<uint64_t>* event_lists[2] = {&detection_events, &flooder.negative_weight_detection_events};
    for (const std::vector<uint64_t>* events : event_lists) {
        for (uint64_t det : *events) {
            if (det >= nodes.size())
                continue;
            pm::DetectorNode& node = nodes[det];
            if (node.region_that_arrived == nullptr)
                continue;
            from_sweep.insert(node.region_that_arrived_top);
        }
    }

    assert(
        from_live_list == from_sweep &&
        "§M2.9.2: the arena live list and the detection-event sweep disagree on the top-level regions");
}
#endif

}  // namespace

void HarvestResult::clear() {
    residual.clear();
    residual_dual_sum.clear();
    committed = pm::MatchingResult();
    committed_pairs_frozen = 0;
    committed_pairs_tree = 0;
    committed_pairs_blossom_cycle = 0;
    committed_boundary = 0;
    num_trees = 0;
    largest_tree_size = 0;
    exposed_root_blossoms = 0;
    dual_sum_at_truncation = 0;
    max_region_dual = 0;
    max_blossom_nesting_depth = 0;
    max_blossom_members = 0;
    max_exposed_blossom_depth = 0;
    max_exposed_blossom_members = 0;
    matched_blossom_shatters = 0;
    harvest_dependent_depth = 0;
    enumerate_ns = 0;
    reduce_ns = 0;
    base_descent_ns = 0;
    shatter_ns = 0;
}

bool HarvestResult::identical_output_to(const HarvestResult& other) const {
    return residual == other.residual && residual_dual_sum == other.residual_dual_sum && committed == other.committed &&
           committed_pairs_frozen == other.committed_pairs_frozen &&
           committed_pairs_tree == other.committed_pairs_tree &&
           committed_pairs_blossom_cycle == other.committed_pairs_blossom_cycle &&
           committed_boundary == other.committed_boundary && num_trees == other.num_trees &&
           largest_tree_size == other.largest_tree_size && exposed_root_blossoms == other.exposed_root_blossoms &&
           dual_sum_at_truncation == other.dual_sum_at_truncation && max_region_dual == other.max_region_dual;
}

std::string HarvestResult::describe_difference(const HarvestResult& other) const {
    std::ostringstream out;
    auto compare_scalar = [&](const char* name, long long lhs, long long rhs) {
        if (lhs != rhs)
            out << " " << name << "(" << lhs << " vs " << rhs << ")";
    };
    if (residual != other.residual) {
        out << " residual([";
        for (size_t i = 0; i < residual.size(); i++)
            out << (i ? "," : "") << residual[i];
        out << "] vs [";
        for (size_t i = 0; i < other.residual.size(); i++)
            out << (i ? "," : "") << other.residual[i];
        out << "])";
    }
    if (residual_dual_sum != other.residual_dual_sum)
        out << " residual_dual_sum";
    compare_scalar("committed.obs_mask", (long long)committed.obs_mask, (long long)other.committed.obs_mask);
    compare_scalar("committed.weight", (long long)committed.weight, (long long)other.committed.weight);
    compare_scalar("committed_pairs_frozen", committed_pairs_frozen, other.committed_pairs_frozen);
    compare_scalar("committed_pairs_tree", committed_pairs_tree, other.committed_pairs_tree);
    compare_scalar("committed_pairs_blossom_cycle", committed_pairs_blossom_cycle, other.committed_pairs_blossom_cycle);
    compare_scalar("committed_boundary", committed_boundary, other.committed_boundary);
    compare_scalar("num_trees", num_trees, other.num_trees);
    compare_scalar("largest_tree_size", largest_tree_size, other.largest_tree_size);
    compare_scalar("exposed_root_blossoms", exposed_root_blossoms, other.exposed_root_blossoms);
    compare_scalar(
        "dual_sum_at_truncation", (long long)dual_sum_at_truncation, (long long)other.dual_sum_at_truncation);
    compare_scalar("max_region_dual", (long long)max_region_dual, (long long)other.max_region_dual);
    std::string text = out.str();
    return text.empty() ? std::string("(no difference)") : text;
}

void HarvestScratch::clear() {
    top_regions.clear();
    live_nodes.clear();
    seen_regions.clear();
    tree_roots.clear();
    seen_roots.clear();
    tree_walk_stack.clear();
    residual_sort_buffer.clear();
    tree_sizes.clear();
}

/// §M2.9.1–§M2.9.3. Direct enumeration, local readout, associative reduction.
///
/// Three sweeps, none of which follows the structure it is reading:
///
///  - **A.** One walk of the region arena's live list. Every region alive at truncation is visited
///    exactly once, so the dual sum and its maximum are flat reductions rather than a descent
///    through blossom nesting; and a top-level region is recognised by `blossom_parent == nullptr`,
///    a register read, so the detection-event sweep and its visited stamp both disappear (§M2.9.2).
///  - **B.** One walk of the node arena's live list. `inner_region == nullptr` identifies a root —
///    one comparator, no `find_root` — and every other node commits its own pair from its own
///    `inner_to_outer_edge`. `children` is never touched, and neither is `AltTreeNode::parent`
///    except under diagnostics.
///  - **C.** Extraction, driven from the representative endpoint of each committed match.
///
/// The commit policy is M1.3's, unchanged and binding: the residual is exactly the exposed defects,
/// one per surviving tree, and everything else is committed. What changed is only how the committed
/// set is enumerated.
template <typename ExtractMatched, typename ExtractExposed>
HarvestResult Harvester::harvest_impl(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    const ExtractMatched& extract_matched,
    const ExtractExposed& extract_exposed) {
    if (use_legacy_enumeration)
        return harvest_impl_legacy(mwpm, detection_events, extract_matched, extract_exposed);

    counters.full_harvests++;
    HarvestResult result;
    scratch.clear();

    auto& flooder = mwpm.flooder;
    auto& nodes = flooder.graph.nodes;
    const pm::DetectorNode* node_base = nodes.data();
    // `process_timeline_until_horizon` has already advanced the clock to exactly the horizon, so
    // every dual read below is the dual at `T`.
    const pm::cumulative_time_int time = flooder.queue.cur_time;

    Stopwatch stage;
    if (collect_diagnostics)
        stage.start();

    // ---- Phase A: enumerate every live region, and classify the top-level ones.
    int frozen_half_pairs = 0;
    for (pm::GraphFillRegion* region : flooder.region_arena.live) {
        counters.regions_visited++;
        // Associative, and exact because the duals are integers, so the order this list happens to
        // be in cannot move the answer (§M2.9.3).
        pm::total_weight_int dual = region->radius.get_distance_at_time(time);
        result.dual_sum_at_truncation += dual;
        result.max_region_dual = std::max(result.max_region_dual, dual);

        if (region->blossom_parent != nullptr)
            continue;  // Nested inside a blossom; its top-level ancestor speaks for it.
        if (region->alt_tree_node != nullptr)
            continue;  // In a surviving tree; phase B owns it.

        if (collect_diagnostics && !region->blossom_children.empty()) {
            int depth = 0;
            int members = 0;
            measure_blossom(*region, depth, members);
            result.max_blossom_nesting_depth = std::max(result.max_blossom_nesting_depth, depth);
            result.max_blossom_members = std::max(result.max_blossom_members, members);
            result.matched_blossom_shatters++;
        }

        if (region->match.region != nullptr) {
            frozen_half_pairs++;
            if (is_extraction_representative(*region))
                scratch.top_regions.push_back(region);
        } else {
            // Not in a tree and not matched to a region: matched to the boundary. An unmatched
            // region outside every tree would mean the primal state is inconsistent.
            assert(region->match.edge.loc_from != nullptr && "top-level region is neither matched nor in a tree");
            result.committed_boundary++;
            scratch.top_regions.push_back(region);
        }
    }
    assert(frozen_half_pairs % 2 == 0 && "a frozen match should contribute both of its endpoints");
    result.committed_pairs_frozen = frozen_half_pairs / 2;

#ifndef NDEBUG
    assert_live_top_regions_match_detection_event_sweep(mwpm, detection_events);
#else
    (void)detection_events;
#endif

    // `largest_tree_size` is the one quantity direct enumeration does not get for free: it needs a
    // root id per node, which is the parent chase §M2.9.1 exists to remove. Computed here, before
    // any node is unlinked, and only when asked for.
    if (collect_diagnostics) {
        for (pm::AltTreeNode* node : mwpm.node_arena.live) {
            auto* root = const_cast<pm::AltTreeNode*>(node->find_root());
            int size = ++scratch.tree_sizes[root];
            result.largest_tree_size = std::max(result.largest_tree_size, size);
        }
    }

    // ---- Phase B: enumerate the live alternating tree nodes and classify each one locally.
    //
    // Committing a pair deletes its node, and the arena removes a deleted object from `live` by
    // swapping the last element into its slot — so the sweep takes a copy first rather than
    // depending on that. The copy is one memcpy of a vector that is a few tens of pointers long.
    scratch.live_nodes.assign(mwpm.node_arena.live.begin(), mwpm.node_arena.live.end());
    for (pm::AltTreeNode* node : scratch.live_nodes) {
        counters.tree_nodes_visited++;
        if (node->inner_region == nullptr) {
            // Root. Its outer region is the residual; deferred so that every pair is committed and
            // frozen before any extraction starts touching detector nodes.
            scratch.tree_roots.push_back(node);
        } else {
            // Committed, tight: the inner/outer pair of a non-root node is already a matched pair
            // of the current primal, joined by a tight edge. This is exactly what M1.3's descent
            // from the root committed, decided here from the node's own three fields.
            const pm::CompressedEdge& edge = node->inner_to_outer_edge;
            result.committed_pairs_tree++;
            if (tight_pairs_out != nullptr) {
                tight_pairs_out->push_back(
                    TightPairRecord{
                        (uint64_t)(edge.loc_from - node_base),
                        (uint64_t)(edge.loc_to - node_base),
                        nested_dual_sum(*edge.loc_from, time),
                        nested_dual_sum(*edge.loc_to, time)});
            }
            pm::GraphFillRegion* inner = node->inner_region;
            pm::GraphFillRegion* outer = node->outer_region;
            inner->add_match(outer, edge);
            flooder.set_region_frozen(*inner);
            flooder.set_region_frozen(*outer);
            inner->alt_tree_node = nullptr;
            outer->alt_tree_node = nullptr;
            scratch.top_regions.push_back(is_extraction_representative(*inner) ? inner : outer);
            if (collect_diagnostics) {
                for (pm::GraphFillRegion* side : {inner, outer}) {
                    if (side->blossom_children.empty())
                        continue;
                    int depth = 0;
                    int members = 0;
                    measure_blossom(*side, depth, members);
                    result.max_blossom_nesting_depth = std::max(result.max_blossom_nesting_depth, depth);
                    result.max_blossom_members = std::max(result.max_blossom_members, members);
                    result.matched_blossom_shatters++;
                }
            }
            mwpm.node_arena.del(node);
        }
    }

    if (collect_diagnostics) {
        result.enumerate_ns = stage.elapsed_ns();
        stage.start();
    }

    // ---- Phase B2: one exposed defect per surviving tree. The base descent of §M2.9.4 is the one
    // traversal that stays: an arbitrary member of a root blossom would break `Y(u) == T`.
    int deepest_base_descent = 0;
    for (pm::AltTreeNode* root : scratch.tree_roots) {
        result.num_trees++;
        pm::GraphFillRegion* exposed = root->outer_region;
        if (!exposed->blossom_children.empty()) {
            result.exposed_root_blossoms++;
            if (collect_diagnostics) {
                int depth = 0;
                int members = 0;
                measure_blossom(*exposed, depth, members);
                result.max_blossom_nesting_depth = std::max(result.max_blossom_nesting_depth, depth);
                result.max_blossom_members = std::max(result.max_blossom_members, members);
                result.max_exposed_blossom_depth = std::max(result.max_exposed_blossom_depth, depth);
                result.max_exposed_blossom_members = std::max(result.max_exposed_blossom_members, members);
                deepest_base_descent = std::max(deepest_base_descent, depth);
            }
        }
        exposed->alt_tree_node = nullptr;
        mwpm.node_arena.del(root);

        counters.base_descents++;
        pm::DetectorNode* base_node = find_exposed_base_node(*exposed, time);
        assert(base_node != nullptr);
        pm::total_weight_int base_dual_sum = nested_dual_sum(*base_node, time);

        uint64_t exposed_defect = 0;
        result.committed_pairs_blossom_cycle += (int)extract_exposed(*exposed, exposed_defect);
        assert(exposed_defect == (uint64_t)(base_node - node_base));

        result.residual.push_back(exposed_defect);
        result.residual_dual_sum.push_back(base_dual_sum);
    }

    if (collect_diagnostics) {
        result.base_descent_ns = stage.elapsed_ns();
        stage.start();
    }

    // ---- Phase C: extract everything that is matched, once per pair, from its representative.
    for (pm::GraphFillRegion* region : scratch.top_regions) {
        counters.extractions++;
        extract_matched(region);
    }

    if (collect_diagnostics) {
        result.shatter_ns = stage.elapsed_ns();
        stage.start();
    }

    // The residual is emitted sorted ascending, with its dual sums kept aligned. In the hardware
    // model this is a prefix sum over the root-asserting nodes; in software it stays a sort, through
    // a scratch buffer so the only allocations per shot are the result's own. The contract that the
    // vector is sorted in `G`'s detector ids is consumed by M3–M6 and does not change (§M2.9.3).
    scratch.residual_sort_buffer.clear();
    scratch.residual_sort_buffer.reserve(result.residual.size());
    for (size_t i = 0; i < result.residual.size(); i++)
        scratch.residual_sort_buffer.emplace_back(result.residual[i], result.residual_dual_sum[i]);
    std::sort(scratch.residual_sort_buffer.begin(), scratch.residual_sort_buffer.end());
    for (size_t i = 0; i < scratch.residual_sort_buffer.size(); i++) {
        result.residual[i] = scratch.residual_sort_buffer[i].first;
        result.residual_dual_sum[i] = scratch.residual_sort_buffer[i].second;
    }
    assert(result.residual.size() == (size_t)result.num_trees);

    if (collect_diagnostics) {
        result.reduce_ns = stage.elapsed_ns();
        // The modelled serial depth of harvest, in dependent operations (§M2.9.6). Enumeration is
        // one step: every region and every node classifies itself from fields it already holds. The
        // committed observable and weight are then associative reductions over the commits, so a
        // balanced tree of depth `ceil(log2(commits))`. Running concurrently with that reduction,
        // and not after it, are the two chases that stay: the base descent through the nesting of
        // an exposed root blossom, and the shatter of a matched blossom, each one level per level
        // of nesting. The residual compaction is a prefix sum over the roots.
        size_t commits = (size_t)result.committed_pairs_frozen + (size_t)result.committed_pairs_tree +
                         (size_t)result.committed_pairs_blossom_cycle + (size_t)result.committed_boundary;
        int reduce_depth = reduction_depth(commits);
        int compact_depth = reduction_depth(result.residual.size());
        int chase_depth = std::max(deepest_base_descent, result.max_blossom_nesting_depth);
        result.harvest_dependent_depth = 1 + std::max(std::max(reduce_depth, compact_depth), chase_depth);
    }

    reset_for_next_shot(mwpm);
    return result;
}

/// §M3.4's production path: phase A and phase C of the harvest above, with phases B and B2 gone.
///
/// The two dropped phases are the tree commits and M1.4's base descent, and both do work only when
/// an alternating tree survives — which is exactly the condition that escalates the shot and throws
/// their output away (§M3.1). On a `COMPLETE` timeline they are not skipped optimistically; there
/// is provably nothing for them to do, and the `assert` below says so.
///
/// What this does *not* drop is the dual sum. It is a flat reduction over the same live-region walk
/// extraction already performs — one integer read and one add per region — and the §M4.2 per-shot
/// certificate is the cheap independent check on the primal that the design keeps precisely because
/// the escalation path has removed every other one.
template <typename ExtractMatched>
HarvestResult Harvester::extract_only_impl(pm::Mwpm& mwpm, const ExtractMatched& extract_matched) {
    counters.extract_only_harvests++;
    HarvestResult result;
    scratch.clear();

    auto& flooder = mwpm.flooder;
    const pm::cumulative_time_int time = flooder.queue.cur_time;

    Stopwatch stage;
    if (collect_diagnostics)
        stage.start();

    int frozen_half_pairs = 0;
    for (pm::GraphFillRegion* region : flooder.region_arena.live) {
        counters.regions_visited++;
        pm::total_weight_int dual = region->radius.get_distance_at_time(time);
        result.dual_sum_at_truncation += dual;
        result.max_region_dual = std::max(result.max_region_dual, dual);

        if (region->blossom_parent != nullptr)
            continue;  // Nested inside a blossom; its top-level ancestor speaks for it.

        // The precondition, restated where it would be violated. A top-level region in a tree on a
        // path that has decided no tree survives means the O(1) status branch and the primal state
        // disagree, which is the failure §M3.3 X9 exists to catch.
        assert(region->alt_tree_node == nullptr && "extract-only path reached a region still in an alternating tree");

        if (collect_diagnostics && !region->blossom_children.empty()) {
            int depth = 0;
            int members = 0;
            measure_blossom(*region, depth, members);
            result.max_blossom_nesting_depth = std::max(result.max_blossom_nesting_depth, depth);
            result.max_blossom_members = std::max(result.max_blossom_members, members);
            result.matched_blossom_shatters++;
        }

        if (region->match.region != nullptr) {
            frozen_half_pairs++;
            if (is_extraction_representative(*region))
                scratch.top_regions.push_back(region);
        } else {
            assert(region->match.edge.loc_from != nullptr && "top-level region is neither matched nor in a tree");
            result.committed_boundary++;
            scratch.top_regions.push_back(region);
        }
    }
    assert(frozen_half_pairs % 2 == 0 && "a frozen match should contribute both of its endpoints");
    result.committed_pairs_frozen = frozen_half_pairs / 2;

    if (collect_diagnostics) {
        result.enumerate_ns = stage.elapsed_ns();
        stage.start();
    }

    for (pm::GraphFillRegion* region : scratch.top_regions) {
        counters.extractions++;
        extract_matched(region);
    }

    if (collect_diagnostics) {
        result.shatter_ns = stage.elapsed_ns();
        // Enumeration is one step; the committed observable and weight are an associative reduction
        // over the commits; the only chase left is the shatter of a matched blossom, one level per
        // level of nesting. There is no base descent and no residual compaction on this path, which
        // is the whole of what §M3.4 removes from the depth.
        size_t commits = (size_t)result.committed_pairs_frozen + (size_t)result.committed_boundary;
        result.harvest_dependent_depth = 1 + std::max(reduction_depth(commits), result.max_blossom_nesting_depth);
    }

    assert(result.residual.empty() && result.num_trees == 0);
    reset_for_next_shot(mwpm);
    return result;
}

/// M1.3's enumeration, verbatim, kept as the H1 oracle. Do not "fix" this to match the new path:
/// its whole job is to be the thing the new path is compared against.
template <typename ExtractMatched, typename ExtractExposed>
HarvestResult Harvester::harvest_impl_legacy(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    const ExtractMatched& extract_matched,
    const ExtractExposed& extract_exposed) {
    HarvestResult result;
    scratch.clear();

    auto& flooder = mwpm.flooder;
    auto& nodes = flooder.graph.nodes;
    const pm::DetectorNode* node_base = nodes.data();
    // `process_timeline_until_horizon` has already advanced the clock to exactly the horizon, so
    // every dual read below is the dual at `T`.
    const pm::cumulative_time_int time = flooder.queue.cur_time;

    // The negative-weight detection events are extracted alongside the shot's own, exactly as the
    // stock decode path does: they are real detection events whose regions have to be resolved.
    const std::vector<uint64_t>* event_lists[2] = {&detection_events, &flooder.negative_weight_detection_events};

    // ---- Phase A: classify. Reads only; nothing is mutated until every region is accounted for.
    for (const std::vector<uint64_t>* events : event_lists) {
        for (uint64_t det : *events) {
            if (det >= nodes.size())
                continue;
            pm::DetectorNode& node = nodes[det];
            if (node.region_that_arrived == nullptr)
                continue;
            pm::GraphFillRegion* top = node.region_that_arrived_top;
            if (!scratch.seen_regions.insert(top).second)
                continue;
            scratch.top_regions.push_back(top);
        }
    }

    int frozen_half_pairs = 0;
    for (pm::GraphFillRegion* top : scratch.top_regions) {
        // Every live region is `top` or nested inside it, so this covers the whole dual solution.
        top->do_op_for_each_descendant_and_self([&](pm::GraphFillRegion* region) {
            pm::total_weight_int dual = region->radius.get_distance_at_time(time);
            result.dual_sum_at_truncation += dual;
            result.max_region_dual = std::max(result.max_region_dual, dual);
        });

        if (top->alt_tree_node != nullptr) {
            auto* root = const_cast<pm::AltTreeNode*>(top->alt_tree_node->find_root());
            if (scratch.seen_roots.insert(root).second)
                scratch.tree_roots.push_back(root);
        } else if (top->match.region != nullptr) {
            frozen_half_pairs++;
        } else {
            // Not in a tree and not matched to a region: matched to the boundary. An unmatched
            // region outside every tree would mean the primal state is inconsistent.
            assert(top->match.edge.loc_from != nullptr && "top-level region is neither matched nor in a tree");
            result.committed_boundary++;
        }
    }
    assert(frozen_half_pairs % 2 == 0 && "a frozen match should contribute both of its endpoints");
    result.committed_pairs_frozen = frozen_half_pairs / 2;

    // ---- Phase B: commit the surviving trees and expose one defect per tree.
    for (pm::AltTreeNode* root : scratch.tree_roots) {
        result.num_trees++;

        int tree_size = 0;
        scratch.tree_walk_stack.clear();
        scratch.tree_walk_stack.push_back(root);
        while (!scratch.tree_walk_stack.empty()) {
            pm::AltTreeNode* tree_node = scratch.tree_walk_stack.back();
            scratch.tree_walk_stack.pop_back();
            tree_size++;
            if (tree_node->inner_region != nullptr) {
                // Committed, tight: the inner/outer pair of a non-root node is already a matched
                // pair of the current primal, joined by a tight edge.
                result.committed_pairs_tree++;
                if (tight_pairs_out != nullptr) {
                    const pm::CompressedEdge& edge = tree_node->inner_to_outer_edge;
                    tight_pairs_out->push_back(
                        TightPairRecord{
                            (uint64_t)(edge.loc_from - node_base),
                            (uint64_t)(edge.loc_to - node_base),
                            nested_dual_sum(*edge.loc_from, time),
                            nested_dual_sum(*edge.loc_to, time)});
                }
            }
            for (const pm::AltTreeEdge& child : tree_node->children)
                scratch.tree_walk_stack.push_back(child.alt_tree_node);
        }
        result.largest_tree_size = std::max(result.largest_tree_size, tree_size);

        pm::GraphFillRegion* exposed = root->outer_region;
        if (!exposed->blossom_children.empty())
            result.exposed_root_blossoms++;

        // Matches and freezes every non-root inner/outer pair, and dismantles the tree. The root
        // is left alone beyond clearing its back-pointer: its outer region is the residual.
        mwpm.shatter_descendants_into_matches_and_freeze(*root);

        pm::DetectorNode* base_node = find_exposed_base_node(*exposed, time);
        assert(base_node != nullptr);
        pm::total_weight_int base_dual_sum = nested_dual_sum(*base_node, time);

        uint64_t exposed_defect = 0;
        result.committed_pairs_blossom_cycle += (int)extract_exposed(*exposed, exposed_defect);
        assert(exposed_defect == (uint64_t)(base_node - node_base));

        result.residual.push_back(exposed_defect);
        result.residual_dual_sum.push_back(base_dual_sum);
    }

    // ---- Phase C: extract everything that is matched. Shattering resets the detector nodes it
    // consumes, so re-reading `region_that_arrived` each time is what stops a pair being extracted
    // twice from its two endpoints. This mirrors the stock extraction helpers.
    for (const std::vector<uint64_t>* events : event_lists) {
        for (uint64_t det : *events) {
            if (det >= nodes.size())
                continue;
            pm::DetectorNode& node = nodes[det];
            if (node.region_that_arrived == nullptr)
                continue;
            extract_matched(node.region_that_arrived_top);
        }
    }

    // The residual is emitted sorted ascending, with its dual sums kept aligned. Sorted through a
    // scratch buffer so the only allocations per shot are the result's own.
    scratch.residual_sort_buffer.clear();
    scratch.residual_sort_buffer.reserve(result.residual.size());
    for (size_t i = 0; i < result.residual.size(); i++)
        scratch.residual_sort_buffer.emplace_back(result.residual[i], result.residual_dual_sum[i]);
    std::sort(scratch.residual_sort_buffer.begin(), scratch.residual_sort_buffer.end());
    for (size_t i = 0; i < scratch.residual_sort_buffer.size(); i++) {
        result.residual[i] = scratch.residual_sort_buffer[i].first;
        result.residual_dual_sum[i] = scratch.residual_sort_buffer[i].second;
    }
    assert(result.residual.size() == (size_t)result.num_trees);

    reset_for_next_shot(mwpm);
    return result;
}

bool any_alternating_tree_survives_by_sweep(const pm::Mwpm& mwpm) {
    for (const pm::GraphFillRegion* region : mwpm.flooder.region_arena.live) {
        if (region->alt_tree_node != nullptr)
            return true;
    }
    return false;
}

void reset_for_next_shot(pm::Mwpm& mwpm) {
    // Every region and alternating tree node has already been handed back to its arena by the
    // extraction above, so — unlike `Mwpm::reset` — the arenas are left alone. Destroying them
    // would free the pools and make the next shot malloc every region again, which costs more than
    // the whole truncated timeline.
    assert(
        mwpm.node_arena.allocated.size() == mwpm.node_arena.available.size() &&
        "harvest left an alternating tree node checked out");
    assert(
        mwpm.flooder.region_arena.allocated.size() == mwpm.flooder.region_arena.available.size() &&
        "harvest left a graph fill region checked out");
    // The same statement read off the intrusive live lists rather than the free lists. They are
    // maintained independently, so agreeing here is a real check on §M2.9.1's enumeration: a node
    // or region that was deleted without being unlinked, or vice versa, shows up as a mismatch on
    // the very next shot rather than as a silently wrong harvest much later.
    assert(mwpm.node_arena.live.empty());
    assert(mwpm.flooder.region_arena.live.empty());

    // Shattering already reset every node the shot owned, via `cleanup_shell_area`, and a node a
    // shrinking region released was zeroed on the way out. What a truncated timeline leaves behind
    // that a completed one does not is the queue: it drains itself when the matching completes,
    // but here it can still hold exempt shrink events. A node with an event still in the queue is
    // exactly a node whose tracker still believes it has one queued, and that belief would stop
    // the next shot from ever enqueueing for that node.
    //
    // So clear those trackers and drop the queue — work proportional to what is left over, rather
    // than to the size of the graph. `Mwpm::reset`'s whole-graph sweep would cost more per shot
    // than the timeline it is cleaning up after.
    for (auto& bucket : mwpm.flooder.queue.bit_buckets) {
        for (auto& event : bucket) {
            if (event.tentative_event_type == pm::LOOK_AT_NODE)
                event.data_look_at_node->node_event_tracker.clear();
        }
    }
    mwpm.flooder.queue.clear();

    // Shrink-event trackers live in graph fill regions, which are default-constructed afresh when
    // the arena hands them out again, so they need no cleaning here.

#ifndef NDEBUG
    // The reasoning above is subtle enough to be worth checking outright in debug builds: the
    // instance must be indistinguishable from a fresh one.
    for (const auto& node : mwpm.flooder.graph.nodes) {
        assert(node.region_that_arrived == nullptr);
        assert(node.region_that_arrived_top == nullptr);
        assert(node.reached_from_source == nullptr);
        assert(node.radius_of_arrival == 0);
        assert(node.wrapped_radius_cached == 0);
        assert(node.observables_crossed_from_source == 0);
        assert(!node.node_event_tracker.has_queued_time);
    }
#endif
}

HarvestResult Harvester::harvest_to_obs(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    HarvestResult result;
    pm::MatchingResult committed;
    result = harvest_impl(
        mwpm,
        detection_events,
        [&](pm::GraphFillRegion* region) {
            committed += mwpm.shatter_blossom_and_extract_matches(region);
        },
        [&](pm::GraphFillRegion& exposed, uint64_t& exposed_defect_out) {
            return shatter_exposed_blossom_and_extract_matches(
                mwpm, exposed, mwpm.flooder.queue.cur_time, exposed_defect_out, committed);
        });
    result.committed = committed;
    return result;
}

HarvestResult Harvester::harvest_to_match_edges(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<pm::CompressedEdge>& match_edges) {
    return harvest_impl(
        mwpm,
        detection_events,
        [&](pm::GraphFillRegion* region) {
            mwpm.shatter_blossom_and_extract_match_edges(region, match_edges);
        },
        [&](pm::GraphFillRegion& exposed, uint64_t& exposed_defect_out) {
            return shatter_exposed_blossom_and_extract_match_edges(
                mwpm, exposed, mwpm.flooder.queue.cur_time, exposed_defect_out, match_edges);
        });
}

HarvestResult Harvester::extract_only_to_obs(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    (void)detection_events;
    // Debug invariant 13a's second half, checked where the decision is acted on rather than only
    // where it is taken: the O(1) branch and the full sweep must agree, on every shot (§M3.3 X9).
    assert(!any_alternating_tree_survives(mwpm) && !any_alternating_tree_survives_by_sweep(mwpm));
    pm::MatchingResult committed;
    HarvestResult result = extract_only_impl(mwpm, [&](pm::GraphFillRegion* region) {
        committed += mwpm.shatter_blossom_and_extract_matches(region);
    });
    result.committed = committed;
    return result;
}

HarvestResult Harvester::extract_only_to_match_edges(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<pm::CompressedEdge>& match_edges) {
    (void)detection_events;
    assert(!any_alternating_tree_survives(mwpm) && !any_alternating_tree_survives_by_sweep(mwpm));
    return extract_only_impl(mwpm, [&](pm::GraphFillRegion* region) {
        mwpm.shatter_blossom_and_extract_match_edges(region, match_edges);
    });
}

HarvestResult harvest_to_obs(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    Harvester harvester;
    return harvester.harvest_to_obs(mwpm, detection_events);
}

HarvestResult harvest_to_match_edges(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<pm::CompressedEdge>& match_edges) {
    Harvester harvester;
    return harvester.harvest_to_match_edges(mwpm, detection_events, match_edges);
}

}  // namespace spec_matching
}  // namespace pm
