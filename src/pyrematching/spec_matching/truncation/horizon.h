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

#ifndef PYREMATCHING_SPEC_MATCHING_TRUNCATION_HORIZON_H
#define PYREMATCHING_SPEC_MATCHING_TRUNCATION_HORIZON_H

#include <cmath>
#include <limits>

#include "pyrematching/sparse_blossom/ints.h"
#include "pyrematching/sparse_blossom/tracker/queued_event_tracker.h"

namespace pm {
namespace spec_matching {

/// The truncation horizon `T`, in the flooder's cumulative time units.
///
/// `pm::NO_HORIZON` is the sentinel meaning "no truncation"; it is the value a `GraphFlooder`
/// carries outside `process_timeline_until_horizon`.
using horizon_int = pm::cumulative_time_int;

/// Converts a horizon (or any distance) expressed in DEM float weight units into the flooder's
/// integer time units.
///
/// `normalising_constant` must be `MatchingGraph::normalising_constant`, i.e. the value returned by
/// `UserGraph::to_matching_graph`, which is already twice
/// `UserGraph::get_edge_weight_normalising_constant(num_distinct_weights)`.
///
/// The `/ 2 ... * 2` dance is not redundant: `UserGraph::iter_discretized_edges` discretises an
/// edge as `2 * round(w * C / 2)` so that every edge weight is even and every collision happens at
/// an integer time. Converting a horizon by the same formula puts `T` on exactly the same lattice
/// as the edge weights, which is what makes `Y(u) == T` an exact equality rather than an
/// approximate one.
///
/// This is the *only* place weights are converted to time units. Store the result; never convert
/// in two places.
inline horizon_int to_time_units(double weight, double normalising_constant) {
    if (std::isinf(weight))
        return weight > 0 ? pm::NO_HORIZON : -pm::NO_HORIZON;
    return 2 * (horizon_int)std::llround(weight * normalising_constant / 2);
}

/// Inverse of `to_time_units`, for reporting weights back in DEM float weight units.
inline double to_weight_units(pm::total_weight_int time_units, double normalising_constant) {
    return (double)time_units / normalising_constant;
}

}  // namespace spec_matching
}  // namespace pm

#endif  // PYREMATCHING_SPEC_MATCHING_TRUNCATION_HORIZON_H
