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

#include "specmatching/spec_matching/escalation/escalate.h"

#include <algorithm>
#include <cassert>

#include "specmatching/sparse_blossom/driver/mwpm_decoding.h"

namespace pm {
namespace spec_matching {

bool escalate_to_stock(
    pm::Mwpm& g_mwpm,
    const std::vector<uint64_t>& dets,
    uint8_t* obs,
    size_t num_observables,
    pm::total_weight_int& weight,
    SpecMatchingProfile* prof) {
    HiResTimer escalation_timer;
    HiResTimer stock_timer;
    if (prof != nullptr)
        escalation_timer.start();

    // Overwrite, not accumulate. `decode_detection_events` XORs its observable bits into whatever
    // is already there (and, on the >64-observable path, adds to `weight` rather than assigning),
    // so zeroing here is what makes "Phase 1's partial result is discarded" true rather than
    // "Phase 1's partial result is XOR-ed with the right answer".
    std::fill(obs, obs + num_observables, (uint8_t)0);
    weight = 0;

    if (prof != nullptr)
        stock_timer.start();
    // The **whole** syndrome, in `G`'s raw detector ids. Not the residual, and not §M2.1's
    // post-preamble effective set: stock runs its own negative-weight preamble (invariant 14).
    pm::decode_detection_events(g_mwpm, dets, obs, weight, false);
    if (prof != nullptr) {
        prof->stock_ns = stock_timer.elapsed_ns();
        prof->escalated = true;
        prof->weight_out = weight;
        // The truncated dual went with the rest of Phase 1's partial result; there is nothing left
        // for the §M4.2 certificate to divide by, and reporting the discarded number would invite
        // exactly the "weight_out / dual_sum" comparison that no longer means anything here.
        prof->dual_sum_at_truncation = 0;
        prof->escalation_ns = escalation_timer.elapsed_ns();
    }
    return true;
}

bool escalate_to_stock_edges(
    pm::Mwpm& g_mwpm, const std::vector<uint64_t>& dets, std::vector<int64_t>& edges, SpecMatchingProfile* prof) {
    HiResTimer escalation_timer;
    HiResTimer stock_timer;
    if (prof != nullptr)
        escalation_timer.start();

    edges.clear();
    if (prof != nullptr)
        stock_timer.start();
    pm::decode_detection_events_to_edges(g_mwpm, dets, edges);
    if (prof != nullptr) {
        prof->stock_ns = stock_timer.elapsed_ns();
        prof->escalated = true;
        prof->dual_sum_at_truncation = 0;
        prof->escalation_ns = escalation_timer.elapsed_ns();
    }
    return true;
}

}  // namespace spec_matching
}  // namespace pm
