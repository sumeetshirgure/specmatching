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

#ifndef PYREMATCHING_TWO_PHASE_ESCALATION_ESCALATE_H
#define PYREMATCHING_TWO_PHASE_ESCALATION_ESCALATE_H

#include <cstdint>
#include <vector>

#include "pyrematching/sparse_blossom/matcher/mwpm.h"
#include "pyrematching/two_phase/perf/two_phase_profile.h"

namespace pm {
namespace two_phase {

/// M3 — stock escalation, the residual fallback.
///
/// When Phase 1 leaves a non-empty residual the partial result is discarded and the **whole shot**
/// is re-decoded with stock exact decode on `G`. The decoder's output is then exact MWPM on every
/// shot, by construction: there is no approximation budget anywhere in the pipeline.
///
/// Three things about this are load-bearing and are asserted by §M3.3's tests rather than left to
/// a reader's good will:
///
///  - **The whole syndrome is re-decoded, not the residual.** Phase 1's committed pairs are optimal
///    *within* the truncated dual, and stock is free to pair a committed defect differently.
///    Feeding stock anything less than the full detection-event set is a correctness bug, not an
///    optimisation (invariant 14, test X3).
///  - **Phase 1's output is discarded, not combined.** There is nothing to XOR and nothing to add.
///  - **Negative-weight DEMs go in raw.** Stock's own preamble applies unchanged, because stock is
///    being called on the raw syndrome. Do *not* pass §M2.1's post-preamble effective set; that set
///    is an input to `H`, never to `G`.
///
/// `g_mwpm` is the caller's long-lived `pm::Mwpm` on `G` — one per decoder, reset between shots
/// exactly as the stock driver does. Do not construct one per escalation.

/// Obs flavour. `obs` must point at `g_mwpm.flooder.graph.num_observables` bytes; it and `weight`
/// are **overwritten** with stock's answer, not accumulated into.
///
/// Returns true, always: the caller has already decided the shot escalates. The return value exists
/// so the call site reads as the branch it is, and so `prof->escalated` and the return cannot drift
/// apart.
bool escalate_to_stock(
    pm::Mwpm& g_mwpm,
    const std::vector<uint64_t>& dets,
    uint8_t* obs,
    size_t num_observables,
    pm::total_weight_int& weight,
    TwoPhaseProfile* prof = nullptr);

/// Edges flavour (M5). `edges` is **cleared** and refilled with stock's answer: pairs of detector
/// ids, `-1` for the boundary, in the inherited `decode_detection_events_to_edges` format.
///
/// Requires `g_mwpm` to carry a search flooder, which is what turns matched detection events into
/// the edges of a correction. Throws otherwise, from the stock call.
bool escalate_to_stock_edges(
    pm::Mwpm& g_mwpm, const std::vector<uint64_t>& dets, std::vector<int64_t>& edges, TwoPhaseProfile* prof = nullptr);

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_ESCALATION_ESCALATE_H
