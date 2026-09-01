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

#ifndef SPECMATCHING_SPEC_MATCHING_CERTIFICATE_MAX_DUAL_H
#define SPECMATCHING_SPEC_MATCHING_CERTIFICATE_MAX_DUAL_H

#include <cstdint>
#include <vector>

#include "specmatching/sparse_blossom/matcher/mwpm.h"
#include "specmatching/spec_matching/truncation/horizon.h"

namespace pm {
namespace spec_matching {

/// M7 — the terminal max-dual certificate.
///
/// **The theorem (§M7.0).** Let `H` be the ball graph: an edge `(u, v)` of weight `d_G(u, v)` for
/// every defect pair with `d_G(u, v) <= 2T`, and a boundary edge of weight `bcost(u)` for every `u`
/// with `bcost(u) <= T`. Run *stock* sparse blossom on `H`, no horizon, to completion. If it
/// completes and its **terminal** dual satisfies `max_u Y(u) <= T`, then extending that dual by
/// `y_S = 0` off `H` is feasible for the full complete-graph matching LP — every omitted edge has
/// `d_G(u, v) > 2T >= Y(u) + Y(v)`, so strictly positive slack — and complementary slackness holds
/// on every edge the matching uses. Primal-feasible plus dual-feasible plus complementary slackness
/// is optimal, so `H`'s matching is a global MWPM: the same pairs and the same weight as stock
/// exact decode on `G`.
///
/// **Three things this file exists to get right.**
///
///  - The certificate is about the **terminal** dual, not the trajectory. A run whose duals spike
///    above `T` mid-flight and relax back below it by completion is exact, so the test is applied
///    once, after the loop reaches `NO_EVENT`. Early termination on a running max is *rejected*
///    (§M7.2): it escalates such shots spuriously and reintroduces the per-event work the whole
///    scheme exists to remove. Do not re-add it.
///  - `H` having no perfect matching is a **valid escalation trigger**, not a failure: if the true
///    optimum had all duals `<= T`, every one of its matched pairs would be a tight edge of weight
///    `<= 2T` and would therefore be present in `H`.
///  - Catching that case is **not sufficient on its own.** `H` can admit a *suboptimal* perfect
///    matching whose terminal `max_u Y(u) > T`; such a run completes and does not throw, and
///    returning its matching would be a silent wrong answer. The dual is therefore inspected on
///    **every** completing run — that is what makes this a certificate rather than an exception
///    handler, and it is the single most important correctness point in M7.
enum class CertificateStatus {
    /// Completed on `H` with `max_u Y(u) <= T`. Globally optimal; extract and return.
    CERTIFIED,
    /// `H` admits no perfect matching, so some dual of the true optimum exceeds `T`. Escalate.
    NO_PERFECT_MATCHING,
    /// Completed on `H`, but the terminal dual escapes the horizon, so `H` may be missing an edge
    /// the true optimum needs. Escalate. This is the case that does *not* announce itself.
    DUAL_EXCEEDS_HORIZON,
};

struct CertificateOutcome {
    CertificateStatus status{CertificateStatus::CERTIFIED};
    /// `max_u Y(u)` at completion — the quantity the certificate tests. Zero when `H` could not
    /// complete, where it is not defined rather than merely unmeasured.
    pm::total_weight_int max_dual{0};
    /// The one terminal scan, timed. Paid only on completing shots, and the only cost this scheme
    /// adds over an unmodified stock-on-`H` solve (§M7.7). Zero unless timing was requested.
    long long dual_scan_ns{0};

    inline bool certified() const {
        return status == CertificateStatus::CERTIFIED;
    }
    /// Every non-certified outcome escalates, and for the two reasons the design distinguishes.
    inline bool escalates() const {
        return status != CertificateStatus::CERTIFIED;
    }
};

/// `max_u Y(u)` over `detection_events` (and the graph's negative-weight detection events, which on
/// `H` is empty by construction), evaluated at the flooder's current time.
///
/// `Y(u)` is `Sum_{S ∋ u} y_S`, read through `nested_dual_sum`: the top region's own dual plus the
/// node's wrapped radius, which is the sum of the frozen sub-region radii below it. Read-only — it
/// must be called *before* extraction starts shattering the nesting it derives the answer from.
pm::total_weight_int max_nested_dual(const pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

/// The same quantity, computed by walking each defect's `blossom_parent` chain and summing `y_S`
/// explicitly — the walk M1 test 2 uses.
///
/// It exists to be a genuinely independent second opinion on `max_nested_dual`, which goes through
/// `region_that_arrived_top->radius + compute_wrapped_radius()` instead. Debug invariant 2 asserts
/// the two agree on every certified shot; a disagreement means the nesting bookkeeping and the dual
/// readout have drifted apart, which nothing else on this path would notice.
pm::total_weight_int max_nested_dual_by_chain_walk(const pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

/// Runs stock sparse blossom on `mwpm` to completion and applies the certificate.
///
/// `horizon` is `T` in the flooder's discretised time units — the same value that filtered `H` at
/// `2 * T_int`, converted once via the shared normalising constant (§0's unit rule). Duals and
/// cumulative time share that metric (M2.0), so the comparison is unit-clean by construction.
///
/// **No horizon is ever set on the flooder.** `GraphFlooder::horizon` holds `pm::NO_HORIZON`
/// throughout and is asserted to, on entry and on exit: that assert is the machine-checkable form
/// of §M7.1's "no hardcoded horizon anywhere", and it is debug invariant 1.
///
/// On `CERTIFIED` the instance is left in its completed state for the caller to extract from. On
/// either escalating outcome it is left standing too — the caller tears it down with
/// `abandon_shot` (`Mwpm::reset`), because it is the caller that knows the shot is being discarded.
CertificateOutcome run_stock_and_certify(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, horizon_int horizon, bool time_dual_scan = false);

}  // namespace spec_matching
}  // namespace pm

#endif  // SPECMATCHING_SPEC_MATCHING_CERTIFICATE_MAX_DUAL_H
