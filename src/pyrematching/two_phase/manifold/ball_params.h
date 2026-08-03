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

#ifndef PYREMATCHING_TWO_PHASE_MANIFOLD_BALL_PARAMS_H
#define PYREMATCHING_TWO_PHASE_MANIFOLD_BALL_PARAMS_H

#include <cstdint>
#include <stdexcept>
#include <string>

namespace pm {
namespace two_phase {

/// Compile-time parameters of the ball tables (§M2.2).
///
/// All lengths are in DEM float weight units; they are converted to the flooder's integer time
/// units exactly once, at table construction, with `to_time_units` and the graph's own
/// `normalising_constant` (the unit rule of §0).
struct BallParams {
    /// The largest horizon the compiled tables support. Decoding at `T > T_max` is a hard error.
    double T_max{2.0};
    /// Ball radius. **Hard requirement `R >= 2 * T_max`**: two defects first interact when
    /// `Y(u) + Y(v) == d_G(u, v)`, and `Y <= T` throughout a truncated run, so a pair beyond `2T`
    /// is unreachable — but a pair *within* `2T` that the ball omitted would silently change the
    /// answer.
    double R{4.0};
    /// Shell granularity for the optional restart ladder (§M2.8). `0` means a single shell covering
    /// the whole ball, which is the default and must always work on its own. Shell `s` covers
    /// `[s * shell_width, (s + 1) * shell_width)`.
    double shell_width{0.0};
    /// Also store the outer node path behind every ball entry, so a correction can be lifted back
    /// to `G`'s edges (M5). Costs a lot of memory; off by default.
    bool need_edge_lift{false};
    /// Recompute every canonical path with a perturbed tie-break and count the entries whose
    /// observable set differs. Diagnostic only — a non-zero count is expected on a surface code and
    /// never blocks compilation.
    bool certify_masks{false};
    uint64_t seed{0};

    bool operator==(const BallParams& rhs) const {
        return T_max == rhs.T_max && R == rhs.R && shell_width == rhs.shell_width &&
               need_edge_lift == rhs.need_edge_lift && certify_masks == rhs.certify_masks && seed == rhs.seed;
    }
    bool operator!=(const BallParams& rhs) const {
        return !(*this == rhs);
    }

    /// Debug-build invariant 6, checked at ball-table construction and again at decoder
    /// construction. An undersized `R` changes the answer without any crash, so it is checked in
    /// release builds too.
    void validate() const {
        if (!(T_max > 0))
            throw std::invalid_argument("BallParams::T_max must be positive, got " + std::to_string(T_max));
        if (!(R >= 2 * T_max))
            throw std::invalid_argument(
                "BallParams requires R >= 2 * T_max (the exactness theorem of §M2.0), got R = " + std::to_string(R) +
                ", T_max = " + std::to_string(T_max));
        if (shell_width < 0)
            throw std::invalid_argument("BallParams::shell_width must be non-negative");
    }
};

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_MANIFOLD_BALL_PARAMS_H
