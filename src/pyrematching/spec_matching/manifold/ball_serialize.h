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

#ifndef PYREMATCHING_SPEC_MATCHING_MANIFOLD_BALL_SERIALIZE_H
#define PYREMATCHING_SPEC_MATCHING_MANIFOLD_BALL_SERIALIZE_H

#include <string>

#include "pyrematching/spec_matching/manifold/ball_tables.h"

namespace pm {
namespace spec_matching {

/// Bumped whenever the on-disk layout changes. Loading an artifact with a different version is a
/// hard error, never a best-effort migration: a silently misread table changes the answer.
constexpr uint32_t BALL_ARTIFACT_VERSION = 2;

/// Writes the tables to `path`. Versioned header, graph hash and `BallParams` first, then the CSR
/// pools verbatim.
void save_ball_tables(const BallTables& tables, const std::string& path);

/// Reads tables written by `save_ball_tables` and validates them against `graph`.
///
/// A version, hash or parameter mismatch throws — the artifact was compiled against a different
/// DEM or different parameters, and using it anyway would silently decode a different problem.
BallTables load_ball_tables(const std::string& path, const pm::MatchingGraph& graph, const BallParams& expected);

}  // namespace spec_matching
}  // namespace pm

#endif  // PYREMATCHING_SPEC_MATCHING_MANIFOLD_BALL_SERIALIZE_H
