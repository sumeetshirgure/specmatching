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

#ifndef PYREMATCHING_TWO_PHASE_DRIVER_TWO_PHASE_PYBIND_H
#define PYREMATCHING_TWO_PHASE_DRIVER_TWO_PHASE_PYBIND_H

#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

namespace py = pybind11;

namespace pm_pybind {

/// §M6.3. Registers `TwoPhaseConfig`, `TwoPhaseDecoder`, `TwoPhaseAggregateStats` and
/// `summarize` on the extension module, following `user_graph.pybind.cc`'s conventions.
void pybind_two_phase(py::module& m);

}  // namespace pm_pybind

#endif  // PYREMATCHING_TWO_PHASE_DRIVER_TWO_PHASE_PYBIND_H
