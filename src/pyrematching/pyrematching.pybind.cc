// Copyright 2022 PyMatching Contributors
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

#include "pybind11/pybind11.h"
#include "pyrematching/rand/rand_gen.pybind.h"
#include "pyrematching/sparse_blossom/driver/namespaced_main.h"
#include "pyrematching/sparse_blossom/driver/user_graph.pybind.h"
#include "pyrematching/two_phase/driver/two_phase.pybind.h"

namespace py = pybind11;
using namespace py::literals;

int pyrematching_main(const std::vector<std::string> &args) {
    std::vector<const char *> argv;
    argv.push_back("pyrematching.main");
    for (const auto &arg : args) {
        argv.push_back(arg.data());
    }
    return pm::main(argv.size(), argv.data());
}

PYBIND11_MODULE(_cpp_pyrematching, m) {
    auto matching_graph = pm_pybind::pybind_user_graph(m);
    pm_pybind::pybind_user_graph_methods(m, matching_graph);
    pm_pybind::pybind_rand_gen_methods(m);
    pm_pybind::pybind_two_phase(m);
    m.def("main", &pyrematching_main, pybind11::kw_only(), pybind11::arg("command_line_args"), R"pbdoc(
Runs the command line tool version of pyrematching with the given arguments.
)pbdoc");
}