cmake_minimum_required(VERSION 3.4...3.18)
project(PyReMatching2)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(src)
set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(SIMD_WIDTH 128)
if (NOT(MSVC))
    if (CMAKE_SYSTEM_PROCESSOR MATCHES x86_64)
         set(ARCH_OPT "-O3" "-mno-avx2")
    else ()
         set(ARCH_OPT "-O3")
    endif ()
else ()
    set(ARCH_OPT "-O2")
    # https://stackoverflow.com/a/8591946 (RTC1 incompatible with O2 flag)
    STRING (REGEX REPLACE "/RTC(su|[1su])" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif ()

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)

FetchContent_MakeAvailable(googletest)

FetchContent_Declare(stim
        GIT_REPOSITORY https://github.com/quantumlib/Stim.git 
        GIT_TAG v1.16.0
      )
FetchContent_MakeAvailable(stim)
if (NOT (MSVC))
    target_compile_options(libstim PRIVATE -fno-strict-aliasing -fPIC ${ARCH_OPT})
else ()
    target_compile_options(libstim PRIVATE -fPIC ${ARCH_OPT})
endif ()

set(SOURCE_FILES_NO_MAIN
        src/pyrematching/sparse_blossom/driver/namespaced_main.cc
        src/pyrematching/sparse_blossom/driver/mwpm_decoding.cc
        src/pyrematching/sparse_blossom/flooder/graph.cc
        src/pyrematching/sparse_blossom/flooder/detector_node.cc
        src/pyrematching/sparse_blossom/flooder_matcher_interop/compressed_edge.cc
        src/pyrematching/sparse_blossom/flooder/graph_fill_region.cc
        src/pyrematching/sparse_blossom/flooder/match.cc
        src/pyrematching/sparse_blossom/flooder/graph_flooder.cc
        src/pyrematching/sparse_blossom/matcher/alternating_tree.cc
        src/pyrematching/sparse_blossom/matcher/mwpm.cc
        src/pyrematching/sparse_blossom/flooder_matcher_interop/region_edge.cc
        src/pyrematching/sparse_blossom/flooder_matcher_interop/mwpm_event.cc
        src/pyrematching/sparse_blossom/tracker/flood_check_event.cc
        src/pyrematching/sparse_blossom/diagram/animation_main.cc
        src/pyrematching/sparse_blossom/diagram/mwpm_diagram.cc
        src/pyrematching/sparse_blossom/search/search_graph.cc
        src/pyrematching/sparse_blossom/search/search_detector_node.cc
        src/pyrematching/sparse_blossom/search/search_flooder.cc
        src/pyrematching/sparse_blossom/driver/user_graph.cc
        src/pyrematching/rand/rand_gen.cc
        src/pyrematching/two_phase/truncation/truncated_timeline.cc
        src/pyrematching/two_phase/truncation/harvest.cc
        src/pyrematching/two_phase/truncation/exposed_blossom.cc
        )

set(TEST_FILES
        src/pyrematching/sparse_blossom/driver/namespaced_main.test.cc
        src/pyrematching/sparse_blossom/driver/mwpm_decoding.test.cc
        src/pyrematching/sparse_blossom/flooder_matcher_interop/varying.test.cc
        src/pyrematching/sparse_blossom/flooder/graph.test.cc
        src/pyrematching/sparse_blossom/flooder/detector_node.test.cc
        src/pyrematching/sparse_blossom/flooder_matcher_interop/compressed_edge.test.cc
        src/pyrematching/sparse_blossom/flooder/graph_fill_region.test.cc
        src/pyrematching/sparse_blossom/flooder/match.test.cc
        src/pyrematching/sparse_blossom/flooder/graph_flooder.test.cc
        src/pyrematching/sparse_blossom/matcher/alternating_tree.test.cc
        src/pyrematching/sparse_blossom/matcher/mwpm.test.cc
        src/pyrematching/sparse_blossom/tracker/flood_check_event.test.cc
        src/pyrematching/sparse_blossom/tracker/radix_heap_queue.test.cc
        src/pyrematching/sparse_blossom/flooder_matcher_interop/mwpm_event.test.cc
        src/pyrematching/sparse_blossom/tracker/queued_event_tracker.test.cc
        src/pyrematching/sparse_blossom/tracker/cyclic.test.cc
        src/pyrematching/sparse_blossom/diagram/mwpm_diagram.test.cc
        src/pyrematching/sparse_blossom/search/search_graph.test.cc
        src/pyrematching/sparse_blossom/search/search_flooder.test.cc
        src/pyrematching/sparse_blossom/driver/user_graph.test.cc
        tests/two_phase/truncation_timeline.test.cc
        tests/two_phase/truncation_harvest.test.cc
        tests/two_phase/truncation_exposed_blossom.test.cc
        tests/two_phase/truncation_fuzz.test.cc
        )

set(PERF_FILES
        src/pyrematching/perf/main.perf.cc
        src/pyrematching/perf/util.perf.cc
        src/pyrematching/sparse_blossom/driver/mwpm_decoding.perf.cc
        src/pyrematching/sparse_blossom/driver/user_graph.perf.cc
        src/pyrematching/sparse_blossom/flooder_matcher_interop/varying.perf.cc
        src/pyrematching/sparse_blossom/tracker/radix_heap_queue.perf.cc
        )

set(PYTHON_API_FILES
        src/pyrematching/sparse_blossom/driver/user_graph.pybind.cc
        src/pyrematching/rand/rand_gen.pybind.cc
        src/pyrematching/pyrematching.pybind.cc
        )

add_executable(pyrematching src/main.cc ${SOURCE_FILES_NO_MAIN})

target_compile_options(pyrematching PRIVATE ${ARCH_OPT})
if (NOT (MSVC))
    target_link_options(pyrematching PRIVATE -pthread ${ARCH_OPT})
endif ()
target_link_libraries(pyrematching libstim)
install(TARGETS pyrematching RUNTIME DESTINATION bin)

enable_testing()

add_executable(pyrematching_tests ${SOURCE_FILES_NO_MAIN} ${TEST_FILES})
# The two-phase tests share a helper header under tests/, which is outside the `src` include root.
target_include_directories(pyrematching_tests PRIVATE ${CMAKE_CURRENT_LIST_DIR})
# -UNDEBUG so the debug-build invariants (asserts, and the M1 horizon gate counters) are live in the
# test binary even when the tree is configured as a Release build.
target_compile_options(pyrematching_tests PRIVATE -fsanitize=address -g -fsanitize=undefined -fno-omit-frame-pointer -coverage -UNDEBUG)
if (NOT (MSVC))
    target_link_options(pyrematching_tests PRIVATE -pthread -fsanitize=address -fsanitize=undefined -coverage)
else ()
    target_link_options(pyrematching_tests PRIVATE -fsanitize=address -fsanitize=undefined -coverage)
endif ()
target_link_libraries(pyrematching_tests GTest::gtest_main GTest::gmock_main libstim)

add_executable(pyrematching_perf ${SOURCE_FILES_NO_MAIN} ${PERF_FILES})
target_compile_options(pyrematching_perf PRIVATE ${ARCH_OPT})
if (NOT (MSVC))
    target_link_options(pyrematching_perf PRIVATE -pthread -O3)
endif ()
target_link_libraries(pyrematching_perf libstim)

# M1 exit artifact: sweeps (d, p, T) and reports the residual/fallback statistics the M1 go/no-go
# read is made from. Built optimised, like the perf harness.
add_executable(two_phase_m1_artifact benchmarks/two_phase/m1_exit_artifact.cc ${SOURCE_FILES_NO_MAIN})
target_compile_options(two_phase_m1_artifact PRIVATE ${ARCH_OPT})
if (NOT (MSVC))
    target_link_options(two_phase_m1_artifact PRIVATE -pthread -O3)
endif ()
target_link_libraries(two_phase_m1_artifact libstim)

add_library(libpyrematching ${SOURCE_FILES_NO_MAIN})
set_target_properties(libpyrematching PROPERTIES PREFIX "")
target_include_directories(libpyrematching PUBLIC src)
target_compile_options(libpyrematching PRIVATE ${ARCH_OPT} -fPIC)
if(NOT(MSVC))
    target_link_options(libpyrematching PRIVATE -pthread -O3)
endif()
target_link_libraries(libpyrematching libstim)
install(TARGETS libpyrematching LIBRARY DESTINATION)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/" DESTINATION "include" FILES_MATCHING PATTERN "*.h" PATTERN "*.inl")

include(GoogleTest)
gtest_discover_tests(pyrematching_tests)

IF(EXISTS ${CMAKE_CURRENT_LIST_DIR}/extern/pybind11/CMakeLists.txt)
    add_subdirectory(extern/pybind11)
    pybind11_add_module(_cpp_pyrematching ${PYTHON_API_FILES} ${SOURCE_FILES_NO_MAIN})
    target_link_libraries(_cpp_pyrematching PRIVATE libstim)
    target_compile_options(_cpp_pyrematching PRIVATE ${ARCH_OPT})
else()
    message("WARNING: Skipped the pybind11 module _cpp_pyrematching because the `pybind11` git submodule isn't present. To fix, run `git submodule update --init --recursive`")
endif()

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
#target_compile_definitions(_cpp_pyrematching
#        PRIVATE VERSION_INFO=${PYREMATCHING_VERSION_INFO})
