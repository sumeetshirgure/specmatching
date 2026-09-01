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

#include "specmatching/spec_matching/manifold/ball_serialize.h"

#include <algorithm>
#include <cstdio>
#include <stdexcept>

namespace pm {
namespace spec_matching {

namespace {

const char BALL_ARTIFACT_MAGIC[8] = {'P', 'Y', 'R', 'M', 'B', 'A', 'L', 'L'};

struct Writer {
    std::FILE* handle;

    template <typename T>
    void pod(const T& value) {
        if (std::fwrite(&value, sizeof(T), 1, handle) != 1)
            throw std::runtime_error("Failed to write ball artifact");
    }

    template <typename T>
    void vec(const std::vector<T>& values) {
        pod<uint64_t>(values.size());
        if (!values.empty() && std::fwrite(values.data(), sizeof(T), values.size(), handle) != values.size())
            throw std::runtime_error("Failed to write ball artifact");
    }
};

struct Reader {
    std::FILE* handle;

    template <typename T>
    T pod() {
        T value;
        if (std::fread(&value, sizeof(T), 1, handle) != 1)
            throw std::runtime_error("Ball artifact is truncated");
        return value;
    }

    template <typename T>
    void vec(std::vector<T>& out) {
        uint64_t size = pod<uint64_t>();
        out.resize(size);
        if (size != 0 && std::fread(out.data(), sizeof(T), size, handle) != size)
            throw std::runtime_error("Ball artifact is truncated");
    }
};

/// Everything in `BallTables` except the pools, in one block so the header stays one `fwrite`.
struct ArtifactHeader {
    char magic[8];
    uint32_t version;
    uint32_t reserved;
    BallParams params;
    double normalising_constant;
    int64_t r_int;
    int64_t t_max_int;
    int64_t shell_width_int;
    uint64_t num_shells;
    uint64_t num_nodes;
    uint64_t num_observables;
    uint64_t graph_hash;
};

void write_stats(Writer& writer, const BallStats& stats) {
    writer.pod(stats.num_nodes);
    writer.pod(stats.total_entries);
    writer.pod(stats.mean_ball_size);
    writer.pod(stats.max_ball_size);
    writer.vec(stats.shell_entry_counts);
    writer.pod(stats.nodes_with_boundary);
    writer.pod(stats.mean_ball_word_len);
    writer.pod(stats.max_ball_word_len);
    writer.pod(stats.ambiguous_mask_pairs);
    writer.pod(stats.bytes_targets);
    writer.pod(stats.bytes_weights);
    writer.pod(stats.bytes_masks);
    writer.pod(stats.bytes_paths);
    writer.pod(stats.bytes_words);
    writer.pod(stats.bytes_boundary);
    writer.pod(stats.bytes_offsets);
    writer.pod(stats.bytes_total);
    writer.pod(stats.compile_wall_seconds);
}

void read_stats(Reader& reader, BallStats& stats) {
    stats.num_nodes = reader.pod<uint64_t>();
    stats.total_entries = reader.pod<uint64_t>();
    stats.mean_ball_size = reader.pod<double>();
    stats.max_ball_size = reader.pod<uint64_t>();
    reader.vec(stats.shell_entry_counts);
    stats.nodes_with_boundary = reader.pod<uint64_t>();
    stats.mean_ball_word_len = reader.pod<double>();
    stats.max_ball_word_len = reader.pod<uint64_t>();
    stats.ambiguous_mask_pairs = reader.pod<uint64_t>();
    stats.bytes_targets = reader.pod<uint64_t>();
    stats.bytes_weights = reader.pod<uint64_t>();
    stats.bytes_masks = reader.pod<uint64_t>();
    stats.bytes_paths = reader.pod<uint64_t>();
    stats.bytes_words = reader.pod<uint64_t>();
    stats.bytes_boundary = reader.pod<uint64_t>();
    stats.bytes_offsets = reader.pod<uint64_t>();
    stats.bytes_total = reader.pod<uint64_t>();
    stats.compile_wall_seconds = reader.pod<double>();
}

/// Closes the artifact handle however the scope is left, so a throw mid-read does not leak it.
struct FileGuard {
    std::FILE* handle;
    ~FileGuard() {
        if (handle != nullptr)
            std::fclose(handle);
    }
};

}  // namespace

void save_ball_tables(const BallTables& tables, const std::string& path) {
    std::FILE* handle = std::fopen(path.c_str(), "wb");
    if (handle == nullptr)
        throw std::runtime_error("Failed to open ball artifact for writing: " + path);
    FileGuard guard{handle};
    Writer writer{handle};

    ArtifactHeader header{};
    std::copy(std::begin(BALL_ARTIFACT_MAGIC), std::end(BALL_ARTIFACT_MAGIC), header.magic);
    header.version = BALL_ARTIFACT_VERSION;
    header.params = tables.params;
    header.normalising_constant = tables.normalising_constant;
    header.r_int = tables.r_int;
    header.t_max_int = tables.t_max_int;
    header.shell_width_int = tables.shell_width_int;
    header.num_shells = tables.num_shells;
    header.num_nodes = tables.num_nodes;
    header.num_observables = tables.num_observables;
    header.graph_hash = tables.graph_hash;
    writer.pod(header);

    writer.vec(tables.ball_offsets);
    writer.vec(tables.ball_shell_offsets);
    writer.vec(tables.ball_target);
    writer.vec(tables.ball_w_int);
    writer.vec(tables.ball_mask_offsets);
    writer.vec(tables.ball_mask_ids);
    writer.vec(tables.ball_path_offsets);
    writer.vec(tables.ball_path_nodes);
    writer.vec(tables.has_bcost);
    writer.vec(tables.bcost_w_int);
    writer.vec(tables.bcost_mask_offsets);
    writer.vec(tables.bcost_mask_ids);
    writer.vec(tables.bcost_path_offsets);
    writer.vec(tables.bcost_path_nodes);
    writer.vec(tables.ball_word_offsets);
    writer.vec(tables.ball_word_base);
    writer.vec(tables.ball_words);
    writer.vec(tables.ball_word_rank);
    writer.vec(tables.ball_entry_by_rank);
    write_stats(writer, tables.stats);
}

BallTables load_ball_tables(const std::string& path, const pm::MatchingGraph& graph, const BallParams& expected) {
    std::FILE* handle = std::fopen(path.c_str(), "rb");
    if (handle == nullptr)
        throw std::runtime_error("Failed to open ball artifact for reading: " + path);
    FileGuard guard{handle};
    Reader reader{handle};

    auto header = reader.pod<ArtifactHeader>();
    if (!std::equal(std::begin(BALL_ARTIFACT_MAGIC), std::end(BALL_ARTIFACT_MAGIC), header.magic))
        throw std::runtime_error("Not a ball artifact: " + path);
    if (header.version != BALL_ARTIFACT_VERSION)
        throw std::runtime_error(
            "Ball artifact " + path + " has version " + std::to_string(header.version) + ", expected " +
            std::to_string(BALL_ARTIFACT_VERSION));
    if (header.params != expected)
        throw std::runtime_error("Ball artifact " + path + " was compiled with different BallParams");
    uint64_t hash = hash_matching_graph(graph);
    if (header.graph_hash != hash)
        throw std::runtime_error("Ball artifact " + path + " was compiled against a different detector graph");
    if (header.num_nodes != graph.nodes.size())
        throw std::runtime_error("Ball artifact " + path + " has the wrong node count");

    BallTables tables;
    tables.params = header.params;
    tables.normalising_constant = header.normalising_constant;
    tables.r_int = header.r_int;
    tables.t_max_int = header.t_max_int;
    tables.shell_width_int = header.shell_width_int;
    tables.num_shells = header.num_shells;
    tables.num_nodes = header.num_nodes;
    tables.num_observables = header.num_observables;
    tables.graph_hash = header.graph_hash;

    reader.vec(tables.ball_offsets);
    reader.vec(tables.ball_shell_offsets);
    reader.vec(tables.ball_target);
    reader.vec(tables.ball_w_int);
    reader.vec(tables.ball_mask_offsets);
    reader.vec(tables.ball_mask_ids);
    reader.vec(tables.ball_path_offsets);
    reader.vec(tables.ball_path_nodes);
    reader.vec(tables.has_bcost);
    reader.vec(tables.bcost_w_int);
    reader.vec(tables.bcost_mask_offsets);
    reader.vec(tables.bcost_mask_ids);
    reader.vec(tables.bcost_path_offsets);
    reader.vec(tables.bcost_path_nodes);
    reader.vec(tables.ball_word_offsets);
    reader.vec(tables.ball_word_base);
    reader.vec(tables.ball_words);
    reader.vec(tables.ball_word_rank);
    reader.vec(tables.ball_entry_by_rank);
    read_stats(reader, tables.stats);
    return tables;
}

}  // namespace spec_matching
}  // namespace pm
