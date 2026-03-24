#pragma once

#include <cstddef>
#include <cstdint>

int find_matching_cluster_cuda(
    const uint64_t *cluster_fps,
    const uint64_t *current_fp,
    int num_clusters,
    size_t cfp_size,
    size_t mol_id_offset,
    size_t fp_size,
    size_t cfp_popcount_index,
    float threshold
);
