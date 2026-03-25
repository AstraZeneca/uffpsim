#pragma once

#include <cstddef>
#include <cstdint>

struct CudaInnerClusteringContext;

CudaInnerClusteringContext* create_cuda_inner_clustering_context(
    size_t cfp_size,
    size_t mol_id_offset,
    size_t fp_size,
    size_t cfp_popcount_index,
    float threshold,
    size_t initial_cluster_capacity,
    size_t initial_chunk_elements
);

void destroy_cuda_inner_clustering_context(CudaInnerClusteringContext* context);

int upload_cfp_chunk_cuda(
    CudaInnerClusteringContext* context,
    const uint64_t *host_cfp_chunk,
    size_t chunk_elements
);

int upsert_cluster_fp_cuda(
    CudaInnerClusteringContext* context,
    const uint64_t *host_cluster_fp,
    size_t cluster_index
);

int find_matching_cluster_in_chunk_cuda(
    CudaInnerClusteringContext* context,
    size_t chunk_fp_index,
    int num_clusters
);

int update_cluster_fp_with_chunk_fp_cuda(
    CudaInnerClusteringContext* context,
    size_t chunk_fp_index,
    size_t cluster_index
);

int download_cluster_fps_cuda(
    CudaInnerClusteringContext* context,
    uint64_t *host_cluster_fps,
    size_t num_clusters
);
