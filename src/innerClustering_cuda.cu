#include <cstdint>
#include <limits>
#include <climits>
#include <stdio.h>
#include <cstring>

#include <cuda_runtime.h>

#include "innerClustering_cuda.h"

struct CudaInnerClusteringContext {
    size_t cfpSize = 0;
    size_t molIdOffset = 0;
    size_t fpSize = 0;
    size_t cfpPopcountIndex = 0;
    float threshold = 0.0f;

    uint64_t *dClusterFps = nullptr;
    uint64_t *dChunkFps = nullptr;
    int *dBestClusterId = nullptr;

    size_t clusterCapacity = 0;
    size_t chunkCapacityElements = 0;
    size_t loadedChunkElements = 0;
};

namespace {

inline bool validate_layout(size_t cfpSize, size_t molIdOffset, size_t fpSize, size_t cfpPopcountIndex) {
    if (cfpSize == 0 || cfpPopcountIndex >= cfpSize || molIdOffset >= cfpSize) {
        return false;
    }
    const size_t fpEnd = molIdOffset + fpSize;
    if (fpEnd > cfpSize) {
        return false;
    }
    // cfpPopcountIndex must lie outside the FP data range; otherwise updateClusterFpKernel
    // would corrupt it with an OR'd FP word.
    if (cfpPopcountIndex >= molIdOffset && cfpPopcountIndex < fpEnd) {
        return false;
    }
    return true;
}

inline int ensure_cluster_capacity(CudaInnerClusteringContext *context, size_t requiredCapacity) {
    if (requiredCapacity <= context->clusterCapacity) {
        return 0;
    }

    size_t newCapacity = context->clusterCapacity == 0 ? 1 : context->clusterCapacity;
    while (newCapacity < requiredCapacity) {
        if (newCapacity > SIZE_MAX / 2) {
            return -2;  // overflow guard: cannot double further
        }
        newCapacity *= 2;
    }

    // Guard against size_t overflow in byte calculation
    if (newCapacity > SIZE_MAX / context->cfpSize / sizeof(uint64_t)) {
        return -2;
    }
    uint64_t *newDevicePtr = nullptr;
    const size_t newBytes = newCapacity * context->cfpSize * sizeof(uint64_t);
    if (cudaMalloc(reinterpret_cast<void **>(&newDevicePtr), newBytes) != cudaSuccess) {
        return -2;
    }

    if (context->dClusterFps != nullptr && context->clusterCapacity > 0) {
        const size_t oldBytes = context->clusterCapacity * context->cfpSize * sizeof(uint64_t);
        if (cudaMemcpy(newDevicePtr, context->dClusterFps, oldBytes, cudaMemcpyDeviceToDevice) != cudaSuccess) {
            cudaFree(newDevicePtr);
            return -2;
        }
        cudaFree(context->dClusterFps);
    }

    context->dClusterFps = newDevicePtr;
    context->clusterCapacity = newCapacity;
    return 0;
}

inline int ensure_chunk_capacity(CudaInnerClusteringContext *context, size_t requiredElements) {
    if (requiredElements <= context->chunkCapacityElements) {
        return 0;
    }

    size_t newCapacity = context->chunkCapacityElements == 0 ? 1 : context->chunkCapacityElements;
    while (newCapacity < requiredElements) {
        if (newCapacity > SIZE_MAX / 2) {
            return -2;  // overflow guard: cannot double further
        }
        newCapacity *= 2;
    }

    uint64_t *newDevicePtr = nullptr;
    const size_t newBytes = newCapacity * sizeof(uint64_t);
    if (cudaMalloc(reinterpret_cast<void **>(&newDevicePtr), newBytes) != cudaSuccess) {
        return -2;
    }

    if (context->dChunkFps != nullptr) {
        cudaFree(context->dChunkFps);
    }

    context->dChunkFps = newDevicePtr;
    context->chunkCapacityElements = newCapacity;
    return 0;
}

__global__ void findMatchingClusterKernel(
    const uint64_t *clusterFps,
    const uint64_t *currentFp,
    int numClusters,
    size_t cfpSize,
    size_t molIdOffset,
    size_t fpSize,
    size_t cfpPopcountIndex,
    float threshold,
    int *bestClusterId
) {
    int clusterIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (clusterIdx >= numClusters) {
        return;
    }

    // A plain global-memory load is not sufficient here because another thread may
    // update bestClusterId concurrently. Use an atomic read so this thread sees a
    // coherent value. We can only skip work when our clusterIdx cannot improve the
    // current best match; otherwise a lower-index match could still exist.
    const int bestSoFar = atomicAdd(bestClusterId, 0);
    if (bestSoFar != INT_MAX && clusterIdx >= bestSoFar) {
        return;
    }

    const uint64_t *clusterFp = clusterFps + (static_cast<size_t>(clusterIdx) * cfpSize);
    uint64_t commonPopcnt = 0;
    for (size_t offset = 0; offset < fpSize; offset++) {
        const int updatedBest = atomicAdd(bestClusterId, 0);
        if (updatedBest != INT_MAX && clusterIdx >= updatedBest) {
            return;
        }
        commonPopcnt += static_cast<uint64_t>(__popcll(clusterFp[molIdOffset + offset] & currentFp[molIdOffset + offset]));
    }

    float distFactor = threshold * static_cast<float>(
        currentFp[cfpPopcountIndex] + clusterFp[cfpPopcountIndex] - commonPopcnt
    );

    if (static_cast<float>(commonPopcnt) >= distFactor) {
        atomicMin(bestClusterId, clusterIdx);
    }
}

__global__ void updateClusterFpKernel(
    uint64_t *clusterFps,
    const uint64_t *currentFp,
    size_t clusterIndex,
    size_t cfpSize,
    size_t molIdOffset,
    size_t fpSize,
    size_t cfpPopcountIndex
) {
    const uint64_t *src = currentFp;
    uint64_t *clusterFp = clusterFps + (clusterIndex * cfpSize);

    if (threadIdx.x == 0) {
        // Use an atomic store so the zero-write is coherent with subsequent atomicAdds
        // even if the launch configuration is ever changed to use multiple blocks.
        atomicExch(reinterpret_cast<unsigned long long *>(&clusterFp[cfpPopcountIndex]),
                   static_cast<unsigned long long>(0));
    }
    __syncthreads();

    for (size_t local = threadIdx.x; local < fpSize; local += blockDim.x) {
        size_t idx = molIdOffset + local;
        uint64_t merged = clusterFp[idx] | src[idx];
        clusterFp[idx] = merged;
        atomicAdd(reinterpret_cast<unsigned long long *>(&clusterFp[cfpPopcountIndex]),
                  static_cast<unsigned long long>(__popcll(merged)));
    }
}

} // namespace

CudaInnerClusteringContext* create_cuda_inner_clustering_context(
    size_t cfp_size,
    size_t mol_id_offset,
    size_t fp_size,
    size_t cfp_popcount_index,
    float threshold,
    size_t initial_cluster_capacity,
    size_t initial_chunk_elements
) {
    if (!validate_layout(cfp_size, mol_id_offset, fp_size, cfp_popcount_index)) {
        fprintf(stderr, "CUDA error: invalid CFP layout in context creation\n");
        return nullptr;
    }

    CudaInnerClusteringContext *context = new CudaInnerClusteringContext();
    context->cfpSize = cfp_size;
    context->molIdOffset = mol_id_offset;
    context->fpSize = fp_size;
    context->cfpPopcountIndex = cfp_popcount_index;
    context->threshold = threshold;

    if (cudaMalloc(reinterpret_cast<void **>(&context->dBestClusterId), sizeof(int)) != cudaSuccess) {
        delete context;
        return nullptr;
    }

    if (ensure_cluster_capacity(context, initial_cluster_capacity > 0 ? initial_cluster_capacity : 1) != 0) {
        cudaFree(context->dBestClusterId);
        delete context;
        return nullptr;
    }

    if (ensure_chunk_capacity(context, initial_chunk_elements > 0 ? initial_chunk_elements : cfp_size) != 0) {
        cudaFree(context->dClusterFps);
        cudaFree(context->dBestClusterId);
        delete context;
        return nullptr;
    }

    return context;
}

void destroy_cuda_inner_clustering_context(CudaInnerClusteringContext* context) {
    if (context == nullptr) {
        return;
    }
    if (context->dClusterFps != nullptr) {
        cudaFree(context->dClusterFps);
    }
    if (context->dChunkFps != nullptr) {
        cudaFree(context->dChunkFps);
    }
    if (context->dBestClusterId != nullptr) {
        cudaFree(context->dBestClusterId);
    }
    delete context;
}

int upload_cfp_chunk_cuda(
    CudaInnerClusteringContext* context,
    const uint64_t *host_cfp_chunk,
    size_t chunk_elements
) {
    if (context == nullptr || host_cfp_chunk == nullptr || chunk_elements == 0) {
        return -2;
    }

    if (ensure_chunk_capacity(context, chunk_elements) != 0) {
        return -2;
    }

    const size_t chunkBytes = chunk_elements * sizeof(uint64_t);
    if (cudaMemcpy(context->dChunkFps, host_cfp_chunk, chunkBytes, cudaMemcpyHostToDevice) != cudaSuccess) {
        return -2;
    }
    context->loadedChunkElements = chunk_elements;
    return 0;
}

int upsert_cluster_fp_cuda(
    CudaInnerClusteringContext* context,
    const uint64_t *host_cluster_fp,
    size_t cluster_index
) {
    if (context == nullptr || host_cluster_fp == nullptr) {
        return -2;
    }

    if (ensure_cluster_capacity(context, cluster_index + 1) != 0) {
        return -2;
    }

    const size_t offset = cluster_index * context->cfpSize;
    const size_t fullBytes = context->cfpSize * sizeof(uint64_t);
    const size_t tailOffset = context->molIdOffset;
    const size_t tailCount = context->cfpSize - context->molIdOffset;
    const size_t tailBytes = tailCount * sizeof(uint64_t);

    if (cudaMemset(context->dClusterFps + offset, 0, fullBytes) != cudaSuccess) {
        return -2;
    }

    if (cudaMemcpy(
            context->dClusterFps + offset + tailOffset,
            host_cluster_fp + tailOffset,
            tailBytes,
            cudaMemcpyHostToDevice
        ) != cudaSuccess) {
        return -2;
    }
    return 0;
}

int find_matching_cluster_in_chunk_cuda(
    CudaInnerClusteringContext* context,
    size_t chunk_fp_index,
    int num_clusters
) {
    if (context == nullptr || num_clusters <= 0) {
        return -1;
    }

    const size_t chunk_offset = chunk_fp_index * context->cfpSize;
    if (chunk_offset + context->cfpSize > context->loadedChunkElements) {
        return -2;
    }

    int initialBest = INT_MAX;
    if (cudaMemcpy(context->dBestClusterId, &initialBest, sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) {
        return -2;
    }

    constexpr int threadsPerBlock = 512;
    const int blocks = (num_clusters + threadsPerBlock - 1) / threadsPerBlock;
    findMatchingClusterKernel<<<blocks, threadsPerBlock>>>(
        context->dClusterFps,
        context->dChunkFps + chunk_offset,
        num_clusters,
        context->cfpSize,
        context->molIdOffset,
        context->fpSize,
        context->cfpPopcountIndex,
        context->threshold,
        context->dBestClusterId
    );

    cudaError_t launchErr = cudaGetLastError();
    if (launchErr != cudaSuccess) {
        fprintf(stderr, "CUDA kernel launch error: %s\n", cudaGetErrorString(launchErr));
        return -2;
    }

    cudaError_t syncErr = cudaDeviceSynchronize();
    if (syncErr != cudaSuccess) {
        fprintf(stderr, "CUDA kernel sync error: %s\n", cudaGetErrorString(syncErr));
        return -2;
    }

    int bestClusterId = INT_MAX;
    if (cudaMemcpy(&bestClusterId, context->dBestClusterId, sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess) {
        return -2;
    }

    if (bestClusterId == INT_MAX) {
        return -1;
    }
    return bestClusterId;
}

int update_cluster_fp_with_chunk_fp_cuda(
    CudaInnerClusteringContext* context,
    size_t chunk_fp_index,
    size_t cluster_index
) {
    if (context == nullptr) {
        return -2;
    }

    const size_t chunk_offset = chunk_fp_index * context->cfpSize;
    if (chunk_offset + context->cfpSize > context->loadedChunkElements) {
        return -2;
    }

    if (cluster_index >= context->clusterCapacity) {
        return -2;
    }

    constexpr int threadsPerBlock = 512;
    updateClusterFpKernel<<<1, threadsPerBlock>>>(
        context->dClusterFps,
        context->dChunkFps + chunk_offset,
        cluster_index,
        context->cfpSize,
        context->molIdOffset,
        context->fpSize,
        context->cfpPopcountIndex
    );

    cudaError_t launchErr = cudaGetLastError();
    if (launchErr != cudaSuccess) {
        fprintf(stderr, "CUDA update kernel launch error: %s\n", cudaGetErrorString(launchErr));
        return -2;
    }

    cudaError_t syncErr = cudaDeviceSynchronize();
    if (syncErr != cudaSuccess) {
        fprintf(stderr, "CUDA update kernel sync error: %s\n", cudaGetErrorString(syncErr));
        return -2;
    }

    return 0;
}

int download_cluster_fps_cuda(
    CudaInnerClusteringContext* context,
    uint64_t *host_cluster_fps,
    size_t num_clusters
) {
    if (context == nullptr || host_cluster_fps == nullptr) {
        return -2;
    }

    if (num_clusters == 0) {
        return 0;
    }

    if (num_clusters > context->clusterCapacity) {
        fprintf(stderr, "CUDA error: download_cluster_fps_cuda: num_clusters (%zu) > clusterCapacity (%zu)\n",
                num_clusters, context->clusterCapacity);
        return -2;
    }

    const size_t bytes = num_clusters * context->cfpSize * sizeof(uint64_t);
    if (cudaMemcpy(host_cluster_fps, context->dClusterFps, bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
        return -2;
    }
    return 0;
}
