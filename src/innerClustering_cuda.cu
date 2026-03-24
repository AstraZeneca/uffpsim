#include <cstdint>
#include <limits>

#include <cuda_runtime.h>

namespace {

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

    const uint64_t *clusterFp = clusterFps + (static_cast<size_t>(clusterIdx) * cfpSize);
    uint64_t commonPopcnt = 0;
    for (size_t offset = 0; offset < fpSize; offset++) {
        commonPopcnt += static_cast<uint64_t>(__popcll(clusterFp[molIdOffset + offset] & currentFp[molIdOffset + offset]));
    }

    float distFactor = threshold * static_cast<float>(
        currentFp[cfpPopcountIndex] + clusterFp[cfpPopcountIndex] - commonPopcnt
    );

    if (static_cast<float>(commonPopcnt) >= distFactor) {
        atomicMin(bestClusterId, clusterIdx);
    }
}

} // namespace

int find_matching_cluster_cuda(
    const uint64_t *cluster_fps,
    const uint64_t *current_fp,
    int num_clusters,
    size_t cfp_size,
    size_t mol_id_offset,
    size_t fp_size,
    size_t cfp_popcount_index,
    float threshold
) {
    if (num_clusters <= 0) {
        return -1;
    }

    const size_t clustersBytes = static_cast<size_t>(num_clusters) * cfp_size * sizeof(uint64_t);
    const size_t currentFpBytes = cfp_size * sizeof(uint64_t);

    uint64_t *dClusterFps = nullptr;
    uint64_t *dCurrentFp = nullptr;
    int *dBestClusterId = nullptr;

    if (cudaMalloc(&dClusterFps, clustersBytes) != cudaSuccess) {
        return -2;
    }
    if (cudaMalloc(&dCurrentFp, currentFpBytes) != cudaSuccess) {
        cudaFree(dClusterFps);
        return -2;
    }
    if (cudaMalloc(&dBestClusterId, sizeof(int)) != cudaSuccess) {
        cudaFree(dClusterFps);
        cudaFree(dCurrentFp);
        return -2;
    }

    if (cudaMemcpy(dClusterFps, cluster_fps, clustersBytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(dCurrentFp, current_fp, currentFpBytes, cudaMemcpyHostToDevice) != cudaSuccess) {
        cudaFree(dClusterFps);
        cudaFree(dCurrentFp);
        cudaFree(dBestClusterId);
        return -2;
    }

    int initialBest = std::numeric_limits<int>::max();
    if (cudaMemcpy(dBestClusterId, &initialBest, sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) {
        cudaFree(dClusterFps);
        cudaFree(dCurrentFp);
        cudaFree(dBestClusterId);
        return -2;
    }

    constexpr int threadsPerBlock = 256;
    const int blocks = (num_clusters + threadsPerBlock - 1) / threadsPerBlock;
    findMatchingClusterKernel<<<blocks, threadsPerBlock>>>(
        dClusterFps,
        dCurrentFp,
        num_clusters,
        cfp_size,
        mol_id_offset,
        fp_size,
        cfp_popcount_index,
        threshold,
        dBestClusterId
    );

    if (cudaGetLastError() != cudaSuccess || cudaDeviceSynchronize() != cudaSuccess) {
        cudaFree(dClusterFps);
        cudaFree(dCurrentFp);
        cudaFree(dBestClusterId);
        return -2;
    }

    int bestClusterId = std::numeric_limits<int>::max();
    if (cudaMemcpy(&bestClusterId, dBestClusterId, sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess) {
        cudaFree(dClusterFps);
        cudaFree(dCurrentFp);
        cudaFree(dBestClusterId);
        return -2;
    }

    cudaFree(dClusterFps);
    cudaFree(dCurrentFp);
    cudaFree(dBestClusterId);

    if (bestClusterId == std::numeric_limits<int>::max()) {
        return -1;
    }
    return bestClusterId;
}
