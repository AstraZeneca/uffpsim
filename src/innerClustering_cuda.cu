#include <cstdint>
#include <limits>
#include <stdio.h>

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

    if (cfpPopcountIndex >= cfpSize || molIdOffset >= cfpSize) {
        return;
    }

    const size_t fpEnd = molIdOffset + fpSize;
    if (fpEnd > cfpSize) {
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

    if (cluster_fps == nullptr || current_fp == nullptr) {
        fprintf(stderr, "CUDA error: null host pointer passed to find_matching_cluster_cuda\n");
        return -2;
    }

    if (cfp_size == 0 || cfp_popcount_index >= cfp_size || mol_id_offset >= cfp_size) {
        fprintf(stderr, "CUDA error: invalid CFP layout (cfp_size=%zu, pop_idx=%zu, mol_offset=%zu)\n",
                cfp_size, cfp_popcount_index, mol_id_offset);
        return -2;
    }

    if (fp_size > (cfp_size - mol_id_offset)) {
        fprintf(stderr, "CUDA error: invalid fp window (cfp_size=%zu, mol_offset=%zu, fp_size=%zu)\n",
                cfp_size, mol_id_offset, fp_size);
        return -2;
    }

    const size_t clustersBytes = static_cast<size_t>(num_clusters) * cfp_size * sizeof(uint64_t);
    const size_t currentFpBytes = cfp_size * sizeof(uint64_t);

    uint64_t *dClusterFps = nullptr;
    uint64_t *dCurrentFp = nullptr;
    int *dBestClusterId = nullptr;

    if (cudaMalloc(reinterpret_cast<void **>(&dClusterFps), clustersBytes) != cudaSuccess) {
        return -2;
    }

    if (cudaMalloc(reinterpret_cast<void **>(&dCurrentFp), currentFpBytes) != cudaSuccess) {
        cudaFree(dClusterFps);
        return -2;
    }

    if (cudaMalloc(reinterpret_cast<void **>(&dBestClusterId), sizeof(int)) != cudaSuccess) {
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

    cudaError_t launchErr = cudaGetLastError();
    if (launchErr != cudaSuccess) {
        cudaFree(dClusterFps);
        cudaFree(dCurrentFp);
        cudaFree(dBestClusterId);
        fprintf(stderr, "CUDA kernel launch error: %s\n", cudaGetErrorString(launchErr));
        return -2;
    }

    cudaError_t syncErr = cudaDeviceSynchronize();
    if (syncErr != cudaSuccess) {
        cudaFree(dClusterFps);
        cudaFree(dCurrentFp);
        cudaFree(dBestClusterId);
        fprintf(stderr, "CUDA kernel sync error: %s\n", cudaGetErrorString(syncErr));
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
