/**
 * @brief This class represents an inner clustering agent which performs inner clustering within the popcount bins.
 * @author Rajendra Kumar
 */

#include <string>
#include <iostream>
#include <algorithm>
#include <map>
#include <tuple>
#include <set>
#include <vector>
#include <cstring>
#include <omp.h>
#include <H5Cpp.h>

#include "innerClustering.h"
#include "fpstore.h"
#ifdef USE_CUDA
#include "innerClustering_cuda.h"
#endif

InnerClusteringAgent::InnerClusteringAgent(FingerprintStore *fpstore, std::string mode, bool parallel) {
    _fpstore = fpstore;
    if (mode == "memory" || mode == "disk" || mode == "cuda")
        _mode = mode;
    else
        throw "clustering mode should be either memory, disk, or cuda";

    _parallel = parallel;
}

bool InnerClusteringAgent::performInnerClusteringAndWrite() {
    std::cout<< "Performing inner clustering ..." << std::endl;
    uint64_t start_time = utils::get_posix_clock_time();
    h5::Group *popcountBinsGroup =  utils::createOrOpenGroup(_fpstore->_root_group, _fpstore->_popCountBinsGroupName);

    // for each popcount bin
    for (auto popcount: _fpstore->_popCountBins) {
        h5::Group *popcountGroup = utils::createOrOpenGroup(popcountBinsGroup, _fpstore->_getPopCountGroupName(popcount));

        // it checks whether inner clustering has been done before with same threshold
        bool to_be_done = true;
        if (popcountGroup->attrExists("innerClusteringThreshold")) {
            float previousThreshold = utils::get_scaler_attribute<float>(popcountGroup, "innerClusteringThreshold");
            if (_fpstore->_innerClusteringThreshold == previousThreshold) to_be_done = false;
        }

        // perform inner clustering
        if (to_be_done) {
            if (_mode == "cuda")
                _doInnerClusteringCUDA(popcount, popcountGroup);
            else
                _doInnerClustering(popcount, popcountGroup);
        }

        delete popcountGroup;
    }

    delete popcountBinsGroup;
    std::cout << " ... inner clustering done. Time taken: " << (utils::get_posix_clock_time() - start_time) / (1000 *60) << " mins." << std::endl;
    py::gil_scoped_acquire acquire;
    return true;
}

void InnerClusteringAgent::_doInnerClustering(int popCount, h5::Group *popcountGroup) {
    hsize_t totalSize = utils::sizeOfFingerPrintDatasetInH5(popcountGroup, "cfpData");
    unsigned int num_fps_per_chunk = 1000; 
    hsize_t chunkSize = _fpstore->_CFPSize * num_fps_per_chunk; // stored in same number of chunks, might be efficient
    uint64_t *cfp = (uint64_t *) malloc(sizeof(uint64_t) * chunkSize); // allocated memory for a cfp chunk, will be filled in each iteration
    std::vector<utils::dt_inner_clusters_fingerprints> clusters; // to store the inner clusters
    uint64_t chunkCounter = 0, remainingSize = totalSize, currentSize = 0, fpIndex=0;
    std::vector<uint64_t> cFpIndexToClusterID; // store the index of cfp in the original dataset as consecutive elements and its cluster id for disk mode
    unsigned long sortCounter = 0;

    while (remainingSize > 0) {// read from file and perform clustering, low level of memory usage
        currentSize = remainingSize < chunkSize ? remainingSize : chunkSize;
        utils::getFingerprintFromIndex(popcountGroup, "cfpData", fpIndex, currentSize, cfp);
        for(size_t i = 0; i < currentSize; i += _fpstore->_CFPSize, fpIndex += _fpstore->_CFPSize) {
            if (cfp[i + _fpstore->_CFPPopCountIndex] > 0) { // for edge case when fp is removed and set to zero.
            
            if (_mode == "memory") // add this fingerprint to clusters where cfp in each cluster remains in memory
                _addFingerprintToCluster(&cfp[i], popCount, clusters);

            if(_mode == "disk") // add this fingerprint to clusters where only cfp index and its cluster id remains in memory
                _addFingerprintToClusterDiskMemory(&cfp[i], fpIndex, popCount, clusters, cFpIndexToClusterID);
            }
        }

        chunkCounter += 1;
        remainingSize = remainingSize < chunkSize ? 0 : totalSize - (chunkCounter * chunkSize);

        // sort clusters by its number of fingerprints every 10000 fingerprints
        if (sortCounter % 10000 == 0) {
            std::sort(clusters.begin(), clusters.end(), [](const utils::dt_inner_clusters_fingerprints& a, const utils::dt_inner_clusters_fingerprints& b) {
                return a.num_fps > b.num_fps;
            });

            /*
            int wave_chunk_size = 10;
            int waves = floor((float)clusters.size()/wave_chunk_size);
            if (clusters.size() > 1000) {
                int counter = 0;
                for (int ci = 0; ci<wave_chunk_size; ci++) {
                    for(int fi =0 ; fi < waves; fi++) {
                        //std::cout<< counter<<" " << (fi*wave_chunk_size)+ci << std::endl;
                        std::swap(clusters[counter], clusters[(fi*wave_chunk_size)+ci]);
                        counter += 1;
                    }
                }
            }*/
        }
        sortCounter += num_fps_per_chunk;

    }

    // write the clusters to the file, it assumes that cfp cluster-wise is present in memory
    if (_mode == "memory")
        _writeClusters(popcountGroup, clusters);

    // write the clusters to the file, the cfp from the original dataset is read and written as cluster-wise cfp
    // using cfp index and its cluster id
    if(_mode == "disk")
        _writeClustersDiskMemory(popcountGroup, clusters, cFpIndexToClusterID);

    std::cout << "PopCount: "<< popCount << "; No. of FPs: " << totalSize/_fpstore->_CFPSize << "; No. of clusters: " << clusters.size() << std::endl;
    
    for(auto &cluster: clusters) utils::free_dt_inner_clusters_fingerprints(cluster);
    free(cfp);
    clusters.resize(0);
    cFpIndexToClusterID.resize(0);
}

void InnerClusteringAgent::_doInnerClusteringCUDA(int popCount, h5::Group *popcountGroup) {
    hsize_t totalSize = utils::sizeOfFingerPrintDatasetInH5(popcountGroup, "cfpData");
    unsigned int num_fps_per_chunk = 1000;
    hsize_t chunkSize = _fpstore->_CFPSize * num_fps_per_chunk;
    uint64_t *cfp = (uint64_t *) malloc(sizeof(uint64_t) * chunkSize);
    std::vector<utils::dt_inner_clusters_fingerprints> clusters;
    uint64_t chunkCounter = 0, remainingSize = totalSize, currentSize = 0, fpIndex=0;
    unsigned long sortCounter = 0;

    while (remainingSize > 0) {
        currentSize = remainingSize < chunkSize ? remainingSize : chunkSize;
        utils::getFingerprintFromIndex(popcountGroup, "cfpData", fpIndex, currentSize, cfp);
        for(size_t i = 0; i < currentSize; i += _fpstore->_CFPSize, fpIndex += _fpstore->_CFPSize) {
            if (cfp[i + _fpstore->_CFPPopCountIndex] > 0) {
                int32_t clusterId = -1;

#ifdef USE_CUDA
                const size_t num_clusters = clusters.size();
                if (num_clusters > 0) {
                    std::vector<uint64_t> contiguous_cluster_fp(num_clusters * _fpstore->_CFPSize);
                    for (size_t cid = 0; cid < num_clusters; cid++) {
                        std::memcpy(
                            contiguous_cluster_fp.data() + (cid * _fpstore->_CFPSize),
                            clusters[cid].clusterFp,
                            sizeof(uint64_t) * _fpstore->_CFPSize
                        );
                    }

                    clusterId = find_matching_cluster_cuda(
                        contiguous_cluster_fp.data(),
                        &cfp[i],
                        (int) num_clusters,
                        _fpstore->_CFPSize,
                        _fpstore->_molIdOffset,
                        _fpstore->_fpSize,
                        _fpstore->_CFPPopCountIndex,
                        _fpstore->_innerClusteringThreshold
                    );
                }
#else
                std::cout <<"Warning: CUDA is not enabled. Running inner clustering without CUDA acceleration." << std::endl;
                uint64_t common_popcnt = 0;
                float dist_factor = 0;
                for(size_t ic = 0; ic < clusters.size(); ic++) {
                    uint64_t *clusterFP = clusters[ic].clusterFp;
                    common_popcnt = bitwise_and_popcount(clusterFP + _fpstore->_molIdOffset, &cfp[i] + _fpstore->_molIdOffset, _fpstore->_fpSize);
                    dist_factor = _fpstore->_innerClusteringThreshold *  (&cfp[i])[_fpstore->_CFPPopCountIndex] + _fpstore->_innerClusteringThreshold * (clusterFP[_fpstore->_CFPPopCountIndex] - common_popcnt);
                    if (common_popcnt >= dist_factor) {
                        clusterId = (int32_t) ic;
                        break;
                    }
                }
#endif

                if(clusterId == -1) {
                    uint64_t *clusterFP = (uint64_t*) malloc(sizeof(uint64_t) * _fpstore->_CFPSize);
                    for (size_t j = _fpstore->_molIdOffset; j<_fpstore->_CFPSize; j++)
                        clusterFP[j] = cfp[i + j];

                    utils::dt_inner_clusters_fingerprints cluster;
                    cluster.popCount = popCount;
                    cluster.clusterFp = clusterFP;
                    cluster.num_clusters = 1;
                    cluster.num_fps = 1;
                    cluster.fp = (uint64_t*) malloc(sizeof(uint64_t) * _fpstore->_CFPSize);
                    std::memcpy(cluster.fp, &cfp[i], sizeof(uint64_t) * _fpstore->_CFPSize);
                    clusters.push_back(cluster);
                } else {
                    int fp_end_index = _fpstore->_molIdOffset + _fpstore->_fpSize;
                    uint64_t *clusterFP = clusters[clusterId].clusterFp;
                    clusterFP[_fpstore->_CFPPopCountIndex] = 0;
                    for (int j = _fpstore->_molIdOffset; j<fp_end_index; j++) {
                        clusterFP[j] |= cfp[i + j];
                        clusterFP[_fpstore->_CFPPopCountIndex] += popcntll(clusterFP[j]);
                    }

                    int previous_fps_full_size = clusters[clusterId].num_fps * _fpstore->_CFPSize;
                    int new_fps_full_size = previous_fps_full_size + _fpstore->_CFPSize;
                    clusters[clusterId].num_fps += 1;
                    clusters[clusterId].fp = (uint64_t*) realloc(clusters[clusterId].fp, sizeof(uint64_t) * new_fps_full_size);
                    std::memcpy(clusters[clusterId].fp + previous_fps_full_size, &cfp[i], sizeof(uint64_t) * _fpstore->_CFPSize);
                }
            }
        }

        chunkCounter += 1;
        remainingSize = remainingSize < chunkSize ? 0 : totalSize - (chunkCounter * chunkSize);

        if (sortCounter % 10000 == 0) {
            std::sort(clusters.begin(), clusters.end(), [](const utils::dt_inner_clusters_fingerprints& a, const utils::dt_inner_clusters_fingerprints& b) {
                return a.num_fps > b.num_fps;
            });
        }
        sortCounter += num_fps_per_chunk;
    }

    _writeClusters(popcountGroup, clusters);

    std::cout << "PopCount: "<< popCount << "; No. of FPs: " << totalSize/_fpstore->_CFPSize << "; No. of clusters: " << clusters.size() << std::endl;

    for(auto &cluster: clusters) utils::free_dt_inner_clusters_fingerprints(cluster);
    free(cfp);
    clusters.resize(0);
}

void InnerClusteringAgent::_addFingerprintToCluster(uint64_t *currentFP, int popCount, std::vector<utils::dt_inner_clusters_fingerprints> &clusters) {
    uint64_t common_popcnt = 0, num_clusters = clusters.size();
    float dist_factor = 0; // in place of real distance, we calculate factor, which avoid division operation
    int fp_end_index = _fpstore->_molIdOffset + _fpstore->_fpSize;
    int32_t clusterId = -1;

    if (num_clusters > 1000 && _parallel){ // parallel version when the number of clusters is large
        #pragma omp parallel default(none) shared(clusterId, clusters, num_clusters, currentFP, _fpstore)
        {
            uint64_t common_popcnt = 0;
            float dist_factor = 0; // in place of real distance, we calculate factor, which avoid division operation
            #pragma omp for 
            for(uint64_t ic = 0; ic < num_clusters; ic++) {
                uint64_t *clusterFP = clusters[ic].clusterFp;
                common_popcnt = bitwise_and_popcount(clusterFP + _fpstore->_molIdOffset, currentFP + _fpstore->_molIdOffset, _fpstore->_fpSize);
                dist_factor = _fpstore->_innerClusteringThreshold *  (currentFP[_fpstore->_CFPPopCountIndex] + clusterFP[_fpstore->_CFPPopCountIndex] - common_popcnt);
                if (common_popcnt >= dist_factor) {
                    #pragma omp atomic write
                    clusterId = ic;
                    #pragma omp cancel for
                }
                #pragma omp cancellation point for
            }
        }

    } else {
        for(uint64_t ic = 0; ic < num_clusters; ic++) {
            uint64_t *clusterFP = clusters[ic].clusterFp;
            common_popcnt = bitwise_and_popcount(clusterFP + _fpstore->_molIdOffset, currentFP + _fpstore->_molIdOffset, _fpstore->_fpSize);
            dist_factor = _fpstore->_innerClusteringThreshold *  (currentFP[_fpstore->_CFPPopCountIndex] + clusterFP[_fpstore->_CFPPopCountIndex] - common_popcnt);
            if (common_popcnt >= dist_factor) {
                clusterId = ic;
                break;
            }
        }
    }

    if(clusterId == -1)   { // means this fp is not in any cluster
        uint64_t *clusterFP = (uint64_t*) malloc(sizeof(uint64_t) * _fpstore->_CFPSize);
        for (size_t j = _fpstore->_molIdOffset; j<_fpstore->_CFPSize; j++)
            clusterFP[j] = currentFP[j];

        utils::dt_inner_clusters_fingerprints cluster;
        cluster.popCount = popCount;
        cluster.clusterFp = clusterFP;
        cluster.num_clusters = 1;
        cluster.num_fps = 1;
        cluster.fp = (uint64_t*) malloc(sizeof(uint64_t) * _fpstore->_CFPSize);
        memcpy(cluster.fp, currentFP, sizeof(uint64_t) * _fpstore->_CFPSize);
        clusters.push_back(cluster);
    } else { // append current fp to existing cluster
        uint64_t *clusterFP = clusters[clusterId].clusterFp;
        clusterFP[_fpstore->_CFPPopCountIndex] = 0;
        for (int j = _fpstore->_molIdOffset; j<fp_end_index; j++) {
            clusterFP[j] |= currentFP[j];
            clusterFP[_fpstore->_CFPPopCountIndex] += popcntll(clusterFP[j]);
        }

        // append current fp
        int previous_fps_full_size = clusters[clusterId].num_fps * _fpstore->_CFPSize;
        int new_fps_full_size = previous_fps_full_size + _fpstore->_CFPSize;
        clusters[clusterId].num_fps += 1;
        clusters[clusterId].fp = (uint64_t*) realloc(clusters[clusterId].fp, sizeof(uint64_t) * new_fps_full_size);
        memcpy(clusters[clusterId].fp + previous_fps_full_size, currentFP, sizeof(uint64_t)*_fpstore->_CFPSize);
    }
}

bool InnerClusteringAgent::_writeClusters(h5::Group *popCountGroup, std::vector<utils::dt_inner_clusters_fingerprints> &inner_clusters){

    // in case if there is no clusters, check if previously, there has been clusters and then remove it.
    // it is for extreme edge case when fp modified and popCount bin has no other fingerprints.
    if (inner_clusters.size() == 0 && popCountGroup->nameExists(_fpstore->_clustersGroupName)) {
        popCountGroup->unlink(_fpstore->_clustersGroupName);
        return true;
    }

    h5::Group *clustersGroup = utils::removeAndOpenGroup(popCountGroup, _fpstore->_clustersGroupName);
    std::vector<int> clusterId;
    unsigned long fpsArrayCounter = 0;
    for (size_t cid = 0; cid < inner_clusters.size(); cid++) {
        clusterId.push_back(cid);
        fpsArrayCounter += inner_clusters[cid].num_fps * _fpstore->_CFPSize;
        inner_clusters[cid].clusterFp[0] = fpsArrayCounter; // add in clusterFp, where fp will end in the main fp-array
        utils::appendFingerPrintDatasetToH5(clustersGroup, _fpstore->_clustersFPGroupName, inner_clusters[cid].clusterFp, _fpstore->_CFPSize, _fpstore->_CFPSize);
        utils::appendFingerPrintDatasetToH5(clustersGroup, _fpstore->_fpArrayInClusterGroupName, inner_clusters[cid].fp, inner_clusters[cid].num_fps * _fpstore->_CFPSize, _fpstore->_CFPSize);
    }

    // write cluster-ids vector to file
    utils::addDataSetToH5(clustersGroup, _fpstore->_ClusterIdDataSetName, H5::PredType::NATIVE_INT, clusterId.data(), {clusterId.size()});

    // write innerClusteringThreshold to popcount-bin, used to check whether clustering is required
    utils::set_scaler_attribute(popCountGroup, "innerClusteringThreshold", h5::PredType::NATIVE_FLOAT, _fpstore->_innerClusteringThreshold);

    delete clustersGroup;
    return true;
}

bool InnerClusteringAgent::reDoInnerClusteringInMemory(float threshold) {
    _fpstore->_innerClusteringThreshold = threshold;

    // perform inner clustering with new threshold;
    performInnerClusteringAndWrite();

    // save new threshold
    utils::set_scaler_attribute(_fpstore->_root_group, "innerClusteringThreshold", h5::PredType::NATIVE_FLOAT, _fpstore->_innerClusteringThreshold);
    py::gil_scoped_acquire acquire;
    return true;
}

void InnerClusteringAgent::performInnerClusteringAfterUpdate(std::set<int> modifiedPopCounts) {
    std::cout<< "Performing inner clustering ..." << std::endl;
    uint64_t start_time = utils::get_posix_clock_time();
    h5::Group *popcountBinsGroup =  utils::createOrOpenGroup(_fpstore->_root_group, _fpstore->_popCountBinsGroupName);

    // iterate over modified popcounts and perform inner clustering
    for (auto popcount: modifiedPopCounts) {
        h5::Group *popcountGroup = utils::createOrOpenGroup(popcountBinsGroup, _fpstore->_getPopCountGroupName(popcount));
        if (_mode == "cuda")
            _doInnerClusteringCUDA(popcount, popcountGroup);
        else
            _doInnerClustering(popcount, popcountGroup);
        delete popcountGroup;
    }

    delete popcountBinsGroup;
    std::cout << " ... inner clustering done. Time taken: " << (utils::get_posix_clock_time() - start_time) / (1000 *60) << " mins." << std::endl;
}

void InnerClusteringAgent::_addFingerprintToClusterDiskMemory(uint64_t *currentFP, uint64_t cfpIndex, int popCount, 
                                                              std::vector<utils::dt_inner_clusters_fingerprints> &clusters, 
                                                              std::vector<uint64_t> &cFpIndexToClusterID) {
    uint64_t common_popcnt = 0, num_clusters = clusters.size();
    float dist_factor = 0; // in place of real distance, we calculate factor, which avoid division operation
    int fp_end_index = _fpstore->_molIdOffset + _fpstore->_fpSize;
    uint32_t clusterId = -1;

    if (num_clusters > 1000 && _parallel){ // parallel version when number of clusters is large
        #pragma omp parallel default(none) shared(clusterId, clusters, num_clusters, currentFP, _fpstore)
        {   
            uint64_t common_popcnt = 0;
            float dist_factor = 0;
            #pragma omp for 
            for(uint64_t ic = 0; ic < num_clusters; ic++) {
                uint64_t *clusterFP = clusters[ic].clusterFp;
                common_popcnt = bitwise_and_popcount(clusterFP + _fpstore->_molIdOffset, currentFP + _fpstore->_molIdOffset, _fpstore->_fpSize);
                dist_factor = _fpstore->_innerClusteringThreshold * (currentFP[_fpstore->_CFPPopCountIndex] + clusterFP[_fpstore->_CFPPopCountIndex] - common_popcnt);
                if (common_popcnt >= dist_factor) {
                    #pragma omp atomic write
                    clusterId = ic;
                    #pragma omp cancel for
                }
                #pragma omp cancellation point for
            }
        }
    } else {
        for(uint64_t ic = 0; ic < num_clusters; ic++) {
            uint64_t *clusterFP = clusters[ic].clusterFp;
            common_popcnt = bitwise_and_popcount(clusterFP + _fpstore->_molIdOffset, currentFP + _fpstore->_molIdOffset, _fpstore->_fpSize);
            dist_factor = _fpstore->_innerClusteringThreshold * (currentFP[_fpstore->_CFPPopCountIndex] + clusterFP[_fpstore->_CFPPopCountIndex] - common_popcnt);
            if (common_popcnt >=  dist_factor) {
                clusterId = ic;
                break;
            }
        }
    }

    if(clusterId == -1)   {  // current-fp does not belong to any cluster, new cluster is created
        uint64_t *clusterFP = (uint64_t*) malloc(sizeof(uint64_t) * _fpstore->_CFPSize);
        for (size_t j = _fpstore->_molIdOffset; j<_fpstore->_CFPSize; j++)
            clusterFP[j] = currentFP[j];

        utils::dt_inner_clusters_fingerprints cluster;
        cluster.popCount = popCount;
        cluster.clusterFp = clusterFP;
        cluster.num_clusters = 1;
        cluster.num_fps = 1;
        cFpIndexToClusterID.emplace_back(cfpIndex);
        cFpIndexToClusterID.emplace_back(num_clusters);
        clusters.push_back(cluster);
    } else { // current-fp belongs to an existing cluster
        uint64_t *clusterFP = clusters[clusterId].clusterFp;
        clusterFP[_fpstore->_CFPPopCountIndex] = 0;
        for (int j = _fpstore->_molIdOffset; j<fp_end_index; j++) {
            clusterFP[j] |= currentFP[j];
            clusterFP[_fpstore->_CFPPopCountIndex] += popcntll(clusterFP[j]);
        }

        // append current fp
        clusters[clusterId].num_fps += 1;
        cFpIndexToClusterID.emplace_back(cfpIndex);
        cFpIndexToClusterID.emplace_back(clusterId);
    }
}

bool InnerClusteringAgent::_writeClustersDiskMemory(h5::Group *popCountGroup, std::vector<utils::dt_inner_clusters_fingerprints> &inner_clusters, std::vector<uint64_t> &cFpIndexToClusterID){

    // in case if there is no clusters, check if previously, there has been clusters and then remove it.
    // it is for extreme edge case when fp modified and popCount bin has no other fingerprints.
    if (inner_clusters.size() == 0 && popCountGroup->nameExists(_fpstore->_clustersGroupName)) {
        popCountGroup->unlink(_fpstore->_clustersGroupName);
        return true;
    }

    h5::Group *clustersGroup = utils::removeAndOpenGroup(popCountGroup, _fpstore->_clustersGroupName);
    std::vector<int> clusterId;
    std::vector<uint64_t> fpRunningIndices; // it will be used to add fp in fpArrayInCluster later, will be used to add next cluster's fp in fpArrayInCluster
    fpRunningIndices.push_back(0); // index of first fingerprint of first cluster that will be always zero
    unsigned long fpsArrayCounter = 0;
    for (size_t cid = 0; cid < inner_clusters.size(); cid++) { // first write cluster-fp data in the file
        clusterId.push_back(cid);
        fpsArrayCounter += inner_clusters[cid].num_fps * _fpstore->_CFPSize;
        fpRunningIndices.push_back(fpsArrayCounter); // initialization, starting index of first fingerprint of next cluster
        inner_clusters[cid].clusterFp[0] = fpsArrayCounter; // add in clusterFp, where fp will end in the main fp-array
        utils::appendFingerPrintDatasetToH5(clustersGroup, _fpstore->_clustersFPGroupName, inner_clusters[cid].clusterFp, _fpstore->_CFPSize, _fpstore->_CFPSize);
    }

    utils::createEmptyFpDataset(clustersGroup, _fpstore->_fpArrayInClusterGroupName, fpsArrayCounter, _fpstore->_CFPSize); // create dataset in file with zeroes values
    uint64_t *tCFpData = (uint64_t *) malloc(sizeof(uint64_t) * _fpstore->_CFPSize);
    for(uint64_t i=0; i < cFpIndexToClusterID.size(); i += 2) { // write data to fpArrayInCluster
        utils::getFingerprintFromIndex(popCountGroup, "cfpData", cFpIndexToClusterID[i], _fpstore->_CFPSize, tCFpData); // get data
        utils::updateFingerprintAtIndex(clustersGroup, _fpstore->_fpArrayInClusterGroupName, fpRunningIndices[cFpIndexToClusterID[i+1]], _fpstore->_CFPSize, tCFpData); //copy data
        fpRunningIndices[cFpIndexToClusterID[i+1]] += _fpstore->_CFPSize;
    }
    free(tCFpData);

    // write cluster-ids vector to file
    utils::addDataSetToH5(clustersGroup, _fpstore->_ClusterIdDataSetName, H5::PredType::NATIVE_INT, clusterId.data(), {clusterId.size()});

    // write innerClusteringThreshold to popcount-bin, used to check whether clustering is required
    utils::set_scaler_attribute(popCountGroup, "innerClusteringThreshold", h5::PredType::NATIVE_FLOAT, _fpstore->_innerClusteringThreshold);

    delete clustersGroup;
    return true;
}