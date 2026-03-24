/**
 * @brief This class represents an inner clustering agent which performs inner clustering within the popcount bins.
 * @author Rajendra Kumar
 */

#include <string>
#include <map>
#include <tuple>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <H5Cpp.h>
#include <pybind11/pybind11.h>

#include "lib.h"

namespace py = pybind11;
namespace h5 = H5;
using namespace py::literals;

class FingerprintStore;


/**
 * @brief The InnerClusteringAgent class is responsible for performing inner clustering on molecular fingerprints and writing the results to an HDF5 file.
 *
 * The class provides methods for adding molecular fingerprints to clusters, performing inner clustering, and writing
 * the results to an HDF5 file.
 * It supports both memory-based and disk-based operations, allowing for efficient processing of large datasets.
 *
 * @param fpstore A pointer to the FingerprintStore object, which provides access to molecular fingerprint data.
 * @param mode A string indicating the mode of operation, either "memory" for memory-based operations or "disk" for disk-based operations.
 */
class InnerClusteringAgent {

    private:
        FingerprintStore *_fpstore = nullptr;
        std::string _mode = "memory";

        /**
         * @brief Adds a molecular fingerprint to a cluster.
         *
         * This method takes a molecular fingerprint represented by a 64-bit integer array, its population count, and a vector of inner clustering results.
         * It then searches for an existing cluster with a similar fingerprint or creates a new cluster if no suitable cluster is found.
         * The molecular fingerprint is added to the appropriate cluster based on its similarity to the existing fingerprints in the cluster.
         *
         * @param currentFP A pointer to the 64-bit integer array representing the molecular fingerprint.
         * @param popCount The population count of the molecular fingerprint.
         * @param clusters A vector of inner clustering results where the molecular fingerprint will be added to a cluster.
         */        
        void _addFingerprintToCluster(uint64_t *currentFP, int popCount, std::vector<utils::dt_inner_clusters_fingerprints> &clusters);

        /**
         * @brief Adds a molecular fingerprint to a cluster in disk memory mode.
         *
         * This method takes a molecular fingerprint represented by a 64-bit integer array, its compact fingerprint index, its population count, and a vector of inner clustering results.
         * It then searches for an existing cluster with a similar fingerprint or creates a new cluster if no suitable cluster is found.
         * The molecular fingerprint is added to the appropriate cluster based on its similarity to the existing fingerprints in the cluster.
         *
         * @param currentFP A pointer to the 64-bit integer array representing the molecular fingerprint.
         * @param cfpIndex The compact fingerprint index of the molecular fingerprint.
         * @param popCount The population count of the molecular fingerprint.
         * @param clusters A vector of inner clustering results where the molecular fingerprint will be added to a cluster.
         * @param cFpIndexToClusterID A vector that maps compact fingerprint indices to cluster IDs.
         */
        void _addFingerprintToClusterDiskMemory(uint64_t *currentFP, uint64_t cfpIndex, int popCount, std::vector<utils::dt_inner_clusters_fingerprints> &clusters, std::vector<uint64_t> &cFpIndexToClusterID);

        /**
         * @brief Performs the inner clustering process for a specific population count.
         *
         * This method takes a population count and a pointer to an HDF5 group representing the population count group.
         * It then iterates through the molecular fingerprints with the specified population count, performs the inner clustering algorithm, 
         * and updates the inner clustering results.
         * The inner clustering results are written in the HDF5 file.
         *
         * @param popCount The population count for which the inner clustering process will be performed.
         * @param popcountGroup A pointer to the HDF5 group representing the population count group.
         */        
        void _doInnerClustering(int popCount, h5::Group *popcountGroup);

        /**
         * @brief Performs CUDA-accelerated inner clustering process for a specific population count.
         *
         * This method is intended to run inner clustering with CUDA acceleration when available.
         *
         * @param popCount The population count for which the inner clustering process will be performed.
         * @param popcountGroup A pointer to the HDF5 group representing the population count group.
         */
        void _doInnerClusteringCUDA(int popCount, h5::Group *popcountGroup);

        /**
         * @brief Writes the inner clustering results to an HDF5 file for a specific population count.
         *
         * This method takes a pointer to an HDF5 group representing the population count group and a vector of inner clustering results.
         * It then writes the inner clustering results to the HDF5 file, creating the necessary datasets and attributes.
         *
         * @param popCountGroup A pointer to the HDF5 group representing the population count group.
         * @param inner_clusters A vector of inner clustering results for the specific population count.
         * @return A boolean value indicating whether the writing process was successful.
         */
        bool _writeClusters(h5::Group *popCountGroup, std::vector<utils::dt_inner_clusters_fingerprints> &inner_clusters);

        /**
         * @brief Writes the inner clustering results to an HDF5 file for a specific population count in disk memory mode.
         *
         * This method takes a pointer to an HDF5 group representing the population count group, a vector of inner clustering results,
         * and a vector of compact fingerprint indices to cluster IDs.
         * It then writes the inner clustering results to the HDF5 file, creating the necessary datasets and attributes.
         *
         * @param popCountGroup A pointer to the HDF5 group representing the population count group.
         * @param inner_clusters A vector of inner clustering results for the specific population count.
         * @param cFpIndexToClusterID A vector that maps compact fingerprint indices to cluster IDs.
         * @return A boolean value indicating whether the writing process was successful.
         */
        bool _writeClustersDiskMemory(h5::Group *popCountGroup, std::vector<utils::dt_inner_clusters_fingerprints> &inner_clusters, std::vector<uint64_t> &cFpIndexToClusterID);
        

    public:
        bool _parallel;
        
        /**
         * @brief The InnerClusteringAgent class is responsible for performing inner clustering on molecular fingerprints and writing the results to an HDF5 file.
         *
         * The class provides methods for adding molecular fingerprints to clusters, performing inner clustering, and writing
         * the results to an HDF5 file.
         * It supports both memory-based and disk-based operations, allowing for efficient processing of large datasets.
         *
         * @param fpstore A pointer to the FingerprintStore object, which provides access to molecular fingerprint data.
         * @param mode A string indicating the mode of operation, either "memory" for memory-based operations or "disk" for disk-based operations.
         */
        InnerClusteringAgent(FingerprintStore *fpstore, std::string mode = "memory", bool parallel = false);

        /**
         * @brief Performs inner clustering on molecular fingerprints and writes the results to an HDF5 file.
         *
         * The method iterates through the molecular fingerprints, performs the inner clustering algorithm, and updates the inner clustering results.
         * The results are then written to an HDF5 file, creating the necessary datasets and attributes.
         *
         * @return A boolean value indicating whether the writing process was successful.
         */
        bool performInnerClusteringAndWrite();

        /**
         * @brief Re-performs inner clustering on molecular fingerprints in memory-based mode and writes the results to an HDF5 file.
         *
         * The method iterates through the molecular fingerprints, performs the inner clustering algorithm, and updates the inner clustering results.
         * The results are then written to an HDF5 file, creating the necessary datasets and attributes.
         *
         * @param threshold A float value representing the threshold for re-performing the inner clustering process.
         * @return A boolean value indicating whether the writing process was successful.
         */
        bool reDoInnerClusteringInMemory(float threshold);

        /**
         * @brief Performs the inner clustering process on molecular fingerprints after updating the results based on modified population counts.
         *
         * The method iterates through the molecular fingerprints with modified population counts, performs the inner clustering algorithm, and updates the inner clustering results.
         * The updated results are then used to perform further processing or analysis.
         *
         * @param modifiedPopCounts A set of integers representing the population counts that have been modified.
         */
        void performInnerClusteringAfterUpdate(std::set<int> modifiedPopCounts);
};