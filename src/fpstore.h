/**
 * @brief This class represents the FingerprintStore class that handles the storage and retrieval of molecular fingerprints.
 * @author Rajendra Kumar
 */

#include <string>
#include <map>
#include <tuple>

#include <H5Cpp.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "molDataTable.h"
#include "lib.h"
#include "popcnt.hpp"

namespace py = pybind11;
namespace h5 = H5;
using namespace py::literals;
class InnerClusteringAgent;

/**
 * @brief A class for storing and managing molecular fingerprint data in an HDF5 file.
 *
 * The FingerprintStore class provides methods for reading, writing, and manipulating molecular fingerprint data in an HDF5 file.
 * It also includes various properties and helper functions to facilitate data access and manipulation.
 *
 */
class FingerprintStore
{
    // constants
    const std::string MagicFingerPrintNumberValue = "uffpsim-search-892167801775531034751912048751";
    const std::string MagicFingerPrintNumberKey = "magic-number";

private:
    // Private member variables
    H5std_string _filename;                      // HDF5 file name
    h5::H5File *_file = nullptr;                 // HDF5 file object
    std::string _file_mode = "r";                // HDF5 file mode (read-only `r`, write `w`, or append `a`)
    std::string _fp_params;                      // Fingerprint parameters in JSON format as string 
    MolDataTable *_molDataTable;                 // MolDataTable object for handling molecular data in the HDF5 file
    InnerClusteringAgent *_innerClusteringAgent; // InnerClusteringAgent object for performing inner clustering operations
    h5::Group *_popCountBinsGroup = nullptr;       // HDF5 group for population count bins
    std::map<int, h5::Group*> _popCountToGroupMap; // map from population count to corresponding HDF5 group for that population count
    std::map<int, h5::Group*> _popCountToClustersGroupMap; // map from population count to corresponding HDF5 group for clusters with that population count

    /**
     * @brief Populates the internal data structures with the molecular fingerprints and other relevant information from the HDF5 file.
     *
     * This method reads the necessary datasets from the HDF5 file and stores them in memory for efficient access.
     * It also initializes the internal variables and structures required for further operations.
     *
     * @return True if the data population was successful, false otherwise.
     */
    bool _populateDataInMemory(bool onlyClusterFPs = false);

    /**
     * @brief Populates the internal data structures with molecular fingerprints for a specific population count.
     *
     * This method reads the molecular fingerprints for molecules with the given population count from the HDF5 file and stores them in memory for efficient access.
     * It also initializes the internal variables and structures required for further operations related to clustering.
     *
     * @param popCount The population count for which molecular fingerprints need to be populated.
     * @param inner_clusters_fingerprints A pointer to the data structure where the molecular fingerprints will be stored.
     * @param onlyClusterFPs A boolean flag indicating whether to populate only the fingerprint data for clustered molecules.
     */
    void _populateClustersInMemory(int popCount, utils::dt_inner_clusters_fingerprints *inner_clusters_fingerprints, bool onlyClusterFPs = false);

    /**
     * @brief Builds a sorted vector of molecular fingerprint data based on the popcount from the given fingerprints.
     *
     * This method takes a vector of molecular fingerprints in the form of tuples (ID, SMILES, Fingerprint) and sorts them based on the popcount.
     * The sorted molecular fingerprint data is then stored in a vector of type utils::dt_mol_fp_data.
     *
     * @param fingerprints A vector of tuples containing molecular fingerprint information (ID, SMILES, Fingerprint).
     * @return A sorted vector of type utils::dt_mol_fp_data containing molecular fingerprint data.
     */
    std::vector<utils::dt_mol_fp_data> _buildSortedMolFpData(const std::vector<std::tuple<std::string, std::string, std::string>> fingerprints);
    
    /**
     * @brief Saves the sorted molecular fingerprint data to the HDF5 file.
     *
     * This method takes a sorted vector of molecular fingerprint data (type utils::dt_mol_fp_data) and saves it to the HDF5 file.
     * The molecular fingerprint data is written to the appropriate datasets within the HDF5 file.
     *
     * @param sortedMolFpData A reference to a sorted vector of type utils::dt_mol_fp_data containing molecular fingerprint data.
     */
    void _saveFpMolData(std::vector<utils::dt_mol_fp_data> &sortedMolFpData);

    /**
     * @brief Saves the modified molecular fingerprint data to the HDF5 file.
     *
     * This method takes a sorted vector of molecular fingerprint data (type utils::dt_mol_fp_data) and a corresponding vector of existing molecular data (type TMolData).
     * It compares the molecular fingerprints in the sortedMolFpData vector with the existing molecular data in the sortedExistingMolsData vector, and saves the modified molecular fingerprint data to the HDF5 file.
     * The modified molecular fingerprint data is written to the appropriate datasets within the HDF5 file.
     *
     * @param sortedMolFpData A reference to a sorted vector of type utils::dt_mol_fp_data containing molecular fingerprint data.
     * @param sortedExistingMolsData A reference to a vector of type TMolData containing existing molecular data.
     */
    void _saveModifiedFpMolData(std::vector<utils::dt_mol_fp_data> &sortedMolFpData, std::vector<TMolData> sortedExistingMolsData);

    /**
     * @brief Writes the root attributes to the HDF5 file.
     *
     * This method writes the root attributes (such as metadata or configuration information) to the HDF5 file.
     * The root attributes are stored in the root group of the HDF5 file.
     * The method ensures that the root attributes are written correctly and efficiently.
     *
     * @return True if the root attributes are written successfully, false otherwise.
     */
    bool _writeRootAttributes();

    /**
     * @brief Populates the root attributes from the HDF5 file.
     *
     * This method reads the root attributes (such as metadata or configuration information) from the HDF5 file.
     * The root attributes are stored in the root group of the HDF5 file.
     * The method ensures that the root attributes are read correctly and efficiently.
     *
     * @return True if the root attributes are populated successfully, false otherwise.
     */
    bool _populateRootAttributes();

public:
    // Public member variables
    const std::string _ClusterIdDataSetName = "clusterIDs";
    const std::string _popCountDatasetName = "popCountArray";
    const std::string _popCountBinsGroupName = "popCountBins";
    const std::string _clustersGroupName = "clusters";
    const std::string _clustersFPGroupName = "clusterFP";
    const std::string _fpArrayInClusterGroupName = "fpArrayInCluster";
    h5::Group *_root_group = nullptr;                                      // HDF5 root group
    std::vector<int> _popCountBins;                                        // Population count bins for molecular fingerprints
    int _molIdMaxLength;                                                   // Maximum length of molecular IDs
    size_t _fpSize = 0;                                                    // Size of fingerprints in 64-bit integer size, real size is (_fpSize*64)
    size_t _CFPSize;                                                       // Size of compact fingerprints in 64-bit integer size (_molIdOffset + _fpSize + 1)
    size_t _molIdOffset;                                                   // Offset for molecular IDs in CF at its beginning
    size_t _CFPPopCountIndex;                                              // Index for compact fingerprint population count in _CFPSize - last element of CFp
    float _innerClusteringThreshold = 0.2;                                 // Inner clustering threshold for clustering molecular fingerprints
    std::string _info = "";                                                // Additional information about the fingerprint store, it could be JSON format string
    utils::dt_inner_clusters_fingerprints *_fp_inner_clusters_by_popcount; // Data structure for storing inner clustering results by population count



    // Constructor and destructor
    /**
     * @brief Constructs a new instance of the FingerprintStore class.
     *
     * This constructor initializes the FingerprintStore object with the provided parameters and sets up the HDF5 file for reading or writing.
     * It also performs initial setup tasks such as setting the maximum length of molecular IDs, fingerprint size, and other configurations.
     *
     * @param filename The name of the HDF5 file to be used for storing molecular fingerprint data.
     * @param molIdMaxLength The maximum length of molecular IDs. Default value is 15.
     * @param mode The mode in which the HDF5 file should be opened. It can be "r" for read-only, "w" for write, or "a" for append. Default value is "r".
     * @param fpSize The size of molecular fingerprints in bits. Default value is 0.
     * @param fp_params The fingerprint parameters in JSON format as a string. Default value is an empty string.
     * @param info Additional information about the fingerprint store, which can be a JSON format string. Default value is an empty string.
     * @param clusterThreshold The inner clustering threshold for clustering molecular fingerprints. Default value is 0.2.
     * @param clusterMode The mode for performing inner clustering. It can be "memory" for in-memory clustering or "file" for file-based clustering. Default value is "memory".
     */
    FingerprintStore(const std::string &filename, int molIdMaxLength = 15, std::string mode = "r", 
                     int fpSize = 0, std::string fp_params = "", std::string info = "", 
                     float clusterThreshold = 0.2, std::string clusterMode = "memory", bool clusterParallel = false);
                     
    ~FingerprintStore();

    /**
     * @brief A static member function that generates a group name for molecular fingerprints based on the population count.
     *
     * The _getPopCountGroupName function takes an integer representing the population count as input and returns a string representing the group name.
     *
     * @param popCount An integer representing the population count.
     * @return A string representing the group name for molecular fingerprints based on the population count.
     */
    static std::string _getPopCountGroupName(int popCount) { return "popCountBin_" + std::to_string(popCount); }
    
    /**
     * @brief Closes the HDF5 file associated with the FingerprintStore object.
     *
     * This method closes the HDF5 file, freeing up any system resources that were being used by the file.
     * After calling this method, any further operations on the HDF5 file will result in an error.
     *
     * It is important to call this method when you are done with the FingerprintStore object to ensure proper cleanup and resource management.
     */
    void close();

    // VARIOUS PROPERTIES

    /**
     * @brief Retrieves the maximum length of molecular IDs from the FingerprintStore object.
     *
     * This method returns the maximum length of molecular IDs that were set during the initialization of the FingerprintStore object.
     * The maximum length of molecular IDs is used to ensure consistency and proper handling of molecular IDs in the HDF5 file.
     *
     * @return The maximum length of molecular IDs as an integer.
     */
    int getMolIdMaxLength() { return _molIdMaxLength; };

    /**
     * @brief Retrieves the size of molecular fingerprints in bits from the FingerprintStore object.
     *
     * This method returns the size of molecular fingerprints in bits that were set during the initialization of the FingerprintStore object.
     * The size of molecular fingerprints is used to ensure consistency and proper handling of molecular fingerprints in the HDF5 file.
     *
     * @return The size of molecular fingerprints in bits as an integer.
     */
    int getFPBitsSize() { return _fpSize * 64; };

    /**
     * @brief Retrieves the inner clustering threshold from the FingerprintStore object.
     *
     * This method returns the inner clustering threshold that was set during the initialization of the FingerprintStore object.
     *
     * @return The inner clustering threshold as a float.
     */
    float getInnerClusteringThreshold() { return _innerClusteringThreshold; };

    /**
     * @brief Retrieves the fingerprint parameters from the FingerprintStore object.
     *
     * This method returns the fingerprint parameters that were set during the initialization of the FingerprintStore object.
     * The fingerprint parameters are typically provided in JSON format and contain information about the fingerprinting algorithm and its parameters.
     *
     * @return The fingerprint parameters as a string in JSON format.
     */
    std::string getFingerprintParameters() { return _fp_params; };

    /**
     * @brief Retrieves the additional information from the FingerprintStore object.
     *
     * This method returns the additional information that was set during the initialization of the FingerprintStore object.
     * The additional information can be used to store any additional details about the data stored in the HDF5 file.
     *
     * @return The additional information as a string.
     */
    std::string getInfo() { return _info; };

    /**
     * @brief Retrieves the total number of molecules stored in the HDF5 file.
     *
     * This method returns the count of molecules in the HDF5 file, which is obtained from the molecular data table.
     *
     * @return The total number of molecules as a 64-bit unsigned integer.
     */
    uint64_t getNumberOfMolecules() { return (uint64_t)_molDataTable->_nrecords; };

   /**
     * @brief Retrieves the population count bins from the HDF5 file.
     *
     * This method returns the population count bins that were set during the initialization of the FingerprintStore object.
     * The population count bins are used to group molecular fingerprints based on their popcount (population count).
     *
     * @return A vector of integers representing the population count bins.
     */ 
    std::vector<int> getPopCountBins() { return _popCountBins; };

    /**
     * @brief Checks if the magic number exists in the HDF5 file associated with the FingerprintStore object.
     *
     * This method reads the root attributes of the HDF5 file and checks if the magic number exists.
     * The magic number is a unique identifier used to verify the integrity of the data stored in the file.
     *
     * @return True if the magic number exists, false otherwise.
     */
    bool magicNumberExists();

    /**
     * @brief Sets the magic number in the HDF5 file associated with the FingerprintStore object.
     *
     * This method writes the magic number to the root attributes of the HDF5 file.
     * The magic number is a unique identifier used to verify the integrity of the data stored in the file.
     *
     * @return True if the magic number is successfully set, false otherwise.
     */
    bool setMagicNumber();

    /**
     * @brief Sets the additional information in the HDF5 file associated with the FingerprintStore object.
     *
     * This method writes the additional information to the root attributes of the HDF5 file.
     * The additional information can be used to store any additional details about the data stored in the file.
     *
     * @param info The additional information to be set.
     * @return True if the additional information is successfully set, false otherwise.
     */
    bool setInfo(std::string info);

    /**
     * @brief Loads the molecular fingerprint data into memory for efficient access.
     *
     * This method reads the necessary molecular fingerprint data from the HDF5 file and stores it in memory for further operations.
     * It is only applicable when the HDF5 file is opened in read-only mode.
     * 
     * It loads all the cluster-fingerprints and fingerprints into memory for full-memory search.
     *
     * @return None.
     */
    void loadDataInMemory(bool onlyClusterFPs = false)
    {
        if (_file_mode == "r")
            _populateDataInMemory(onlyClusterFPs);
    }

    /**
     * @brief Frees the memory allocated for storing molecular fingerprint data.
     *
     * This method releases the memory allocated for storing molecular fingerprint data in memory.
     * It is important to call this method when the molecular fingerprint data is no longer needed to free up system resources.
     *
     * @return None.
     */    
    void freeMemory();

    /**
     * @brief Appends new molecular fingerprint data to the HDF5 file.
     *
     * This method takes a vector of molecular fingerprint data in the form of tuples (ID, SMILES, Fingerprint) and appends 
     * them to the HDF5 file.
     * The new molecular fingerprint data is written to the appropriate datasets within the HDF5 file.
     *
     * @param fingerprints A vector of tuples containing molecular fingerprint information (ID, SMILES, Fingerprint).
     * @return None.
     */    
    void AppendFingerprints(const std::vector<std::tuple<std::string, std::string, std::string>> fingerprints);

    /**
     * @brief Updates the molecular fingerprint data in the HDF5 file.
     *
     * This method takes a vector of molecular fingerprint data in the form of tuples (ID, SMILES, Fingerprint) as input and
     * updates the corresponding datasets within the HDF5 file.
     * The method compares the existing molecular fingerprints with the new fingerprints and updates the HDF5 file accordingly.
     *
     * @param fingerprints A vector of tuples containing molecular fingerprint information (ID, SMILES, Fingerprint).
     * @return None.
     */
    void updateFingerprints(const std::vector<std::tuple<std::string, std::string, std::string>> fingerprints);

    /**
     * @brief Performs inner clustering on the molecular fingerprint data and writes the results to the HDF5 file.
     *
     * This method performs inner clustering on the molecular fingerprint data using the specified inner clustering threshold.
     * The inner clustering algorithm groups similar molecular fingerprints together based on their popcount (population count).
     * The results of the inner clustering, such as cluster IDs and compact fingerprints, are written to the appropriate datasets 
     * within the HDF5 file.
     *
     * The method ensures that the inner clustering results are accurate and efficiently stored in the HDF5 file.
     *
     * @return None.
     */    
    void performInnerClusteringAndWrite();

    /**
     * @brief Re-Performs inner clustering on the molecular fingerprint data using the specified threshold.
     *
     *
     * @param threshold The inner clustering threshold for clustering molecular fingerprints.
     * @return None.
     */    
    void reDoInnerClusteringInMemory(float threshold);

    /**
     * @brief Builds a mapping from molecular IDs to their corresponding indices in the molecular data table.
     *
     * This method creates an internal map that associates each molecular ID with its index in the molecular data table.
     * The mapping is used for efficient lookup and retrieval of molecular data based on molecular IDs.
     *
     * The method relies on the `buildMolIdToIndexMap` function of the `MolDataTable` class to perform the mapping.
     *
     * @return None.
     */
    void buildMolIdToIndexMap() { _molDataTable->buildMolIdToIndexMap(); }

    /**
     * @brief Retrieves the SMILES string associated with a given molecular ID from the HDF5 file.
     *
     * This method takes a molecular ID as input and searches for the corresponding molecular data in the HDF5 file.
     * It then retrieves the SMILES string associated with the given molecular ID.
     *
     * @param id The molecular ID for which the SMILES string needs to be retrieved.
     * @return The SMILES string associated with the given molecular ID.
     * @throws std::runtime_error If the molecular ID is not found in the HDF5 file.
     */
    std::string getSmilesFromID(std::string id);

    uint64_t* getFPsForCluster(int popCount, int fpStartIndex, int fpEndIndex);

    void initH5GroupsMappingForPopCountBins();
};