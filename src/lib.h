/**
 * @brief These are utility functions used in this package.
 * @author Rajendra Kumar
 */

#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include <H5Cpp.h>

namespace h5 = H5;

namespace utils {

    /**
     * This struct is used to store the cluster fingerprints and the fingerprints.
     * It can be used to store all the clusters or a single cluster.
     * 
     * Members:
     * - popCount: An integer representing the population count.
     * - num_clusters: An unsigned long integer representing the number of clusters.
     * - clusterFp: A pointer to an array of uint64_t representing the cluster fingerprints.
     * - num_fps: An unsigned long integer representing the number of fingerprints.
     * - fp: A pointer to an array of uint64_t representing the fingerprints.
     * 
     * Note: Ensure that memory for clusterFp and fp are allocated using malloc/realloc.
     */
    typedef struct {
        int popCount;
        unsigned long num_clusters = 0;
        uint64_t *clusterFp = nullptr;
        unsigned long num_fps = 0;
        uint64_t *fp = nullptr;
    } dt_inner_clusters_fingerprints;

    /**
     * This struct is used to store the cluster fingerprints and the fingerprints along with a score.
     * It can be used to store all the clusters or a single cluster with their corresponding scores.
     * 
     * Members:
     * - inner_clusters_fingerprints: A dt_inner_clusters_fingerprints struct representing the cluster fingerprints and fingerprints.
     * - score: A float representing the score associated with the cluster fingerprints and fingerprints.
     */
    typedef struct {
        utils::dt_inner_clusters_fingerprints inner_clusters_fingerprints;
        float score;
    } dt_inner_clusters_fingerprints_maxscore;

    /**
     * This struct represents a molecule's fingerprint data.
     * It contains the population count, a pointer to the fingerprint array, the molecule's ID, and the molecule's SMILES string.
     *
     * Members:
     * - popCount: An integer representing the population count of the fingerprint.
     * - fp: A pointer to an array of uint64_t representing the fingerprint.
     * - mol_id: A character array of size 64 to store the molecule's ID.
     * - smiles: A character array of size 1024 to store the molecule's SMILES string. The default value is an empty string.
     */
    typedef struct {
        int popCount;
        uint64_t *fp;
        char mol_id[64];
        char smiles[1024] = "";
    } dt_mol_fp_data;


    /*
     * This struct is used to store a single result for a query.

     * Members:
     * - id: A pointer to a std::string representing the id of the hit.
     * - score: A float representing the score of the hit.
     */
    typedef struct {
       std::string *id;
       float score;
    } dt_result;

    /*
     * This struct is used to store query data for batch search.

     * Members:
     * - cfp: A pointer to a uint64_t representing the query fingerprint.
     * - done: A boolean representing whether the search is completed.
     * - max_coeff: A float representing the maximum similarity score of the query.
     * - results_size: An integer representing the number of results.
     * - results: A pointer to an array of dt_result representing the results.
     */
    typedef struct {
        uint64_t *cfp = NULL;
        bool done = false;
        float max_coeff = 0;
        int results_size = 0;
        dt_result *results = NULL;
    } dt_batch_query_data;

    /**
     * This struct is used to store queries data for its popcount-bin for batch search.
     * 
     * Members:
     * - qdata: A pointer to an array of dt_batch_query_data representing the array of queries data.
     * - qsize: An integer representing the number of queries.
     * - popCount: An unsigned long integer representing the population count of the queries.
     * - filteredPopCountBinsWithMaxScore: A std::vector of dt_inner_clusters_fingerprints_maxscore
     */
    typedef struct {
        dt_batch_query_data *qdata = NULL;
        int qsize = 0;
        uint64_t popCount = 0;
        std::vector<utils::dt_inner_clusters_fingerprints_maxscore> filteredPopCountBinsWithMaxScore;
    } dt_batch_data;

    /***********************  utils.cpp ******************/

    /**
     * This function frees the memory allocated for the cluster fingerprints and fingerprints in the given
     * dt_inner_clusters_fingerprints struct.
     * 
     * @param inner_clusters_fingerprints: A reference to the dt_inner_clusters_fingerprints struct.
     * 
     * Note: This function should be called to free the memory allocated for clusterFp and fp before the
     * struct goes out of scope.
     */
    void free_dt_inner_clusters_fingerprints(dt_inner_clusters_fingerprints &inner_clusters_fingerprints);

    /**
     * This function generates a random string of the specified length.
     *
     * @param length: An unsigned integer representing the desired length of the random string.
     *
     * @return: A std::string containing the generated random string.
     *
     * Note: The generated random string will consist of alphanumeric characters (both uppercase and lowercase).
     *       The function uses the current system time as a seed for random number generation.
     */
    std::string randomString(size_t length);

    /**
     * This function retrieves the current time in seconds using the POSIX clock.
     *
     * @return: A double value representing the current time in seconds since the POSIX epoch (00:00:00 UTC, January 1, 1970).
     *
     * Note: The function uses the clock_gettime() function with the CLOCK_REALTIME clock to retrieve the current time.
     *       The returned time is accurate to the nearest microsecond.
     */
    double get_posix_clock_time();

    /**
     * This function calculates the bitwise AND popcount of two vectors of uint64_t.
     *
     * @param a: A vector of uint64_t representing the first input for the bitwise AND operation.
     * @param b: A vector of uint64_t representing the second input for the bitwise AND operation.
     *
     * @return: An integer representing the popcount (number of set bits) of the result of the bitwise AND operation on 'a' and 'b'.
     *
     * Note: The function assumes that both input vectors 'a' and 'b' have the same size.
     *       The bitwise AND operation is performed on each corresponding pair of elements from 'a' and 'b'.
     *       The popcount of the result is then returned.
     */
    int bitwiseAndPopcount(std::vector<uint64_t> a, std::vector<uint64_t> b);

    /**
     * Converts a bit-string representation of a molecule's fingerprint into a compact fingerprint array with its id and pop-count.
     *
     * @param bit_string: A const reference to the bit string representation of the molecule's fingerprint.
     * @param mol_id: A std::string containing the molecule's ID.
     * @param cfp_arr: A pointer to the array where the compact fingerprint will be stored.
     * @param id_offset: An integer representing the offset for the molecule's ID i.e. maximum number of characters in id.
     * @param total_size: An integer representing the total size of the compact fingerprint array.
     *
     * The function converts the bit string representation into a compact fingerprint array by performing
     * the following steps:
     * 1. Converts the molecule's ID into a compact representation and stores it in the compact fingerprint array.
     * 2. Converts the molecule's fingerprint into a compact representation and stores it in the compact fingerprint array.
     *
     * The compact fingerprint array is stored in the memory pointed to by 'cfp_arr'. The molecule's ID is stored
     * starting from zeroth postion to the 'id_offset' position in the array, and the fingerprint is stored after the ID.
     * Last element of the compact fingerprint array is its population count.
     *
     * Note: The caller is responsible for ensuring that the 'cfp_arr' array has enough space to store the compact fingerprint.
     */
    void BitStrToCompactFPArray(const std::string &bit_string, std::string mol_id, uint64_t *cfp_arr, int id_offset, int total_size);

    /**
     * Retrieves the molecule's ID from a compact fingerprint array.
     *
     * @param cfp_arr: A pointer to the compact fingerprint array.
     * @param max_mol_id_size: An integer representing the maximum size of the molecule's ID.
     *
     * @return: A std::string containing the molecule's ID.
     *
     * Note: The molecule's ID is stored in the compact fingerprint array starting from the zeroth position.
     *       The function extracts the molecule's ID from the compact fingerprint array and returns it as a string.
     */
    std::string getMolIdFromCompactFPArray(uint64_t *cfp_arr, int max_mol_id_size);
    
    /**
     * This function converts a bit string representation of a molecule's fingerprint into a compact fingerprint array.
     *
     * @param bit_string: A reference to the bit string representation of the molecule's fingerprint.
     *
     * @return: A vector of uint64_t containing the compact fingerprint array.
     *
     * Note: The compact fingerprint array is a more efficient representation of the fingerprint, suitable for further processing.
     *       The function assumes that the bit string representation is valid and follows a specific format.
     */
    std::vector<uint64_t> getCompactFingerPrintArray(std::string &bit_string);


    /***********************  h5_utils.cpp ******************/

    /**
     * This function creates a new group with the given name within the parent group, or opens an existing group
     * with the given name if it already exists.
     * 
     * @param parentGroup: A pointer to the parent group within which the new group or existing group should be created/opened.
     * @param name: A const reference to the name of the new group or existing group.
     * 
     * @return: A pointer to the newly created or opened group.
     * 
     * Note: The caller is responsible for closing the returned group using the appropriate H5::Group::close() method.
     */
    h5::Group* createOrOpenGroup(h5::Group *parentGroup, const std::string& name);

    /**
     * This function removes an existing group with the given name from the parent group, and then opens a new group with the same name.
     * 
     * @param parentGroup: A pointer to the parent group from which the existing group should be removed and the new group should be opened.
     * @param name: A const reference to the name of the existing group to be removed and the new group to be opened.
     * 
     * @return: A pointer to the newly opened group.
     * 
     * Note: The caller is responsible for closing the returned group using the appropriate H5::Group::close() method.
     */
    h5::Group* removeAndOpenGroup(h5::Group *parentGroup, const std::string& name);

    /**
     * This template function sets the value of a scalar attribute with the given name in the given group.
     * 
     * @param group: A pointer to the group in which the attribute should be set.
     * @param name: A const reference to the name of the attribute to be set.
     * @param dtype: The H5::PredType representing the data type of the attribute.
     * @param value: The value to be set for the attribute.
     * 
     */
    template<typename T>
    bool set_scaler_attribute(h5::Group *group, const std::string& name, h5::PredType dtype, const T value);
    
    /**
     * This template function retrieves the value of a scalar attribute with the given name from the given group.
     * 
     * @param group: A pointer to the group from which the attribute should be retrieved.
     * @param name: A const reference to the name of the attribute to be retrieved.
     * @param dtype: The H5::PredType representing the data type of the attribute.
     * 
     * @return: The value of the attribute. Note: If the attribute does not exist 0 is returned.
     */
    template<typename T>
    T get_scaler_attribute(h5::Group *group, const std::string& name);

    /**
     * This function sets the value of a string attribute with the given name in the given group.
     * 
     * @param group: A pointer to the group in which the attribute should be set.
     * @param name: A const reference to the name of the attribute to be set.
     * @param value: A const reference to the value to be set for the attribute.
     * 
     * Note: The caller is responsible for ensuring that the attribute exists and that the data type matches the specified function signature.
     */
    bool set_string_attribute(h5::Group *group, const std::string name, const std::string value);

    /**
     * This function retrieves the value of a string attribute with the given name from the given group.
     * 
     * @param group: A pointer to the group from which the attribute should be retrieved.
     * @param name: A const reference to the name of the attribute to be retrieved.
     * 
     * @return: The value of the attribute. Note: If the attribute does not exist an "not-found" is returned.
     */
    std::string get_string_attribute(h5::Group *group, const std::string name);


    /**
     * This function adds a dataset with the given name, data type, and data to the given group in an HDF5 file.
     * 
     * @param group: A pointer to the group in which the dataset should be added.
     * @param setName: A const reference to the name of the dataset to be added.
     * @param dtype: The H5::PredType representing the data type of the dataset.
     * @param data: A pointer to the data to be added to the dataset.
     * @param dims: A vector of integers representing the dimensions of the dataset.
     * 
     * Note: The caller is responsible for ensuring that the dataset name is unique within the group.
     */
    void addDataSetToH5(h5::Group *group, const std::string setName, h5::PredType dtype, const void *data, std::vector<size_t> dims);

        /**
     * This template function reads a dataset with the given name, data type, and stores the data in the output vector.
     *
     * @param group: A pointer to the group from which the dataset should be read.
     * @param setName: A const reference to the name of the dataset to be read.
     * @param dtype: The H5::PredType representing the data type of the dataset.
     * @param output: A reference to the vector where the read data will be stored.
     *
     * @return: A boolean value indicating the success of the read operation.
     *          - true: If the dataset is successfully read and the data is stored in the output vector.
     *          - false: If an error occurs during the read operation.
     *
     * Note: The caller is responsible for ensuring that the dataset exists in the group and that the data type matches the specified template type.
     */
    template<typename T> bool readDataSetFromH5(h5::Group *group, const std::string setName, h5::PredType dtype, std::vector<T> &output);

    /**
     * This function appends a dataset of fingerprints to an HDF5 file.
     *
     * @param group: A pointer to the HDF5 group where the dataset will be created.
     * @param setName: A const reference to the name of the dataset to be created.
     * @param cfp: A pointer to the array of fingerprints to be stored in the dataset.
     * @param size: number of elements in cfp array (size of fingerprints * number of fingerprints).
     * @param CFPSize: The size of each fingerprint in the array.
     *
     * The function creates a dataset with the given name and stores the fingerprints in the dataset.
     * If the dataset already exists, the new fingerprints are appended to the existing dataset.
     *
     * Note: The caller is responsible for ensuring that the dataset name is unique within the group.
     *       The memory pointed to by 'cfp' should be valid for the duration of the function call.
     */
    void appendFingerPrintDatasetToH5(h5::Group *group, const std::string setName, uint64_t *cfp, hsize_t size, hsize_t CFPSize);

    /**
     * This function creates an empty dataset of fingerprints with the given name in an HDF5 file.
     *
     * @param group: A pointer to the HDF5 group where the dataset will be created.
     * @param setName: A const reference to the name of the dataset to be created.
     * @param size: An hsize_t value representing the initial size of the dataset.
     * @param CFPSize: An hsize_t value representing the size of each fingerprint in the dataset.
     *
     * The function creates a dataset with the given name and allocates memory for the specified number of fingerprints.
     * The dataset is initially empty, and the caller can then update fingerprints to it using the `updateFingerprintAtIndex` function.
     *
     * Note: The caller is responsible for ensuring that the dataset name is unique within the group.
     *       The function does not modify the HDF5 file or its contents.
     */
    void createEmptyFpDataset(h5::Group *group, const std::string setName, hsize_t size, hsize_t CFPSize);

    /**
     * This function calculates and returns the size of a fingerprint dataset in an HDF5 file.
     *
     * @param group: A pointer to the HDF5 group where the dataset is located.
     * @param setName: A const reference to the name of the dataset.
     *
     * @return: An hsize_t value representing the size of the fingerprint dataset.
     *          - If the dataset exists and is of the correct type, the function returns the size of the dataset.
     *          - If the dataset does not exist or is of an incorrect type, the function returns 0.
     *
     * Note: The caller is responsible for ensuring that the dataset name is unique within the group.
     *       The function does not modify the HDF5 file or its contents.
     */
    hsize_t sizeOfFingerPrintDatasetInH5(h5::Group *group, const std::string setName);

    /**
     * This function retrieves a fingerprint from a dataset at the specified index in an HDF5 file.
     *
     * @param group: A pointer to the HDF5 group where the dataset is located.
     * @param setName: A const reference to the name of the dataset.
     * @param index: An hsize_t value representing the index of the fingerprint to retrieve.
     * @param CFPSize: An hsize_t value representing the size of each fingerprint in the dataset.
     * @param cfp: A pointer to an array of uint64_t where the retrieved fingerprint will be stored.
     *
     * @return: A boolean value indicating the success of the operation.
     *          - true: If the fingerprint is successfully retrieved and stored in the 'cfp' array.
     *          - false: If an error occurs during the retrieval operation, such as the dataset not existing,
     *                   an invalid index, or an inability to read the fingerprint data.
     *
     * Note: The caller is responsible for ensuring that the dataset name is unique within the group,
     *       the 'cfp' array has enough space to store the fingerprint, and the HDF5 file is open and valid.
     *       The function does not modify the HDF5 file or its contents.
     */
    bool getFingerprintFromIndex(h5::Group *group, const std::string setName, hsize_t index, hsize_t CFPSize, uint64_t *cfp);

    /**
     * This function updates a fingerprint at the specified index in an HDF5 dataset.
     *
     * @param group: A pointer to the HDF5 group where the dataset is located.
     * @param setName: A const reference to the name of the dataset.
     * @param index: An hsize_t value representing the index of the fingerprint to update.
     * @param CFPSize: An hsize_t value representing the size of each fingerprint in the dataset.
     * @param cfp: A pointer to an array of uint64_t containing the new/modified fingerprint to be stored.
     *
     * @return: A boolean value indicating the success of the operation.
     *          - true: If the fingerprint is successfully updated at the specified index.
     *          - false: If an error occurs during the update operation, such as the dataset not existing,
     *                   an invalid index, or an inability to write the fingerprint data.
     *
     * Note: The caller is responsible for ensuring that the dataset name is unique within the group,
     *       the 'cfp' array contains a valid fingerprint, and the HDF5 file is open and valid.
     *       The function modifies the HDF5 file or its contents.
     */
    bool updateFingerprintAtIndex(h5::Group *group, const std::string setName, hsize_t index, hsize_t CFPSize, uint64_t *cfp);

} // namespace utils