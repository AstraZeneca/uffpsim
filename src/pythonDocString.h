/**
 * @brief Contains the doc-strings that will be exposed on python side
 * @author Rajendra Kumar
 */

const char* FpStoreSetInfoDoc = R"%(Sets the additional information in the HDF5 file associated with the FingerprintStore object.
This method writes the additional information to the root attributes of the HDF5 file.
The additional information can be used to store any additional details about the data stored in the file.

Parameters
----------
info : str
    The additional information to be set.

Return
-------
True
)%";

const char* FpStoreAppendFingerprintsDoc = R"%(Appends new molecular fingerprint data to the HDF5 file.

This method takes a vector of molecular fingerprint data in the form of tuples (ID, SMILES, Fingerprint) and appends 
them to the HDF5 file.

The new molecular fingerprint data is written to the appropriate datasets within the HDF5 file.

Parameters
----------
fingerprints : list of tuples
     A list of tuples containing molecular fingerprint information (ID, SMILES, Fingerprint) where fingerprint is a string of 0 and 1.
)%";

const char* FpStorePerformInnerClusteringAndWriteDoc = R"%(Performs inner clustering on the molecular fingerprint data and writes the results to the HDF5 file.

This method performs inner clustering on the molecular fingerprint data using the specified inner clustering threshold.
The inner clustering algorithm groups similar molecular fingerprints together based on their popcount (population count).
The results of the inner clustering, such as cluster IDs and compact fingerprints, are written to the appropriate datasets 
within the HDF5 file.

The method ensures that the inner clustering results are accurate and efficiently stored in the HDF5 file.

Return
-------
None
)%";

const char* FpStoreReDoInnerClusteringInMemory = R"%(Re-Performs inner clustering on the molecular fingerprint data using the specified threshold.

Parameters
----------
threshold : Float
    The inner clustering threshold for clustering molecular fingerprints.

)%";

const char* FpStoreMagicNumberExistsDoc = R"%(Checks if the magic number exists in the HDF5 file.

This method reads the root attributes of the HDF5 file and checks if the magic number exists.
The magic number is a unique identifier used to verify the integrity of the data stored in the file.

Return
-------
bool (True/False) If the magic number exists, it returns True; otherwise, it returns False.)%";

const char* FpStoreSetMagicNumberDoc = R"%(Sets the magic number in the HDF5 file.

This method writes the magic number to the root attributes of the HDF5 file.
The magic number is a unique identifier used to verify the integrity of the data stored in the file.

Return
-------
bool (True/False) if the magic number is successfully set, false otherwise.
)%";

const char* FpStorebuildMolIdToIndexMap = R"%(Build the mapping of Molecule ID to index in molecule data table.

This method builds a map that associates each molecule ID with its corresponding index in the molecule data table.
This is required to access SMILES based on the given molecule ID.

)%";

const char* FpStoreGetSmilesFromIdDoc = R"%(Retrieves the SMILES string associated with a given molecular ID from the HDF5 file.

This method takes a molecular ID as input and searches for the corresponding molecular data in the HDF5 file.\
It then retrieves the SMILES string associated with the given molecular ID.

.. note: if mapping of molecule-id to index in molecule data table is not performed, it is performed in the background. To perform this mapping,
         use :func:`build_mol_id_to_index_map()`.

Parameters
----------
id : str
    The molecular ID for which the SMILES string needs to be retrieved.

Returns
-------
str The SMILES string associated with the given molecular ID.

Raises
------
    std::runtime_error If the molecular ID is not found in the HDF5 file.
)%";

const char* FpStoreUpdateFingerprintsDoc = R"%(Updates the molecular fingerprint data in the HDF5 file.

This method takes a vector of molecular fingerprint data in the form of tuples (ID, SMILES, Fingerprint) as input and
updates the corresponding datasets within the HDF5 file.

The method compares the existing molecular fingerprints with the new fingerprints and updates the HDF5 file accordingly.

Parameters
----------
fingerprints : list of tuples
    A vector of tuples containing molecular fingerprint information (ID, SMILES, Fingerprint).

Returns
-------
None
)%";

const char* FpStoreLoadDataInMemoryDoc = R"%(Loads the molecular fingerprint data into memory for efficient access.

This method reads the necessary molecular fingerprint data from the HDF5 file and stores it in memory for further operations.
It is only applicable when the HDF5 file is opened in read-only mode.

It loads all the cluster-fingerprints and fingerprints into memory for full-memory search.

)%";

const char* FpStoreCloseDoc = R"%(Closes the HDF5 file associated with the FingerprintStore object.

This method closes the HDF5 file, freeing up any system resources that were being used by the file.
After calling this method, any further operations on the HDF5 file will result in an error.

It is important to call this method when you are done with the FingerprintStore object to ensure proper cleanup and resource management.
)%";

const char* FpSearchEngineCloseDoc = R"%(The close method closes the fingerprint store file and releases any associated resources.

.. note:: The close method should be called when the FPSearchEngine object is no longer needed to free up system resources.

)%";

const char* FpSearchEngineSearchDoc = R"%(The search method performs a single-query molecular similarity search using the Tanimoto coefficient.

Parameters
----------
query : str
    The query molecular fingerprint string.
threshold : float
    The similarity threshold for filtering the results.
limits : int
    The maximum number of results to return.

Returns
-------
A vector of tuples, where each tuple contains the molecule ID and similarity score.

)%";

const char* FpSearchEngineBatchSearchDoc = R"%(The batchSearch method performs a batch molecular similarity search using the Tanimoto coefficient.

Parameters
----------
queries : list of str
    A list of query fingerprint strings.
threshold: float
    The similarity threshold for filtering the results.
limits: int
    The maximum number of results to return for each query.

Returns
-------
A list of lists, where each inner list contains tuples, and each tuple contains the molecule ID and similarity score for a query.

)%";