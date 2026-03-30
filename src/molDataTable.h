/**
 * @brief This class represents a molecular data table stored in HDF5 file.
 * @author Rajendra Kumar
 */

#include <vector>
#include <unordered_map>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdlib.h>
#include <H5Cpp.h>

namespace h5 = H5;

#define MolDataTableFieldsCount (hsize_t) 4
#define SMILES_SIZE_IN_MOL_DATA_TABLE (size_t) 1024
#define ID_SIZE_IN_MOL_DATA_TABLE (size_t) 64

typedef struct {
    char id[ID_SIZE_IN_MOL_DATA_TABLE] = {""};
    int popcount = 0;
    hsize_t fpIndex = 0;
    char smiles[SMILES_SIZE_IN_MOL_DATA_TABLE] = {""};
} TMolData;

typedef struct {
    char id[ID_SIZE_IN_MOL_DATA_TABLE] = {""};
    uint64_t index = 0;
} TMolIdIndex;



/**
 * The MolDataTable class is responsible for managing and interacting with a molecule data table in an HDF5 file.
 * It provides methods for creating, appending, retrieving, and updating data in the table.
 *
 * The class uses the HDF5 library for file I/O operations and provides error handling and exception handling.
 *
 * The class has private member variables for managing the HDF5 file, table, and data types, as well as public member functions for interacting with the table.
 *
 * The class also includes private helper methods for setting up data types, creating the table, and converting molecule IDs to data paths for indexing.
 *
 * Note: The specific implementation details and requirements of the class may vary depending on the requirements of the application.
 */
class MolDataTable {
    

    private:    
        bool _table_exist = true;                        // flag to check if table `MoleculeDataTable` exists in HDF5 or not
        std::string _file_mode = "r";                   // HDF5 file mode inherited from FingerprintStore
        h5::Group *_root_group;                          // HDF5 root group pointer
        hid_t _loc_id;                                   // location ID of the root group
        hid_t _table_id;                                 // ID of the table `MoleculeDataTable`
        hsize_t NFields = 4;                             // number of fields in `MoleculeDataTable`
        TMolData _temp_mol_data[1];                      // temporary storage for TMolData
        size_t _row_total_size_bytes = sizeof(TMolData); // total size of TMolData in bytes
        size_t _row_offset_bytes[MolDataTableFieldsCount] = { // offset of each field in TMolData in bytes
            HOFFSET(TMolData, id),
            HOFFSET(TMolData, popcount),
            HOFFSET(TMolData, fpIndex),
            HOFFSET(TMolData, smiles),
        };
        size_t _row_fields_size_bytes[MolDataTableFieldsCount] = { // size of each field in bytes
            sizeof(_temp_mol_data[0].id),
            sizeof(_temp_mol_data[0].popcount),
            sizeof(_temp_mol_data[0].fpIndex),
            sizeof(_temp_mol_data[0].smiles),
        };

        // HDF5 field types and names for `MoleculeDataTable`
        const char *_field_names[MolDataTableFieldsCount] = {"ID", "popcount", "fpIndex", "SMILES"};

        hid_t _id_type;                              // HDF5 type for ID field
        hid_t _smiles_type;                          // HDF5 type for SMILES field
        hid_t _field_types[MolDataTableFieldsCount]; // HDF5 types for all four fields in table
        hsize_t _chunk_size = 1000;
        int _compress = 1;

        std::unordered_map<std::string, uint64_t> _molIdToIndexMap; // molecule ID to molDataTable index mapping
        bool _isMolIdToIndexMapBuilt = false; // flag to check if the mapping is built

        /**
         * This private method is responsible for setting up the data types for the fields in the `MoleculeDataTable`.
         * It initializes the HDF5 types for the ID, popcount, fpIndex, and SMILES fields, and stores them in the `_field_types` array.
         *
         * Note: This method should be called before creating the table or appending data to it.
         */
        void _setupDataTypes();

        /**
         * This private method is responsible for creating the `MoleculeDataTable` in the HDF5 file.
         * It initializes the table with the specified number of records, sets up the data types for the fields, and writes the initial data to the table.
         *
         * @param moldata: A pointer to an array of `TMolData` structures containing the initial data to be written to the table.
         * @param nrecord: The number of records to be written to the table.
         *
         * @return: An HDF5 error code (herr_t) indicating the success or failure of the operation.
         *
         * Note: This method should be called after setting up the data types using the `_setupDataTypes` method.
         */
        herr_t _createTable(TMolData *moldata, int nrecord);

        /**
         * This private method is responsible for converting a molecule ID into a data path for indexing.
         * It takes a molecule ID as input and returns a vector of strings representing the data path.
         *
         * @param mol_id: A pointer to a C-style string containing the molecule ID.
         *
         * @return: A vector of strings representing the data path for indexing.
         *
         */
        std::vector<std::string> _molIdToIndexDataPath(const char *mol_id);

    public:
        /**
         * The constructor of the MolDataTable class initializes the class with the provided HDF5 group and table existence flag.
         * It sets up the necessary member variables and initializes the HDF5 file, table, and data types.
         *
         * @param root_group: A pointer to the HDF5 group where the molecule data table will be created or accessed.
         * @param table_exist: A boolean flag indicating whether the molecule data table already exists in the HDF5 file.
         *
         * Note: The constructor performs error handling and exception handling to ensure that the initialization process is successful.
         */
        MolDataTable(h5::Group *root_group, bool table_exist, const std::string &file_mode = "r");

        // Public properties
        hsize_t _nrecords = 0;
        static const char *table_name;
        static const char *mol_id_index_table_name;

    private:
        bool _canWriteSerializedMolIdToIndexMap() const;
        bool _loadMolIdToIndexMapFromH5();
        bool _writeMolIdToIndexMapToH5();
        void _invalidateSerializedMolIdToIndexMap();

    public:
        
        /**
         * The append method is responsible for appending new data records to the molecule data table in the HDF5 file.
         * It takes a pointer to an array of `TMolData` structures containing the data to be appended and the number of
         * records to be appended as input parameters.
         *
         * The method performs error handling and exception handling to ensure that the data is appended successfully.
         *
         * @param moldata: A pointer to an array of `TMolData` structures containing the data to be appended.
         * @param nrecord: The number of records to be appended to the table.
         *
         * @return: A boolean value indicating whether the data was successfully appended to the table.
         *
         * Note: The append method assumes that the table and data types have already been initialized and set up.
         */
        bool append(TMolData *moldata, int nrecord);

        /**
         * @brief Builds a mapping from molecular IDs to their corresponding indices in the molecular data table.
         *
         * This method creates an internal map that associates each molecular ID with its index in the molecular data table.
         * The mapping is used for efficient lookup and retrieval of molecular data based on molecular IDs.
         *
         * @return None.
         */
        void buildMolIdToIndexMap();

        /**
         * The getIndexOfMolDataFromID method is responsible for retrieving the index of mol_id in the `MoleculeDataTable`
         * in the HDF5 file.
         * It takes a molecule ID as input and returns a vector of unsigned 64-bit integer with one element containing index.
         * If mol_id is not found in the table, an empty vector is returned.
         *
         * @param id: A string representing the molecule ID for which the index(es) are to be retrieved.
         *
         * @return: A vector of unsigned 64-bit integers representing the index of the mol_id in the table.
         *
         */
        std::vector<uint64_t> getIndexOfMolDataFromID(std::string id);

        /**
         * The getMolDataFromId method is responsible for retrieving the molecule data record with the specified 
         * molecule ID from the `MoleculeDataTable` in the HDF5 file.
         * It takes a molecule ID as input and a pointer to a `TMolData` structure where the retrieved data will be stored.
         *
         *
         * @param id: A string representing the molecule ID for which the data record is to be retrieved.
         * @param molData: A pointer to a `TMolData` structure where the retrieved data will be stored.
         *
         */
        void getMolDataFromId(std::string id, TMolData *molData);

        /**
         * The replaceMolDataForId method is responsible for replacing the molecule data record with the specified molecule 
         * ID in the `MoleculeDataTable` in the HDF5 file.
         * It takes a molecule ID as input and a pointer to a `TMolData` structure containing the new data to be stored.
         *
         *
         * @param id: A string representing the molecule ID for which the data record is to be replaced.
         * @param molData: A pointer to a `TMolData` structure containing the new data to be stored.
         */
        void replaceMolDataForId(std::string id, TMolData *molData);
};
