/**
 * @brief This class represents a molecular data table stored in HDF5 file.
 * @author Rajendra Kumar
 */

#include<iostream>
#include<cstring>

#include "molDataTable.h"
#include "lib.h"

const char* MolDataTable::table_name = "MoleculeDataTable";
const char* MolDataTable::mol_id_index_table_name = "MolIdIndexTable";

MolDataTable::MolDataTable(h5::Group *root_group, bool table_exist, const std::string &file_mode)
{
    _root_group = root_group;
    _file_mode = file_mode;
    _loc_id = root_group->getLocId();
    _setupDataTypes(); // set-up data types

    _table_exist = table_exist;

    if(_table_exist) {
        hsize_t nfields;
        H5TBget_table_info(_loc_id, table_name, &nfields, &_nrecords);
    }
}

void MolDataTable::_setupDataTypes() {
    _id_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(_id_type, ID_SIZE_IN_MOL_DATA_TABLE);
    _smiles_type =  H5Tcopy(H5T_C_S1);
    H5Tset_size(_smiles_type, SMILES_SIZE_IN_MOL_DATA_TABLE);
    _field_types[0] = _id_type;
    _field_types[1] =  H5T_NATIVE_INT;
    _field_types[2] = H5T_NATIVE_HSIZE;
    _field_types[3] = _smiles_type;
}

herr_t MolDataTable::_createTable(TMolData *moldata, int nrecord)
{
    // make the table
    herr_t ret = H5TBmake_table(table_name, _loc_id, table_name, MolDataTableFieldsCount, nrecord, _row_total_size_bytes,
                                _field_names, _row_offset_bytes, _field_types, _chunk_size, NULL, _compress, moldata);
    if (ret >= 0) {
        _table_exist = true;
    }

    return ret;
}

std::vector<std::string> MolDataTable::_molIdToIndexDataPath(const char *mol_id) {
    int id_length = strlen(mol_id);
    std::vector<std::string> pathVector;
    int pathLength = 6;
    int remainingLength;

    // split mol_id into sub-paths, each path will be group in HDF5 file
    for (remainingLength = 0; remainingLength < (id_length-pathLength); remainingLength += pathLength) {
        char subpath[pathLength+1];
        strncpy(subpath, mol_id+remainingLength, pathLength);
        subpath[pathLength] = '\0';
        std::string p = subpath;
        pathVector.push_back(subpath);
    }

    // append the last sub-path, which will be dataset name in HDF5 file
    if (pathVector.size() > 0) {
        size_t dsetNameLength = id_length - remainingLength;
        char dsetName[dsetNameLength + 1];
        strncpy(dsetName, mol_id+remainingLength, dsetNameLength);
        dsetName[dsetNameLength] = '\0';
        pathVector.push_back(dsetName);
    }

    return pathVector;
};

bool MolDataTable::append(TMolData *moldata, int nrecord)
{
    herr_t ret;

    // first check if table already exists, if it doesn't, create it and append records
    if (_table_exist) {
        ret = H5TBappend_records(_loc_id, table_name, nrecord, _row_total_size_bytes, _row_offset_bytes, _row_fields_size_bytes, moldata);
    } else {
        ret = _createTable(moldata, nrecord);
    }

    if (ret >= 0) {
        _nrecords += nrecord;
        _molIdToIndexMap.clear();

        // invalidate the serialized molId to index mapping in the file,
        // so that it can be re-created with new data when needed
        if (_isMolIdToIndexMapBuilt) {
            _isMolIdToIndexMapBuilt = false;
            _invalidateSerializedMolIdToIndexMap();
        }
    }

    return ret < 0 ? false : true;
}

bool MolDataTable::_canWriteSerializedMolIdToIndexMap() const {
    return (_file_mode == "w") || (_file_mode == "a");
}

bool MolDataTable::_loadMolIdToIndexMapFromH5() {
    if (_nrecords == 0) {
        _molIdToIndexMap.clear();
        return true;
    }

    if (H5Lexists(_loc_id, mol_id_index_table_name, H5P_DEFAULT) <= 0) {
        return false;
    }

    hsize_t nfields = 0;
    hsize_t nrecords = 0;
    if (H5TBget_table_info(_loc_id, mol_id_index_table_name, &nfields, &nrecords) < 0) {
        return false;
    }

    if (nfields != 2 || nrecords != _nrecords) {
        return false;
    }

    TMolIdIndex *indexRows = new TMolIdIndex[nrecords];
    size_t row_total_size = sizeof(TMolIdIndex);
    size_t row_offsets[2] = {HOFFSET(TMolIdIndex, id), HOFFSET(TMolIdIndex, index)};
    size_t row_sizes[2] = {sizeof(indexRows[0].id), sizeof(indexRows[0].index)};

    herr_t ret = H5TBread_records(_loc_id, mol_id_index_table_name, 0, nrecords, row_total_size, row_offsets, row_sizes, indexRows);
    if (ret < 0) {
        delete[] indexRows;
        return false;
    }

    _molIdToIndexMap.clear();
    _molIdToIndexMap.reserve(nrecords);
    for (hsize_t i = 0; i < nrecords; ++i) {
        _molIdToIndexMap[indexRows[i].id] = indexRows[i].index;
    }

    delete[] indexRows;
    return true;
}

bool MolDataTable::_writeMolIdToIndexMapToH5() {
    if (!_canWriteSerializedMolIdToIndexMap()) {
        return false;
    }

    if (_molIdToIndexMap.size() == 0) {
        return true;
    }

    if (H5Lexists(_loc_id, mol_id_index_table_name, H5P_DEFAULT) > 0) {
        H5Ldelete(_loc_id, mol_id_index_table_name, H5P_DEFAULT);
    }

    hsize_t nrecords = _molIdToIndexMap.size();
    TMolIdIndex *indexRows = new TMolIdIndex[nrecords];

    hsize_t i = 0;
    for (const auto &it : _molIdToIndexMap) {
        memset(indexRows[i].id, 0, sizeof(indexRows[i].id));
        strncpy(indexRows[i].id, it.first.c_str(), sizeof(indexRows[i].id) - 1);
        indexRows[i].index = it.second;
        i++;
    }

    hid_t id_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(id_type, ID_SIZE_IN_MOL_DATA_TABLE);
    hid_t field_types[2] = {id_type, H5T_NATIVE_UINT64};
    const char *field_names[2] = {"ID", "index"};
    size_t row_total_size = sizeof(TMolIdIndex);
    size_t row_offsets[2] = {HOFFSET(TMolIdIndex, id), HOFFSET(TMolIdIndex, index)};
    size_t row_sizes[2] = {sizeof(indexRows[0].id), sizeof(indexRows[0].index)};

    herr_t ret = H5TBmake_table(mol_id_index_table_name, _loc_id, mol_id_index_table_name,
                                2, nrecords, row_total_size, field_names, row_offsets,
                                field_types, _chunk_size, NULL, _compress, indexRows);

    H5Tclose(id_type);
    delete[] indexRows;

    return ret >= 0;
}

void MolDataTable::_invalidateSerializedMolIdToIndexMap() {
    if (!_canWriteSerializedMolIdToIndexMap()) {
        return;
    }

    if (H5Lexists(_loc_id, mol_id_index_table_name, H5P_DEFAULT) > 0) {
        H5Ldelete(_loc_id, mol_id_index_table_name, H5P_DEFAULT);
    }
}

void MolDataTable::buildMolIdToIndexMap() {
    if (_isMolIdToIndexMapBuilt) {
        return;
    }

    if (_loadMolIdToIndexMapFromH5()) {
        _isMolIdToIndexMapBuilt = true;
        return;
    }

    std::cout << " Building map for molecule Id to index in molecule data table..." << std::endl;
    _molIdToIndexMap.reserve(_nrecords); // reserve space for the map
    hsize_t chunk_size = 1000000;
    TMolData *molData = new TMolData[chunk_size];
    for(hsize_t i = 0; i < _nrecords; i+=chunk_size) {
        hsize_t num_to_read = std::min(_nrecords - i, chunk_size);
        H5TBread_records(_loc_id, table_name, i, num_to_read, _row_total_size_bytes, _row_offset_bytes, _row_fields_size_bytes, molData);
        for (hsize_t j = 0; j < num_to_read; j++) {
            std::string id = molData[j].id;
            _molIdToIndexMap[id] = i + j;
        }
    }
    delete[] molData;

    _isMolIdToIndexMapBuilt = true;
    _writeMolIdToIndexMapToH5();
    std::cout << " ... finished building map." << std::endl;
}

std::vector<uint64_t> MolDataTable::getIndexOfMolDataFromID(std::string id) {
    buildMolIdToIndexMap();

    std::vector<uint64_t> index;
    auto it = _molIdToIndexMap.find(id);
    if (it != _molIdToIndexMap.end()) {
        index.push_back(it->second);
    } else {
        std::cerr << "Molecule id '"<<id<<"' not found in molecule data table" << std::endl;
    }

    // if id not found in the table, index will be empty
    return index;
};

void MolDataTable::getMolDataFromId(std::string id, TMolData *molData) {
    std::vector<uint64_t> index = getIndexOfMolDataFromID(id); // get index for this mol_id from `molDataTableIndex`
    if (index.size() > 0) { // if index found, read record from `molDataTable`
        H5TBread_records(_loc_id, table_name, index[0], 1, _row_total_size_bytes, _row_offset_bytes, _row_fields_size_bytes, molData);
    }
};

void MolDataTable::replaceMolDataForId(std::string id, TMolData *molData) {
    std::vector<uint64_t> index = getIndexOfMolDataFromID(id); // get index for this mol_id from `molDataTableIndex`
    if (index.size() > 0) { // if index found, read record from `molDataTable`
        H5TBwrite_records(_loc_id, table_name, index[0], 1, _row_total_size_bytes, _row_offset_bytes, _row_fields_size_bytes, molData);
        _molIdToIndexMap.clear();
        _isMolIdToIndexMapBuilt = false;
        _invalidateSerializedMolIdToIndexMap();
    }
}