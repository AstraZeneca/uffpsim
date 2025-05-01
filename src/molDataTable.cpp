/**
 * @brief This class represents a molecular data table stored in HDF5 file.
 * @author Rajendra Kumar
 */

#include<iostream>
#include<cstring>

#include "molDataTable.h"
#include "lib.h"

const char* MolDataTable::table_name = "MoleculeDataTable";

MolDataTable::MolDataTable(h5::Group *root_group, bool table_exist)
{
    _root_group = root_group;
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

    return ret < 0 ? false : true;
}

void MolDataTable::buildMolIdToIndexMap() {
    if (_isMolIdToIndexMapBuilt) {
        std::cout << " MolIdToIndexMap already built, skipping..." << std::endl;
        return;
    }

    std::cout << " Building map for molecule Id to index in molecule data table..." << std::endl;
    TMolData molData;
    for(hsize_t i = 0; i < _nrecords; i++) {
        H5TBread_records(_loc_id, table_name, i, 1, _row_total_size_bytes, _row_offset_bytes, _row_fields_size_bytes, &molData);
        std::string id = molData.id;
        _molIdToIndexMap[id] = i;
    }

    _isMolIdToIndexMapBuilt = true;
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
    }
}