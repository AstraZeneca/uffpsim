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

    // Exploitation of b-tree algorithm in HDF5 for faster retrieval of records by mol_id
    // where mol_id is splitted into sub-paths and stored as parent-child groups.
    // Each group has two groups as `datasets` and `paths`. `paths` group contains all subsequent sub-paths.
    // `datasets` group contains actual data i.e. index in MolDataTable for this mol_id.
    if (ret >= 0){
        h5::Group *molDataIndexGroup = utils::createOrOpenGroup(_root_group, "molDataTableIndex");
        for(int n=0; n < nrecord; n++) {
            const uint64_t index[1] = { (uint64_t) _nrecords};
            std::vector<size_t> dims = {1};
            
            std::vector<std::string> pathVector = _molIdToIndexDataPath(moldata[n].id);
            if (pathVector.size() == 0) { // if only one sub-path it is mol_id and put it in datasets
                h5::Group *group = utils::createOrOpenGroup(molDataIndexGroup, "datasets");
                utils::addDataSetToH5(group, moldata[n].id, h5::PredType::NATIVE_UINT64, index, dims );
                group->close();
                delete group;
            } else { // if multiple sub-paths, create parent-child groups and put it in paths, datasets
                std::vector<h5::Group*> pathGroups;
                h5::Group *group = nullptr;
                for(size_t g = 0; g < pathVector.size()-1; g++) { // create parent-child groups through `paths` group
                    if (g == 0) {
                        group = utils::createOrOpenGroup(molDataIndexGroup, "paths");
                        pathGroups.push_back(group);
                        group = utils::createOrOpenGroup(group, pathVector[g]);
                    } else {
                        group = utils::createOrOpenGroup(pathGroups[g-1], "paths");
                        pathGroups.push_back(group);
                        group = utils::createOrOpenGroup(group, pathVector[g]);
                    }
                    pathGroups.push_back(group);
                }
                
                if (group != nullptr) {// last of the path, put it in datasets
                    h5::Group *datasetGroup = utils::createOrOpenGroup(group, "datasets");
                    utils::addDataSetToH5(datasetGroup, pathVector[pathVector.size()-1], h5::PredType::NATIVE_UINT64, index, dims );
                    datasetGroup->close();
                    delete datasetGroup;
                }

                for(auto &grp: pathGroups) {
                    grp->close();
                    delete grp;
                }
                pathGroups.resize(0);
            }
            
            _nrecords += 1;
        }
        molDataIndexGroup->close();
        delete molDataIndexGroup;
    }

    return ret < 0 ? false : true;
}

std::vector<uint64_t> MolDataTable::getIndexOfMolDataFromID(std::string id) {
    h5::Group *molDataIndexGroup = utils::createOrOpenGroup(_root_group, "molDataTableIndex");
    std::vector<std::string> pathVector = _molIdToIndexDataPath(id.c_str()); // first split mol_id into sub-paths
    std::vector<uint64_t> index;

    if (pathVector.size() == 0) { // if only one sub-path it is mol_id and get it from datasets
        h5::Group *group = utils::createOrOpenGroup(molDataIndexGroup, "datasets");
        
        if(group->nameExists(id)) {
            utils::readDataSetFromH5(group, id, H5::PredType::NATIVE_UINT64, index);
        }
        delete group;
    } else { // if multiple sub-paths, get it from sub-paths, and end datasets
        std::vector<h5::Group*> pathGroups;
        h5::Group *group = nullptr;
        for(size_t g = 0; g < pathVector.size()-1; g++) { // get parent-child groups through `paths` group
            if (g == 0) {
                group = utils::createOrOpenGroup(molDataIndexGroup, "paths");
                pathGroups.push_back(group);
                group = utils::createOrOpenGroup(group, pathVector[g]);
            } else {
                group = utils::createOrOpenGroup(pathGroups[g-1], "paths");
                pathGroups.push_back(group);
                group = utils::createOrOpenGroup(group, pathVector[g]);
            }
            pathGroups.push_back(group);
        }

        if (group->nameExists("datasets")) { // reached last of the path, get it from datasets
            h5::Group *datasetGroup = utils::createOrOpenGroup(group, "datasets");
            std::string datasetName = pathVector[pathVector.size()-1];
            if(datasetGroup->nameExists(datasetName)) {
                utils::readDataSetFromH5(datasetGroup, datasetName, H5::PredType::NATIVE_UINT64, index);
            }
            delete datasetGroup;
        }


        for(auto &grp: pathGroups) {
            delete grp;
        }
        pathGroups.resize(0);
    }
    
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