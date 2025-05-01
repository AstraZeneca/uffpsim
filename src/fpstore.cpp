/**
 * @brief This class represents the FingerprintStore class that handles the storage and retrieval of molecular fingerprints.
 * @author Rajendra Kumar
 */

#include <string>
#include <iostream>
#include <map>
#include <tuple>
#include <set>
#include <vector>
#include <H5Cpp.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>

#include "fpstore.h"
#include "innerClustering.h"

namespace py = pybind11;
namespace h5 = H5;

FingerprintStore::~FingerprintStore() {
    close();
}

void FingerprintStore::close() {

    if (_root_group != nullptr) {
        _root_group->close();
        delete _root_group;
        _root_group = nullptr;
    }

    if (_file != nullptr) {
        int types = H5F_OBJ_DATASET | H5F_OBJ_GROUP | H5F_OBJ_DATATYPE | H5F_OBJ_ATTR;
        ssize_t num_open_objects = _file->getObjCount(types); // number of open objects in the file
        if (num_open_objects > 0) { // number of open objects > 0, close them
            hid_t *obj_ids = (hid_t *) malloc(sizeof(hid_t) * num_open_objects);
            _file->getObjIDs(types, num_open_objects, obj_ids); // get all open object IDs
            for(int i = 0; i < num_open_objects; i++) { // close each open object
                _file->closeObjId(obj_ids[i]);
            }
            free(obj_ids);
        }
        _file->close();
        delete _file;
        _file = nullptr;
    }
}

void FingerprintStore::freeMemory() {
    for(int pidx = 0; pidx < _popCountBins.size(); pidx++) {
        free(_fp_inner_clusters_by_popcount[pidx].clusterFp);
        free(_fp_inner_clusters_by_popcount[pidx].fp);
    }
    free(_fp_inner_clusters_by_popcount);
}

FingerprintStore::FingerprintStore(const std::string& filename, int molIdMaxLength, std::string mode, int fpSize, 
                                   std::string fp_params, std::string info, float clusterThreshold, std::string clusterMode, bool clusterParallel) {
    _filename = H5std_string(filename);
    _file_mode = mode;

    unsigned int mode_flag = H5F_ACC_RDONLY;
    if (mode == "w") {
        mode_flag = H5F_ACC_TRUNC;
        _molIdMaxLength = molIdMaxLength;
        _innerClusteringThreshold = clusterThreshold;
        _fp_params = fp_params;
        _molIdOffset = (_molIdMaxLength + 1)/8;
        _fpSize = fpSize/64;
        _CFPSize = _molIdOffset + _fpSize + 1;
        _CFPPopCountIndex = _molIdOffset + _fpSize;
        _info = info;
    } else if (mode == "a") {   
        mode_flag = H5F_ACC_RDWR;
    }

    _file = new h5::H5File(_filename, mode_flag);
    _root_group = new h5::Group(_file->openGroup("/"));
    _molDataTable = new MolDataTable(_root_group, _root_group->nameExists(MolDataTable::table_name));

    if (((mode == "a") || (mode == "r")) && !magicNumberExists()) {
        throw std::runtime_error("The input h5 file is not created by uffpsim!");
    }

    if (mode == "w") {
        setMagicNumber();
        _writeRootAttributes(); // write all root attributes to the hdf5 file
    }

    if (mode != "w") _populateRootAttributes(); // read and load all root attributes from the hdf5 file

    if ((mode == "a") || (mode == "w")) {
        _innerClusteringAgent = new InnerClusteringAgent(this, clusterMode, clusterParallel);
    }
    
}

bool FingerprintStore::magicNumberExists() {
    std::string magic_number = utils::get_string_attribute(_root_group, MagicFingerPrintNumberKey);
    if (magic_number != "not-found") {
        if (magic_number == MagicFingerPrintNumberValue)
            return true;
    }
    return false;
}

bool FingerprintStore::setMagicNumber() {
    return utils::set_string_attribute(_root_group, MagicFingerPrintNumberKey, MagicFingerPrintNumberValue);
}

bool FingerprintStore::_writeRootAttributes() {
    utils::set_scaler_attribute(_root_group, "molIdMaxLength", h5::PredType::NATIVE_INT, _molIdMaxLength);
    utils::set_scaler_attribute(_root_group, "fpSize", h5::PredType::NATIVE_INT, _fpSize);
    utils::set_scaler_attribute(_root_group, "CFPSize", h5::PredType::NATIVE_INT, _CFPSize);
    utils::set_scaler_attribute(_root_group, "CFPPopCountIndex", h5::PredType::NATIVE_INT, _CFPPopCountIndex);
    utils::set_scaler_attribute(_root_group, "molIdOffset", h5::PredType::NATIVE_INT, _molIdOffset);
    utils::set_scaler_attribute(_root_group, "innerClusteringThreshold", h5::PredType::NATIVE_FLOAT, _innerClusteringThreshold);
    utils::set_string_attribute(_root_group, "fp_params", _fp_params);
    utils::set_string_attribute(_root_group, "info", _info);
    return true;
}

bool FingerprintStore::setInfo(std::string info) {
    if (_file_mode == "w" || _file_mode == "a") {
        _info = info;
        return utils::set_string_attribute(_root_group, "info", _info);
    } else {
        std::cerr << "Warning: Cannot set info attribute in read-only mode!" << std::endl;
    }
    return false;
}

bool FingerprintStore::_populateRootAttributes() {
    if (_popCountBins.size() == 0 && _root_group->nameExists(_popCountDatasetName)) {// read popcount bins
        utils::readDataSetFromH5(_root_group, _popCountDatasetName, h5::PredType::NATIVE_INT, _popCountBins);
    }
    _molIdMaxLength = utils::get_scaler_attribute<int>(_root_group, "molIdMaxLength");
    _fpSize = utils::get_scaler_attribute<int>(_root_group, "fpSize");
    _CFPSize = utils::get_scaler_attribute<int>(_root_group, "CFPSize");
    _CFPPopCountIndex = utils::get_scaler_attribute<int>(_root_group, "CFPPopCountIndex");
    _molIdOffset = utils::get_scaler_attribute<int>(_root_group, "molIdOffset");
    _innerClusteringThreshold = utils::get_scaler_attribute<float>(_root_group, "innerClusteringThreshold");
    _fp_params = utils::get_string_attribute(_root_group, "fp_params");
    _info = utils::get_string_attribute(_root_group, "info");
    return true;
}

bool FingerprintStore::_populateDataInMemory() {
    if(_root_group->nameExists(_popCountBinsGroupName)) {
        h5::Group *popcountBinsGroup =  utils::createOrOpenGroup(_root_group, _popCountBinsGroupName);
        _fp_inner_clusters_by_popcount = new utils::dt_inner_clusters_fingerprints[_popCountBins.size()];
        for (size_t i=0; i < _popCountBins.size(); i++) { // for each popcount bin
            if (popcountBinsGroup->exists(_getPopCountGroupName(_popCountBins[i]))) {
                utils::dt_inner_clusters_fingerprints inner_clusters_fingerprints;
                _populateClustersInMemory(_popCountBins[i], &inner_clusters_fingerprints); // populate all the clusters
                _fp_inner_clusters_by_popcount[i] = inner_clusters_fingerprints;
            } else {
                throw std::runtime_error("The following path in h5 file does not exist: " + _getPopCountGroupName(_popCountBins[i])
                                        + "!\n  There are inconsistencies in h5 file. Please rebuild the file!");
            }
        }
        delete popcountBinsGroup;

    } else {
        throw std::runtime_error("The following path in h5 file does not exist: " + _popCountBinsGroupName
                                            + "!\n  There are inconsistencies in h5 file. Please rebuild the file!");
    }
    return true;
}

void FingerprintStore::_populateClustersInMemory(int popCount, utils::dt_inner_clusters_fingerprints *inner_clusters_fingerprints){
    h5::Group *popCountBinsGroup =  utils::createOrOpenGroup(_root_group, _popCountBinsGroupName);
    h5::Group *popCountGroup = utils::createOrOpenGroup(popCountBinsGroup, _getPopCountGroupName(popCount));

    // it is for extreme edge case when fp modified and popCount bin has no other fingerprints, so no clusters.
    if (!popCountGroup->nameExists(_clustersGroupName)) {
        return;
    }

    h5::Group *clustersGroup = new h5::Group(popCountGroup->openGroup(_clustersGroupName));
    std::vector<int> clusterIDs;
    
    if (utils::readDataSetFromH5(clustersGroup, _ClusterIdDataSetName, h5::PredType::NATIVE_UINT64, clusterIDs)) {
        inner_clusters_fingerprints->popCount = popCount;
        inner_clusters_fingerprints->num_clusters = clusterIDs.size();
        
        if (clustersGroup->nameExists(_clustersFPGroupName) && clustersGroup->nameExists(_fpArrayInClusterGroupName)) {
            hsize_t lengthOfFPsInClusterFP = utils::sizeOfFingerPrintDatasetInH5(clustersGroup, _clustersFPGroupName);
            inner_clusters_fingerprints->clusterFp = (uint64_t*) malloc(sizeof(uint64_t) * lengthOfFPsInClusterFP);
            utils::getFingerprintFromIndex(clustersGroup, _clustersFPGroupName, 0, lengthOfFPsInClusterFP, inner_clusters_fingerprints->clusterFp);

            hsize_t lengthOfFPsInClusters = utils::sizeOfFingerPrintDatasetInH5(clustersGroup, _fpArrayInClusterGroupName);
            inner_clusters_fingerprints->fp = (uint64_t*) malloc(sizeof(uint64_t) * lengthOfFPsInClusters);
            utils::getFingerprintFromIndex(clustersGroup, _fpArrayInClusterGroupName, 0, lengthOfFPsInClusters, inner_clusters_fingerprints->fp);

            inner_clusters_fingerprints->num_fps = lengthOfFPsInClusters / _CFPSize;
        }

    } else {
        throw std::runtime_error("The following path in h5 file does not exist: " + _ClusterIdDataSetName 
                                    + "!\n  There are inconsistencies in h5 file. Please rebuild the file!");
    }
    delete clustersGroup;
    delete popCountGroup;
    delete popCountBinsGroup;
}

void FingerprintStore::AppendFingerprints(const std::vector<std::tuple<std::string, std::string, std::string>> fingerprints) {
    if (fingerprints.size() == 0) {
        std::cerr << "Warning: No fingerprints to append!" << std::endl;
        return;
    }

    // transformed mol and fp data from input
    std::vector<utils::dt_mol_fp_data> mol_fp_data = _buildSortedMolFpData(fingerprints);

    // write fingerprints and mol_data
    _saveFpMolData(mol_fp_data);
    
    // clear memories
    for(auto &molFp: mol_fp_data) free(molFp.fp);
    mol_fp_data.resize(0);
}

std::vector<utils::dt_mol_fp_data> FingerprintStore::_buildSortedMolFpData(const std::vector<std::tuple<std::string, std::string, std::string>> fingerprints) {
    // transformed mol and fp data from input
    std::vector<utils::dt_mol_fp_data> mol_fp_data;
    int32_t num_fps = fingerprints.size();
    
    // initialization
    mol_fp_data.reserve(num_fps);
    for(int32_t ifp=0; ifp<num_fps; ifp++) {
        utils::dt_mol_fp_data mol_data;
        mol_data.fp = (uint64_t*) malloc(sizeof(uint64_t) * _CFPSize); // will hold the compacted fingerprint
        mol_fp_data.push_back(mol_data);
    }

    for (int32_t ifp=0; ifp<num_fps; ifp++) {
        std::string mol_id = std::get<0>(fingerprints[ifp]);
        std::string fp_string = std::get<1>(fingerprints[ifp]);
        std::string smiles = std::get<2>(fingerprints[ifp]);

        // id is copied here as it is modified later on
        strcpy(mol_fp_data[ifp].mol_id, mol_id.c_str());

        // only store smiles when within limit, avoid any memory-error caused by it
        if (smiles.size() < 1024) strcpy(mol_fp_data[ifp].smiles, smiles.c_str());

        // compact the fingerprint 
        utils::BitStrToCompactFPArray(fp_string, mol_id, mol_fp_data[ifp].fp, _molIdOffset, _CFPSize);

        // assign popCount
        mol_fp_data[ifp].popCount = (int) mol_fp_data[ifp].fp[_CFPPopCountIndex];
    }

    // sort mol_fp_data according to popcount
    std::sort(mol_fp_data.begin(), mol_fp_data.end(), [this](const utils::dt_mol_fp_data a, const utils::dt_mol_fp_data b) {
        return a.fp[_CFPPopCountIndex] < b.fp[_CFPPopCountIndex];
    });

    return mol_fp_data;
}

void FingerprintStore::_saveFpMolData(std::vector<utils::dt_mol_fp_data> &sortedMolFpData) {
    std::vector<int> localPopCountBins; // list of popcount acquired in current input run
    int previous_popcount = -1;
    int current_popcount = previous_popcount;
    uint64_t *cfp = NULL;  // pointer to cumulative 1D cfp array
    hsize_t size = 0;      // total size of cumulative cfp array
    hsize_t index = 0;    // index in cumulative cfp array for current mol
    int counter = 0;      // current mol data index for assignment
    TMolData *molDataRecords = (TMolData*) malloc(sizeof(TMolData)*sortedMolFpData.size()); // mol data for table
    h5::Group *popcountBinsGroup =  utils::createOrOpenGroup(_root_group, _popCountBinsGroupName);
    h5::Group *popcountGroup = nullptr;

    // open first popcount group and try to get index for already saved data
    popcountGroup = utils::createOrOpenGroup(popcountBinsGroup, _getPopCountGroupName(sortedMolFpData[0].fp[_CFPPopCountIndex]));
    index = utils::sizeOfFingerPrintDatasetInH5(popcountGroup, "cfpData");

    for (const auto &data : sortedMolFpData) {
        current_popcount = data.fp[_CFPPopCountIndex];

        // write fingerprints to HDF5 file when completed for the last popcount bin
        if (current_popcount != previous_popcount && previous_popcount != -1) {
            utils::appendFingerPrintDatasetToH5(popcountGroup, "cfpData", cfp, size, _CFPSize);
            localPopCountBins.push_back(previous_popcount);
            
            // reset data
            delete popcountGroup;
            free(cfp);
            cfp = NULL;
            size = 0;

            // get saved size of fingerprints from file for current popcount bin, used to get index of individual mol
            popcountGroup = utils::createOrOpenGroup(popcountBinsGroup, _getPopCountGroupName(current_popcount));
            index = utils::sizeOfFingerPrintDatasetInH5(popcountGroup, "cfpData");
        }

        // copy incoming fp to outgoing fp
        size = size + (hsize_t) _CFPSize;
        cfp = (uint64_t*) realloc(cfp, sizeof(uint64_t) * size);
        uint64_t *cfp_current = cfp + (size - _CFPSize);
        memcpy(cfp_current, data.fp, sizeof(uint64_t)*_CFPSize);

        // assign mol_data values for table
        strcpy(molDataRecords[counter].id, data.mol_id);
        strcpy(molDataRecords[counter].smiles, data.smiles);
        molDataRecords[counter].fpIndex = index;
        molDataRecords[counter].popcount = data.popCount;

        previous_popcount = current_popcount;
        index += _CFPSize;
        counter += 1;
    }

    // last remaining popcount group
    utils::appendFingerPrintDatasetToH5(popcountGroup, "cfpData", cfp, size, _CFPSize);
    localPopCountBins.push_back(current_popcount);
    delete popcountGroup;
    delete popcountBinsGroup;
    free(cfp);

    for (auto pc : localPopCountBins) {
        if (std::find(_popCountBins.begin(), _popCountBins.end(), pc) == _popCountBins.end() ) _popCountBins.push_back(pc);
    }
    std::sort(_popCountBins.begin(), _popCountBins.end(), [](int a, int b) { return a < b; });
    utils::addDataSetToH5(_root_group, _popCountDatasetName, H5::PredType::NATIVE_INT, _popCountBins.data(), {_popCountBins.size()});

    // update moldata table and free memory
    _molDataTable->append(molDataRecords, counter);
    free(molDataRecords);
}

std::string FingerprintStore::getSmilesFromID(std::string id) {
    TMolData molData;
    _molDataTable->getMolDataFromId(id, &molData);
    if (molData.popcount > 0) {
        return molData.smiles;
    } else {
        return "";
    }
}

void FingerprintStore::updateFingerprints(const std::vector<std::tuple<std::string, std::string, std::string>> fingerprints) {
    // transformed mol and fp data from input
    std::vector<utils::dt_mol_fp_data> all_mol_fp_data = _buildSortedMolFpData(fingerprints);

    // separate new and modified data
    std::vector<utils::dt_mol_fp_data> new_mol_fp_data, modified_mol_fp_data;
    std::vector<TMolData> existingMolsData;
    std::set<int> modifiedPopCounts;
    new_mol_fp_data.reserve(fingerprints.size());
    modified_mol_fp_data.reserve(fingerprints.size());

    for (const auto &molFp: all_mol_fp_data) {
        // now separate out new and modified fingerprint
        TMolData existingMolData; 
        _molDataTable->getMolDataFromId(molFp.mol_id, &existingMolData);
        if (existingMolData.popcount == 0) { // it is a new data
            new_mol_fp_data.push_back(molFp);
            modifiedPopCounts.insert(molFp.popCount);
        } else {
            modified_mol_fp_data.push_back(molFp);
            existingMolsData.push_back(existingMolData);
            modifiedPopCounts.insert(molFp.popCount);
            modifiedPopCounts.insert(existingMolData.popcount);
        }
    }

    // save modified FP and Mol Data
    if (modified_mol_fp_data.size() > 0) _saveModifiedFpMolData(modified_mol_fp_data, existingMolsData);

    // save new FP and Mol data
    if (new_mol_fp_data.size() > 0) _saveFpMolData(new_mol_fp_data);
    
    // redo clustering in pop-count bin fp is modified
    _innerClusteringAgent->performInnerClusteringAfterUpdate(modifiedPopCounts);
    
    // clear memories
    for(auto &molFp: all_mol_fp_data) free(molFp.fp);
    all_mol_fp_data.resize(0);
    new_mol_fp_data.resize(0);
    modified_mol_fp_data.resize(0);
}

void FingerprintStore::_saveModifiedFpMolData(std::vector<utils::dt_mol_fp_data> &sortedMolFpData, std::vector<TMolData> sortedExistingMolsData) {
    h5::Group *popcountBinsGroup =  utils::createOrOpenGroup(_root_group, _popCountBinsGroupName); // open top level popcount group

    for (size_t i=0; i< sortedMolFpData.size(); i++) {
        uint64_t existingCfp[_CFPSize];
        h5::Group *existingPopcountGroup = utils::createOrOpenGroup(popcountBinsGroup, _getPopCountGroupName(sortedExistingMolsData[i].popcount));
        bool success = utils::getFingerprintFromIndex(existingPopcountGroup, "cfpData", sortedExistingMolsData[i].fpIndex, _CFPSize, existingCfp);
        if (success) { // if existing fingerprint is found
            if (sortedMolFpData[i].popCount == sortedExistingMolsData[i].popcount) { // if modified and existing have same popcount
                bool new_fp = false;
                for (size_t f = _molIdOffset; f < (_CFPSize - 1); f++) { // check if incoming fingerprint is different from existing
                    if (sortedMolFpData[i].fp[f] != existingCfp[f]) { // if different, fingerprint need to be updated
                        new_fp = true;
                        break;
                    }
                }
                if(new_fp) { // update fingerprint. do not need to update mol-data table because popcount is unchanged
                    utils::updateFingerprintAtIndex(existingPopcountGroup, "cfpData", sortedExistingMolsData[i].fpIndex, _CFPSize, existingCfp);
                }
            } else { // if popcount is changed
                // first make existing fingerprint unusable by all its value to zero and update in HDF5 file
                for (size_t f = _molIdOffset; f < _CFPSize; f++) {
                    existingCfp[f] = 0;
                }
                utils::updateFingerprintAtIndex(existingPopcountGroup, "cfpData", sortedExistingMolsData[i].fpIndex, _CFPSize, existingCfp);

                // now add the new fingerprint in new popcount-bin in HDF5 file
                h5::Group *newPopcountGroup = utils::createOrOpenGroup(popcountBinsGroup, _getPopCountGroupName(sortedMolFpData[i].popCount));
                hsize_t index = utils::sizeOfFingerPrintDatasetInH5(newPopcountGroup, "cfpData");
                utils::appendFingerPrintDatasetToH5(newPopcountGroup, "cfpData", sortedMolFpData[i].fp, _CFPSize, _CFPSize);

                // now update the data in mol-data-table
                TMolData newMolData;
                strcpy(newMolData.id, sortedMolFpData[i].mol_id);
                strcpy(newMolData.smiles, sortedMolFpData[i].smiles);
                newMolData.fpIndex = index;
                newMolData.popcount = sortedMolFpData[i].popCount;
                _molDataTable->replaceMolDataForId(newMolData.id, &newMolData);

                delete newPopcountGroup;
            }
        }
        delete existingPopcountGroup;
    }
    delete popcountBinsGroup;
}

void FingerprintStore::performInnerClusteringAndWrite() { 
    if (_innerClusteringAgent != nullptr)
        _innerClusteringAgent->performInnerClusteringAndWrite();
    else
        std::cout << "Inner clustering not allowed!"<<std::endl;
};

void FingerprintStore::reDoInnerClusteringInMemory(float threshold) {
    if (_innerClusteringAgent != nullptr)
        _innerClusteringAgent->reDoInnerClusteringInMemory(threshold);
    else
        std::cout << "Inner clustering not allowed!"<<std::endl;
};