/**
 * @brief This class represents a search engine for molecular fingerprints.
 * @author Rajendra Kumar
 */

#include <algorithm>
#include <iostream>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>

#include "searchEngine.h"

namespace py = pybind11;

FPSearchEngine::FPSearchEngine(const std::string& filename, std::string mode) {
    _fpStore = new FingerprintStore(filename);
    _molIdMaxLength = _fpStore->_molIdMaxLength;
    _fpSize = _fpStore->_fpSize;
    _CFPSize = _fpStore->_CFPSize;
    _molIdOffset =  _fpStore->_molIdOffset;
    _CFPPopCountIndex = _fpStore->_CFPPopCountIndex;
    _fpEndIndex = _molIdOffset + _fpSize;
    _mode = mode;

    if (mode == "disk") {
        _normal_search = std::bind(&FPSearchEngine::_normal_search_disk, this, std::placeholders::_1, std::placeholders::_2, 
                                    std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
        _fpStore->loadDataInMemory(true); // load only cluster fps in memory for disk-based search
    } else {
        _normal_search = std::bind(&FPSearchEngine::_normal_search_memory, this, std::placeholders::_1, std::placeholders::_2, 
                                    std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
        _fpStore->loadDataInMemory(); // load all fps in memory for memory-based search
    }
    

    for (int i = 0; i <= _fpSize*128; i++) {
        _div_lookup_table[i+1] = 1.0/(i+1);
    }
    py::gil_scoped_acquire acquire;
       
}

FPSearchEngine::~FPSearchEngine() {
    if (_fpStore!= nullptr)
        delete _fpStore;
}

void FPSearchEngine::close() {
    _fpStore->freeMemory();
    _fpStore->close();
    delete _fpStore;
    _fpStore = nullptr;
}

uint64_t* FPSearchEngine::prepareQuery(const std::string& query) {
    uint64_t *queryCFp = new uint64_t[_fpStore->_CFPSize];
    utils::BitStrToCompactFPArray(query, "_", queryCFp, _fpStore->_molIdOffset, _fpStore->_CFPSize);
    return queryCFp;
}

std::vector<utils::dt_inner_clusters_fingerprints_maxscore> FPSearchEngine::filterPopcountBins(uint64_t queryPopcount, float threshold) {
    std::vector<utils::dt_inner_clusters_fingerprints_maxscore> filteredPopCountBinsWithMaxScore;
    for (size_t i=0; i < _fpStore->_popCountBins.size(); i++) {
        int popCount = _fpStore->_popCountBins[i];
        float maxScore = (float) std::min(popCount, (int) queryPopcount) / std::max(popCount, (int)queryPopcount);
        if (maxScore >= threshold ) {
            utils::dt_inner_clusters_fingerprints *inner_clusters_fingerprints = _fpStore->_fp_inner_clusters_by_popcount + i;
            filteredPopCountBinsWithMaxScore.push_back({*inner_clusters_fingerprints, maxScore});
        }
    }

    // Sort by score
    std::sort(filteredPopCountBinsWithMaxScore.begin(), filteredPopCountBinsWithMaxScore.end(), [](const utils::dt_inner_clusters_fingerprints_maxscore &a, const utils::dt_inner_clusters_fingerprints_maxscore &b) {
        return a.score > b.score;
    });
    return filteredPopCountBinsWithMaxScore;
}

std::vector<std::tuple<std::string, float>> FPSearchEngine::search(const std::string& query, float threshold, int limits) {
    uint64_t *queryCFp = prepareQuery(query);
    std::vector<utils::dt_inner_clusters_fingerprints_maxscore> filteredPopCountBinsWithMaxScore = filterPopcountBins(queryCFp[_fpStore->_CFPPopCountIndex], threshold);
    std::vector<std::tuple<std::string, float>> results;

    _normal_search(filteredPopCountBinsWithMaxScore, queryCFp, threshold, limits, results);

    //sort results by score
    std::sort(results.begin(), results.end(), [](const std::tuple<std::string, float>& a, const std::tuple<std::string, float>& b) {
        return std::get<1>(a) > std::get<1>(b);
    });

    delete[] queryCFp;
    py::gil_scoped_acquire acquire;
    return results;
}

void FPSearchEngine::_normal_search_memory(std::vector<utils::dt_inner_clusters_fingerprints_maxscore> popCountBinsWithMaxScore, 
                                    uint64_t *queryCFp, float threshold, int limits,
                                    std::vector<std::tuple<std::string, float>> &results) {
    uint64_t commonPopCountThreshold = 0;
    float coeff, max_coeff = 0;
    uint64_t common_popcnt = 0;
    
    for( auto inner_clusters_fingerprints_maxScore : popCountBinsWithMaxScore ) {
        float maxScore = inner_clusters_fingerprints_maxScore.score;

        // check if hits count is equal to the required limit
        if (max_coeff >= maxScore) {
            int hits = 0;
            for (auto &r : results) {
                if (std::get<1>(r) >= maxScore) {
                    hits++;
                }
            }
            if (hits >= limits) break;
        }

        utils::dt_inner_clusters_fingerprints inner_clusters_fingerprints = inner_clusters_fingerprints_maxScore.inner_clusters_fingerprints;
        commonPopCountThreshold = (uint64_t) ceil(threshold * std::max(inner_clusters_fingerprints.popCount, (int)queryCFp[_CFPPopCountIndex])); // get the common popcount threshold for this bin

        uint64_t *clusterFp_ptr = inner_clusters_fingerprints.clusterFp;
        uint64_t *fp_ptr = inner_clusters_fingerprints.fp;
        uint64_t inner_start = 0;
        for(size_t cid=0; cid < inner_clusters_fingerprints.num_clusters; cid++, clusterFp_ptr += _CFPSize) {
            common_popcnt = bitwise_and_popcount(clusterFp_ptr+_molIdOffset, queryCFp+_molIdOffset, _fpSize);
            /*for (auto j = _molIdOffset; j < _fpEndIndex; j++) {
                common_popcnt += popcntll(clusterFp_ptr[j] & queryCFp[j]);
            }*/

            if (common_popcnt >= commonPopCountThreshold) {
                uint64_t inner_end = clusterFp_ptr[0];
                for (auto i = inner_start; i < inner_end; i+=_CFPSize, fp_ptr += _CFPSize) {

                    common_popcnt = bitwise_and_popcount(fp_ptr+_molIdOffset, queryCFp+_molIdOffset, _fpSize);;
                    /*for (auto j = _molIdOffset; j < _fpEndIndex; j++) {
                        common_popcnt += popcntll(fp_ptr[j] & queryCFp[j]);
                    }*/

                    if (common_popcnt >= commonPopCountThreshold) {
                        coeff =  TanimotoCoeff(common_popcnt, queryCFp[_CFPPopCountIndex], fp_ptr[_CFPPopCountIndex], _div_lookup_table);
                        if (coeff >= threshold) {
                            results.push_back(std::make_tuple(utils::getMolIdFromCompactFPArray(fp_ptr, _molIdMaxLength), coeff));
                            if (coeff > max_coeff) max_coeff = coeff;
                        }
                    }
                }
            } else {
                fp_ptr += clusterFp_ptr[0] - inner_start;
            }
            inner_start = clusterFp_ptr[0];
        }
    }
}

void FPSearchEngine::_normal_search_disk(std::vector<utils::dt_inner_clusters_fingerprints_maxscore> popCountBinsWithMaxScore, 
                                    uint64_t *queryCFp, float threshold, int limits,
                                    std::vector<std::tuple<std::string, float>> &results) {
    uint64_t commonPopCountThreshold = 0;
    float coeff, max_coeff = 0;
    uint64_t common_popcnt = 0;
    
    for( auto inner_clusters_fingerprints_maxScore : popCountBinsWithMaxScore ) {
        float maxScore = inner_clusters_fingerprints_maxScore.score;

        // check if hits count is equal to the required limit
        if (max_coeff >= maxScore) {
            int hits = 0;
            for (auto &r : results) {
                if (std::get<1>(r) >= maxScore) {
                    hits++;
                }
            }
            if (hits >= limits) break;
        }

        utils::dt_inner_clusters_fingerprints inner_clusters_fingerprints = inner_clusters_fingerprints_maxScore.inner_clusters_fingerprints;
        commonPopCountThreshold = (uint64_t) ceil(threshold * std::max(inner_clusters_fingerprints.popCount, (int)queryCFp[_CFPPopCountIndex])); // get the common popcount threshold for this bin

        uint64_t *clusterFp_ptr = inner_clusters_fingerprints.clusterFp;
        uint64_t inner_start = 0;
        for(size_t cid=0; cid < inner_clusters_fingerprints.num_clusters; cid++, clusterFp_ptr += _CFPSize) {
            common_popcnt = bitwise_and_popcount(clusterFp_ptr+_molIdOffset, queryCFp+_molIdOffset, _fpSize);
            /*for (auto j = _molIdOffset; j < _fpEndIndex; j++) {
                common_popcnt += popcntll(clusterFp_ptr[j] & queryCFp[j]);
            }*/

            if (common_popcnt >= commonPopCountThreshold) {
                uint64_t inner_end = clusterFp_ptr[0];
                uint64_t *fp_ptr = _fpStore->getFPsForCluster(inner_clusters_fingerprints.popCount, inner_start, inner_end); // read fps for this cluster from disk
                for (auto i = inner_start; i < inner_end; i+=_CFPSize, fp_ptr += _CFPSize) {

                    common_popcnt = bitwise_and_popcount(fp_ptr+_molIdOffset, queryCFp+_molIdOffset, _fpSize);;
                    /*for (auto j = _molIdOffset; j < _fpEndIndex; j++) {
                        common_popcnt += popcntll(fp_ptr[j] & queryCFp[j]);
                    }*/

                    if (common_popcnt >= commonPopCountThreshold) {
                        coeff =  TanimotoCoeff(common_popcnt, queryCFp[_CFPPopCountIndex], fp_ptr[_CFPPopCountIndex], _div_lookup_table);
                        if (coeff >= threshold) {
                            results.push_back(std::make_tuple(utils::getMolIdFromCompactFPArray(fp_ptr, _molIdMaxLength), coeff));
                            if (coeff > max_coeff) max_coeff = coeff;
                        }
                    }
                }
                free(fp_ptr - (inner_end - inner_start)); // free memory allocated for fps read from disk
            }
            inner_start = clusterFp_ptr[0];
        }
    }
}

std::vector<std::vector<std::tuple<std::string, float>>> FPSearchEngine::batchSearch(const std::vector<std::string>& queries, float threshold, int limits) {
    // build query Compact Fingerprints
    std::vector<uint64_t*> queriesCFp;
    queriesCFp.reserve(queries.size());
    for (size_t i=0; i < queries.size(); i++) {
        uint64_t *queryCFp = prepareQuery(queries[i]);
        queryCFp[0] = i;
        queriesCFp.push_back(queryCFp);
    }

    //sort queries by popcount
    std::sort(queriesCFp.begin(), queriesCFp.end(), [this](uint64_t *a, uint64_t *b) {
        return a[_CFPPopCountIndex] > b[_CFPPopCountIndex];
    });

    //prepare batch data by sorting queries by popcount
    uint64_t previous_popcount = 0;
    std::vector<utils::dt_batch_data> batch_data;
    for (size_t i=0; i < queriesCFp.size(); i++) {
        if (queriesCFp[i][_CFPPopCountIndex] != previous_popcount) { // create new batch data
            utils::dt_batch_data new_batch_data;         
            std::vector<utils::dt_inner_clusters_fingerprints_maxscore> filteredPopCountBinsWithMaxScore = filterPopcountBins(queriesCFp[i][_CFPPopCountIndex], threshold);
            new_batch_data.filteredPopCountBinsWithMaxScore = filteredPopCountBinsWithMaxScore;
            new_batch_data.popCount = queriesCFp[i][_CFPPopCountIndex];
            batch_data.push_back(new_batch_data);
        }
        // create new query data
        utils::dt_batch_query_data new_query_data;
        new_query_data.cfp = queriesCFp[i];

        // add new query data to batch data
        utils::dt_batch_query_data *qdata = batch_data.at(batch_data.size()-1).qdata;
        int qsize = batch_data.at(batch_data.size()-1).qsize;
        qdata = (utils::dt_batch_query_data*) realloc(qdata, sizeof(utils::dt_batch_query_data) * (qsize + 1)); // memory allocation
        qdata[qsize] = new_query_data;
        batch_data.at(batch_data.size()-1).qdata = qdata;
        batch_data.at(batch_data.size()-1).qsize = qsize + 1;

        previous_popcount = queriesCFp[i][_CFPPopCountIndex];
    }

    std::vector<std::vector<std::tuple<std::string, float>>> finalResults;
    finalResults.resize(queries.size());
    for(auto bdata : batch_data) {
        // perform batch search
        _batch_search(bdata, threshold, limits);

        // process and assign results
        for(auto q=0; q < bdata.qsize; q++) {
            // prepare results
            std::vector<std::tuple<std::string, float>> results;
            results.reserve(bdata.qdata[q].results_size);
            for(auto r=0; r < bdata.qdata[q].results_size; r++) {
                results.push_back(std::make_tuple(*bdata.qdata[q].results[r].id, bdata.qdata[q].results[r].score));
            }
            
            // sort results by score
            std::sort(results.begin(), results.end(), [](const std::tuple<std::string, float>& a, const std::tuple<std::string, float>& b) {
                return std::get<1>(a) > std::get<1>(b);
            });

            // assign results to finalResults
            finalResults[bdata.qdata[q].cfp[0]] = results;
        }
    }

    //release memory for batch_data
    for(auto bdata : batch_data) {
        for(auto q=0; q < bdata.qsize; q++) {
            for(auto r=0; r < bdata.qdata[q].results_size; r++)
                delete bdata.qdata[q].results[r].id; // id created through new operator
            free(bdata.qdata[q].results); // by realloc
        }
        free(bdata.qdata); // by realloc
    }

    // remove memory for queryCFp
    for (size_t i=0; i < queriesCFp.size(); i++) {
        delete[] queriesCFp[i];
    }

    py::gil_scoped_acquire acquire;
    return finalResults;
}

void FPSearchEngine::_batch_search(utils::dt_batch_data &batch_data, float threshold, int limits) {
    uint64_t commonPopCountThreshold = 0;
    uint64_t queriesPopcount = batch_data.popCount;
    float coeff = 0;
    uint64_t max_common_popcnt = 0;
    uint64_t common_popcnt = 0;
    bool all_done = false;
    utils::dt_batch_query_data *query_data = batch_data.qdata;
    
    for( auto &inner_clusters_fingerprints_maxScore : batch_data.filteredPopCountBinsWithMaxScore) {
        float maxScore = inner_clusters_fingerprints_maxScore.score;

        // check if all queries are done within requested limits
        all_done = true;
        query_data = batch_data.qdata;
        for(auto q=0; q < batch_data.qsize; q++, query_data++) {
            if (!query_data->done && query_data->max_coeff >= maxScore) {
                int hits = 0;
                for(auto r=0; r < query_data->results_size; r++) {
                    if (query_data->results[r].score >= maxScore) hits++;
                }
                if (hits >= limits) query_data->done = true;
            }
            all_done &= query_data->done;
        }
        if (all_done) break;

        utils::dt_inner_clusters_fingerprints inner_clusters_fingerprints = inner_clusters_fingerprints_maxScore.inner_clusters_fingerprints;
        commonPopCountThreshold = (uint64_t) ceil(threshold * std::max(inner_clusters_fingerprints.popCount, (int)queriesPopcount)); // get the common popcount threshold for this bin

        uint64_t *clusterFp_ptr = inner_clusters_fingerprints.clusterFp;
        uint64_t *fp_ptr = inner_clusters_fingerprints.fp;
        uint64_t inner_start = 0;
        for(size_t cid=0; cid < inner_clusters_fingerprints.num_clusters; cid++, clusterFp_ptr += _CFPSize) { //searching through clusters

            // check if any query need to be search in current cluster
            max_common_popcnt = 0;
            query_data = batch_data.qdata;
            for(auto q=0; q < batch_data.qsize; q++, query_data++) {
                if (!query_data->done) {
                    common_popcnt = bitwise_and_popcount(clusterFp_ptr + _molIdOffset, query_data->cfp + _molIdOffset, _fpSize);
                    //for (auto j = _molIdOffset; j < _fpEndIndex; j++) {
                    //    common_popcnt += popcntll(clusterFp_ptr[j] & query_data->cfp[j]);
                    //}
                    if (common_popcnt > max_common_popcnt) max_common_popcnt = common_popcnt;
                }
            }

            // queries are searched in current cluster
            if (max_common_popcnt >= commonPopCountThreshold) { // potentially hit could be found in current cluster for at least one query
                uint64_t inner_end = clusterFp_ptr[0];
                for (auto i = inner_start; i < inner_end; i+=_CFPSize, fp_ptr += _CFPSize) {

                    query_data = batch_data.qdata;
                    for(auto q=0; q < batch_data.qsize; q++, query_data++) { // similarity of each query with each fingerprint of current cluster
                        if (!query_data->done) {
                            common_popcnt = bitwise_and_popcount(fp_ptr + _molIdOffset, query_data->cfp + _molIdOffset, _fpSize);
                            //for (auto j = _molIdOffset; j < _fpEndIndex; j++) {
                            //    common_popcnt += popcntll(fp_ptr[j] & query_data->cfp[j]);
                            //}

                            if (common_popcnt >= commonPopCountThreshold) { // potential hit found
                                coeff =  TanimotoCoeff(common_popcnt, query_data->cfp[_CFPPopCountIndex], fp_ptr[_CFPPopCountIndex], _div_lookup_table);

                                if (coeff >= threshold) { // exact hit found, add to results
                                    utils::dt_result *results = query_data->results;
                                    results = (utils::dt_result*) realloc(results, sizeof(utils::dt_result) * (query_data->results_size + 1));
                                    std::string *id = new std::string(utils::getMolIdFromCompactFPArray(fp_ptr, _molIdMaxLength));
                                    results[query_data->results_size] = {id, coeff};
                                    query_data->results = results;
                                    query_data->results_size += 1;
                                    if (coeff > query_data->max_coeff) query_data->max_coeff = coeff;
                                }
                            }
                        }
                    }
                }
            } else {
                fp_ptr += clusterFp_ptr[0] - inner_start;
            }
            inner_start = clusterFp_ptr[0];
        }
    }
}
