/**
 * @brief These are general utility functions used in this package.
 * @author Rajendra Kumar
 */

#include <iostream>
#include <string>
#include <cstring>
#include <cstdint>
#include <cstdlib>
#include <time.h>

#include "lib.h"
#include "popcnt.hpp"

int utils::bitwiseAndPopcount(std::vector<uint64_t> a, std::vector<uint64_t> b) {
    int count = -1;
    if (a.size() == b.size() && (a.size()*64 % 512) == 0) {
        count = bitwise_and_popcount(a.data(), b.data(), a.size());
    }
    return count;
}

void utils::BitStrToCompactFPArray(const std::string &bit_string, std::string mol_id, uint64_t *cfp_arr, int id_offset, int total_size) {
    mol_id.reserve(8*id_offset); // reserve full space
    memcpy(cfp_arr, mol_id.data(), sizeof(uint64_t) * id_offset); // copy mol_id to cfp_arr

    uint64_t popcount = 0;
    for (size_t i =0; i < bit_string.length(); i += 64) {
        int index = (i/64) + id_offset;
        uint64_t *fp_int = &cfp_arr[index];
        *fp_int = std::stoull(bit_string.substr(i, 64), 0, 2);
        popcount += popcntll(*fp_int);
    }

    cfp_arr[total_size-1] = popcount;
}

std::string utils::getMolIdFromCompactFPArray(uint64_t *cfp_arr, int max_mol_id_size) {
    char mol_id[max_mol_id_size+1];
    memcpy(mol_id, cfp_arr, max_mol_id_size + 1);
    return mol_id;
}

std::vector<uint64_t> utils::getCompactFingerPrintArray(std::string &bit_string) {
    std::vector<uint64_t> cfpArray;
    int cfp_size = (bit_string.size()/64)+1;
    cfpArray.resize(cfp_size);

    int popcount = 0;
    for (size_t i =0; i < bit_string.length(); i += 64) {
        uint64_t cfp_unit = std::stoull(bit_string.substr(i, 64), 0, 2);
        cfpArray[i/64] = cfp_unit;
        popcount += popcntll(cfp_unit);
    }
    cfpArray[cfp_size-1] = popcount;
    return cfpArray;
}

double utils::get_posix_clock_time ()
{
    struct timespec ts;

    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0)
        return (ts.tv_sec * 1000000000 + (double)ts.tv_nsec)/1000000.0;
    else
        return 0;
}

void utils::free_dt_inner_clusters_fingerprints(dt_inner_clusters_fingerprints &inner_clusters_fingerprints) {
    if(inner_clusters_fingerprints.clusterFp != nullptr)
        free(inner_clusters_fingerprints.clusterFp);
    if(inner_clusters_fingerprints.fp != nullptr)
        free(inner_clusters_fingerprints.fp);
};
