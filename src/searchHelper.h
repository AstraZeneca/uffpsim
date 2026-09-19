#pragma once
#include <vector>
#include <memory>
#include <cstdlib>
#include <cstdint>
#include <array>
#include <algorithm>
#include <cmath>

#include "lib.h"

// Custom stateless deleter for std::unique_ptr
struct FreeDeleter {
    void operator()(void* ptr) const { std::free(ptr); }
};

// Define a struct to hold the managed memory and the index mapping together
struct PopcntCache {
    std::unique_ptr<uint64_t[], FreeDeleter> memory_owner;
    std::vector<uint64_t*> view;
};


// The allocation function
PopcntCache create_popcnt_clusters_cache(const std::vector<utils::dt_inner_clusters_fingerprints_maxscore>& popCountBinsWithMaxScore) {
    // 1. Calculate total elements needed
    size_t total_elements = 0;
    for (const auto& bin : popCountBinsWithMaxScore) {
        total_elements += bin.inner_clusters_fingerprints.num_clusters;
    }

    // 2. Allocate zeroed memory via calloc
    std::unique_ptr<uint64_t[], FreeDeleter> flat_cache(
        static_cast<uint64_t*>(std::calloc(total_elements, sizeof(uint64_t)))
    );

    // 3. Build the 2D row indexer
    std::vector<uint64_t*> index_view;
    index_view.reserve(popCountBinsWithMaxScore.size());

    uint64_t* current_row = flat_cache.get();
    for (size_t i = 0; i < popCountBinsWithMaxScore.size(); ++i) {
        index_view.push_back(current_row);
        current_row += popCountBinsWithMaxScore[i].inner_clusters_fingerprints.num_clusters;
    }

    // Returns via NRVO (Zero-copy move)
    return { std::move(flat_cache), std::move(index_view) };
}

// Define a lightweight structure to hold the stack data and its valid count
struct ThresholdBuffer {
    std::array<float, 11> data; // Maximum possible elements (10 steps + 1 extra threshold)
    size_t size = 0;            // Keeps track of how many elements were actually written
};

// Inline function for header placement and zero runtime overhead
inline ThresholdBuffer generate_threshold_steps(float threshold) {
    constexpr float kStep = 0.1f;
    ThresholdBuffer buffer;

    // 1. Calculate the starting bound
    const int min_tenth = std::max(0, static_cast<int>(std::ceil(threshold * 10.0f)));
    
    // 2. Unrollable loop filling the stack array
    for (int tenth = 9; tenth >= min_tenth; --tenth) {
        buffer.data[buffer.size++] = static_cast<float>(tenth) * kStep;
    }
    
    // 3. Conditional boundary check
    if (buffer.size > 0 && buffer.data[buffer.size - 1] > threshold) {
        buffer.data[buffer.size++] = threshold;
    }

    return buffer; // Optimized by compiler via NRVO (no-copy stack transfer)
}


// Define a struct to hold the managed memory and the 2D view mapping together
struct ClustersDoneCache {
    std::unique_ptr<uint8_t[], FreeDeleter> memory_owner;
    std::vector<uint8_t*> view;
};

// Inline function for easy header placement and zero runtime overhead
inline ClustersDoneCache create_clusters_done_cache(const std::vector<utils::dt_inner_clusters_fingerprints_maxscore>& popCountBinsWithMaxScore) {
    // 1. Calculate total elements needed (1 byte per element)
    size_t total_elements = 0;
    for (const auto& bin : popCountBinsWithMaxScore) {
        total_elements += bin.inner_clusters_fingerprints.num_clusters;
    }

    // 2. Allocate zeroed memory via calloc (guarantees all elements are 0)
    std::unique_ptr<uint8_t[], FreeDeleter> flat_cache(
        static_cast<uint8_t*>(std::calloc(total_elements, sizeof(uint8_t)))
    );

    // 3. Build the 2D row indexer for standard [i][j] lookup
    std::vector<uint8_t*> index_view;
    index_view.reserve(popCountBinsWithMaxScore.size());

    uint8_t* current_row = flat_cache.get();
    for (size_t i = 0; i < popCountBinsWithMaxScore.size(); ++i) {
        index_view.push_back(current_row);
        current_row += popCountBinsWithMaxScore[i].inner_clusters_fingerprints.num_clusters;
    }

    // Returns via NRVO (Zero-copy move)
    return { std::move(flat_cache), std::move(index_view) };
}
