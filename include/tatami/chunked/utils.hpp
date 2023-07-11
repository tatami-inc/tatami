#ifndef TATAMI_CHUNK_UTILS_HPP
#define TATAMI_CHUNK_UTILS_HPP

#include <memory>
#include "OracleChunkCache.hpp"
#include "LruChunkCache.hpp"

/**
 * @file utils.hp
 * @brief Utilities for chunked extraction.
 */

namespace tatami {

/**
 * @brief Typical options for chunked extraction.
 *
 * This is relevant to any matrix representation that stores its data in regular rectangular chunks.
 * We define a "chunk set" as the set of chunks required to read a row/column (or a contiguous block/indexed subset thereof) during a `tatami::DenseExtractor::fetch()` or `tatami::SparseExtractor::Fetch()` call.
 * We aim to cache one or more complete chunk sets; this means that we can re-use the cached chunks when adjacent rows/columns are requested, rather than re-reading them (e.g., from disk).
 */
struct TypicalChunkCacheOptions {
    /**
     * Size of the in-memory cache in bytes.
     * Larger caches improve access speed at the cost of memory usage.
     * Small values may be ignored if `require_minimum_cache` is `true`.
     */
    size_t maximum_cache_size = 100000000;

    /**
     * Whether to automatically enforce a minimum size for the cache, regardless of `maximum_cache_size`.
     * This minimum is chosen to ensure that a single chunk set can be retained in memory,
     * so that the same chunks are not repeatedly re-read when iterating over consecutive rows/columns of the matrix.
     */
    bool require_minimum_cache = true;
};

/**
 * @brief Workspace for typical chunk extraction.
 *
 * Implements a workspace to initialize the chunk caches (i.e., the `LruChunkCache` and `OracleChunkCache`) for extraction from a chunked matrix representation.
 * This is intended to be part of an `Extractor` class to allow extraction to switch between different caches (e.g., when `ExtractorBase::set_oracle()` is called).
 * It also handles the calculation of various cache size statistics.
 *
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam ChunkSet_ Class that contains a single chunk set.
 */
template<typename Index_, class ChunkSet_>
struct TypicalChunkCacheWorkspace {
    /**
     * Default constructor.
     */
    TypicalChunkCacheWorkspace() = default;

    /**
     * @param primary_length Length of the primary dimension of each chunk set.
     * The primary dimension contains the elements to be extracted from each cached chunk set.
     * For example, if we were iterating through rows of a matrix, `primary_length` would be the number of rows spanned by each chunk set.
     * @param secondary_length Length of the secondary dimension of each chunk set.
     * This is the dimension that is not the primary, e.g., if we were iterating through matrix rows, the `secondary_length` would be the number of columns spanned by each chunk set.
     * @param cache_size_in_elements Total size of the cache in terms of the number of elements.
     * This is usually derived from `TypicalChunkCacheOptions::maximum_cache_size` and the size of each element.
     * @param require_minimum_cache Whether to enforce a minimum size of the cache for efficient extraction, see `TypicalChunkCacheOptions` for details.
     */
    TypicalChunkCacheWorkspace(Index_ primary_length, Index_ secondary_length, size_t cache_size_in_elements, bool require_minimum_cache) {
        chunk_set_size_in_elements = static_cast<size_t>(primary_length) * static_cast<size_t>(secondary_length);
        num_chunk_sets_in_cache = static_cast<double>(cache_size_in_elements) / chunk_set_size_in_elements;

        if (num_chunk_sets_in_cache == 0 && require_minimum_cache) {
            num_chunk_sets_in_cache = 1;
        }

        if (num_chunk_sets_in_cache) {
            lru_cache.reset(new LruChunkCache<Index_, ChunkSet_>(num_chunk_sets_in_cache));
        }
    }

public:
    /**
     * Length of the primary dimension of each chunk set.
     */
    Index_ primary_length;

    /**
     * Size of each chunk set, in terms of the number of elements.
     */
    size_t chunk_set_size_in_elements;

    /**
     * Number of chunk sets that can fit in the cache.
     */
    size_t num_chunk_sets_in_cache;

    /**
     * Cache of least recently used chunk sets.
     * This is only allocated if `num_chunk_sets_in_cache` is positive and `oracle_cache` is NULL.
     */
    std::unique_ptr<LruChunkCache<Index_, ChunkSet_> > lru_cache;

    /**
     * Cache of to-be-used chunk sets, based on an `Oracle`'s predictions.
     * This may be NULL, see `set_oracle()` for more details.
     */
    std::unique_ptr<OracleChunkCache<Index_, Index_, ChunkSet_> > oracle_cache; 

public:
    /**
     * Set up the oracle cache.
     * This will only have an effect if the number of chunk sets in the cache is greater than 1,
     * otherwise the predictions have no effect on the choice of chunk set to retain.
     * Callers should check for a non-NULL `oracle_cache` before attempting to use it,
     * otherwise they should use the `lru_cache` (provided `num_chunk_sets_in_cache > 0`).
     *
     * @param o Oracle to provide predictions for subsequent accesses on the primary dimension.
     */
    void set_oracle(std::unique_ptr<Oracle<Index_> > o) {
        // The oracle won't have any effect if fewer than one chunk set can be cached.
        if (num_chunk_sets_in_cache > 1) {
            size_t max_predictions = static_cast<size_t>(num_chunk_sets_in_cache) * primary_length * 2; // double the cache size, basically.
            oracle_cache.reset(new OracleChunkCache<Index_, Index_, ChunkSet_>(std::move(o), max_predictions, num_chunk_sets_in_cache));
            lru_cache.reset();
        }
    }
};

}

#endif
