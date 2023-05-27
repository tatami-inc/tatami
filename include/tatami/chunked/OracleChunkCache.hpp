#ifndef TATAMI_ORACLE_CHUNK_CACHE_HPP
#define TATAMI_ORACLE_CHUNK_CACHE_HPP

#include <unordered_map>
#include <vector>
#include "../utils/Oracles.hpp"

/**
 * @file OracleChunkCache.hpp
 * @brief Create a oracle-aware cache for matrix chunks.
 */

namespace tatami {

/**
 * @brief Oracle-aware cache for chunks.
 *
 * @tparam Id_ Type of chunk identifier, typically integer.
 * @tparam Index_ Type of row/column index produced by the `Oracle`.
 * @tparam ChunkContents_ Class that contains the contents of a single chunk.
 *
 * Implement an oracle-aware cache for chunks.
 * This is typically used for `Matrix` representations where the data is costly to load (e.g., from file) and an oracle is provided to predict future accesses.
 * In such cases, the appropriate chunks of data can be loaded and cached such that all future requests will just fetch the cached data.
 */
template<typename Id_, typename Index_, class ChunkContents_> 
class OracleChunkCache {
    OracleStream<Index_> prediction_stream;
    std::vector<std::pair<Id_, Index_> > predictions_made;
    size_t predictions_fulfilled = 0;
    size_t max_predictions;

    size_t max_chunks;
    std::unordered_map<Id_, Index_> cache_exists, next_cache_exists;
    std::vector<ChunkContents_> cache_data, next_cache_data;
    std::vector<std::pair<Id_, Index_> > chunks_in_need; 

public:
    /**
     * @param oracle Pointer to an `Oracle` to be used for predictions.
     * @param per_iteration Maximum number of predictions to make per iteration.
     * @param num_chunks Maximum number of chunks to store.
     */
    OracleChunkCache(std::unique_ptr<Oracle<Index_> > oracle, size_t per_iteration, size_t num_chunks) :
        prediction_stream(std::move(oracle)), 
        max_predictions(per_iteration),
        max_chunks(num_chunks), 
        cache_data(max_chunks),
        next_cache_data(max_chunks)
    {
        chunks_in_need.reserve(max_chunks);
    } 

    /**
     * @cond
     */
    // For testing only.
    OracleChunkCache() = default;
    /**
     * @endcond
     */

public:
    /**
     * Fetch the next chunk according to the stream of predictions provided by the `Oracle`.
     *
     * @tparam Ifunction_ Function to identify the chunk containing each predicted row/column.
     * @tparam Sfunction_ Function to swap two chunks' contents.
     * @tparam Rfunction_ Function to determine whether a chunk object is a non-mock instance.
     * @tparam Afunction_ Function to allocate memory to a mock chunk for storage of chunk contents.
     * @tparam Pfunction_ Function to populate zero, one or more chunks with their contents.
     *
     * @param identify Function that accepts an `i`, an `Index_` containing the predicted row/column index.
     * This should return a pair containing (1) the identifier of the chunk containing `i`, and (2) the index of row/column `i` inside that chunk.
     * @param swap Function that accepts two `ChunkContents_&` and swaps their contents.
     * @param ready Function that accepts a `const ChunkContents_&` and returns a boolean indicating whether it has already been allocated.
     * This should return `true` for objects that have been used in `allocate()`, and `false` otherwise.
     * @param allocate Function that accepts a single default-initialized `ChunkContents_` object,
     * and allocates sufficient memory to it in order to hold a chunk's contents when used in `populate()`.
     * @param populate Function that accepts:
     *
     * - `chunks_in_needs`, a `const std::vector<std::pair<Id_, Index_> >&` specifying the identity of the chunks to be filled in the first `Id_` element of each pair.
     *   The second `Index_` element specifies the index of the `ChunkContents_` in `chunk_data` in which to store the contents of each chunk.
     * - `chunk_data`, a `std::vector<ChunkContents>&` containing the cached chunk contents.
     *
     * This function should iterate over the `chunks_in_need` and populates the corresponding entries in `chunk_data`.
     *
     * @return Pair containing (1) a pointer to a chunk's contents and (2) the index of the next predicted row/column inside the retrieved chunk.
     */
    template<class Ifunction_, class Sfunction_, class Rfunction_, class Afunction_, class Pfunction_>
    std::pair<const ChunkContents_*, Index_> next_chunk(Ifunction_ identify, Sfunction_ swap, Rfunction_ ready, Afunction_ allocate, Pfunction_ populate) {
        if (predictions_made.size() > predictions_fulfilled) {
            const auto& chosen = predictions_made[predictions_fulfilled++];
            return std::make_pair(cache_data.data() + chosen.first, chosen.second);
        }

        next_cache_exists.clear();
        chunks_in_need.clear();
        size_t used = 0;

        predictions_made.clear();
        predictions_made.reserve(max_predictions);
        for (size_t p = 0; p < max_predictions; ++p) {
            Index_ current;
            if (!prediction_stream.next(current)) {
                break;
            }

            auto chunk_id = identify(current);
            auto curchunk = chunk_id.first;
            auto curindex = chunk_id.second;

            auto it = next_cache_exists.find(curchunk);
            if (it == next_cache_exists.end()) {
                if (used == max_chunks) {
                    prediction_stream.back();
                    break;
                }

                next_cache_exists[curchunk] = used;
                predictions_made.emplace_back(used, curindex);

                auto it2 = cache_exists.find(curchunk);
                if (it2 != cache_exists.end()) {
                    swap(next_cache_data[used], cache_data[it2->second]);
                } else {
                    chunks_in_need.emplace_back(curchunk, used);
                }

                ++used;
            } else {
                predictions_made.emplace_back(it->second, curindex);
            }
        }

        /**
         * Doing a linear scan across chunks to find the allocated but unused cache
         * elements. This is the simplest and safest approach; trying to keep
         * track of the unused caches would require a scan to fill a map/set
         * anyway, and deleting an iterator from cache_exists only works if
         * cache_exists actually contains all cache elements, which it might
         * not be if we reached max_predictions without filling up the cache.
         *
         * In any case, a linear scan should be pretty good for consecutive
         * access; the first cache elements would be the oldest, so there
         * wouldn't be any wasted iterations to find available cache elements.
         */
        size_t search = 0;
        for (const auto& c : chunks_in_need) {
            if (!ready(next_cache_data[c.second])) {
                while (search < cache_data.size() && !ready(cache_data[search])) {
                    ++search;
                }
                if (search < cache_data.size()) {
                    swap(next_cache_data[c.second], cache_data[search]);
                    ++search;
                } else {
                    // This should be called no more than 'max_chunks' times across the lifetime of this Cache object. 
                    // At any given point in time, allocated chunks will be interspersed between 'cache_data' and 'next_cache_data', 
                    // so either 'next_cache_data[c.second]' is already ready, or we'll find a ready chunk from the linear scan through 'cache_data'.
                    allocate(next_cache_data[c.second]);
                }
            }
        }

        populate(chunks_in_need, next_cache_data);

        cache_data.swap(next_cache_data);
        cache_exists.swap(next_cache_exists);
        predictions_fulfilled = 1; // well, because we just used one.
        const auto& chosen = predictions_made.front(); // assuming at least one prediction was made, otherwise, why was this function even called?
        return std::make_pair(cache_data.data() + chosen.first, chosen.second);
    }
};

}

#endif
