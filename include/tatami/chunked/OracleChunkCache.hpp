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
public:
    /**
     * Contents of each chunk, as returned by `ChunkManager_::mock()`.
     */
    typedef decltype(std::declval<ChunkManager_>().create()) ChunkContents;

private:
    std::unordered_map<Id_, Index_> cache_exists, next_cache_exists;
    std::vector<ChunkContents_> cache_data, next_cache_data;
    std::vector<std::pair<Id, Index_> > chunks_in_need; 

    OracleStream<Index_> prediction_stream;
    std::vector<std::pair<Id_, Index_> > predictions_made;
    size_t predictions_fulfilled = 0;
    size_t max_predictions;
    size_t max_chunks;

public:
    /**
     * @param oracle Pointer to an `Oracle` to be used for predictions.
     * @param per_iteration Maximum number of predictions to make per iteration.
     * @param num_chunks Maximum number of chunks to store.
     */
    OracleChunkCache(std::unique_ptr<Oracle<Index_> > oracle, size_t per_iteration, size_t num_chunks) : 
        prediction_stream(std::move(o)), max_chunks(num_chunks), max_predictions(per_iteration) {}

private:
    /**
     * Fetch the next chunk according to the stream of predictions provided by the `Oracle`.
     *
     * @param identify Function that accepts `(Index_ i, Id_& id, Index_& internal)`, where `i` is the predicted row/column index.
     * On output, `id` is set to the identifier of the chunk containing `i`,
     * and `internal` is set to the index of row/column `i` inside that chunk.
     * @param populate Function that accepts `(const std::vector<std::pair<Id_, Index_> >& chunks_in_need, std::vector<ChunkContents>& chunk_data)`.
     * This should iterate over the `chunks_in_need` and populates the corresponding entries in `chunk_data`.
     * Each entry of `chunks_in_need` contains the chunk identifier and the index of the `ChunkContents` in `chunk_data`;
     * it is expected that each chunk's contents is written to the specified entry of `chunk_data`.
     * @param mock Function that accepts no arguments and returns a mock `ChunkContents` object containing no data and with minimal allocated memory.
     * This will be used as a placeholder for some internal book-keeping.
     * @param allocate Function that accepts a single `ChunkContents` object created by `mock()`,
     * and allocates sufficient memory to it in order to hold a chunk's contents when used in `populate()`.
     *
     * @return Pair containing (1) a pointer to a chunk's contents and (2) the index of the next predicted row/column inside the retrieved chunk.
     */
    template<class Ifunction_, class Pfunction_, class Mfunction_, class Afunction_>
    std::pair<const ChunkContents_*, Index_> next_chunk(Ifunction_ identify, Pfunction_ populate, Mfunction_ mock, Afunction_ allocate) {
        if (predictions_made.size() > predictions_fulfilled) {
            const auto& chosen = predictions_made[predictions_fulfilled++];
            return std::make_pair(cache_data.data() + chosen.first, chosen.second);
        }

        bool starting = cache_data.empty();
        if (starting) {
            cache_data.resize(max_chunks);
            next_cache_data.resize(max_chunks, mock()); // using a placeholder value for availability checks.
            chunks_in_need.reserve(max_chunks);
        } else {
            next_cache_exists.clear();
            chunks_in_need.clear();
        }

        size_t used = 0;

        predictions_made.clear();
        predictions_made.reserve(max_predictions);
        for (size_t p = 0; p < max_predictions; ++p) {
            Index_ current;
            if (!prediction_stream.next(current)) {
                break;
            }

            Id_ curchunk;
            Index_ curindex;
            identify(current, curchunk, curindex);

            auto it = next_cache_exists.find(curchunk);
            if (it == next_cache_exists.end()) {
                next_cache_exists[curchunk] = used;
                predictions_made.emplace_back(used, curindex);

                auto it2 = cache_exists.find(curchunk);
                if (it2 != cache_exists.end()) {
                    next_cache_data[used].swap(cache_data[it2->second]);
                } else {
                    chunks_in_need.emplace_back(curchunk, used);
                }

                ++used;
                if (used == max_chunks) {
                    break;
                }
            } else {
                predictions_made.emplace_back(it->second, curindex);
            }
        }

        /**
         * Doing a linear scan across chunks to find the currently unused cache
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
            // If we're lucky enough to already be sitting on an allocation
            // (e.g., if the last call to this function didn't consume all
            // cache elements), we just use it directly.
            if (next_cache_data[c.second].empty()) { 
                while (cache_data[search].empty()) { 
                    ++search;
                }

#ifdef DEBUG
                // This should never reach the end, as there should always be
                // 'num_chunks' non-empty elements across cache_data and
                // next_cache_data. Nonetheless, we add an assertion here.
                if (search >= cache_data.size()) {
                    throw std::runtime_error("!!! internal cache management error for oracles !!!");
                }
#endif

                next_cache_data[c.second].swap(cache_data[search]);
                ++search;
            }

            allocate(next_cache_data[c.second]); // eventually no-op when all available caches are of the right size.
        }

        populate(chunks_in_need, cache_data);

        cache_data.swap(next_cache_data);
        cache_exists.swap(next_cache_exists);
        predictions_fulfilled = 1; // well, because we just used one.
        const auto& chosen = predictions_made.front(); // assuming at least one prediction was made, otherwise, why was this function even called?
        return std::make_pair(cache_data.data() + chosen.first, chosen.second);
    }
};

}

#endif
