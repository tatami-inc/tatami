#ifndef TATAMI_LRU_CHUNK_CACHE_HPP
#define TATAMI_LRU_CHUNK_CACHE_HPP

#include <unordered_map>
#include <list>

/**
 * @file LruChunkCache.hpp
 * @brief Create a LRU cache for matrix chunks.
 */

namespace tatami {

/**
 * @tparam Id_ Type of chunk identifier, typically integer.
 * @tparam ChunkContents_ Class containing the contents of a single chunk.
 *
 * @brief LRU cache for chunks.
 *
 * Implement a least-recently-used cache for chunks.
 * This is typically used for `Matrix` representations where the data is costly to load (e.g., from file) and no oracle is provided to predict future accesses.
 * In such cases, chunks of data can be loaded and cached such that any possible future access to an already-loaded chunk will just fetch it from cache.
 */
template<typename Id_, class ChunkContents_> 
class LruChunkCache {
private:
    typedef std::pair<ChunkContents_, Id_> Element;
    std::list<Element> cache_data;
    std::unordered_map<Id_, typename std::list<Element>::iterator> cache_exists;
    size_t max_chunks; 

public:
    /**
     * @param m Maximum number of chunks to store.
     */
    LruChunkCache(size_t m) : max_chunks(m) {}

public:
    /**
     * @tparam Cfunction_ Function to create a new `ChunkContents_` object.
     * @tparam Pfunction_ Function to populate a `ChunkContents_` object with the contents of a chunk.
     *
     * @param id Identifier for the chunk.
     * @param create Function that accepts no arguments and returns a `ChunkContents_` object.
     * @param populate Function that accepts a chunk ID and a reference to a `ChunkContents_` object,
     * and populates the latter with the contents of the former.
     * 
     * @return Reference to a chunk.
     * If the chunk already exists in the cache, it is returned directly.
     * If the chunk does not exist and there is still space in the cache, a new chunk is created and populated with the contents of chunk `id`.
     * If the chunk does not exist and there is no space in the cache, the least recently used chunk is evicted and its `ChunkContents_` is populated with the contents of chunk `id`.
     */
    template<class Cfunction_, class Pfunction_>
    const ChunkContents_& find_chunk(Id_ id, Cfunction_ create, Pfunction_ populate) {
        if (max_chunks == 1) {
            // Minor optimization if there's just one chunk, in which case we can
            // skip the search and the list splice if there's a hit.
            if (!cache_data.empty()) {
                const auto& solo = cache_data.front();
                if (solo.second == id) {
                    return solo.first;
                }
            }
        } else {
            auto it = cache_exists.find(id);
            if (it != cache_exists.end()) {
                auto chosen = it->second;
                cache_data.splice(cache_data.end(), cache_data, chosen); // move to end.
                return chosen->first;
            } 
        }

        typename std::list<Element>::iterator location;
        if (cache_data.size() < max_chunks) {
            cache_data.emplace_back(create(), id);
            location = std::prev(cache_data.end());
        } else {
            location = cache_data.begin();
            cache_exists.erase(location->second);
            location->second = id;
            cache_data.splice(cache_data.end(), cache_data, location); // move to end.
        }
        cache_exists[id] = location;

        auto& chunk = location->first;
        populate(id, chunk);
        return chunk;
    }
};

}

#endif
