#ifndef TATAMI_CUSTOM_CHUNK_MANAGER_HPP
#define TATAMI_CUSTOM_CHUNK_MANAGER_HPP

#include <vector>
#include <algorithm>

namespace tatami {

template<typename Chunk_>
class CustomChunkManager {
public:
    std::vector<Chunk_> chunks;
    size_t chunk_nrow;
    size_t chunk_ncol;
    size_t num_chunks_per_row;
    size_t num_chunks_per_column;
    bool row_major;

private:
    template<bool accrow_>
    size_t get_primary_chunkdim() const {
        if constexpr(accrow_) {
            return chunk_nrow;
        } else {
            return chunk_ncol;
        }
    }

    template<bool accrow_>
    size_t get_secondary_chunkdim() const {
        if constexpr(accrow_) {
            return chunk_ncol;
        } else {
            return chunk_nrow;
        }
    }

    template<bool accrow_>
    size_t get_primary_num_chunks() const {
        if constexpr(accrow_) {
            return num_chunks_per_row;
        } else {
            return num_chunks_per_column;
        }
    }

public:
    static constexpr bool sparse_chunk = Chunk_::sparse;

    struct SparseCache {
        SparseCache(size_t primary_dim) : cache_indices(primary_dim), cache_values(primary_dim) {}

        std::vector<std::vector<typename Chunk_::index_type> > cache_indices;
        std::vector<std::vector<typename Chunk_::value_type> > cache_values;

        std::vector<size_t> buffer_indptrs;
        std::vector<typename Chunk_::index_type> buffer_indices;
        std::vector<typename Chunk_::value_type> buffer_values;
    };

    struct DenseCache {
        DenseCache(size_t length, size_t primary_dim) : cache(length * primary_dim) {}

        std::vector<typename Chunk_::value_type> cache;

        std::vector<typename Chunk_::value_type> buffer;
    };

    typedef typename std::conditional<sparse_chunk, SparseCache, DenseCache>::type Cache;

    template<bool accrow_, bool exact_>
    Cache create_chunk_cache(size_t length) const {
        if constexpr(sparse_chunk) {
            return SparseCache(exact_ ? 1 : get_primary_chunkdim<accrow_>());
        } else {
            return DenseCache(exact_ ? 1 : length, get_primary_chunkdim<accrow_>());
        }
    }

public:
    template<bool accrow_, bool exact_>
    void extract_block(size_t i, size_t start, size_t length, Cache& cache) {
        if (!length) {
            return;
        }

        size_t primary_num_chunks = get_primary_num_chunks<accrow_>();
        size_t primary_chunkdim = get_primary_chunkdim<accrow_>();
        size_t secondary_chunkdim = get_secondary_chunkdim<accrow_>();

        size_t primary_chunk_index = (i / primary_chunkdim);
        size_t shift = (row_major ? num_chunks_per_row : num_chunks_per_column);
        size_t increment = (row_major == accrow_ ? 1 : shift);
        size_t offset = primary_chunk_index * shift;

        size_t start_chunk_index = start / secondary_chunkdim;
        size_t end_chunk_index = (start + length + secondary_chunkdim - 1) / secondary_chunkdim; // i.e., ceiling of the integer division.

        size_t primary_start_pos = 0, primary_end_pos = primary_chunkdim;
        if constexpr(exact_) {
            primary_start_pos = i - primary_chunk_index * primary_chunkdim;
            primary_end_pos = primary_start_pos + 1;
        }

        size_t dense_cache_offset = 0;
        size_t secondary_start_pos = start_chunk_index * secondary_chunkdim;

        for (size_t c = start_chunk_index; c < end_chunk_index; ++c) {
            const auto& chunk = chunks[offset];
            size_t from = (c == start_chunk_index ? start - secondary_start_pos : 0);
            size_t to = (c + 1 == end_chunk_index ? start + length - secondary_start_pos : secondary_chunkdim);

            if constexpr(sparse_chunk) {
                chunk.inflate(cache.buffer_values, cache.buffer_indices, cache.buffer_indptrs);
                typename Chunk_::index_type secondary_offset = secondary_start_pos;

                if (chunk.row_major == accrow_) {
                    for (size_t p = primary_start_pos; p < primary_end_pos; ++p) {
                        auto start = cache.buffer_indptrs[p], end = cache.buffer_indptrs[p + 1];
                        if (start == end) {
                            continue;
                        }

                        if (from) {
                            start = std::lower_bound(cache.buffer_indices.begin() + start, cache.buffer_indices.begin() + end, from) - cache.buffer_indices.begin();
                        }
                        if (to != secondary_chunkdim) {
                            end = std::lower_bound(cache.buffer_indices.begin() + start, cache.buffer_indices.begin() + end, to) - cache.buffer_indices.begin();
                        }

                        auto pout = (exact_ ? 0 : p);
                        cache.cache_values[pout].insert(cache.cache_values[p].end(), cache.buffer_values.begin() + start, cache.buffer_values.begin() + end);
                        for (size_t i = start; i < end; ++i) {
                            cache.cache_indices[pout].push_back(cache.buffer_indices[i] + secondary_offset);
                        }
                    }

                } else {
                    for (size_t s = from; s < to; ++s) {
                        auto start = cache.buffer_indptrs[s], end = cache.buffer_indptrs[s + 1];
                        if constexpr(exact_) {
                            start = std::lower_bound(cache.buffer_indices.begin() + start, cache.buffer_indices.begin() + end, primary_start_pos) - cache.buffer_indices.begin();
                            if (start != end && cache.buffer_indices[start] == primary_start_pos) {
                                cache.cache_values[0].push_back(cache.buffer_values[start]);
                                cache.cache_indices[0].push_back(s + secondary_offset);
                            }

                        } else {
                            for (size_t i = start; i < end; ++i) {
                                auto p = cache.buffer_indices[i];
                                cache.cache_values[p].push_back(cache.buffer_values[i]);
                                cache.cache_indices[p].push_back(s + secondary_offset);
                            }
                        }
                    }
                }

            } else {
                chunk.inflate(cache.buffer);
                auto bptr = cache.buffer.data();
                auto cptr = cache.cache.data() + dense_cache_offset;

                if (chunk.row_major == accrow_) {
                    for (size_t p = primary_start_pos; p < primary_end_pos; ++p) {
                        std::copy(bptr + from, bptr + to, cptr);
                        bptr += secondary_chunkdim;
                        cptr += length;
                    }

                } else {
                    for (size_t p = primary_start_pos; p < primary_end_pos; ++p) {
                        auto copy_bptr = bptr;
                        auto copy_cptr = cptr;
                        for (size_t s = from; s < to; ++s) {
                            *copy_cptr = *copy_bptr;
                            ++copy_cptr;
                            copy_bptr += primary_chunkdim;
                        }
                        ++bptr;
                        cptr += length;
                    }
                }

                dense_cache_offset += (to - from);
            }

            offset += increment;
            secondary_start_pos += secondary_chunkdim; 
        }
    }

public:
    template<bool accrow_, bool exact_, typename Index_>
    void extract_dense_block(size_t i, const std::vector<Index_>& indices, Cache& cache) {
        if (indices.empty()) {
            return;
        }

        size_t primary_num_chunks = get_primary_num_chunks<accrow_>();
        size_t primary_chunkdim = get_primary_chunkdim<accrow_>();
        size_t secondary_chunkdim = get_secondary_chunkdim<accrow_>();

        size_t primary_chunk_index = (i / primary_chunkdim);
        size_t shift = (row_major ? num_chunks_per_row : num_chunks_per_column);
        size_t increment = (row_major == accrow_ ? 1 : shift);
        size_t offset = primary_chunk_index * shift;

        size_t start_chunk_index = indices.front() / secondary_chunkdim;
        std::vector<Index_> collected;
        auto iIt = indices.begin();

        size_t primary_start_pos = 0, primary_end_pos = primary_chunkdim;
        if constexpr(exact_) {
            primary_start_pos = i - primary_chunk_index * primary_chunkdim;
            primary_end_pos = primary_start_pos + 1;
        }

        size_t dense_cache_offset = 0;
        size_t secondary_start_pos = start_chunk_index * secondary_chunkdim;

        for (size_t c = start_chunk_index; c < primary_num_chunks && iIt != indices.end(); ++c) {
            Index_ secondary_end_pos = secondary_start_pos + secondary_chunkdim;
            collected.clear();
            while (iIt != indices.end() && *iIt < secondary_end_pos) {
                collected.push_back(*iIt - secondary_start_pos);
                ++iIt;
            }

            if (!collected.empty()) {
                const auto& chunk = chunks[offset];

                if constexpr(sparse_chunk) {
                    chunk.inflate(cache.buffer_values, cache.buffer_indices, cache.buffer_indptrs);
                    typename Chunk_::index_type secondary_offset = secondary_start_pos;

                    if (chunk.row_major == accrow_) {
                        for (size_t p = primary_start_pos; p < primary_end_pos; ++p) {
                            auto start = cache.buffer_indptrs[p], end = cache.buffer_indptrs[p + 1];
                            if (start == end) {
                                continue;
                            }

                            if (collected.front()) {
                                start = std::lower_bound(cache.buffer_indices.begin() + start, cache.buffer_indices.begin() + end, collected.front()) - cache.buffer_indices.begin();
                            }

                            auto pout = (exact_ ? 0 : p);
                            auto cIt = collected.begin();
                            for (size_t i = start; i < end; ++i) {
                                Index_ target = cache.buffer_indices[i];
                                while (cIt != collected.end() && *cIt < target) {
                                    ++cIt;
                                }
                                if (cIt == collected.end()) {
                                    break;
                                }
                                if (*cIt == target) {
                                    cache.cache_values[pout].push_back(cache.buffer_values[i]);
                                    cache.cache_indices[pout].push_back(cache.buffer_indices[i] + secondary_offset);
                                    ++cIt;
                                }
                            }
                        }

                    } else {
                        for (auto s : collected) {
                            auto start = cache.buffer_indptrs[s], end = cache.buffer_indptrs[s + 1];
                            if constexpr(exact_) {
                                start = std::lower_bound(cache.buffer_indices.begin() + start, cache.buffer_indices.begin() + end, primary_start_pos) - cache.buffer_indices.begin();
                                if (start != end && cache.buffer_indices[start] == primary_start_pos) {
                                    cache.cache_values[0].push_back(cache.buffer_values[start]);
                                    cache.cache_indices[0].push_back(s + secondary_offset);
                                }

                            } else {
                                for (size_t i = start; i < end; ++i) {
                                    auto p = cache.buffer_indices[i];
                                    cache.cache_values[p].push_back(cache.buffer_values[i]);
                                    cache.cache_indices[p].push_back(s + secondary_offset);
                                }
                            }
                        }
                    }

                } else {
                    chunk.inflate(cache.buffer);
                    auto bptr = cache.buffer.data();
                    auto cptr = cache.cache.data() + dense_cache_offset;

                    if (chunk.row_major == accrow_) {
                        for (size_t p = primary_start_pos; p < primary_end_pos; ++p) {
                            auto copy_cptr = cptr;
                            for (auto x : collected) {
                                *copy_cptr = bptr[x];
                                ++copy_cptr;
                            }
                            bptr += secondary_chunkdim;
                            cptr += indices.size();
                        }

                    } else {
                        for (size_t p = primary_start_pos; p < primary_end_pos; ++p) {
                            auto copy_cptr = cptr;
                            for (auto x : collected) {
                                *copy_cptr = bptr[x * primary_chunkdim];
                                ++copy_cptr;
                            }
                            ++bptr;
                            cptr += indices.size();
                        }
                    }

                    dense_cache_offset += collected.size();
                }
            }

            offset += increment;
            secondary_start_pos += secondary_chunkdim;
        }
    }
};

}

#endif
