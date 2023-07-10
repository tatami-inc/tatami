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

    template<bool accrow_>
    std::pair<size_t, size_t> get_chunk_offset_and_stride(size_t primary_chunk_index) const {
        size_t offset, increment;
        if (row_major) {
            if constexpr(accrow_) {
                offset = primary_chunk_index * num_chunks_per_row;
                increment = 1;
            } else {
                offset = primary_chunk_index;
                increment = num_chunks_per_row;
            }
        } else {
            if constexpr(accrow_) {
                offset = primary_chunk_index;
                increment = num_chunks_per_column;
            } else {
                offset = primary_chunk_index * num_chunks_per_column;
                increment = 1;
            }
        }
        return std::make_pair(offset, increment);
    }

    static size_t get_primary_upper_bound(size_t primary_chunk_index, size_t primary_maxdim, size_t primary_chunkdim) {
        return std::min(primary_chunkdim, primary_maxdim - primary_chunk_index * primary_chunkdim);
    }

public:
    static constexpr bool sparse_chunk = Chunk_::sparse;
    typedef typename Chunk_::index_type Chunkdex;

    struct SparseCache {
        SparseCache(size_t primary_dim) : cache_indices(primary_dim), cache_values(primary_dim) {}

        std::vector<std::vector<Chunkdex> > cache_indices;
        std::vector<std::vector<typename Chunk_::value_type> > cache_values;

        std::vector<size_t> buffer_indptrs;
        std::vector<Chunkdex> buffer_indices;
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
            return DenseCache(length, exact_ ? 1 : get_primary_chunkdim<accrow_>());
        }
    }

public:
    template<bool accrow_, bool exact_>
    void extract(size_t primary_chunk_index, size_t primary_chunk_offset, size_t primary_maxdim, size_t start, size_t length, Cache& cache) const {
        if (!length) {
            return;
        }

        size_t primary_num_chunks = get_primary_num_chunks<accrow_>();
        size_t primary_chunkdim = get_primary_chunkdim<accrow_>();
        size_t secondary_chunkdim = get_secondary_chunkdim<accrow_>();

        size_t primary_end_pos = (exact_ ? 0 : get_primary_upper_bound(primary_chunk_index, primary_maxdim, primary_chunkdim));

        auto chunk_stride = get_chunk_offset_and_stride<accrow_>(primary_chunk_index);
        auto offset = chunk_stride.first;
        auto increment = chunk_stride.second;

        size_t start_chunk_index = start / secondary_chunkdim;
        size_t end_chunk_index = (start + length + secondary_chunkdim - 1) / secondary_chunkdim; // i.e., ceiling of the integer division.
        offset += increment * start_chunk_index;

        size_t dense_cache_offset = 0;
        size_t secondary_start_pos = start_chunk_index * secondary_chunkdim;

        for (size_t c = start_chunk_index; c < end_chunk_index; ++c) {
            const auto& chunk = chunks[offset];
            size_t from = (c == start_chunk_index ? start - secondary_start_pos : 0);
            size_t to = (c + 1 == end_chunk_index ? start + length - secondary_start_pos : secondary_chunkdim);

            if constexpr(sparse_chunk) {
                chunk.inflate(cache.buffer_values, cache.buffer_indices, cache.buffer_indptrs);
                Chunkdex secondary_offset = secondary_start_pos;

                if (chunk.row_major == accrow_) {
                    auto refine_start_and_end = [&](size_t start, size_t end) -> void {
                        if (from) {
                            auto it = cache.buffer_indices.begin();
                            start = std::lower_bound(it + start, it + end, static_cast<Chunkdex>(from)) - it;
                        }
                        if (to != secondary_chunkdim) {
                            auto it = cache.buffer_indices.begin();
                            end = std::lower_bound(it + start, it + end, static_cast<Chunkdex>(to)) - it;
                        }
                    };

                    if constexpr(exact_) {
                        auto start = cache.buffer_indptrs[primary_chunk_offset], end = cache.buffer_indptrs[primary_chunk_offset + 1];
                        if (start < end) {
                            refine_start_and_end(start, end);
                            cache.cache_values[0].insert(cache.cache_values[0].end(), cache.buffer_values.begin() + start, cache.buffer_values.begin() + end);
                            for (size_t i = start; i < end; ++i) {
                                cache.cache_indices[0].push_back(cache.buffer_indices[i] + secondary_offset);
                            }
                        }

                    } else {
                        for (size_t p = 0; p < primary_end_pos; ++p) {
                            auto start = cache.buffer_indptrs[p], end = cache.buffer_indptrs[p + 1];
                            if (start < end) {
                                refine_start_and_end(start, end);
                                cache.cache_values[p].insert(cache.cache_values[p].end(), cache.buffer_values.begin() + start, cache.buffer_values.begin() + end);
                                for (size_t i = start; i < end; ++i) {
                                    cache.cache_indices[p].push_back(cache.buffer_indices[i] + secondary_offset);
                                }
                            }
                        }
                    }

                } else {
                    if constexpr(exact_) {
                        Chunkdex target = primary_chunk_offset;
                        for (size_t s = from; s < to; ++s) {
                            auto start = cache.buffer_indptrs[s], end = cache.buffer_indptrs[s + 1];
                            auto it = cache.buffer_indices.begin();
                            start = std::lower_bound(it + start, it + end, target) - it;
                            if (start != end && cache.buffer_indices[start] == target) {
                                cache.cache_values[0].push_back(cache.buffer_values[start]);
                                cache.cache_indices[0].push_back(s + secondary_offset);
                            }
                        }

                    } else {
                        for (size_t s = from; s < to; ++s) {
                            auto start = cache.buffer_indptrs[s], end = cache.buffer_indptrs[s + 1];
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
                    if constexpr(exact_) {
                        bptr += primary_chunk_offset * secondary_chunkdim;
                        std::copy(bptr + from, bptr + to, cptr);

                    } else {
                        for (size_t p = 0; p < primary_end_pos; ++p) {
                            std::copy(bptr + from, bptr + to, cptr);
                            bptr += secondary_chunkdim;
                            cptr += length;
                        }
                    }

                } else {
                    if constexpr(exact_) {
                        bptr += from * primary_chunkdim + primary_chunk_offset;
                        for (size_t s = from; s < to; ++s) {
                            *cptr = *bptr;
                            ++cptr;
                            bptr += primary_chunkdim;
                        }

                    } else {
                        bptr += from * primary_chunkdim;
                        for (size_t p = 0; p < primary_end_pos; ++p) {
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
                }

                dense_cache_offset += (to - from);
            }

            offset += increment;
            secondary_start_pos += secondary_chunkdim; 
        }
    }

public:
    template<bool accrow_, bool exact_, typename Index_>
    void extract(size_t primary_chunk_index, size_t primary_chunk_offset, size_t primary_maxdim, const std::vector<Index_>& indices, Cache& cache) const {
        if (indices.empty()) {
            return;
        }

        size_t primary_num_chunks = get_primary_num_chunks<accrow_>();
        size_t primary_chunkdim = get_primary_chunkdim<accrow_>();
        size_t secondary_chunkdim = get_secondary_chunkdim<accrow_>();

        size_t primary_end_pos = (exact_ ? 0 : get_primary_upper_bound(primary_chunk_index, primary_maxdim, primary_chunkdim));

        auto chunk_stride = get_chunk_offset_and_stride<accrow_>(primary_chunk_index);
        auto offset = chunk_stride.first;
        auto increment = chunk_stride.second;

        size_t start_chunk_index = indices.front() / secondary_chunkdim;
        std::vector<Index_> collected;
        auto iIt = indices.begin();
        offset += start_chunk_index * increment;

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
                    Chunkdex secondary_offset = secondary_start_pos;

                    auto collect_sparse = [&](size_t start, size_t end, size_t pout) -> void {
                        if (collected.front()) {
                            auto it = cache.buffer_indices.begin();
                            start = std::lower_bound(it + start, it + end, static_cast<Chunkdex>(collected.front())) - it;
                        }

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
                    };

                    if (chunk.row_major == accrow_) {
                        if constexpr(exact_) {
                            auto start = cache.buffer_indptrs[primary_chunk_offset], end = cache.buffer_indptrs[primary_chunk_offset + 1];
                            if (start < end) {
                                collect_sparse(start, end, 0);
                            }

                        } else {
                            for (size_t p = 0; p < primary_end_pos; ++p) {
                                auto start = cache.buffer_indptrs[p], end = cache.buffer_indptrs[p + 1];
                                if (start < end) {
                                    collect_sparse(start, end, p);
                                }
                            }
                        }

                    } else {
                        if constexpr(exact_) {
                            Chunkdex target = primary_chunk_offset;
                            for (auto s : collected) {
                                auto start = cache.buffer_indptrs[s], end = cache.buffer_indptrs[s + 1];
                                auto it = cache.buffer_indices.begin();
                                start = std::lower_bound(it + start, it + end, target) - it;
                                if (start != end && cache.buffer_indices[start] == target) {
                                    cache.cache_values[0].push_back(cache.buffer_values[start]);
                                    cache.cache_indices[0].push_back(s + secondary_offset);
                                }
                            }

                        } else {
                            for (auto s : collected) {
                                auto start = cache.buffer_indptrs[s], end = cache.buffer_indptrs[s + 1];
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
                        if constexpr(exact_) {
                            bptr += primary_chunk_offset * secondary_chunkdim;
                            for (auto x : collected) {
                                *cptr = bptr[x];
                                ++cptr;
                            }

                        } else {
                            for (size_t p = 0; p < primary_end_pos; ++p) {
                                auto copy_cptr = cptr;
                                for (auto x : collected) {
                                    *copy_cptr = bptr[x];
                                    ++copy_cptr;
                                }
                                bptr += secondary_chunkdim;
                                cptr += indices.size();
                            }
                        }

                    } else {
                        if constexpr(exact_) {
                            bptr += primary_chunk_offset;
                            for (auto x : collected) {
                                *cptr = bptr[x * primary_chunkdim];
                                ++cptr;
                            }

                        } else {
                            for (size_t p = 0; p < primary_end_pos; ++p) {
                                auto copy_cptr = cptr;
                                for (auto x : collected) {
                                    *copy_cptr = bptr[x * primary_chunkdim];
                                    ++copy_cptr;
                                }
                                ++bptr;
                                cptr += indices.size();
                            }
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
