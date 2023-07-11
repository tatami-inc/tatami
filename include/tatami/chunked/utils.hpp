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
 * @brief Wrap a dense chunk for a `CustomChunkedMatrix`.
 *
 * Implements a simple wrapper around a dense chunk for use inside a `CustomChunkedMatrix`.
 * Each chunk should hold a 2-dimensional array of numeric values.
 * This wrapper is considered to be "simple", as extraction of any data involves inflation of the entire chunk.
 *
 * The `Chunk_` class should provide the following:
 *
 * - A `static constexpr bool row_major`, specifying whether the array is row-major.
 * - A `typedef value_type`, specifying the type of the value in the array.
 * - A `nrow() const` method, defining the number of rows in the array.
 * - A `ncol() const` method, defining the number of columns in the array.
 * - A `void inflate(std::vector<value_type>& buffer) const` method that fills `buffer` with the contents of the array.
 *   This should be filled in row-major format if `row_major = true` and in column-major format otherwise.
 *
 * @tparam Chunk_ Class to represent the chunk.
 */
template<class Chunk_>
struct SimpleDenseChunkWrapper {
    /**
     * Type of the value stored in this chunk.
     */
    typedef typename Chunk_::value_type value_type;

    /**
     * Workspace for chunk extraction.
     * This can be used in multiple `fetch()` calls, possibly across different chunks.
     */
    typedef std::vector<value_type> Workspace;

public:
    /**
     * Default constructor.
     */
    SimpleDenseChunkWrapper() = default;

    /**
     * @param c Chunk to be wrapped.
     */
    SimpleDenseChunkWrapper(Chunk_ c) : chunk(std::move(c)) {}

private:
    Chunk_ chunk;

    template<bool accrow_>
    auto get_primary_chunkdim() const {
        if constexpr(accrow_) {
            return chunk.nrow();
        } else {
            return chunk.ncol();
        }
    }

    template<bool accrow_>
    auto get_secondary_chunkdim() const {
        if constexpr(accrow_) {
            return chunk.ncol();
        } else {
            return chunk.nrow();
        }
    }

public:
    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the chunk.
     * @tparam Output_ Numeric type for the output.
     *
     * @param primary_start Index of the first element on the primary dimension to be extracted.
     * If `accrow_ = true`, this is the first row, otherwise it is the first column.
     * @param primary_length Number of elements on the primary dimension to be extracted.
     * If `accrow_ = true`, this is the number of rows, otherwise it is the number of columns.
     * This is guaranteed to be positive.
     * @param secondary_start Index of the first element on the secondary dimension to be extracted.
     * If `accrow_ = true`, this is the first column, otherwise it is the first row.
     * @param secondary_length Number of elements on the secondary dimension to be extracted.
     * If `accrow_ = true`, this is the number of columns, otherwise it is the number of rows.
     * This is guaranteed to be positive.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output Pointer to an output array of length no less than `primary_length * stride`.
     * @param stride Stride separating corresponding values from consecutive elements on the primary dimension.
     *
     * This method extracts the specified values from the chunk into `output`.
     * For example, if `accrow_ = true`, we would extract rows `[primary_start, primary_start + length)` and columns `[secondary_start, secondary_start + secondary_length)`.
     * For a primary dimension index `p` and secondary dimension index `s`, the value from the chunk should be stored in `output[(p - primary_start) * stride + (s - secondary_start)]`.
     * This layout allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomChunkedMatrix`.
     */
    template<bool accrow_, typename Index_, typename Output_>
    void extract(Index_ primary_start, Index_ primary_length, Index_ secondary_start, Index_ secondary_length, Workspace& work, Output_* output, size_t stride) const {
        chunk.inflate(work);
        size_t primary_chunkdim = get_primary_chunkdim<accrow_>(); // use size_t to avoid integer overflow with Index_.
        size_t secondary_chunkdim = get_secondary_chunkdim<accrow_>();

        if constexpr(Chunk_::row_major == accrow_) {
            auto srcptr = work.data() + primary_start * secondary_chunkdim + secondary_start;
            for (size_t p = 0; p < primary_length; ++p) {
                std::copy(srcptr, srcptr + secondary_length, output);
                srcptr += secondary_chunkdim;
                output += stride;
            }

        } else {
            auto srcptr = work.data() + secondary_start * primary_chunkdim + primary_start;
            for (size_t p = 0; p < primary_length; ++p) {
                auto copy_srcptr = srcptr;
                auto copy_output = output;
                for (size_t s = 0; s < secondary_length; ++s) {
                    *copy_output = *copy_srcptr;
                    ++copy_output;
                    copy_srcptr += primary_chunkdim;
                }
                ++srcptr;
                output += stride;
            }
        }
    }

    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the chunk.
     * @tparam Output_ Numeric type for the output.
     *
     * @param primary_start Index of the first element on the primary dimension to be extracted.
     * If `accrow_ = true`, this is the first row, otherwise it is the first column.
     * @param primary_length Number of elements on the primary dimension to be extracted.
     * If `accrow_ = true`, this is the number of rows, otherwise it is the number of columns.
     * This is guaranteed to be positive.
     * @param secondary_indices Indices of the elements on the secondary dimension to be extracted.
     * If `accrow_ = true`, these are column indices, otherwise these are row indices.
     * This is guaranteed to be non-empty.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output Pointer to an output array of length no less than `primary_length * stride`.
     * @param stride Stride separating corresponding values from consecutive elements on the primary dimension.
     *
     * This method extracts the specified values from the chunk into `output`.
     * For example, if `accrow_ = true`, we would extract rows `[primary_start, primary_start + length)` and columns in `secondary_indices`.
     * For a primary dimension index `p` and secondary dimension index `secondary_indices[i]`, the value from the chunk should be stored in `output[(p - primary_start) * stride + i]`.
     * This layout allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomChunkedMatrix`.
     */
    template<bool accrow_, typename Index_, typename Output_>
    void extract(Index_ primary_start, Index_ primary_length, const std::vector<Index_>& secondary_indices, Workspace& work, Output_* output, size_t stride) const {
        chunk.inflate(work);
        size_t primary_chunkdim = get_primary_chunkdim<accrow_>(); // use size_t to avoid integer overflow with Index_.
        size_t secondary_chunkdim = get_secondary_chunkdim<accrow_>(); 

        if constexpr(Chunk_::row_major == accrow_) {
            auto srcptr = work.data() + primary_start * secondary_chunkdim;
            for (size_t p = 0; p < primary_length; ++p) {
                auto copy_output = output;
                for (auto x : secondary_indices) {
                    *copy_output = srcptr[x];
                    ++copy_output;
                }
                srcptr += secondary_chunkdim;
                output += stride;
            }

        } else {
            auto srcptr = work.data() + primary_start;
            for (size_t p = 0; p < primary_length; ++p) {
                auto copy_output = output;
                for (auto x : secondary_indices) {
                    *copy_output = srcptr[x * primary_chunkdim];
                    ++copy_output;
                }
                ++srcptr;
                output += stride;
            }
        }
    }
};

/**
 * @brief Wrap a sparse chunk for a `CustomChunkedMatrix`.
 *
 * Implements a simple wrapper around a sparse chunk for use inside a `CustomChunkedMatrix`.
 * Each chunk should hold a 2-dimensional compressed sparse submatrix of numeric values.
 * This wrapper is considered to be "simple", as extraction of any data involves inflation of the entire chunk.
 *
 * The `Chunk_` class should provide the following:
 *
 * - A `static constexpr bool row_major`, specifying whether the submatrix is in the compressed sparse row layout.
 * - A `typedef value_type`, specifying the type of the value in the submatrix.
 * - A `typedef index_type`, specifying the type of the index in the submatrix.
 * - A `nrow() const` method, defining the number of rows in the submatrix.
 * - A `ncol() const` method, defining the number of columns in the submatrix.
 * - A `void inflate(std::vector<value_type>& values, std::vector<index_type>& indices, std::vector<size_t>& pointers) const` method that fills `values`, `indices` and `pointers` with the standard compressed sparse data.
 *   This should be a CSR matrix if `row_major = true` and a CSC matrix otherwise.
 *
 * @tparam Chunk_ Class to represent the chunk.
 */
template<class Chunk_>
struct SimpleSparseChunkWrapper {
    /**
     * Type of the value stored in this chunk.
     */
    typedef typename Chunk_::value_type value_type;

    /**
     * Type of the index stored in this chunk.
     */
    typedef typename Chunk_::index_type index_type;

    /**
     * @brief Workspace for chunk extraction.
     *
     * This can be used in multiple `fetch()` calls, possibly across different chunks.
     */
    struct Workspace {
        /**
         * @cond
         */
        std::vector<value_type> values;
        std::vector<index_type> indices;
        std::vector<size_t> indptrs;
        /**
         * @endcond
         */
    };

public:
    /**
     * Default constructor.
     */
    SimpleSparseChunkWrapper() = default;

    /**
     * @param c Chunk to be wrapped.
     */
    SimpleSparseChunkWrapper(Chunk_ c) : chunk(std::move(c)) {}

private:
    Chunk_ chunk;

    template<bool accrow_>
    size_t get_primary_chunkdim() const {
        if constexpr(accrow_) {
            return chunk.nrow();
        } else {
            return chunk.ncol();
        }
    }

    template<bool accrow_>
    size_t get_secondary_chunkdim() const {
        if constexpr(accrow_) {
            return chunk.ncol();
        } else {
            return chunk.nrow();
        }
    }

    template<typename Index_>
    static void refine_start_and_end(size_t& start, size_t& end, Index_ desired_start, Index_ desired_end, Index_ max_end, const std::vector<typename Chunk_::index_type>& indices) {
        if (desired_start) {
            auto it = indices.begin();
            start = std::lower_bound(it + start, it + end, static_cast<typename Chunk_::index_type>(desired_start)) - it;
        }

        if (desired_end != max_end) {
            if (desired_end == desired_start + 1) {
                if (start != end && indices[start] == desired_start) {
                    end = start + 1;
                } else {
                    end = start;
                }
            } else {
                auto it = indices.begin();
                end = std::lower_bound(it + start, it + end, static_cast<typename Chunk_::index_type>(desired_end)) - it;
            }
        }
    }

public:
    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the chunk.
     * @tparam OutputValue_ Numeric type for the output values.
     * @tparam OutputIndex_ Integer type for the output indices.
     *
     * @param primary_start Index of the first element on the primary dimension to be extracted.
     * If `accrow_ = true`, this is the first row, otherwise it is the first column.
     * @param primary_length Number of elements on the primary dimension to be extracted.
     * If `accrow_ = true`, this is the number of rows, otherwise it is the number of columns.
     * This is guaranteed to be positive.
     * @param secondary_start Index of the first element on the secondary dimension to be extracted.
     * If `accrow_ = true`, this is the first column, otherwise it is the first row.
     * @param secondary_length Number of elements on the secondary dimension to be extracted.
     * If `accrow_ = true`, this is the number of columns, otherwise it is the number of rows.
     * This is guaranteed to be positive.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output_values Vector of vectors in which to store the output values.
     * The outer vector is of length no less than `primary_length`; each inner vector corresponds to an element of the primary dimension, starting at `primary_start`.
     * @param[out] output_indices Vector of vectors in which to store the output indices.
     * The outer vector is of length no less than `primary_length`; each inner vector corresponds to an element of the primary dimension, starting at `primary_start`.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * This method extracts the specified values from the chunk into `output`.
     * For example, if `accrow_ = true`, we would extract rows `[primary_start, primary_start + length)` and columns `[secondary_start, secondary_start + secondary_length)`.
     * For a primary dimension index `p` and secondary dimension index `s`, the value from the chunk should be appended to `output_values[p - primary_start]`.
     * Similarly, the secondary index should be increased by `shift` and then appended to `output_indices[p - primary_start]`.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomChunkedMatrix`.
     */
    template<bool accrow_, typename Index_, typename OutputValue_, typename OutputIndex_>
    void extract(
        Index_ primary_start,
        Index_ primary_length,
        Index_ secondary_start,
        Index_ secondary_length,
        Workspace& work,
        std::vector<std::vector<OutputValue_> >& output_values,
        std::vector<std::vector<OutputIndex_> >& output_indices,
        OutputIndex_ shift)
    const {
        chunk.inflate(work.values, work.indices, work.indptrs);
        Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
        Index_ secondary_chunkdim = get_secondary_chunkdim<accrow_>();
        Index_ primary_end = primary_start + primary_length;
        Index_ secondary_end = secondary_start + secondary_length;

        if constexpr(Chunk_::row_major == accrow_) {
            for (Index_ p = primary_start; p < primary_end; ++p) {
                auto start = work.indptrs[p], end = work.indptrs[p + 1];

                if (start < end) {
                    refine_start_and_end(start, end, secondary_start, secondary_end, secondary_chunkdim, work.indices);

                    auto& current_values = output_values[p - primary_start];
                    current_values.insert(current_values.end(), work.values.begin() + start, work.values.begin() + end);

                    auto& current_indices = output_indices[p - primary_start];
                    for (size_t i = start; i < end; ++i) {
                        current_indices.push_back(work.indices[i] + shift);
                    }
                }
            }

        } else {
            for (size_t s = secondary_start; s < secondary_end; ++s) {
                auto start = work.indptrs[s], end = work.indptrs[s + 1];
                refine_start_and_end(start, end, primary_start, primary_end, primary_chunkdim, work.indices);

                for (size_t i = start; i < end; ++i) {
                    auto p = work.indices[i] - primary_start;
                    output_values[p].push_back(work.values[i]);
                    output_indices[p].push_back(s + shift);
                }
            }
        }
    }

    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the chunk.
     * @tparam OutputValue_ Numeric type for the output values.
     * @tparam OutputIndex_ Integer type for the output indices.
     *
     * @param primary_start Index of the first element on the primary dimension to be extracted.
     * If `accrow_ = true`, this is the first row, otherwise it is the first column.
     * @param primary_length Number of elements on the primary dimension to be extracted.
     * If `accrow_ = true`, this is the number of rows, otherwise it is the number of columns.
     * This is guaranteed to be positive.
     * @param secondary_indices Indices of the elements on the secondary dimension to be extracted.
     * If `accrow_ = true`, these are column indices, otherwise these are row indices.
     * This is guaranteed to be non-empty.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output_values Vector of vectors in which to store the output values.
     * The outer vector is of length no less than `primary_length`; each inner vector corresponds to an element of the primary dimension, starting at `primary_start`.
     * @param[out] output_indices Vector of vectors in which to store the output indices.
     * The outer vector is of length no less than `primary_length`; each inner vector corresponds to an element of the primary dimension, starting at `primary_start`.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * This method extracts the specified values from the chunk into `output`.
     * For example, if `accrow_ = true`, we would extract rows `[primary_start, primary_start + length)` and columns in `secondary_indices`.
     * For a primary dimension index `p` and secondary dimension index `secondary_indices[i]`, the value from the chunk should be appended to `output_values[p - primary_start]`.
     * Similarly, the secondary index should be increased by `shift` and then appended to `output_indices[p - primary_start]`.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomChunkedMatrix`.
     */
    template<bool accrow_, typename Index_, typename OutputValue_, typename OutputIndex_>
    void extract(
        Index_ primary_start,
        Index_ primary_length,
        const std::vector<Index_>& secondary_indices,
        Workspace& work,
        std::vector<std::vector<OutputValue_> >& output_values,
        std::vector<std::vector<OutputIndex_> >& output_indices,
        OutputIndex_ shift)
    const {
        chunk.inflate(work.values, work.indices, work.indptrs);
        Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
        Index_ primary_end = primary_start + primary_length;

        if constexpr(Chunk_::row_major == accrow_) {
            for (Index_ p = primary_start; p < primary_end; ++p) {
                auto start = work.indptrs[p], end = work.indptrs[p + 1];

                if (start < end) {
                    if (secondary_indices.front()) {
                        auto it = work.indices.begin();
                        start = std::lower_bound(it + start, it + end, static_cast<typename Chunk_::index_type>(secondary_indices.front())) - it;
                    }

                    auto sIt = secondary_indices.begin();
                    auto& current_values = output_values[p - primary_start];
                    auto& current_indices = output_indices[p - primary_start];

                    for (size_t i = start; i < end; ++i) {
                        Index_ target = work.indices[i];
                        while (sIt != secondary_indices.end() && *sIt < target) {
                            ++sIt;
                        }
                        if (sIt == secondary_indices.end()) {
                            break;
                        }
                        if (*sIt == target) {
                            current_values.push_back(work.values[i]);
                            current_indices.push_back(work.indices[i] + shift);
                            ++sIt;
                        }
                    }
                }
            }

        } else {
            for (auto s : secondary_indices) {
                auto start = work.indptrs[s], end = work.indptrs[s + 1];
                refine_start_and_end(start, end, primary_start, primary_end, primary_chunkdim, work.indices);

                for (size_t i = start; i < end; ++i) {
                    auto p = work.indices[i] - primary_start;
                    output_values[p].push_back(work.values[i]);
                    output_indices[p].push_back(s + shift);
                }
            }
        }
    }
};

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
        num_chunk_sets_in_cache = (chunk_set_size_in_elements ? cache_size_in_elements / chunk_set_size_in_elements : 1);

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
