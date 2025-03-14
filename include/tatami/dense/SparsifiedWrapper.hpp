#ifndef TATAMI_SPARSIFIED_WRAPPER_HPP
#define TATAMI_SPARSIFIED_WRAPPER_HPP

#include "../base/SparseRange.hpp"
#include "../base/Options.hpp"
#include "../base/Extractor.hpp"

#include <algorithm>

/**
 * @file SparsifiedWrapper.hpp
 *
 * @brief Wrapper classes for sparse extraction from a dense `tatami::Matrix`.
 */

namespace tatami {

/**
 * @brief Wrap a full dense extractor in the sparse interface. 
 *
 * This can be used to quickly implement the full sparse extraction methods for a dense `Matrix` subclass.
 * The dense extraction is performed as usual and every value is treated as a structural non-zero;
 * this wrapper just adds the associated indices to satisfy the sparse extraction interface.
 *
 * @tparam oracle_ Whether an oracle is involved.
 * @tparam Value_ Type of the matrix values.
 * @tparam Index_ Type of the row/column indices.
 */
template<bool oracle_, typename Value_, typename Index_>
class FullSparsifiedWrapper final : public SparseExtractor<oracle_, Value_, Index_> {
public:
    /**
     * @param dense Instance of a dense extractor that retrieves the full extent of the non-target dimension.
     * If `oracle_ = true`, this should be an instance of a `MyopicDenseExtractor` subclass;
     * otherwise it should be an `OracularDenseExtractor` instance.
     * @param extent Extent of the row/column extracted by `d`.
     * @param opt Options for extraction.
     */
    FullSparsifiedWrapper(std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense, Index_ extent, const Options& opt) :
        my_dense(std::move(dense)), 
        my_extent(extent), 
        my_needs_value(opt.sparse_extract_value), 
        my_needs_index(opt.sparse_extract_index) 
    {}

    /**
     * @param i Index of the element to extract on the target dimension, ignored if `oracle_ = true`.
     * @param[in, out] value_buffer See `MyopicSparseExtractor::fetch()` or `OracularSparseExtractor::fetch()`.
     * @param[in, out] index_buffer See `MyopicSparseExtractor::fetch()` or `OracularSparseExtractor::fetch()`.
     * @return Sparse output, see `MyopicSparseExtractor::fetch()` or `OracularSparseExtractor::fetch()`.
     */
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        SparseRange<Value_, Index_> output(my_extent);
        if (my_needs_value) {
            output.value = my_dense->fetch(i, value_buffer);
        }
        if (my_needs_index) {
            std::iota(index_buffer, index_buffer + my_extent, static_cast<Index_>(0));
            output.index = index_buffer;
        }
        return output;
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > my_dense;
    Index_ my_extent;
    bool my_needs_value;
    bool my_needs_index;
};

/**
 * @brief Wrap a block dense extractor in the sparse interface. 
 *
 * This can be used to quickly implement the block sparse extraction methods for a dense `Matrix` subclass.
 * The dense extraction is performed as usual and every value is treated as a structural non-zero;
 * this wrapper just adds the associated indices to satisfy the sparse extraction interface.
 *
 * @tparam oracle_ Whether an oracle is involved.
 * @tparam Value_ Type of the matrix values.
 * @tparam Index_ Type of the row/column indices.
 */
template<bool oracle_, typename Value_, typename Index_>
class BlockSparsifiedWrapper final : public SparseExtractor<oracle_, Value_, Index_> {
public:
    /**
     * @param dense Instance of a dense extractor for a contiguous block of the non-target dimension.
     * If `oracle_ = true`, this should be an instance of a `MyopicDenseExtractor` subclass;
     * otherwise it should be an `OracularDenseExtractor` instance.
     * @param block_start Start of the block extracted by `dense`.
     * Should be the same as that used to construct `dense`.
     * @param block_length Length of the block extracted by `dense`.
     * Should be the same as that used to construct `dense`.
     * @param opt Options for extraction.
     */
    BlockSparsifiedWrapper(std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense, Index_ block_start, Index_ block_length, const Options& opt) :
        my_dense(std::move(dense)), 
        my_block_start(block_start), 
        my_block_length(block_length), 
        my_needs_value(opt.sparse_extract_value), 
        my_needs_index(opt.sparse_extract_index) 
    {}

    /**
     * @param i Index of the element to extract on the target dimension, ignored if `oracle_ = true`.
     * @param[in, out] value_buffer See `MyopicSparseExtractor::fetch()` or `OracularSparseExtractor::fetch()`.
     * @param[in, out] index_buffer See `MyopicSparseExtractor::fetch()` or `OracularSparseExtractor::fetch()`.
     * @return Sparse output, see `MyopicSparseExtractor::fetch()` or `OracularSparseExtractor::fetch()`.
     */
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        SparseRange<Value_, Index_> output(my_block_length, NULL, NULL);
        if (my_needs_value) {
            output.value = my_dense->fetch(i, value_buffer);
        }
        if (my_needs_index) {
            std::iota(index_buffer, index_buffer + my_block_length, my_block_start);
            output.index = index_buffer;
        }
        return output;
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > my_dense;
    Index_ my_block_start, my_block_length;
    bool my_needs_value;
    bool my_needs_index;
};

/**
 * @brief Wrap an indexed dense extractor in the sparse interface. 
 *
 * This can be used to quickly implement the indexed sparse extraction methods for a dense `Matrix` subclass.
 * The dense extraction is performed as usual and every value is treated as a structural non-zero;
 * this wrapper just adds the associated indices to satisfy the sparse extraction interface.
 *
 * @tparam oracle_ Whether an oracle is involved.
 * @tparam Value_ Type of the matrix values.
 * @tparam Index_ Type of the row/column indices.
 */
template<bool oracle_, typename Value_, typename Index_>
class IndexSparsifiedWrapper final : public SparseExtractor<oracle_, Value_, Index_> {
public:
    /**
     * @param dense Instance of a dense extractor for an indexed subset of the non-target dimension.
     * If `oracle_ = true`, this should be an instance of a `MyopicDenseExtractor` subclass;
     * otherwise it should be an `OracularDenseExtractor` instance.
     * @param indices_ptr Pointer to a vector of sorted and unique indices for the non-target dimension.
     * Should be the same as that used to construct `dense`.
     * @param opt Options for extraction.
     */
    IndexSparsifiedWrapper(std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense, VectorPtr<Index_> indices_ptr, const Options& opt) :
        my_dense(std::move(dense)), 
        my_indices_ptr(std::move(indices_ptr)), 
        my_needs_value(opt.sparse_extract_value), 
        my_needs_index(opt.sparse_extract_index) 
    {}

    /**
     * @param i Index of the element to extract on the target dimension, ignored if `oracle_ = true`.
     * @param[in, out] value_buffer See `MyopicSparseExtractor::fetch()` or `OracularSparseExtractor::fetch()`.
     * @param[in, out] index_buffer See `MyopicSparseExtractor::fetch()` or `OracularSparseExtractor::fetch()`.
     * @return Sparse output, see `MyopicSparseExtractor::fetch()` or `OracularSparseExtractor::fetch()`.
     */
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        const auto& ix = *my_indices_ptr;
        SparseRange<Value_, Index_> output(ix.size());
        if (my_needs_value) {
            output.value = my_dense->fetch(i, value_buffer);
        }
        if (my_needs_index) {
            std::copy(ix.begin(), ix.end(), index_buffer);
            output.index = index_buffer;
        }
        return output;
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > my_dense;
    VectorPtr<Index_> my_indices_ptr;
    bool my_needs_value;
    bool my_needs_index;
};

}

#endif
