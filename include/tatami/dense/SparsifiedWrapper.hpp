#ifndef TATAMI_SPARSIFIED_WRAPPER_HPP
#define TATAMI_SPARSIFIED_WRAPPER_HPP

#include "../base/SparseRange.hpp"
#include "../base/Options.hpp"
#include "../base/Extractor.hpp"

#include <algorithm>

/**
 * @file SparsifiedWrapper.hpp
 *
 * @brief Wrapper class for sparse extraction from a dense `Matrix`.
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
struct FullSparsifiedWrapper : public SparseExtractor<oracle_, Value_, Index_> {
    /**
     * @param d Instance of a dense extractor for the full extent of the row/column.
     * If `oracle_ = true`, this should be an instance of a `MyopicDenseExtractor` subclass;
     * otherwise it should be an `OracleDenseExtractor` instance.
     * @param ex Extent of the row/column extracted by `d`.
     * @param opt Options for extraction.
     */
    FullSparsifiedWrapper(std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > d, Index_ ex, const Options& opt) :
        raw(std::move(d)), 
        extent(ex), 
        needs_value(opt.sparse_extract_value), 
        needs_index(opt.sparse_extract_index) 
    {}

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        SparseRange<Value_, Index_> output(extent, NULL, NULL);
        if (needs_value) {
            output.value = raw->fetch(i, vbuffer);
        }
        if (needs_index) {
            std::iota(ibuffer, ibuffer + extent, static_cast<Index_>(0));
            output.index = ibuffer;
        }
        return output;
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > raw;
    Index_ extent;
    bool needs_value;
    bool needs_index;
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
struct BlockSparsifiedWrapper : public SparseExtractor<oracle_, Value_, Index_> {
    /**
     * @param d Instance of a dense extractor for a block of the row/column.
     * If `oracle_ = true`, this should be an instance of a `MyopicDenseExtractor` subclass;
     * otherwise it should be an `OracleDenseExtractor` instance.
     * @param bs Start of the block extracted by `d`.
     * Should be the same as that used to construct `d`.
     * @param bl Length of the block extracted by `d`.
     * Should be the same as that used to construct `d`.
     * @param opt Options for extraction.
     */
    BlockSparsifiedWrapper(std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > d, Index_ bs, Index_ bl, const Options& opt) :
        raw(std::move(d)), 
        block_start(bs), 
        block_length(bl), 
        needs_value(opt.sparse_extract_value), 
        needs_index(opt.sparse_extract_index) 
    {}

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        SparseRange<Value_, Index_> output(block_length, NULL, NULL);
        if (needs_value) {
            output.value = raw->fetch(i, vbuffer);
        }
        if (needs_index) {
            std::iota(ibuffer, ibuffer + block_length, block_start);
            output.index = ibuffer;
        }
        return output;
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > raw;
    Index_ block_start, block_length;
    bool needs_value;
    bool needs_index;
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
struct IndexSparsifiedWrapper : public SparseExtractor<oracle_, Value_, Index_> {
    /**
     * @param d Instance of a dense extractor for a block of the row/column.
     * If `oracle_ = true`, this should be an instance of a `MyopicDenseExtractor` subclass;
     * otherwise it should be an `OracleDenseExtractor` instance.
     * @param ip Pointer to a vector of sorted and unique row/column indices to extract.
     * Should be the same as that used to construct `d`.
     * @param opt Options for extraction.
     */
    IndexSparsifiedWrapper(std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > d, VectorPtr<Index_> ip, const Options& opt) :
        raw(std::move(d)), 
        indices_ptr(std::move(ip)), 
        needs_value(opt.sparse_extract_value), 
        needs_index(opt.sparse_extract_index) 
    {}

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        const auto& ix = *indices_ptr;
        SparseRange<Value_, Index_> output(ix.size(), NULL, NULL);
        if (needs_value) {
            output.value = raw->fetch(i, vbuffer);
        }
        if (needs_index) {
            std::copy(ix.begin(), ix.end(), ibuffer);
            output.index = ibuffer;
        }
        return output;
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > raw;
    VectorPtr<Index_> indices_ptr;
    bool needs_value;
    bool needs_index;
};

}

#endif
