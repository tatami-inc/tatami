#ifndef TATAMI_VIRTUAL_DENSE_MATRIX_H
#define TATAMI_VIRTUAL_DENSE_MATRIX_H

#include "../base/Matrix.hpp"
#include "../base/utils.hpp"
#include <algorithm>
#include <numeric>

/**
 * @file VirtualDenseMatrix.hpp
 *
 * @brief Virtual class for a dense matrix of some numeric type.
 */

namespace tatami {

/**
 * @brief Virtual class for a dense matrix with a defined type.
 * 
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * This virtual class provides default methods for sparse extraction that just wrap the dense methods.
 * By inheriting from this class, implementers of dense matrices can skip the implementation of irrelevant methods for `Matrix::sparse_row()`, `Matrix::sparse_column()`, etc.,
 */
template <typename Value_, typename Index_ = int>
class VirtualDenseMatrix : public Matrix<Value_, Index_> {
protected:
    VirtualDenseMatrix() = default;

public:
    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

    bool sparse() const { return false; }

    double sparse_proportion() const { return 0; }

private:
    template<DimensionSelectionType selection_>
    struct SparseWrapper : public Extractor<selection_, true, Value_, Index_> {
        SparseWrapper(std::unique_ptr<Extractor<selection_, false, Value_, Index_> > base, bool nv, bool ni) : 
            internal(std::move(base)), needs_value(nv), needs_index(ni) 
        {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = internal->full_length;
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = internal->block_start;
                this->block_length = internal->block_length;
            } else if constexpr(selection_ == DimensionSelectionType::INDEX) {
                this->index_length = internal->index_length;
            }
        }

    public:
        const Index_* index_start() const {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                return internal->index_start();
            } else {
                return NULL;
            }
        }

        void set_oracle(std::unique_ptr<Oracle<Index_> > o) {
            internal->set_oracle(std::move(o));
        }

    public:
        SparseRange<Value_, Index_> fetch(Index_ position, Value_* vbuffer, Index_* ibuffer) {
            const Value_* vout = (needs_value ? internal->fetch(position, vbuffer) : NULL);
            if (needs_index) {
                if constexpr(selection_ == DimensionSelectionType::FULL) {
                    std::iota(ibuffer, ibuffer + this->full_length, static_cast<Index_>(0));
                } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                    std::iota(ibuffer, ibuffer + this->block_length, static_cast<Index_>(this->block_start));
                } else {
                    auto ptr = internal->index_start();
                    std::copy(ptr, ptr + this->index_length, ibuffer);
                }
            } else {
                ibuffer = NULL;
            }
            return SparseRange<Value_, Index_>(extracted_length<selection_, Index_>(*this), vout, ibuffer);
        }

    protected:
        std::unique_ptr<Extractor<selection_, false, Value_, Index_> > internal;
        bool needs_value = false;
        bool needs_index = false;
    };

    typedef SparseWrapper<DimensionSelectionType::FULL> FullSparseWrapper;
    typedef SparseWrapper<DimensionSelectionType::BLOCK> BlockSparseWrapper;
    typedef SparseWrapper<DimensionSelectionType::INDEX> IndexSparseWrapper;

public:
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_row(const Options& opt) const {
        auto ptr = new FullSparseWrapper(this->dense_row(opt), opt.sparse_extract_value, opt.sparse_extract_index);
        return std::unique_ptr<FullSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        auto ptr = new BlockSparseWrapper(this->dense_row(block_start, block_length, opt), opt.sparse_extract_value, opt.sparse_extract_index);
        return std::unique_ptr<BlockSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& opt) const {
        auto ptr = new IndexSparseWrapper(this->dense_row(std::move(indices), opt), opt.sparse_extract_value, opt.sparse_extract_index);
        return std::unique_ptr<IndexSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_column(const Options& opt) const {
        auto ptr = new FullSparseWrapper(this->dense_column(opt), opt.sparse_extract_value, opt.sparse_extract_index);
        return std::unique_ptr<FullSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        auto ptr = new BlockSparseWrapper(this->dense_column(block_start, block_length, opt), opt.sparse_extract_value, opt.sparse_extract_index);
        return std::unique_ptr<BlockSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options& opt) const {
        auto ptr = new IndexSparseWrapper(this->dense_column(std::move(indices), opt), opt.sparse_extract_value, opt.sparse_extract_index);
        return std::unique_ptr<IndexSparseExtractor<Value_, Index_> >(ptr);
    }
};

}

#endif
