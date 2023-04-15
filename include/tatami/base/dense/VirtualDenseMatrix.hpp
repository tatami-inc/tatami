#ifndef TATAMI_VIRTUAL_DENSE_MATRIX_H
#define TATAMI_VIRTUAL_DENSE_MATRIX_H

#include "../ExtractFormat.hpp"
#include "../SparseRange.hpp"
#include "../Options.hpp"
#include "../Matrix.hpp"
#include "../utils.hpp"
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
 * By inheriting from this class, implementers of dense matrices can skip the implementation of irrelevant methods for `Matrix::sparse_row()`, `Matrix::sparse_column_workspace()`, etc.,
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

private:
    struct SparseWrapper : public SparseFormat<Value_, Index_> {
        SparseWrapper(std::unique_ptr<DenseFormat<Value_, Index_> > base) : internal(std::move(base)) {
            this->extracted_limit = internal->extracted_limit;
            this->extracted_length = internal->extracted_length;
            this->extracted_block = internal->extracted_block;
        }

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            const Value_* vout = (needs_value ? internal->fetch(i, vbuffer) : NULL);
            if (needs_index) {
                switch (this->extracted_limit) {
                    case DimensionLimitType::NONE: 
                        std::iota(ibuffer, ibuffer + this->extracted_length, static_cast<Index_>(0));
                        break;
                    case DimensionLimitType::BLOCK: 
                        std::iota(ibuffer, ibuffer + this->extracted_length, static_cast<Index_>(this->extracted_block));
                        break;
                    case DimensionLimitType::INDEX:
                        {
                            auto ptr = internal->extracted_index();
                            std::copy(ptr, ptr + this->extracted_length, ibuffer);
                        }
                        break;
                }
            } else {
                ibuffer = NULL;
            }
            return SparseRange<Value_, Index_>(this->extracted_length, vout, ibuffer);
        }

        const Index_* extracted_index() const {
            return internal->extracted_index();
        }

        bool needs_value = false;
        bool needs_index = false;
        std::unique_ptr<DenseFormat<Value_, Index_> > internal;
    };

    template<bool accrow_>
    std::unique_ptr<SparseFormat<Value_, Index_> > populate(DimensionLimit<Index_> row_limits, DimensionLimit<Index_> column_limits, ExtractOptions options) const {
        bool needs_index = options.sparse_extract_index;
        bool needs_value = options.sparse_extract_value;

        auto ptr = new SparseWrapper(new_extractor<accrow_, false>(this, std::move(row_limits), std::move(column_limits), std::move(options)));
        std::unique_ptr<SparseFormat<Value_, Index_> > output(ptr);

        ptr->needs_index = needs_index;
        ptr->needs_value = needs_value;
        return output;
    }

public:
    std::unique_ptr<SparseFormat<Value_, Index_> > sparse_row(DimensionLimit<Index_> row_limits, DimensionLimit<Index_> column_limits, ExtractOptions options) const {
        return populate<true>(std::move(row_limits), std::move(column_limits), std::move(options));
    }

    std::unique_ptr<SparseFormat<Value_, Index_> > sparse_column(DimensionLimit<Index_> row_limits, DimensionLimit<Index_> column_limits, ExtractOptions options) const {
        return populate<false>(std::move(row_limits), std::move(column_limits), std::move(options));
    }
};

}

#endif
