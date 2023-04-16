#ifndef TATAMI_DENSE_MATRIX_H
#define TATAMI_DENSE_MATRIX_H

#include "VirtualDenseMatrix.hpp"
#include "../utils.hpp"

#include <vector>
#include <algorithm>
#include <stdexcept>

/**
 * @file DenseMatrix.hpp
 *
 * @brief Dense matrix representation.
 *
 * `typedef`s are provided for the usual row- and column-major formats.
 */

namespace tatami {

/**
 * @brief Dense matrix representation.
 *
 * @tparam row_ Whether this is a row-major representation.
 * If `false`, a column-major representation is assumed instead.
 * @tparam Value_ Type of the matrix values.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Storage_ Vector class used to store the matrix values internally.
 * This does not necessarily have to contain `Value_`, as long as the type is convertible to `Value_`.
 * Methods should be available for `size()`, `begin()`, `end()` and `[]`.
 * If a method is available for `data()` that returns a `const Value_*`, it will also be used.
 */
template<bool row_, typename Value_, typename Index_ = int, class Storage_ = std::vector<Value_> >
class DenseMatrix : public VirtualDenseMatrix<Value_, Index_> {
public: 
    /**
     * @param nr Number of rows.
     * @param nc Number of columns.
     * @param source Vector of values, or length equal to the product of `nr` and `nc`.
     */
    DenseMatrix(Index_ nr, Index_ nc, const Storage_& source) : nrows(nr), ncols(nc), values(source) {
        if (nrows * ncols != values.size()) {
            throw std::runtime_error("length of 'values' should be equal to product of 'nrows' and 'ncols'");
        }
        return;
    }

    /**
     * @param nr Number of rows.
     * @param nc Number of columns.
     * @param source Vector of values, or length equal to the product of `nr` and `nc`.
     */
    DenseMatrix(Index_ nr, Index_ nc, Storage_&& source) : nrows(nr), ncols(nc), values(source) {
        if (nrows * ncols != values.size()) {
            throw std::runtime_error("length of 'values' should be equal to product of 'nrows' and 'ncols'");
        }
        return;
    }

private: 
    Index_ nrows, ncols;
    Storage_ values;

public:
    Index_ nrow() const { return nrows; }

    Index_ ncol() const { return ncols; }

    bool prefer_rows() const { return row_; }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

private:
    template<bool accrow_>
    struct DenseBase : public DenseExtractor<Value_, Index_> {
        DenseBase(const DenseMatrix* p, ExtractionOptions<Index_> non_target) : parent(p) {
            this->extracted_limit = non_target.limit.type;
            switch (this->extracted_limit) {
                case DimensionLimitType::NONE:
                    this->extracted_length = (accrow_ ? p->ncols : p->nrows);
                    break;
                case DimensionLimitType::BLOCK:
                    this->extracted_length = non_target.limit.block_length;
                    this->extracted_block = non_target.limit.block_start;
                    break;
                case DimensionLimitType::INDEX:
                    if (non_target.limit.index_start) {
                        this->extracted_length = non_target.limit.index_length;
                        this->extracted_index_internal = non_target.limit.index_start;
                    } else {
                        this->extracted_length = non_target.limit.indices.size();
                        this->extracted_indices = std::move(non_target.limit.indices);
                    }
                    break;
            }
        }

    public:
        const Index_* extracted_index() const { return extracted_index_0(); }
        
        const Value_* fetch(Index_ position, Value_* buffer) {
            if constexpr(row_ == accrow_) {
                switch (this->extracted_limit) {
                    case DimensionLimitType::NONE:
                        return parent->primary<accrow_>(position, buffer, 0, this->extracted_length);
                        break;
                    case DimensionLimitType::BLOCK:
                        return parent->primary<accrow_>(position, buffer, this->extracted_block, this->extracted_block + this->extracted_length);
                        break;
                    case DimensionLimitType::INDEX:
                        return parent->primary<accrow_>(position, buffer, this->extracted_index_0(), this->extracted_length);
                        break;
                }
            } else {
                switch (this->extracted_limit) {
                    case DimensionLimitType::NONE:
                        parent->secondary<accrow_>(position, buffer, 0, this->extracted_length);
                        break;
                    case DimensionLimitType::BLOCK:
                        parent->secondary<accrow_>(position, buffer, this->extracted_block, this->extracted_block + this->extracted_length);
                        break;
                    case DimensionLimitType::INDEX:
                        parent->secondary<accrow_>(position, buffer, this->extracted_index_0(), this->extracted_length);
                        break;
                }
                return buffer;
            }
        }

    private:
        const Index_* extracted_index_internal = NULL;
        std::vector<Index_> extracted_indices;
        const Index_* extracted_index_0() const { return (extracted_index_internal ? extracted_index_internal : extracted_indices.data()); }

        const DenseMatrix* parent;
    };

private:
    template<bool accrow_> 
    Index_ other_dimension() const {
        if constexpr(row_) {
            return ncols;
        } else {
            return nrows;
        }
    }

    template<bool accrow_> 
    const Value_* primary(Index_ x, Value_* buffer, Index_ start, Index_ end) const {
        Index_ shift = x * other_dimension<accrow_>();
        if constexpr(has_data<Value_, Storage_>::value) {
            return values.data() + shift + start;
        } else {
            std::copy(values.begin() + shift + start, values.begin() + shift + end, buffer);
            return buffer;
        }
    }

    template<bool accrow_> 
    void secondary(Index_ x, Value_* buffer, Index_ start, Index_ end) const {
        Index_ dim_secondary = other_dimension<accrow_>();
        auto it = values.begin() + x + start * dim_secondary;
        for (Index_ i = start; i < end; ++i, ++buffer, it += dim_secondary) {
            *buffer = *it; 
        }
        return;
    }

    template<bool accrow_> 
    const Value_* primary(Index_ x, Value_* buffer, const Index_* indices, Index_ length) const {
        auto offset = x * other_dimension<accrow_>();
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = values[indices[i] + offset];
        }
        return buffer;
    }

    template<bool accrow_> 
    void secondary(Index_ x, Value_* buffer, const Index_* indices, Index_ length) const {
        Index_ dim_secondary = other_dimension<accrow_>();        
        for (Index_ i = 0; i < length; ++i, ++buffer) {
            *buffer = values[indices[i] * dim_secondary + x]; 
        }
        return;
    }

public:
    std::unique_ptr<DenseExtractor<Value_, Index_> > dense_row(IterationOptions<Index_>, ExtractionOptions<Index_> eopt) const {
        return std::unique_ptr<DenseExtractor<Value_, Index_> >(new DenseBase<true>(this, std::move(eopt)));
    }

    std::unique_ptr<DenseExtractor<Value_, Index_> > dense_column(IterationOptions<Index_>, ExtractionOptions<Index_> eopt) const {
        return std::unique_ptr<DenseExtractor<Value_, Index_> >(new DenseBase<false>(this, std::move(eopt)));
    }
};

/**
 * Column-major matrix.
 * See `tatami::DenseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_ = int, class Storage_ = std::vector<Value_> >
using DenseColumnMatrix = DenseMatrix<false, Value_, Index_, Storage_>;

/**
 * Row-major matrix.
 * See `tatami::DenseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_ = int, class Storage_ = std::vector<Value_> >
using DenseRowMatrix = DenseMatrix<true, Value_, Index_, Storage_>;

}

#endif
