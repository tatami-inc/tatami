#ifndef TATAMI_DENSE_MATRIX_H
#define TATAMI_DENSE_MATRIX_H

#include "VirtualDenseMatrix.hpp"
#include "../base/utils.hpp"

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
        check_dimensions(nr, nc, values.size());
        return;
    }

    /**
     * @param nr Number of rows.
     * @param nc Number of columns.
     * @param source Vector of values, or length equal to the product of `nr` and `nc`.
     */
    DenseMatrix(Index_ nr, Index_ nc, Storage_&& source) : nrows(nr), ncols(nc), values(source) {
        check_dimensions(nr, nc, values.size());
        return;
    }

private: 
    Index_ nrows, ncols;
    Storage_ values;

    static void check_dimensions(size_t nr, size_t nc, size_t expected) { // cast to size_t is deliberate to avoid overflow on Index_ on product.
        if (nr * nc != expected) {
            throw std::runtime_error("length of 'values' should be equal to product of 'nrows' and 'ncols'");
        }
    }

public:
    Index_ nrow() const { return nrows; }

    Index_ ncol() const { return ncols; }

    bool prefer_rows() const { return row_; }

    bool uses_oracle(bool) const { return false; }

    double prefer_rows_proportion() const { return static_cast<double>(row_); }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

private:
    template<bool accrow_, DimensionSelectionType selection_>
    struct DenseBase : public Extractor<selection_, false, Value_, Index_> {
        DenseBase(const DenseMatrix* p) : parent(p) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = (accrow_ ? p->ncol() : p->nrow());
            }
        }

        DenseBase(const DenseMatrix* p, Index_ bs, Index_ bl) : parent(p) {
            if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = bs;
                this->block_length = bl;
            }
        }

        DenseBase(const DenseMatrix* p, std::vector<Index_> idx) : parent(p) {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                indices = std::move(idx);
                this->index_length = indices.size();
            }
        }

    public:
        const Index_* index_start() const {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                return indices.data();
            } else {
                return NULL;
            }
        }

        void set_oracle(std::unique_ptr<Oracle<Index_> >) {
            return;
        }

    public:
        const Value_* fetch(Index_ position, Value_* buffer) {
            if constexpr(row_ == accrow_) {
                if constexpr(selection_ == DimensionSelectionType::FULL) {
                    return parent->primary<accrow_>(position, buffer, static_cast<Index_>(0), this->full_length);
                } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                    return parent->primary<accrow_>(position, buffer, this->block_start, this->block_start + this->block_length);
                } else {
                    return parent->primary<accrow_>(position, buffer, indices.data(), this->index_length);
                }
            } else {
                if constexpr(selection_ == DimensionSelectionType::FULL) {
                    parent->secondary<accrow_>(position, buffer, static_cast<Index_>(0), this->full_length);
                } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                    parent->secondary<accrow_>(position, buffer, this->block_start, this->block_start + this->block_length);
                } else {
                    parent->secondary<accrow_>(position, buffer, indices.data(), this->index_length);
                }
                return buffer;
            }
        }

    private:
        const DenseMatrix* parent;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type indices;
    };

private:
    template<bool accrow_> 
    size_t other_dimension() const { // deliberate cast to avoid integer overflow on Index_ when multiplying to compute offsets.
        if constexpr(row_) {
            return ncols;
        } else {
            return nrows;
        }
    }

    template<bool accrow_> 
    const Value_* primary(Index_ x, Value_* buffer, Index_ start, Index_ end) const {
        size_t shift = x * other_dimension<accrow_>();
        if constexpr(has_data<Value_, Storage_>::value) {
            return values.data() + shift + start;
        } else {
            std::copy(values.begin() + shift + start, values.begin() + shift + end, buffer);
            return buffer;
        }
    }

    template<bool accrow_> 
    void secondary(Index_ x, Value_* buffer, Index_ start, Index_ end) const {
        size_t dim_secondary = other_dimension<accrow_>();
        auto it = values.begin() + x + start * dim_secondary;
        for (Index_ i = start; i < end; ++i, ++buffer, it += dim_secondary) {
            *buffer = *it; 
        }
        return;
    }

    template<bool accrow_> 
    const Value_* primary(Index_ x, Value_* buffer, const Index_* index_start, Index_ index_length) const {
        size_t offset = x * other_dimension<accrow_>();
        for (Index_ i = 0; i < index_length; ++i) {
            buffer[i] = values[index_start[i] + offset];
        }
        return buffer;
    }

    template<bool accrow_> 
    void secondary(Index_ x, Value_* buffer, const Index_* index_start, Index_ index_length) const {
        size_t dim_secondary = other_dimension<accrow_>();        
        for (Index_ i = 0; i < index_length; ++i, ++buffer) {
            *buffer = values[index_start[i] * dim_secondary + x]; 
        }
        return;
    }

public:
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        auto ptr = new DenseBase<true, DimensionSelectionType::FULL>(this);
        return std::unique_ptr<FullDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        auto ptr = new DenseBase<true, DimensionSelectionType::BLOCK>(this, block_start, block_length);
        return std::unique_ptr<BlockDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        auto ptr = new DenseBase<true, DimensionSelectionType::INDEX>(this, std::move(indices));
        return std::unique_ptr<IndexDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        auto ptr = new DenseBase<false, DimensionSelectionType::FULL>(this);
        return std::unique_ptr<FullDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        auto ptr = new DenseBase<false, DimensionSelectionType::BLOCK>(this, block_start, block_length);
        return std::unique_ptr<BlockDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        auto ptr = new DenseBase<false, DimensionSelectionType::INDEX>(this, std::move(indices));
        return std::unique_ptr<IndexDenseExtractor<Value_, Index_> >(ptr);
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
