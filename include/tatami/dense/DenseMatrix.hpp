#ifndef TATAMI_DENSE_MATRIX_H
#define TATAMI_DENSE_MATRIX_H

#include "../base/Matrix.hpp"
#include "SparsifiedWrapper.hpp"
#include "../utils/has_data.hpp"
#include "../utils/PseudoOracularExtractor.hpp"

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
 * @cond
 */
namespace DenseMatrix_internals {

template<typename Value_, typename Index_, class Storage_>
struct PrimaryMyopicFullDense : public MyopicDenseExtractor<Value_, Index_> {
    PrimaryMyopicFullDense(const Storage_& store, Index_ sec) : storage(store), secondary(sec) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        size_t offset = i * static_cast<size_t>(secondary); // cast to size_t to avoid overflow of 'Index_'.
        if constexpr(has_data<Value_, Storage_>::value) {
            return storage.data() + offset;
        } else {
            auto it = storage.begin() + offset;
            std::copy(it, it + secondary, buffer);
            return buffer;
        }
    }

    Index_ number() const {
        return secondary;
    }

private:
    const Storage_& storage;
    Index_ secondary;

public:
    Index_ sparsify_full_length() const {
        return secondary;
    }
};

template<typename Value_, typename Index_, class Storage_>
struct PrimaryMyopicBlockDense : public MyopicDenseExtractor<Value_, Index_> {
    PrimaryMyopicBlockDense(const Storage_& store, size_t sec, Index_ bs, Index_ bl) : 
        storage(store), secondary(sec), block_start(bs), block_length(bl) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        size_t offset = i * static_cast<size_t>(secondary) + block_start; // cast to avoid overflow.
        if constexpr(has_data<Value_, Storage_>::value) {
            return storage.data() + offset;
        } else {
            auto it = storage.begin() + offset;
            std::copy(it, it + block_length, buffer);
            return buffer;
        }
    }

    Index_ number() const {
        return block_length;
    }

private:
    const Storage_& storage;
    size_t secondary;
    Index_ block_start, block_length;

public:
    Index_ sparsify_block_start() const {
        return block_start;
    }

    Index_ sparsify_block_length() const {
        return block_length;
    }
};

template<typename Value_, typename Index_, class Storage_>
struct PrimaryMyopicIndexDense : public MyopicDenseExtractor<Value_, Index_> {
    PrimaryMyopicIndexDense(const Storage_& store, size_t sec, std::vector<Index_> idx) : 
        storage(store), secondary(sec), indices(std::move(idx)) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto copy = buffer;
        size_t offset = i * static_cast<size_t>(secondary); // cast to avoid overflow.
        for (auto x : indices) {
            *copy = storage[offset + x];
            ++copy;
        }
        return buffer;
    }

    Index_ number() const {
        return indices.size();
    }

private:
    const Storage_& storage;
    size_t secondary;
    std::vector<Index_> indices;

public:
    const std::vector<Index_>& sparsify_indices() const {
        return indices;
    }
};

template<typename Value_, typename Index_, class Storage_>
struct SecondaryMyopicFullDense : public MyopicDenseExtractor<Value_, Index_> {
    SecondaryMyopicFullDense(const Storage_& store, Index_ sec, Index_ prim) : 
        storage(store), secondary(sec), primary(prim) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        size_t offset = i; // use size_t to avoid overflow.
        auto copy = buffer;
        for (Index_ x = 0; x < primary; ++x, ++copy, offset += secondary) {
            *copy = storage[offset];
        }
        return buffer;
    }

    Index_ number() const {
        return primary;
    }

private:
    const Storage_& storage;
    Index_ secondary, primary;

public:
    Index_ sparsify_full_length() const {
        return primary;
    }
};

template<typename Value_, typename Index_, class Storage_>
struct SecondaryMyopicBlockDense : public MyopicDenseExtractor<Value_, Index_> {
    SecondaryMyopicBlockDense(const Storage_& store, Index_ sec, Index_ bs, Index_ bl) : 
        storage(store), secondary(sec), block_start(bs), block_length(bl) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        size_t offset = i + block_start * static_cast<size_t>(secondary); // cast to avoid overflow.
        auto copy = buffer;
        for (Index_ x = 0; x < block_length; ++x, ++copy, offset += secondary) {
            *copy = storage[offset];
        }
        return buffer;
    }

    Index_ number() const {
        return block_length;
    }

private:
    const Storage_& storage;
    Index_ secondary;
    Index_ block_start, block_length;

public:
    Index_ sparsify_block_start() const {
        return block_start;
    }

    Index_ sparsify_block_length() const {
        return block_length;
    }
};

template<typename Value_, typename Index_, class Storage_>
struct SecondaryMyopicIndexDense : public MyopicDenseExtractor<Value_, Index_> {
    SecondaryMyopicIndexDense(const Storage_& store, Index_ sec, std::vector<Index_> idx) : 
        storage(store), secondary(sec), indices(std::move(idx)) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto copy = buffer;
        for (auto x : indices) {
            *copy = storage[static_cast<size_t>(secondary) * x + i]; // cast to avoid overflow.
            ++copy;
        }
        return buffer;
    }

    Index_ number() const {
        return indices.size();
    }

private:
    const Storage_& storage;
    Index_ secondary;
    std::vector<Index_> indices;

public:
    const std::vector<Index_>& sparsify_indices() const {
        return indices;
    }
};

}
/**
 * @endcond
 */

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
class DenseMatrix : public Matrix<Value_, Index_> {
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

    bool sparse() const { return false; }

    double sparse_proportion() const { return 0; }

    double prefer_rows_proportion() const { return static_cast<double>(row_); }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

private:
    Index_ primary() const {
        if constexpr(row_) {
            return nrows;
        } else {
            return ncols;
        }
    }

    Index_ secondary() const {
        if constexpr(row_) {
            return ncols;
        } else {
            return nrows;
        }
    }

private:
    auto dense_row_internal(const Options&) const {
        if constexpr(row_) {
            return std::make_unique<DenseMatrix_internals::PrimaryMyopicFullDense<Value_, Index_, Storage_> >(values, secondary());
        } else {
            return std::make_unique<DenseMatrix_internals::SecondaryMyopicFullDense<Value_, Index_, Storage_> >(values, secondary(), primary()); 
        }
    }

    auto dense_row_internal(Index_ block_start, Index_ block_length, const Options&) const {
        if constexpr(row_) {
            return std::make_unique<DenseMatrix_internals::PrimaryMyopicBlockDense<Value_, Index_, Storage_> >(values, secondary(), block_start, block_length);
        } else {
            return std::make_unique<DenseMatrix_internals::SecondaryMyopicBlockDense<Value_, Index_, Storage_> >(values, secondary(), block_start, block_length);
        }
    }

    auto dense_row_internal(std::vector<Index_> indices, const Options&) const {
        if constexpr(row_) {
            return std::make_unique<DenseMatrix_internals::PrimaryMyopicIndexDense<Value_, Index_, Storage_> >(values, secondary(), std::move(indices));
        } else {
            return std::make_unique<DenseMatrix_internals::SecondaryMyopicIndexDense<Value_, Index_, Storage_> >(values, secondary(), std::move(indices));
        }
    }

    auto dense_column_internal(const Options&) const {
        if constexpr(row_) {
            return std::make_unique<DenseMatrix_internals::SecondaryMyopicFullDense<Value_, Index_, Storage_> >(values, secondary(), primary());
        } else {
            return std::make_unique<DenseMatrix_internals::PrimaryMyopicFullDense<Value_, Index_, Storage_> >(values, secondary());
        }
    }

    auto dense_column_internal(Index_ block_start, Index_ block_length, const Options&) const {
        if constexpr(row_) {
            return std::make_unique<DenseMatrix_internals::SecondaryMyopicBlockDense<Value_, Index_, Storage_> >(values, secondary(), block_start, block_length);
        } else {
            return std::make_unique<DenseMatrix_internals::PrimaryMyopicBlockDense<Value_, Index_, Storage_> >(values, secondary(), block_start, block_length);
        }
    }

    auto dense_column_internal(std::vector<Index_> indices, const Options&) const {
        if constexpr(row_) {
            return std::make_unique<DenseMatrix_internals::SecondaryMyopicIndexDense<Value_, Index_, Storage_> >(values, secondary(), std::move(indices));
        } else {
            return std::make_unique<DenseMatrix_internals::PrimaryMyopicIndexDense<Value_, Index_, Storage_> >(values, secondary(), std::move(indices));
        }
    }

public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        return dense_row_internal(opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense_row_internal(block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        return dense_row_internal(std::move(indices), opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        return dense_column_internal(opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense_column_internal(block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        return dense_column_internal(std::move(indices), opt);
    }

public:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(const Options& opt) const {
        auto ptr = dense_row_internal(opt);
        return std::make_unique<MyopicSparsifiedWrapper<DimensionSelectionType::FULL, Value_, Index_, typename decltype(ptr)::element_type> >(std::move(*ptr), opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        auto ptr = dense_row_internal(block_start, block_length, opt);
        return std::make_unique<MyopicSparsifiedWrapper<DimensionSelectionType::BLOCK, Value_, Index_, typename decltype(ptr)::element_type> >(std::move(*ptr), opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& opt) const {
        auto ptr = dense_row_internal(std::move(indices), opt);
        return std::make_unique<MyopicSparsifiedWrapper<DimensionSelectionType::INDEX, Value_, Index_, typename decltype(ptr)::element_type> >(std::move(*ptr), opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(const Options& opt) const {
        auto ptr = dense_column_internal(opt);
        return std::make_unique<MyopicSparsifiedWrapper<DimensionSelectionType::FULL, Value_, Index_, typename decltype(ptr)::element_type> >(std::move(*ptr), opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        auto ptr = dense_column_internal(block_start, block_length, opt);
        return std::make_unique<MyopicSparsifiedWrapper<DimensionSelectionType::BLOCK, Value_, Index_, typename decltype(ptr)::element_type> >(std::move(*ptr), opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options& opt) const {
        auto ptr = dense_column_internal(std::move(indices), opt);
        return std::make_unique<MyopicSparsifiedWrapper<DimensionSelectionType::INDEX, Value_, Index_, typename decltype(ptr)::element_type> >(std::move(*ptr), opt);
    }

public:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense_row(opt));
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_end, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense_row(block_start, block_end, opt));
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense_row(std::move(indices), opt));
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense_column(opt));
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_end, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense_column(block_start, block_end, opt));
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense_column(std::move(indices), opt));
    }

public:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse_row(opt));
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_end, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse_row(block_start, block_end, opt));
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse_row(std::move(indices), opt));
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse_column(opt));
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_end, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse_column(block_start, block_end, opt));
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse_column(std::move(indices), opt));
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
