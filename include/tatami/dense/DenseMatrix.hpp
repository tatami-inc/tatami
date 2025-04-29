#ifndef TATAMI_DENSE_MATRIX_H
#define TATAMI_DENSE_MATRIX_H

#include "../base/Matrix.hpp"
#include "SparsifiedWrapper.hpp"
#include "../utils/has_data.hpp"
#include "../utils/integer_comparisons.hpp"
#include "../utils/PseudoOracularExtractor.hpp"

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <utility>

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
class PrimaryMyopicFullDense final : public MyopicDenseExtractor<Value_, Index_> {
public:
    typedef decltype(std::declval<Storage_>().size()) Size;

    PrimaryMyopicFullDense(const Storage_& storage, Size secondary) : my_storage(storage), my_secondary(secondary) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        Size offset = static_cast<Size>(i) * my_secondary; // cast to the container's size to avoid overflow.
        if constexpr(has_data<Value_, Storage_>::value) {
            return my_storage.data() + offset;
        } else {
            std::copy_n(my_storage.begin() + offset, my_secondary, buffer);
            return buffer;
        }
    }

private:
    const Storage_& my_storage;
    Size my_secondary;
};

template<typename Value_, typename Index_, class Storage_>
class PrimaryMyopicBlockDense final : public MyopicDenseExtractor<Value_, Index_> {
public:
    typedef decltype(std::declval<Storage_>().size()) Size;

    PrimaryMyopicBlockDense(const Storage_& storage, Size secondary, Index_ block_start, Index_ block_length) : 
        my_storage(storage), my_secondary(secondary), my_block_start(block_start), my_block_length(block_length) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        Size offset = static_cast<Size>(i) * my_secondary + my_block_start; // cast to container size to avoid overflow.
        if constexpr(has_data<Value_, Storage_>::value) {
            return my_storage.data() + offset;
        } else {
            std::copy_n(my_storage.begin() + offset, my_block_length, buffer);
            return buffer;
        }
    }

private:
    const Storage_& my_storage;
    Size my_secondary;
    Size my_block_start, my_block_length;
};

template<typename Value_, typename Index_, class Storage_>
class PrimaryMyopicIndexDense final : public MyopicDenseExtractor<Value_, Index_> {
public:
    typedef decltype(std::declval<Storage_>().size()) Size;

    PrimaryMyopicIndexDense(const Storage_& storage, Size secondary, VectorPtr<Index_> indices_ptr) : 
        my_storage(storage), my_secondary(secondary), my_indices_ptr(std::move(indices_ptr)) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        Size offset = static_cast<Size>(i) * my_secondary; // cast to container size to avoid overflow.
        const auto& indices = *my_indices_ptr;
        for (decltype(indices.size()) x = 0, end = indices.size(); x < end; ++x) {
            buffer[x] = my_storage[offset + static_cast<Size>(indices[x])]; // more casting for overflow protection.
        }
        return buffer;
    }

private:
    const Storage_& my_storage;
    Size my_secondary;
    VectorPtr<Index_> my_indices_ptr;
};

template<typename Value_, typename Index_, class Storage_>
class SecondaryMyopicFullDense final : public MyopicDenseExtractor<Value_, Index_> {
public:
    typedef decltype(std::declval<Storage_>().size()) Size;

    SecondaryMyopicFullDense(const Storage_& storage, Index_ secondary, Index_ primary) : 
        my_storage(storage), my_secondary(secondary), my_primary(primary) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        for (decltype(my_primary) x = 0; x < my_primary; ++x) {
            buffer[x] = my_storage[x * my_secondary + static_cast<Size>(i)]; // cast to container size to avoid overflow.
        }
        return buffer;
    }

private:
    const Storage_& my_storage;
    Size my_secondary;
    Size my_primary;
};

template<typename Value_, typename Index_, class Storage_>
class SecondaryMyopicBlockDense final : public MyopicDenseExtractor<Value_, Index_> {
public:
    typedef decltype(std::declval<Storage_>().size()) Size;

    SecondaryMyopicBlockDense(const Storage_& storage, Index_ secondary, Index_ block_start, Index_ block_length) : 
        my_storage(storage), my_secondary(secondary), my_block_start(block_start), my_block_length(block_length) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        Size offset = my_block_start * my_secondary + static_cast<Size>(i); // cast to avoid overflow.
        for (decltype(my_block_length) x = 0; x < my_block_length; ++x) {
            buffer[x] = my_storage[x * my_secondary + offset]; // everything is already Size to avoid overflow.
        }
        return buffer;
    }

private:
    const Storage_& my_storage;
    Size my_secondary;
    Size my_block_start;
    Size my_block_length;
};

template<typename Value_, typename Index_, class Storage_>
class SecondaryMyopicIndexDense final : public MyopicDenseExtractor<Value_, Index_> {
public:
    typedef decltype(std::declval<Storage_>().size()) Size;

    SecondaryMyopicIndexDense(const Storage_& storage, Index_ secondary, VectorPtr<Index_> indices_ptr) : 
        my_storage(storage), my_secondary(secondary), my_indices_ptr(std::move(indices_ptr)) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        Size offset = i;
        const auto& indices = *my_indices_ptr;
        for (decltype(indices.size()) x = 0, end = indices.size(); x < end; ++x) {
            buffer[x] = my_storage[static_cast<Size>(indices[x]) * my_secondary + offset]; // casting to avoid overflow.
        }
        return buffer;
    }

private:
    const Storage_& my_storage;
    Size my_secondary;
    VectorPtr<Index_> my_indices_ptr;
};

}
/**
 * @endcond
 */

/**
 * @brief Dense matrix representation.
 *
 * @tparam Value_ Type of the matrix values.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Storage_ Vector class used to store the matrix values internally.
 * This does not necessarily have to contain `Value_`, as long as the type is convertible to `Value_`.
 * Methods should be available for `size()`, `begin()`, `end()` and `[]`.
 * If a method is available for `data()` that returns a `const Value_*`, it will also be used.
 */
template<typename Value_, typename Index_, class Storage_>
class DenseMatrix : public Matrix<Value_, Index_> {
public: 
    /**
     * @param nrow Number of rows.
     * @param ncol Number of columns.
     * @param values Vector of values of length equal to the product of `nrow` and `ncol`.
     * @param row_major Whether `values` stores the matrix contents in a row-major representation.
     * If `false`, a column-major representation is assumed instead.
     */
    DenseMatrix(Index_ nrow, Index_ ncol, Storage_ values, bool row_major) : my_nrow(nrow), my_ncol(ncol), my_values(std::move(values)), my_row_major(row_major) {
        auto nvalues = my_values.size();
        if (my_nrow == 0) {
            if (nvalues != 0) {
                throw std::runtime_error("length of 'values' should be equal to product of 'nrows' and 'ncols'");
            }
        } else if (!safe_non_negative_equal(nvalues / my_nrow, my_ncol) || nvalues % my_nrow != 0) { // using division instead of 'my_nrow * my_ncol == nvalues' as the product might overflow.
            throw std::runtime_error("length of 'values' should be equal to product of 'nrows' and 'ncols'");
        }
    }

private: 
    Index_ my_nrow, my_ncol;
    Storage_ my_values;
    bool my_row_major;

public:
    Index_ nrow() const { return my_nrow; }

    Index_ ncol() const { return my_ncol; }

    bool prefer_rows() const { return my_row_major; }

    bool uses_oracle(bool) const { return false; }

    bool is_sparse() const { return false; }

    double is_sparse_proportion() const { return 0; }

    double prefer_rows_proportion() const { return static_cast<double>(my_row_major); }

    using Matrix<Value_, Index_>::dense;

    using Matrix<Value_, Index_>::sparse;

private:
    Index_ primary() const {
        if (my_row_major) {
            return my_nrow;
        } else {
            return my_ncol;
        }
    }

    Index_ secondary() const {
        if (my_row_major) {
            return my_ncol;
        } else {
            return my_nrow;
        }
    }

    /*****************************
     ******* Dense myopic ********
     *****************************/
public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, const Options&) const {
        if (my_row_major == row) {
            return std::make_unique<DenseMatrix_internals::PrimaryMyopicFullDense<Value_, Index_, Storage_> >(my_values, secondary());
        } else {
            return std::make_unique<DenseMatrix_internals::SecondaryMyopicFullDense<Value_, Index_, Storage_> >(my_values, secondary(), primary()); 
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, Index_ block_start, Index_ block_length, const Options&) const {
        if (my_row_major == row) { 
            return std::make_unique<DenseMatrix_internals::PrimaryMyopicBlockDense<Value_, Index_, Storage_> >(my_values, secondary(), block_start, block_length);
        } else {
            return std::make_unique<DenseMatrix_internals::SecondaryMyopicBlockDense<Value_, Index_, Storage_> >(my_values, secondary(), block_start, block_length);
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, VectorPtr<Index_> indices_ptr, const Options&) const {
        if (my_row_major == row) {
            return std::make_unique<DenseMatrix_internals::PrimaryMyopicIndexDense<Value_, Index_, Storage_> >(my_values, secondary(), std::move(indices_ptr));
        } else {
            return std::make_unique<DenseMatrix_internals::SecondaryMyopicIndexDense<Value_, Index_, Storage_> >(my_values, secondary(), std::move(indices_ptr));
        }
    }

    /******************************
     ******* Sparse myopic ********
     ******************************/
public:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const Options& opt) const {
        return std::make_unique<FullSparsifiedWrapper<false, Value_, Index_> >(dense(row, opt), (row ? my_ncol : my_nrow), opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<BlockSparsifiedWrapper<false, Value_, Index_> >(dense(row, block_start, block_length, opt), block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        auto ptr = dense(row, indices_ptr, opt);
        return std::make_unique<IndexSparsifiedWrapper<false, Value_, Index_> >(std::move(ptr), std::move(indices_ptr), opt);
    }

    /*******************************
     ******* Dense oracular ********
     *******************************/
public:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense(row, opt));
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_end, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense(row, block_start, block_end, opt));
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense(row, std::move(indices_ptr), opt));
    }

    /********************************
     ******* Sparse oracular ********
     ********************************/
public:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse(row, opt));
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_end, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse(row, block_start, block_end, opt));
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse(row, std::move(indices_ptr), opt));
    }
};

/**
 * @brief Dense column-major matrix.
 *
 * See `tatami::DenseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_, class Storage_ = std::vector<Value_> >
class DenseColumnMatrix final : public DenseMatrix<Value_, Index_, Storage_> {
public:
    /**
     * @param nrow Number of rows.
     * @param ncol Number of columns.
     * @param values Vector of values of length equal to the product of `nr` and `nc`, storing the matrix in column-major format.
     */
    DenseColumnMatrix(Index_ nrow, Index_ ncol, Storage_ values) : DenseMatrix<Value_, Index_, Storage_>(nrow, ncol, std::move(values), false) {}
};

/**
 * @brief Dense row-major matrix.
 *
 * See `tatami::DenseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_, class Storage_ = std::vector<Value_> >
class DenseRowMatrix final : public DenseMatrix<Value_, Index_, Storage_> {
public:
    /**
     * @param nrow Number of rows.
     * @param ncol Number of columns.
     * @param values Vector of values of length equal to the product of `nr` and `nc`, storing the matrix in row-major format.
     */
    DenseRowMatrix(Index_ nrow, Index_ ncol, Storage_ values) : DenseMatrix<Value_, Index_, Storage_>(nrow, ncol, std::move(values), true) {}
};

}

#endif
