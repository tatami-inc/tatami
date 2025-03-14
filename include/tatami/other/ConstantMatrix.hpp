#ifndef TATAMI_CONSTANT_MATRIX_HPP
#define TATAMI_CONSTANT_MATRIX_HPP

#include "../base/Matrix.hpp"
#include "../dense/SparsifiedWrapper.hpp"
#include "../utils/new_extractor.hpp"

#include <algorithm>

/**
 * @file ConstantMatrix.hpp
 * @brief Matrix of constant values.
 */

namespace tatami {

/**
 * @cond
 */
namespace ConstantMatrix_internal {

template<bool oracle_, typename Value_, typename Index_>
class DenseFiller final : public DenseExtractor<oracle_, Value_, Index_> {
public:
    DenseFiller(Index_ length, Value_ value) : my_length(length), my_value(value) {}

    const Value_* fetch(Index_, Value_* buffer) {
        std::fill_n(buffer, my_length, my_value);
        return buffer;
    }
private:
    Index_ my_length;
    Value_ my_value;
};

template<bool oracle_, typename Value_, typename Index_>
class SparseFiller final : public SparseExtractor<oracle_, Value_, Index_> {
public:
    SparseFiller(const Options& opt) : my_needs_value(opt.sparse_extract_value), my_needs_index(opt.sparse_extract_index) {}

    SparseRange<Value_, Index_> fetch(Index_, Value_* value_buffer, Index_* index_buffer) {
        return SparseRange<Value_, Index_>(0, (my_needs_value ? value_buffer : NULL), (my_needs_index ? index_buffer : NULL));
    }
private:
    bool my_needs_value;
    bool my_needs_index;
};

}
/**
 * @endcond
 */

/**
 * @brief Matrix containing a constant value.
 *
 * This is based on the `ConstantArray` class from the **DelayedArray** Bioconductor package.
 * It is occasionally useful as a placeholder, e.g., to populate a matrix with missing values when no data is available.
 *
 * @tparam Value_ Type of the data.
 * @tparam Index_ Type of the row/column indices.
 */
template<typename Value_, typename Index_>
class ConstantMatrix final : public Matrix<Value_, Index_> {
public:
    /**
     * @param nrow Number of rows.
     * @param ncol Number of columns.
     * @param value The constant value for every data element of this matrix.
     */
    ConstantMatrix(Index_ nrow, Index_ ncol, Value_ value) : my_nrow(nrow), my_ncol(ncol), my_value(value) {}

private:
    Index_ my_nrow, my_ncol;
    Value_ my_value;

public:
    Index_ nrow() const {
        return my_nrow;
    } 

    Index_ ncol() const {
        return my_ncol;
    } 

    bool is_sparse() const {
        return my_value == 0;
    }

    double is_sparse_proportion() const {
        return static_cast<double>(my_value == 0);
    }

    bool prefer_rows() const {
        return true; // pretty much arbitrary here.
    }

    double prefer_rows_proportion() const {
        return 1; // pretty much arbitrary here.
    }

    bool uses_oracle(bool) const {
        return false;
    }

    /********************
     *** Myopic dense ***
     ********************/
private:
    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_internal(bool row, MaybeOracle<oracle_, Index_>, const Options&) const {
        return std::make_unique<ConstantMatrix_internal::DenseFiller<oracle_, Value_, Index_> >(row ? my_ncol : my_nrow, my_value);
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_internal(bool, MaybeOracle<oracle_, Index_>, [[maybe_unused]] Index_ block_start, Index_ block_length, const Options&) const {
        return std::make_unique<ConstantMatrix_internal::DenseFiller<oracle_, Value_, Index_> >(block_length, my_value);
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_internal(bool, MaybeOracle<oracle_, Index_>, VectorPtr<Index_> indices_ptr, const Options&) const {
        return std::make_unique<ConstantMatrix_internal::DenseFiller<oracle_, Value_, Index_> >(indices_ptr->size(), my_value);
    }

public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, const Options& opt) const {
        return dense_internal<false>(row, false, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense_internal<false>(row, false, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return dense_internal<false>(row, false, std::move(indices_ptr), opt);
    }

    /**********************
     *** Oracular dense ***
     **********************/
public:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return dense_internal<true>(row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense_internal<true>(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return dense_internal<true>(row, std::move(oracle), std::move(indices_ptr), opt);
    }

    /*********************
     *** Myopic sparse ***
     *********************/
private:
    template<bool oracle_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_internal(bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) const {
        if (my_value == 0) {
            return std::make_unique<ConstantMatrix_internal::SparseFiller<oracle_, Value_, Index_> >(opt);
        } else {
            return std::make_unique<FullSparsifiedWrapper<oracle_, Value_, Index_> >(dense_internal<oracle_>(row, std::move(oracle), opt), (row ? my_ncol : my_nrow), opt);
        }
    }

    template<bool oracle_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        if (my_value == 0) {
            return std::make_unique<ConstantMatrix_internal::SparseFiller<oracle_, Value_, Index_> >(opt);
        } else {
            return std::make_unique<BlockSparsifiedWrapper<oracle_, Value_, Index_> >(dense_internal<oracle_>(row, std::move(oracle), block_start, block_length, opt), block_start, block_length, opt);
        }
    }

    template<bool oracle_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_internal(bool, MaybeOracle<oracle_, Index_>, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        if (my_value == 0) {
            return std::make_unique<ConstantMatrix_internal::SparseFiller<oracle_, Value_, Index_> >(opt);
        } else {
            auto host = std::make_unique<ConstantMatrix_internal::DenseFiller<oracle_, Value_, Index_> >(indices_ptr->size(), my_value);
            return std::make_unique<IndexSparsifiedWrapper<oracle_, Value_, Index_> >(std::move(host), std::move(indices_ptr), opt);
        }
    }

public:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const Options& opt) const {
        return sparse_internal<false>(row, false, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse_internal<false>(row, false, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return sparse_internal<false>(row, false, std::move(indices_ptr), opt);
    }

    /**********************
     *** Oracular sparse ***
     **********************/
public:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return sparse_internal<true>(row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse_internal<true>(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return sparse_internal<true>(row, std::move(oracle), std::move(indices_ptr), opt);
    }
};

}

#endif
