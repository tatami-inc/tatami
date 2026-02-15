#ifndef TATAMI_FORCED_DENSE_HPP
#define TATAMI_FORCED_DENSE_HPP

#include "../base/Matrix.hpp"

/**
 * @file ForcedDense.hpp
 * @brief Forced dense representation.
 */

namespace tatami {

/**
 * @brief Forced dense representation.
 *
 * @tparam Value_ Type of the matrix values.
 * @tparam Index_ Type of the row/column indices.
 *
 * Override `Matrix::is_sparse()` to always return false.
 * Similarly, `Matrix::is_sparse_proportion()` will always return zero.
 * This is occasionally helpful to indicate that downstream applications should ignore structural sparsity in the input `Matrix`.
 *
 * To illustrate, consider a structurally sparse matrix like a `CompressedSparseMatrix`.
 * When performing calculations on such a matrix, we can take advantage of sparsity and only process its non-zero elements for greater efficiency.
 * However, if the `CompressedSparseMatrix` has a high density of non-zero values, the sparsity-aware algorithm may be less efficient than its dense counterpart,
 * e.g., due to non-contiguous and unpredictable memory look-ups based on the index of the non-zero element.
 * In such cases, it may be better to ignore the sparsity and treat it as a dense matrix, which can be achieved by wrapping the `CompressedSparseMatrix` in a `ForcedDense` layer.
 *
 * Aside from `is_sparse()` and `is_sparse_proportion()`, calls to any other method are simply forwarded to the corresponding method of the input `Matrix`. 
 * In particular, calls to `ForcedDense::sparse()` may return a `SparseRange` where the number of non-zero elements is not equal to the non-target extraction length.
 */
template<typename Value_, typename Index_>
class ForcedDense final : public Matrix<Value_, Index_> {
public: 
    /**
     * @param matrix Matrix to be wrapped.
     */
    ForcedDense(std::shared_ptr<const Matrix<Value_, Index_> > matrix) : my_matrix(std::move(matrix)) {}

private: 
    std::shared_ptr<const Matrix<Value_, Index_> > my_matrix;

public:
    Index_ nrow() const { return my_matrix->nrow(); }

    Index_ ncol() const { return my_matrix->ncol(); }

    bool prefer_rows() const { return my_matrix->prefer_rows(); }

    bool uses_oracle(const bool row) const { return my_matrix->uses_oracle(row); }

    bool is_sparse() const { return false; }

    double is_sparse_proportion() const { return 0; }

    double prefer_rows_proportion() const { return my_matrix->prefer_rows_proportion(); }

    using Matrix<Value_, Index_>::dense;

    using Matrix<Value_, Index_>::sparse;

    /*****************************
     ******* Dense myopic ********
     *****************************/
public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(const bool row, const Options& opt) const {
        return my_matrix->dense(row, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(const bool row, const Index_ block_start, const Index_ block_length, const Options& opt) const {
        return my_matrix->dense(row, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(const bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return my_matrix->dense(row, std::move(indices_ptr), opt); 
    }

    /******************************
     ******* Sparse myopic ********
     ******************************/
public:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(const bool row, const Options& opt) const {
        return my_matrix->sparse(row, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(const bool row, const Index_ block_start, const Index_ block_length, const Options& opt) const {
        return my_matrix->sparse(row, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(const bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return my_matrix->sparse(row, indices_ptr, opt);
    }

    /*******************************
     ******* Dense oracular ********
     *******************************/
public:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Options& opt
    ) const {
        return my_matrix->dense(row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt)
    const {
        return my_matrix->dense(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) const {
        return my_matrix->dense(row, std::move(oracle), std::move(indices_ptr), opt);
    }

    /********************************
     ******* Sparse oracular ********
     ********************************/
public:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Options& opt
    ) const {
        return my_matrix->sparse(row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt)
    const {
        return my_matrix->sparse(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) const {
        return my_matrix->sparse(row, std::move(oracle), std::move(indices_ptr), opt);
    }
};

}

#endif
