#ifndef TATAMI_MATRIX_H
#define TATAMI_MATRIX_H

#include "Extractor.hpp"
#include "Options.hpp"
#include <algorithm>
#include <numeric>
#include <memory>

/**
 * @file Matrix.hpp
 *
 * @brief Virtual class for a matrix of some numeric type.
 */

namespace tatami {

/**
 * @brief Virtual class for a matrix with a defined type.
 * 
 * @tparam Value Data value type, should be numeric.
 * @tparam Index Row/column index type, should be integer.
 */
template <typename Value_, typename Index_ = int>
class Matrix {
public:
    /**
     * @cond
     */
    Matrix() = default;
    virtual ~Matrix() = default;
    /**
     * @endcond
     */

    /** 
     * Type of data to be returned by getters.
     */
    typedef Value_ value_type;

    /** 
     * Type of index to be returned by the sparse getters.
     */
    typedef Index_ index_type;

    /*******************************
     **** Basic virtual methods ****
     *******************************/
public:
    /**
     * @return Number of rows.
     */
    virtual Index_ nrow() const = 0;

    /**
     * @return Number of columns.
     */
    virtual Index_ ncol() const = 0;

    /**
     * @return Boolean indicating whether this matrix is sparse.
     *
     * This can be used to choose between `sparse_row()` and `dense_row()` when iterating over the rows (similarly so for the columns).
     */
    virtual bool sparse() const = 0;

    /**
     * @return Approximate proportion of the matrix that is sparse.
     *
     * This is defined as the proportion of matrix elements that lie within sparse submatrices.
     * It is intended for use in `Matrix` representations that consist of combinations of multiple submatrices (e.g., `DelayedBind`),
     * allowing them to derive a suitable value for `sparse()` based on whether most of its submatrices are sparse.
     * (A more granular approach would be to report the density of structural non-zero elements, but this may not be known by all representations at construction time.)
     */
    virtual double sparse_proportion() const = 0;

    /**
     * @return The preferred dimension for extracting values.
     * If `true`, row-wise extraction is preferred; if `false`, column-wise extraction is preferred.
     */
    virtual bool prefer_rows() const = 0;

    /**
     * @return Approximate proportion of the matrix that prefers row-level access.
     *
     * This is defined as the proportion of matrix elements that lie within submatrices that prefer row-level access.
     * It is useful for determining the return value of `prefer_rows()` in combined matrices consisting of both row- and column-preferred submatrices.
     * In such cases, the net preference can be determined based on the combined size of the submatrices for each preference.
     * (A more granular approach would be to report the iteration cost on each dimension, but this is difficult to estimate.)
     */
    virtual double prefer_rows_proportion() const = 0;

    /**
     * @param row Row access if `true`, column access otherwise.
     * @return Whether this matrix's `tatami::Extractor` classes make use of oracle predictions for row (if `row = true`) or column access (otherwise).
     *
     * The output of this method indicates whether callers should construct an oracle for use in `ExtractorBase::set_oracle()`.
     * If `false`, callers should not bother to pass an oracle as it will be ignored.
     */
    virtual bool uses_oracle(bool row) const = 0;

    /**************************************
     **** Dense access virtual methods ****
     **************************************/
public:
    /**
     * @param opt Options for extraction.
     * @return A `FullDenseExtractor` object for dense access to full rows.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const = 0;

    /**
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @param opt Options for extraction.
     * @return A `BlockDenseExtractor` object for dense access to a block of each row.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const = 0;

    /**
     * @param indices Vector of sorted and unique column indices.
     * @param opt Options for extraction.
     * @return A `IndexDenseExtractor` object for dense access to a indexed subset of each row.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const = 0;

    /**
     * @param opt Options for extraction.
     * @return A `FullDenseExtractor` object for dense access to full columns.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const = 0;

    /**
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @param opt Options for extraction.
     * @return A `BlockDenseExtractor` object for dense access to a block of each column.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const = 0;

    /**
     * @param indices Vector of sorted and unique column indices.
     * @param opt Options for extraction.
     * @return A `IndexDenseExtractor` object for dense access to a indexed subset of each column.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const = 0;

    /***************************************
     **** Sparse access virtual methods ****
     ***************************************/
public:
    /**
     * @param opt Options for extraction.
     * @return A `FullSparseExtractor` object for sparse access to full rows.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_row(const Options& opt) const = 0;

    /**
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @param opt Options for extraction.
     * @return A `BlockSparseExtractor` object for sparse access to a block of each row.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options& opt) const = 0;

    /**
     * @param indices Vector of sorted and unique column indices.
     * @param opt Options for extraction.
     * @return A `IndexSparseExtractor` object for sparse access to a indexed subset of each row.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& opt) const = 0;

    /**
     * @param opt Options for extraction.
     * @return A `FullSparseExtractor` object for sparse access to full columns.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_column(const Options& opt) const = 0;

    /**
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @param opt Options for extraction.
     * @return A `BlockSparseExtractor` object for sparse access to a block of each column.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options& opt) const = 0;

    /**
     * @param indices Vector of sorted and unique column indices.
     * @param opt Options for extraction.
     * @return A `IndexSparseExtractor` object for sparse access to a indexed subset of each column.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options& opt) const = 0;

    /***************************************
     **** Dense access default overload ****
     ***************************************/
public:
    /**
     * @return A `FullDenseExtractor` object for dense access to full rows, using default options.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row() const {
        return dense_row(Options());
    }

    /**
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @return A `BlockDenseExtractor` object for dense access to a block of each row, using default options.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length) const {
        return dense_row(block_start, block_length, Options());
    }

    /**
     * @param indices Vector of sorted and unique column indices.
     * @return A `IndexDenseExtractor` object for dense access to a indexed subset of each row, using default options.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices) const {
        return dense_row(std::move(indices), Options());
    }

    /**
     * @return A `FullDenseExtractor` object for dense access to full columns, using default options.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column() const {
        return dense_column(Options());
    }

    /**
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @return A `BlockDenseExtractor` object for dense access to a block of each column, using default options.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length) const {
        return dense_column(block_start, block_length, Options());
    }

    /**
     * @param indices Vector of sorted and unique column indices.
     * @return A `IndexDenseExtractor` object for dense access to a indexed subset of each column, using default options.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices) const {
        return dense_column(std::move(indices), Options()); 
    }

    /*****************************************
     **** Sparse access default overloads ****
     *****************************************/
public:
    /**
     * @return A `FullSparseExtractor` object for sparse access to full rows, using default options.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_row() const {
        return sparse_row(Options());
    }

    /**
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @return A `BlockSparseExtractor` object for sparse access to a block of each row, using default options.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length) const {
        return sparse_row(block_start, block_length, Options());
    }

    /**
     * @param indices Vector of sorted and unique column indices.
     * @return A `IndexSparseExtractor` object for sparse access to a indexed subset of each row, using default options.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices) const {
        return sparse_row(std::move(indices), Options());
    }

    /**
     * @return A `FullSparseExtractor` object for sparse access to full columns, using default options.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_column() const {
        return sparse_column(Options());
    }

    /**
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @return A `BlockSparseExtractor` object for sparse access to a block of each column, using default options.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length) const {
        return sparse_column(block_start, block_length, Options());
    }

    /**
     * @param indices Vector of sorted and unique column indices.
     * @return A `IndexSparseExtractor` object for sparse access to a indexed subset of each column, using default options.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices) const {
        return sparse_column(std::move(indices), Options());
    }
};

/**
 * A convenient shorthand for the most common use case of double-precision matrices.
 */
using NumericMatrix = Matrix<double, int>;

}

#endif
