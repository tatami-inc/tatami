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
 * @tparam Index Row/column index type, should be integer.
 *
 * Pointer to a vector, typically containing unique and sorted indices.
 * We use a shared pointer so that we can cheaply re-use the same sequence of indices for multiple `Matrix` objects.
 */
template<typename Index_>
using VectorPtr = std::shared_ptr<const std::vector<Index_> >;

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

    /******************************
     **** Myopic dense methods ****
     ******************************/
public:
    /**
     * @param row Whether to create a row-wise extractor.
     * @param opt Options for extraction.
     * @return Object for extracting the full extent of each row (if `row = true`) or column (otherwise) in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, const Options& opt) const = 0;

    /**
     * @param row Whether to create a row-wise extractor.
     * @param block_start Index of the column (if `row = true`) or row (otherwise) at the start of the block.
     * @param block_length Number of columns (if `row = true`) or rows (otherwise) in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each row (if `row = true`) or column (otherwise) in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, Index_ block_start, Index_ block_length, const Options& opt) const = 0;

    /**
     * @param row Whether to create a row-wise extractor.
     * @param indices_ptr Pointer to a vector of sorted and unique column indices (if `row = true`) or row indices (otherwise).
     * This should be non-NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row (if `row = true`) or column (otherwise) in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const = 0;

public: // ==== Convenience methods ====
    /**
     * @param opt Options for extraction.
     * @return Object for extracting the full extent of each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        return dense(true, opt);
    }

    /**
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense(true, block_start, block_length, opt);
    }

    /**
     * @param indices_ptr Pointer to a vector of sorted and unique column indices.
     * This should not be NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return dense(true, std::move(indices_ptr), opt);
    }

    /**
     * @param indices Vector of sorted and unique column indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        return dense_row(std::make_shared<std::vector<Index_> >(std::move(indices)), opt);
    }

    /**
     * @param opt Options for extraction.
     * @return Object for extracting the full extent of each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        return dense(false, opt);
    }

    /**
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense(false, block_start, block_length, opt);
    }

    /**
     * @param indices_ptr Pointer to a vector of sorted and unique row indices.
     * This should not be NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return dense(false, std::move(indices_ptr), opt);
    }

    /**
     * @param indices Vector of sorted and unique row indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        return dense_column(std::make_shared<std::vector<Index_> >(std::move(indices)), opt);
    }

public: // ==== Default option overloads ====
    /**
     * Overload of `dense_row()` that uses the default options.
     * @return Object for extracting the full extent of each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row() const {
        return dense_row(Options());
    }

    /**
     * Overload of `dense_row()` that uses the default options.
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @return Object for extracting a contiguous block from each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length) const {
        return dense_row(block_start, block_length, Options());
    }

    /**
     * Overload of `dense_row()` that uses the default options.
     * @param indices_ptr Pointer to a vector of sorted and unique column indices.
     * This should be non-NULL.
     * @return Object for extracting an indexed subset from each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(VectorPtr<Index_> indices_ptr) const {
        return dense_row(std::move(indices_ptr), Options());
    }

    /**
     * Overload of `dense_row()` that uses the default options.
     * @param indices Vector of sorted and unique column indices.
     * @return Object for extracting an indexed subset from each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices) const {
        return dense_row(std::move(indices), Options());
    }

    /**
     * Overload of `dense_column()` that uses the default options.
     * @return Object for extracting the full extent of each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column() const {
        return dense_column(Options());
    }

    /**
     * Overload of `dense_column()` that uses the default options.
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @return Object for extracting a contiguous block from each column in dense form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length) const {
        return dense_column(block_start, block_length, Options());
    }

    /**
     * Overload of `dense_column()` that uses the default options.
     * @param indices_ptr Pointer to a vector of sorted and unique row indices.
     * This should be non-NULL.
     * @return Obejct for extracting an indexed subset from each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(VectorPtr<Index_> indices_ptr) const {
        return dense_column(std::move(indices_ptr), Options()); 
    }

    /**
     * Overload of `dense_column()` that uses the default options.
     * @param indices Vector of sorted and unique row indices.
     * @return Obejct for extracting an indexed subset from each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices) const {
        return dense_column(std::move(indices), Options()); 
    }

    /******************************
     **** Myopic sparse access ****
     ******************************/
public:
    /**
     * @param row Whether to create a row-wise extractor.
     * @param opt Options for extraction.
     * @return Object for extracting the full extent of each row (if `row = true`) or columns (otherwise) in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const Options& opt) const = 0;

    /**
     * @param row Whether to create a row-wise extractor.
     * @param block_start Index of the column (if `row = true`) or row (otherwise) at the start of the block.
     * @param block_length Number of columns (if `row = true`) or rows (otherwise) in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each row (if `row = true`) or column (otherwise) in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_length, const Options& opt) const = 0;

    /**
     * @param row Whether to create a row-wise extractor.
     * @param indices_ptr Pointer to a vector of sorted and unique column indices (if `row = true`) or row indices (otherwise).
     * This should be non-NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row (if `row = true`) or column (otherwise) in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const = 0;

public: // ==== Convenience methods ====
    /**
     * @param opt Options for extraction.
     * @return Object for extracting the full extent of each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(const Options& opt) const {
        return sparse(true, opt);
    }

    /**
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse(true, block_start, block_length, opt);
    }

    /**
     * @param indices_ptr Pointer to a vector of sorted and unique column indices.
     * This should not be NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return sparse(true, std::move(indices_ptr), opt);
    }

    /**
     * @param indices Vector of sorted and unique column indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& opt) const {
        return sparse_row(std::make_shared<std::vector<Index_> >(std::move(indices)), opt);
    }

    /**
     * @param opt Options for extraction.
     * @return Object for extracting the full extent of each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(const Options& opt) const {
        return sparse(false, opt);
    }

    /**
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse(false, block_start, block_length, opt);
    }

    /**
     * @param indices_ptr Pointer to a vector of sorted and unique row indices.
     * This should not be NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return sparse(false, std::move(indices_ptr), opt);
    }

    /**
     * @param indices_ptr Vector of sorted and unique row indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options& opt) const {
        return sparse_column(std::make_shared<std::vector<Index_> >(std::move(indices)), opt);
    }

public: // ==== Default option overloads ====
    /**
     * Overload of `sparse_row()` that uses the default options.
     * @return Object for extracting the full extent of each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row() const {
        return sparse_row(Options());
    }

    /**
     * Overload of `sparse_row()` that uses the default options.
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @return Object for extracting a contiguous block from each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length) const {
        return sparse_row(block_start, block_length, Options());
    }

    /**
     * Overload of `sparse_row()` that uses the default options.
     * @param indices_ptr Pointer to a vector of sorted and unique column indices.
     * This should not be NULL.
     * @return Object for extracting an indexed subset from each row in sparse form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(VectorPtr<Index_> indices_ptr) const {
        return sparse_row(std::move(indices_ptr), Options());
    }

    /**
     * Overload of `sparse_row()` that uses the default options.
     * @param indices Vector of sorted and unique column indices.
     * @return Object for extracting an indexed subset from each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices) const {
        return sparse_row(std::move(indices), Options());
    }

    /**
     * Overload of `sparse_column()` that uses the default options.
     * @return Object for extracting the full extent of each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column() const {
        return sparse_column(Options());
    }

    /**
     * Overload of `sparse_column()` that uses the default options.
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @return Object for extracting a contiguous block from each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length) const {
        return sparse_column(block_start, block_length, Options());
    }

    /**
     * Overload of `sparse_column()` that uses the default options.
     * @param indices_ptr Pointer to a vector of sorted and unique row indices.
     * This should not be NULL.
     * @return Object for extracting an indexed subset from each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(VectorPtr<Index_> indices_ptr) const {
        return sparse_column(std::move(indices_ptr), Options());
    }

    /**
     * Overload of `sparse_column()` that uses the default options.
     * @param indices Vector of sorted and unique row indices.
     * @return Object for extracting an indexed subset from each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices) const {
        return sparse_column(std::move(indices), Options());
    }

    /*******************************
     **** Oracular dense access ****
     *******************************/
public:
    /**
     * @param row Whether to create a row-wise extractor.
     * @param oracle An oracle supplying predictions of the next requested row (if `row = true`) or column (otherwise).
     * @param opt Options for extraction.
     * @return Object for extracting the full extent of each row (if `row = true)` or columns (otherwise) in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const = 0;

    /**
     * @param row Whether to create a row-wise extractor.
     * @param oracle An oracle supplying predictions of the next requested row (if `row = true`) or column (otherwise).
     * @param block_start Index of the column (if `row = true`) or row (otherwise) at the start of the block.
     * @param block_length Number of columns (if `row = true`) or rows (otherwise) in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each row (if `row = true`) or column (otherwise) in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const = 0;

    /**
     * @param row Whether to create a row-wise extractor.
     * @param oracle An oracle supplying predictions of the next requested row (if `row = true`) or column (otherwise).
     * @param indices_ptr Pointer to a vector of sorted and unique column indices (if `row = true`) or row indices (otherwise).
     * This should not be NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row (if `row = true`) or column (otherwise) in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const = 0;

public: // ==== Convenience methods ====
    /**
     * @param oracle An oracle supplying row predictions.
     * @param opt Options for extraction.
     * @return Object for extracting the full extent of each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return dense(true, std::move(oracle), opt);
    }

    /**
     * @param oracle An oracle supplying row predictions.
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each row.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense(true, std::move(oracle), block_start, block_length, opt);
    }

    /**
     * @param oracle An oracle supplying row predictions.
     * @param indices_ptr Pointer to a vector of sorted and unique column indices.
     * This should not be NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return dense(true, std::move(oracle), std::move(indices_ptr), opt);
    }

    /**
     * @param oracle An oracle supplying row predictions.
     * @param indices Vector of sorted and unique column indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return dense_row(std::move(oracle), std::make_shared<std::vector<Index_> >(std::move(indices)), opt);
    } 

    /**
     * @param oracle An oracle supplying column predictions.
     * @param opt Options for extraction.
     * @return Object for extracting the full extent of each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return dense(false, std::move(oracle), opt);
    }

    /**
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @param oracle An oracle supplying column predictions.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense(false, std::move(oracle), block_start, block_length, opt);
    }

    /**
     * @param oracle An oracle supplying column predictions.
     * @param indices_ptr Pointer to a vector of sorted and unique row indices.
     * This should not be NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return dense(false, std::move(oracle), std::move(indices_ptr), opt);
    }

    /**
     * @param oracle An oracle supplying column predictions.
     * @param indices Vector of sorted and unique row indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return dense_column(std::move(oracle), std::make_shared<std::vector<Index_> >(std::move(indices)), opt);
    }

public: // ==== Default option overloads ====
    /**
     * Overload of `dense_row()` that uses the default options.
     * @param oracle An oracle supplying row predictions.
     * @return Object for extracting the full extent of each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle) const {
        return dense_row(std::move(oracle), Options());
    }

    /**
     * Overload of `dense_row()` that uses the default options.
     * @param oracle An oracle supplying row predictions.
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @param oracle An oracle supplying row predictions.
     * @return Object for extracting a contiguous block from each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length) const {
        return dense_row(std::move(oracle), block_start, block_length, Options());
    }

    /**
     * Overload of `dense_row()` that uses the default options.
     * @param indices_ptr Pointer to a vector of sorted and unique column indices.
     * This should not be NULL.
     * @param oracle An oracle supplying row predictions.
     * @return Object for extracting an indexed subset from each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr) const {
        return dense_row(std::move(oracle), std::move(indices_ptr), Options());
    }

    /**
     * Overload of `dense_row()` that uses the default options.
     * @param indices Vector of sorted and unique column indices.
     * @param oracle An oracle supplying row predictions.
     * @return Object for extracting an indexed subset from each row in dense form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices) const {
        return dense_row(std::move(oracle), std::move(indices), Options());
    }

    /**
     * Overload of `dense_column()` that uses the default options.
     * @param oracle An oracle supplying column predictions.
     * @return Object for extracting the full extent of each column in dense form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle) const {
        return dense_column(std::move(oracle), Options());
    }

    /**
     * Overload of `dense_column()` that uses the default options.
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @param oracle An oracle supplying column predictions.
     * @return Object for extracting a contiguous block from each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length) const {
        return dense_column(std::move(oracle), block_start, block_length, Options());
    }

    /**
     * Overload of `dense_column()` that uses the default options.
     * @param indices_ptr Pointer to a vector of sorted and unique row indices.
     * This should not be NULL.
     * @param oracle An oracle supplying column predictions.
     * @return Object for extracting an indexed subset from each column in dense form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr) const {
        return dense_column(std::move(oracle), std::move(indices_ptr), Options()); 
    }

    /**
     * Overload of `dense_column()` that uses the default options.
     * @param indices Vector of sorted and unique row indices.
     * @param oracle An oracle supplying column predictions.
     * @return Object for extracting an indexed subset from each column in dense form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices) const {
        return dense_column(std::move(oracle), std::move(indices), Options()); 
    }

    /*********************************
     **** Oracular sparse methods ****
     *********************************/
public:
    /**
     * @param row Whether to create a row-wise extractor.
     * @param oracle An oracle supplying predictions of the next requested row (if `row = true`) or column (otherwise).
     * @param opt Options for extraction.
     * @return Object for extracting the full extent of each row (if `row = true)` or columns (otherwise) in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const = 0;

    /**
     * @param row Whether to create a row-wise extractor.
     * @param oracle An oracle supplying predictions of the next requested row (if `row = true`) or column (otherwise).
     * @param block_start Index of the column (if `row = true`) or row (otherwise) at the start of the block.
     * @param block_length Number of columns (if `row = true`) or rows (otherwise) in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each row (if `row = true`) or column (otherwise) in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const = 0;

    /**
     * @param row Whether to create a row-wise extractor.
     * @param oracle An oracle supplying predictions of the next requested row (if `row = true`) or column (otherwise).
     * @param indices_ptr Pointer to a vector of sorted and unique column indices (if `row = true`) or row indices (otherwise).
     * This should not be NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row (if `row = true`) or column (otherwise) in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const = 0;

public: // ==== Convenience methods ====
    /**
     * @param oracle An oracle supplying row predictions.
     * @param opt Options for extraction.
     * @return Object for extracting the full extent of each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return sparse(true, std::move(oracle), opt);
    }

    /**
     * @param oracle An oracle supplying row predictions.
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse(true, std::move(oracle), block_start, block_length, opt);
    }

    /**
     * @param oracle An oracle supplying row predictions.
     * @param indices_ptr Pointer to a vector of sorted and unique column indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return sparse(true, std::move(oracle), std::move(indices_ptr), opt);
    }

    /**
     * @param oracle An oracle supplying row predictions.
     * @param indices Vector of sorted and unique column indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return sparse_row(std::move(oracle), std::make_shared<std::vector<Index_> >(std::move(indices)), opt);
    }

    /**
     * @param oracle An oracle supplying column predictions.
     * @param opt Options for extraction.
     * @return Object for extracting the full extent of each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return sparse(false, std::move(oracle), opt);
    }

    /**
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @param oracle An oracle supplying column predictions.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse(false, std::move(oracle), block_start, block_length, opt);
    }

    /**
     * @param oracle An oracle supplying column predictions.
     * @param indices_ptr Pointer to a vector of sorted and unique row indices.
     * This should not be NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return sparse(false, std::move(oracle), std::move(indices_ptr), opt);
    }

    /**
     * @param oracle An oracle supplying column predictions.
     * @param indices Vector of sorted and unique row indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return sparse_column(std::move(oracle), std::make_shared<std::vector<Index_> >(std::move(indices)), opt);
    }

public: // ==== Default option overloads ====
    /**
     * Overload of `sparse_row()` that uses the default options.
     * @param oracle An oracle supplying row predictions.
     * @return Object for extracting the full extent of each row in sparse form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle) const {
        return sparse_row(std::move(oracle), Options());
    }

    /**
     * Overload of `sparse_row()` that uses the default options.
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @param oracle An oracle supplying row predictions.
     * @return Object for extracting a contiguous block from each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length) const {
        return sparse_row(std::move(oracle), block_start, block_length, Options());
    }

    /**
     * Overload of `sparse_row()` that uses the default options.
     * @param oracle An oracle supplying row predictions.
     * @param indices_ptr Pointer to a vector of sorted and unique column indices.
     * @return Object for extracting an indexed subset from each row in sparse form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr) const {
        return sparse_row(std::move(oracle), std::move(indices_ptr), Options());
    }

    /**
     * Overload of `sparse_row()` that uses the default options.
     * @param oracle An oracle supplying row predictions.
     * @param indices Vector of sorted and unique column indices.
     * @return Object for extracting an indexed subset from each row in sparse form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices) const {
        return sparse_row(std::move(oracle), std::move(indices), Options());
    }

    /**
     * Overload of `sparse_column()` that uses the default options.
     * @param oracle An oracle supplying column predictions.
     * @return Object for extracting the full extent of each column in sparse form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle) const {
        return sparse_column(std::move(oracle), Options());
    }

    /**
     * Overload of `sparse_column()` that uses the default options.
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @param oracle An oracle supplying column predictions.
     * @return Object for extracting a contiguous block from each column in sparse form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length) const {
        return sparse_column(std::move(oracle), block_start, block_length, Options());
    }

    /**
     * Overload of `sparse_column()` that uses the default options.
     * @param indices_ptr Pointer to a vector of sorted and unique row indices.
     * @param oracle An oracle supplying column predictions.
     * @return Object for extracting an indexed subset from each column in sparse form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr) const {
        return sparse_column(std::move(oracle), std::move(indices_ptr), Options());
    }

    /**
     * Overload of `sparse_column()` that uses the default options.
     * @param indices Vector of sorted and unique row indices.
     * @param oracle An oracle supplying column predictions.
     * @return Object for extracting an indexed subset from each column in sparse form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices) const {
        return sparse_column(std::move(oracle), std::move(indices), Options());
    }
};

/**
 * A convenient shorthand for the most common use case of double-precision matrices.
 */
using NumericMatrix = Matrix<double, int>;

}

#endif
