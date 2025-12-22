#ifndef TATAMI_MATRIX_H
#define TATAMI_MATRIX_H

#include "Extractor.hpp"
#include "Options.hpp"
#include "Oracle.hpp"
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
 * @brief Virtual class for a matrix. 
 * 
 * @tparam Value Data value type, should be numeric.
 * @tparam Index Row/column index type, should be integer.
 *
 * Interface for a matrix in the **tatami** library.
 * This declares methods to iterate through the matrix by row or column, extracting data in either dense or sparse form.
 * Check out `DenseMatrix` and `CompressedSparseMatrix` for examples of concrete subclasses.
 *
 * To access the matrix data, the `Matrix` methods first return an instance of an extractor class like `MyopicDenseExtractor`, which can then be used to retrieve the matrix contents.
 * Creation of the extractor depends on a few parameters:
 *
 * - Choice of the "target" dimension.
 *   The "target" dimension is defined as the one that is iterated over/indexed into, while the "non-target" dimension is the other dimension.
 *   For example, if we were iterating row-wise through a matrix, the rows would constitute the target dimension, while the columns would be the non-target dimension.
 *   An element of the target dimension is obtained by indexing into that dimension, e.g., if the rows are the target dimension, then any particular row is an element of the target dimension.
 * - Whether to restrict the non-target dimension (see `tatami::DimensionSelectionType`).
 *   We can choose to extract the full extent of the non-target dimension, a contiguous block, or an indexed subset.
 *   For example, if we were iterating row-wise through a matrix, we might only be interested in a subset of columns.
 * - Whether the order of accesses on the target dimension are known.
 *   If so, we can potentially improve efficiency by supplying an `Oracle` during extractor construction.
 *   For example, if we know we will iterate through the matrix by consecutive rows, we could supply a `ConsecutiveOracle` to allow `Matrix` implementations to optimize for this access pattern.
 * - Whether to obtain the contents for a target dimension element in dense or sparse form.
 *   The dense form is simply a contiguous 1-dimensional array of matrix `Value_`s.
 *   The sparse form is a `SparseRange` describing the structural non-zeros for that dimension element.
 * 
 * `Matrix` subclasses should describe whether they are dense/sparse and if they prefer row or column access.
 * This allows users to choose the best method of extracting data from the matrix.
 */
template <typename Value_, typename Index_>
class Matrix {
public:
    /**
     * @cond
     */
    Matrix() = default;
    Matrix(const Matrix&) = default;
    Matrix& operator=(const Matrix&) = default;
    Matrix(Matrix&&) = default;
    Matrix& operator=(Matrix&&) = default;
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
     *
     * It is expected that the number of rows can be represented by a `std::size_t` without overflow, even if `Index_` is of a larger size.
     * This is implied by the fact that `MyopicDenseExtractor::fetch()` and friends accept a pointer for their buffer arguments;
     * the length of the array referenced by these pointers (possibly equal to the number of rows) must fit in a `std::size_t`.
     */
    virtual Index_ nrow() const = 0;

    /**
     * @return Number of columns.
     *
     * It is expected that the number of columns can be represented by a `std::size_t` without overflow, even if `Index_` is of a larger size.
     * This is implied by the fact that `MyopicDenseExtractor::fetch()` and friends accept a pointer for their buffer arguments;
     * the length of the array referenced by these pointers (possibly equal to the number of rows) must fit in a `std::size_t`.
     */
    virtual Index_ ncol() const = 0;

    /**
     * @return Boolean indicating whether this matrix is sparse.
     *
     * This can be used to choose between dense and sparse outputs.
     */
    virtual bool is_sparse() const = 0;

    /**
     * @cond
     */
    // Back-compatibility only.
    bool sparse() const {
        return is_sparse();
    }

    bool sparse_proportion() const {
        return is_sparse_proportion();
    }
    /**
     * @endcond
     */

    /**
     * @return Approximate proportion of the matrix that is sparse.
     *
     * This is defined as the proportion of matrix elements that lie within sparse submatrices.
     * It is intended for use in `Matrix` representations that consist of combinations of multiple submatrices (e.g., `DelayedBind`),
     * allowing them to derive a suitable value for `is_sparse()` based on whether most of its submatrices are sparse.
     * (A more granular approach would be to report the density of structural non-zero elements, but this may not be known by all representations at construction time.)
     */
    virtual double is_sparse_proportion() const = 0;

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
     * Create an extractor that retrieves the full extent of the non-target dimension in dense form.
     * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
     * @param opt Options for extraction.
     * @return Object for extracting each row (if `row = true`) or column (otherwise) in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, const Options& opt) const = 0;

    /**
     * Create an extractor that retrieves a contiguous block of the non-target dimension in dense form.
     * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
     * @param block_start Index of the column (if `row = true`) or row (otherwise) at the start of the block.
     * @param block_length Number of columns (if `row = true`) or rows (otherwise) in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each row (if `row = true`) or column (otherwise) in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, Index_ block_start, Index_ block_length, const Options& opt) const = 0;

    /**
     * Create an extractor that retrieves an indexed subset of the non-target dimension in dense form.
     * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
     * @param indices_ptr Pointer to a vector of sorted and unique column indices (if `row = true`) or row indices (otherwise).
     * This should be non-NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row (if `row = true`) or column (otherwise) in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const = 0;

public: // ==== Convenience methods ====
    /**
     * Create a row-wise extractor that retrieves all columns in dense form.
     * @param opt Options for extraction.
     * @return Object for extracting each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        return dense(true, opt);
    }

    /**
     * Create a row-wise extractor that retrieves a contiguous block of columns in dense form.
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
     * Create a row-wise extractor that retrieves an indexed subset of columns in dense form.
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
     * Create a row-wise extractor that retrieves an indexed subset of columns in dense form.
     * @param indices Vector of sorted and unique column indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        return dense_row(std::make_shared<std::vector<Index_> >(std::move(indices)), opt);
    }

    /**
     * Create a column-wise extractor that retrieves all rows in dense form.
     * @param opt Options for extraction.
     * @return Object for extracting each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        return dense(false, opt);
    }

    /**
     * Create a column-wise extractor that retrieves a contiguous block of rows in dense form.
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
     * Create a column-wise extractor that retrieves an indexed subset of rows in dense form.
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
     * Create a column-wise extractor that retrieves an indexed subset of rows in dense form.
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
     * @return Object for extracting each row in dense form.
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
     * @return Object for extracting each column in dense form.
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
     * Create an extractor that retrieves the full extent of the non-target dimension in sparse form.
     * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
     * @param opt Options for extraction.
     * @return Object for extracting each row (if `row = true`) or columns (otherwise) in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const Options& opt) const = 0;

    /**
     * Create an extractor that retrieves a contiguous block of the non-target dimension in sparse form.
     * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
     * @param block_start Index of the column (if `row = true`) or row (otherwise) at the start of the block.
     * @param block_length Number of columns (if `row = true`) or rows (otherwise) in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each row (if `row = true`) or column (otherwise) in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_length, const Options& opt) const = 0;

    /**
     * Create an extractor that retrieves an indexed subset of the non-target dimension in sparse form.
     * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
     * @param indices_ptr Pointer to a vector of sorted and unique column indices (if `row = true`) or row indices (otherwise).
     * This should be non-NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row (if `row = true`) or column (otherwise) in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const = 0;

public: // ==== Convenience methods ====
    /**
     * Create a row-wise extractor that retrieves all columns in sparse form.
     * @param opt Options for extraction.
     * @return Object for extracting each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(const Options& opt) const {
        return sparse(true, opt);
    }

    /**
     * Create a row-wise extractor that retrieves a contiguous block of columns in sparse form.
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
     * Create a row-wise extractor that retrieves an indexed subset of columns in sparse form.
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
     * Create a row-wise extractor that retrieves an indexed subset of columns in sparse form.
     * @param indices Vector of sorted and unique column indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& opt) const {
        return sparse_row(std::make_shared<std::vector<Index_> >(std::move(indices)), opt);
    }

    /**
     * Create a column-wise extractor that retrieves all rows in sparse form.
     * @param opt Options for extraction.
     * @return Object for extracting each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(const Options& opt) const {
        return sparse(false, opt);
    }

    /**
     * Create a column-wise extractor that retrieves a contiguous block of rows in sparse form.
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
     * Create a column-wise extractor that retrieves an indexed subset of rows in sparse form.
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
     * Create a column-wise extractor that retrieves an indexed subset of rows in sparse form.
     * @param indices Vector of sorted and unique row indices.
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
     * @return Object for extracting each row in sparse form.
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
     * @return Object for extracting each column in sparse form.
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
     * Create an oracle-aware extractor that retrieves the full extent of the non-target dimension in dense form.
     * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
     * @param oracle An oracle supplying predictions of the next requested row (if `row = true`) or column (otherwise).
     * @param opt Options for extraction.
     * @return Object for extracting each row (if `row = true)` or columns (otherwise) in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const = 0;

    /**
     * Create an oracle-aware extractor that retrieves a contiguous block of the non-target dimension in dense form.
     * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
     * @param oracle An oracle supplying predictions of the next requested row (if `row = true`) or column (otherwise).
     * @param block_start Index of the column (if `row = true`) or row (otherwise) at the start of the block.
     * @param block_length Number of columns (if `row = true`) or rows (otherwise) in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each row (if `row = true`) or column (otherwise) in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const = 0;

    /**
     * Create an oracle-aware extractor that retrieves an indexed subset of the non-target dimension in dense form.
     * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
     * @param oracle An oracle supplying predictions of the next requested row (if `row = true`) or column (otherwise).
     * @param indices_ptr Pointer to a vector of sorted and unique column indices (if `row = true`) or row indices (otherwise).
     * This should not be NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row (if `row = true`) or column (otherwise) in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const = 0;

public: // ==== Convenience methods ====
    /**
     * Create an oracle-aware row-wise extractor that retrieves all columns in dense form.
     * @param oracle An oracle supplying row predictions.
     * @param opt Options for extraction.
     * @return Object for extracting each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return dense(true, std::move(oracle), opt);
    }

    /**
     * Create an oracle-aware row-wise extractor that retrieves a contiguous block of columns in dense form.
     * @param oracle An oracle supplying row predictions.
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each row.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense(true, std::move(oracle), block_start, block_length, opt);
    }

    /**
     * Create an oracle-aware row-wise extractor that retrieves an indexed subset of columns in dense form.
     * @param oracle An oracle supplying row predictions.
     * @param indices_ptr Pointer to a vector of sorted and unique column indices.
     * This should not be NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return dense(true, std::move(oracle), std::move(indices_ptr), opt);
    }

    /**
     * Create an oracle-aware row-wise extractor that retrieves an indexed subset of columns in dense form.
     * @param oracle An oracle supplying row predictions.
     * @param indices Vector of sorted and unique column indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<const Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return dense_row(std::move(oracle), std::make_shared<std::vector<Index_> >(std::move(indices)), opt);
    } 

    /**
     * Create an oracle-aware column-wise extractor that retrieves all rows in dense form.
     * @param oracle An oracle supplying column predictions.
     * @param opt Options for extraction.
     * @return Object for extracting each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return dense(false, std::move(oracle), opt);
    }

    /**
     * Create an oracle-aware column-wise extractor that retrieves a contiguous block of rows in dense form.
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @param oracle An oracle supplying column predictions.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense(false, std::move(oracle), block_start, block_length, opt);
    }

    /**
     * Create an oracle-aware column-wise extractor that retrieves an indexed subset of rows in dense form.
     * @param oracle An oracle supplying column predictions.
     * @param indices_ptr Pointer to a vector of sorted and unique row indices.
     * This should not be NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return dense(false, std::move(oracle), std::move(indices_ptr), opt);
    }

    /**
     * Create an oracle-aware column-wise extractor that retrieves an indexed subset of rows in dense form.
     * @param oracle An oracle supplying column predictions.
     * @param indices Vector of sorted and unique row indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each column in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<const Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return dense_column(std::move(oracle), std::make_shared<std::vector<Index_> >(std::move(indices)), opt);
    }

public: // ==== Default option overloads ====
    /**
     * Overload of `dense_row()` that uses the default options.
     * @param oracle An oracle supplying row predictions.
     * @return Object for extracting each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<const Oracle<Index_> > oracle) const {
        return dense_row(std::move(oracle), Options());
    }

    /**
     * Overload of `dense_row()` that uses the default options.
     * @param oracle An oracle supplying row predictions.
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @return Object for extracting a contiguous block from each row in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length) const {
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
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr) const {
        return dense_row(std::move(oracle), std::move(indices_ptr), Options());
    }

    /**
     * Overload of `dense_row()` that uses the default options.
     * @param indices Vector of sorted and unique column indices.
     * @param oracle An oracle supplying row predictions.
     * @return Object for extracting an indexed subset from each row in dense form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<const Oracle<Index_> > oracle, std::vector<Index_> indices) const {
        return dense_row(std::move(oracle), std::move(indices), Options());
    }

    /**
     * Overload of `dense_column()` that uses the default options.
     * @param oracle An oracle supplying column predictions.
     * @return Object for extracting each column in dense form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<const Oracle<Index_> > oracle) const {
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
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length) const {
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
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr) const {
        return dense_column(std::move(oracle), std::move(indices_ptr), Options()); 
    }

    /**
     * Overload of `dense_column()` that uses the default options.
     * @param indices Vector of sorted and unique row indices.
     * @param oracle An oracle supplying column predictions.
     * @return Object for extracting an indexed subset from each column in dense form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<const Oracle<Index_> > oracle, std::vector<Index_> indices) const {
        return dense_column(std::move(oracle), std::move(indices), Options()); 
    }

    /*********************************
     **** Oracular sparse methods ****
     *********************************/
public:
    /**
     * Create an oracle-aware extractor that retrieves the full extent of the non-target dimension in sparse form.
     * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
     * @param oracle An oracle supplying predictions of the next requested row (if `row = true`) or column (otherwise).
     * @param opt Options for extraction.
     * @return Object for extracting each row (if `row = true)` or columns (otherwise) in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const = 0;

    /**
     * Create an oracle-aware extractor that retrieves a contiguous block of the non-target dimension in sparse form.
     * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
     * @param oracle An oracle supplying predictions of the next requested row (if `row = true`) or column (otherwise).
     * @param block_start Index of the column (if `row = true`) or row (otherwise) at the start of the block.
     * @param block_length Number of columns (if `row = true`) or rows (otherwise) in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each row (if `row = true`) or column (otherwise) in dense form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const = 0;

    /**
     * Create an oracle-aware extractor that retrieves an indexed subset of the non-target dimension in sparse form.
     * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
     * @param oracle An oracle supplying predictions of the next requested row (if `row = true`) or column (otherwise).
     * @param indices_ptr Pointer to a vector of sorted and unique column indices (if `row = true`) or row indices (otherwise).
     * This should not be NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row (if `row = true`) or column (otherwise) in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    virtual std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const = 0;

public: // ==== Convenience methods ====
    /**
     * Create an oracle-aware row-wise extractor that retrieves all columns in sparse form.
     * @param oracle An oracle supplying row predictions.
     * @param opt Options for extraction.
     * @return Object for extracting each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return sparse(true, std::move(oracle), opt);
    }

    /**
     * Create an oracle-aware row-wise extractor that retrieves a contiguous block of columns in sparse form.
     * @param oracle An oracle supplying row predictions.
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse(true, std::move(oracle), block_start, block_length, opt);
    }

    /**
     * Create an oracle-aware row-wise extractor that retrieves an indexed subset of columns in sparse form.
     * @param oracle An oracle supplying row predictions.
     * @param indices_ptr Pointer to a vector of sorted and unique column indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return sparse(true, std::move(oracle), std::move(indices_ptr), opt);
    }

    /**
     * Create an oracle-aware row-wise extractor that retrieves an indexed subset of columns in sparse form.
     * @param oracle An oracle supplying row predictions.
     * @param indices Vector of sorted and unique column indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each row in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<const Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return sparse_row(std::move(oracle), std::make_shared<std::vector<Index_> >(std::move(indices)), opt);
    }

    /**
     * Create an oracle-aware column-wise extractor that retrieves all rows in sparse form.
     * @param oracle An oracle supplying column predictions.
     * @param opt Options for extraction.
     * @return Object for extracting each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return sparse(false, std::move(oracle), opt);
    }

    /**
     * Create an oracle-aware column-wise extractor that retrieves a contiguous block of rows in sparse form.
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @param oracle An oracle supplying column predictions.
     * @param opt Options for extraction.
     * @return Object for extracting a contiguous block from each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse(false, std::move(oracle), block_start, block_length, opt);
    }

    /**
     * Create an oracle-aware column-wise extractor that retrieves an indexed subset of rows in sparse form.
     * @param oracle An oracle supplying column predictions.
     * @param indices_ptr Pointer to a vector of sorted and unique row indices.
     * This should not be NULL.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return sparse(false, std::move(oracle), std::move(indices_ptr), opt);
    }

    /**
     * Create an oracle-aware column-wise extractor that retrieves an indexed subset of rows in sparse form.
     * @param oracle An oracle supplying column predictions.
     * @param indices Vector of sorted and unique row indices.
     * @param opt Options for extraction.
     * @return Object for extracting an indexed subset from each column in sparse form.
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<const Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return sparse_column(std::move(oracle), std::make_shared<std::vector<Index_> >(std::move(indices)), opt);
    }

public: // ==== Default option overloads ====
    /**
     * Overload of `sparse_row()` that uses the default options.
     * @param oracle An oracle supplying row predictions.
     * @return Object for extracting each row in sparse form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<const Oracle<Index_> > oracle) const {
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
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length) const {
        return sparse_row(std::move(oracle), block_start, block_length, Options());
    }

    /**
     * Overload of `sparse_row()` that uses the default options.
     * @param oracle An oracle supplying row predictions.
     * @param indices_ptr Pointer to a vector of sorted and unique column indices.
     * @return Object for extracting an indexed subset from each row in sparse form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr) const {
        return sparse_row(std::move(oracle), std::move(indices_ptr), Options());
    }

    /**
     * Overload of `sparse_row()` that uses the default options.
     * @param oracle An oracle supplying row predictions.
     * @param indices Vector of sorted and unique column indices.
     * @return Object for extracting an indexed subset from each row in sparse form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<const Oracle<Index_> > oracle, std::vector<Index_> indices) const {
        return sparse_row(std::move(oracle), std::move(indices), Options());
    }

    /**
     * Overload of `sparse_column()` that uses the default options.
     * @param oracle An oracle supplying column predictions.
     * @return Object for extracting each column in sparse form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<const Oracle<Index_> > oracle) const {
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
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length) const {
        return sparse_column(std::move(oracle), block_start, block_length, Options());
    }

    /**
     * Overload of `sparse_column()` that uses the default options.
     * @param indices_ptr Pointer to a vector of sorted and unique row indices.
     * @param oracle An oracle supplying column predictions.
     * @return Object for extracting an indexed subset from each column in sparse form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr) const {
        return sparse_column(std::move(oracle), std::move(indices_ptr), Options());
    }

    /**
     * Overload of `sparse_column()` that uses the default options.
     * @param indices Vector of sorted and unique row indices.
     * @param oracle An oracle supplying column predictions.
     * @return Object for extracting an indexed subset from each column in sparse form. 
     * This should not outlive the parent `Matrix` from which it was created.
     */
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<const Oracle<Index_> > oracle, std::vector<Index_> indices) const {
        return sparse_column(std::move(oracle), std::move(indices), Options());
    }
};

/**
 * A convenient shorthand for the most common use case of double-precision matrices.
 */
using NumericMatrix = Matrix<double, int>;

}

#endif
