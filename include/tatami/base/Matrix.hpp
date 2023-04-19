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
protected:
    /**
     * @cond
     */
    Matrix() = default;
    virtual ~Matrix() = default;
    /**
     * @endcond
     */

public:
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
     * @return Is this matrix sparse?
     * Defaults to `false` if no specialized method is provided in derived classes.
     */
    virtual bool sparse() const = 0;

    /**
     * @return The preferred dimension for extracting values.
     * If `true`, row-wise extraction is preferred; if `false`, column-wise extraction is preferred.
     */
    virtual bool prefer_rows() const = 0;

    /**
     * @return A `pair` containing the number of matrix elements that prefer row-level access (`first`) or column-level access (`second`).
     *
     * This method is useful for determining the return value of `prefer_rows()` in combined matrices consisting of both row- and column-preferred submatrices.
     * In such cases, the net preference can be determined based on the combined size of the submatrices for each preference.
     *
     * For simpler matrices, the return value contains the total size of the matrix in one of the `double`s and zero in the other.
     */
    virtual std::pair<double, double> dimension_preference () const {
        double size = static_cast<double>(nrow()) * static_cast<double>(ncol());
        if (prefer_rows()) {
            return std::make_pair(size, 0.0);
        } else {
            return std::make_pair(0.0, size);
        }
    }

    /**************************************
     **** Dense access virtual methods ****
     **************************************/
public:
    /**
     * @param opt Options for extraction.
     * @return A `FullDenseExtractor` object for dense access to full rows.
     */
    virtual std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row(const Options<Index_>& opt) const = 0;

    /**
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @param opt Options for extraction.
     * @return A `BlockDenseExtractor` object for dense access to a block of each row.
     */
    virtual std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_end, const Options<Index_>& opt) const = 0;

    /**
     * @param index_start Pointer to an array containing sorted and unique column indices.
     * @param index_length Length of the array pointed to by `index_start`.
     * @param opt Options for extraction.
     * @return A `IndexDenseExtractor` object for dense access to a indexed subset of each row.
     */
    virtual std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(const Index_* index_start, size_t index_length, const Options<Index_>& opt) const = 0;

    /**
     * @param opt Options for extraction.
     * @return A `FullDenseExtractor` object for dense access to full columns.
     */
    virtual std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column(const Options<Index_>& opt) const = 0;

    /**
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @param opt Options for extraction.
     * @return A `BlockDenseExtractor` object for dense access to a block of each column.
     */
    virtual std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_end, const Options<Index_>& opt) const = 0;

    /**
     * @param index_start Pointer to an array containing sorted and unique row indices.
     * @param index_length Length of the array pointed to by `index_start`.
     * @param opt Options for extraction.
     * @return A `IndexDenseExtractor` object for dense access to a indexed subset of each column.
     */
    virtual std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(const Index_* index_start, size_t index_length, const Options<Index_>& opt) const = 0;

    /***************************************
     **** Sparse access virtual methods ****
     ***************************************/
public:
    /**
     * @param opt Options for extraction.
     * @return A `FullSparseExtractor` object for sparse access to full rows.
     */
    virtual std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_row(const Options<Index_>& opt) const = 0;

    /**
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @param opt Options for extraction.
     * @return A `BlockSparseExtractor` object for sparse access to a block of each row.
     */
    virtual std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_end, const Options<Index_>& opt) const = 0;

    /**
     * @param index_start Pointer to an array containing sorted and unique column indices.
     * @param index_length Length of the array pointed to by `index_start`.
     * @param opt Options for extraction.
     * @return A `IndexSparseExtractor` object for sparse access to a indexed subset of each row.
     */
    virtual std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_row(const Index_* index_start, size_t index_length, const Options<Index_>& opt) const = 0;

    /**
     * @param opt Options for extraction.
     * @return A `FullSparseExtractor` object for sparse access to full columns.
     */
    virtual std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_column(const Options<Index_>& opt) const = 0;

    /**
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @param opt Options for extraction.
     * @return A `BlockSparseExtractor` object for sparse access to a block of each column.
     */
    virtual std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_end, const Options<Index_>& opt) const = 0;

    /**
     * @param index_start Pointer to an array containing sorted and unique row indices.
     * @param index_length Length of the array pointed to by `index_start`.
     * @param opt Options for extraction.
     * @return A `IndexSparseExtractor` object for sparse access to a indexed subset of each column.
     */
    virtual std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_column(const Index_* index_start, size_t index_length, const Options<Index_>& opt) const = 0;

    /***************************************
     **** Dense access default overload ****
     ***************************************/
public:
    /**
     * @return A `FullDenseExtractor` object for dense access to full rows, using default options.
     */
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row() const {
        return dense_row(Options<Index_>());
    }

    /**
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @return A `BlockDenseExtractor` object for dense access to a block of each row, using default options.
     */
    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_end) const {
        return dense_row(block_start, block_end, Options<Index_>());
    }

    /**
     * @param index_start Pointer to an array containing sorted and unique column indices.
     * @param index_length Length of the array pointed to by `index_start`.
     * @return A `IndexDenseExtractor` object for dense access to a indexed subset of each row, using default options.
     */
    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(const Index_* index_start, size_t index_length) const {
        return dense_row(index_start, index_length, Options<Index_>());
    }

    /**
     * @return A `FullDenseExtractor` object for dense access to full columns, using default options.
     */
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column() const {
        return dense_column(Options<Index_>());
    }

    /**
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @return A `BlockDenseExtractor` object for dense access to a block of each column, using default options.
     */
    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_end) const {
        return dense_column(block_start, block_end, Options<Index_>());
    }

    /**
     * @param index_start Pointer to an array containing sorted and unique row indices.
     * @param index_length Length of the array pointed to by `index_start`.
     * @return A `IndexDenseExtractor` object for dense access to a indexed subset of each column, using default options.
     */
    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(const Index_* index_start, size_t index_length) const {
        return dense_column(index_start, index_length, Options<Index_>()); 
    }

    /*****************************************
     **** Sparse access default overloads ****
     *****************************************/
public:
    /**
     * @return A `FullSparseExtractor` object for sparse access to full rows, using default options.
     */
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_row() const {
        return sparse_row(Options<Index_>());
    }

    /**
     * @param block_start Index of the column at the start of the block.
     * @param block_length Number of columns in the block.
     * @return A `BlockSparseExtractor` object for sparse access to a block of each row, using default options.
     */
    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_end) const {
        return sparse_row(block_start, block_end, Options<Index_>());
    }

    /**
     * @param index_start Pointer to an array containing sorted and unique column indices.
     * @param index_length Length of the array pointed to by `index_start`.
     * @return A `IndexSparseExtractor` object for sparse access to a indexed subset of each row, using default options.
     */
    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_row(const Index_* index_start, size_t index_length) const {
        return sparse_row(index_start, index_length, Options<Index_>());
    }

    /**
     * @return A `FullSparseExtractor` object for sparse access to full columns, using default options.
     */
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_column() const {
        return sparse_column(Options<Index_>());
    }

    /**
     * @param block_start Index of the row at the start of the block.
     * @param block_length Number of rows in the block.
     * @return A `BlockSparseExtractor` object for sparse access to a block of each column, using default options.
     */
    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_end) const {
        return sparse_column(block_start, block_end, Options<Index_>());
    }

    /**
     * @param index_start Pointer to an array containing sorted and unique row indices.
     * @param index_length Length of the array pointed to by `index_start`.
     * @return A `IndexSparseExtractor` object for sparse access to a indexed subset of each column, using default options.
     */
    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_column(const Index_* index_start, size_t index_length) const {
        return sparse_column(index_start, index_length, Options<Index_>());
    }
};

/**
 * A convenient shorthand for the most common use case of double-precision matrices.
 */
using NumericMatrix = Matrix<double, int>;

}

#endif
