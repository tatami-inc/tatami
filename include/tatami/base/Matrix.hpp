#ifndef TATAMI_MATRIX_H
#define TATAMI_MATRIX_H

#include "ExtractFormat.hpp"
#include "Options.hpp"
#include <algorithm>
#include <numeric>

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
     * @param row_limits Dimension limits on the rows.
     * @param column_limits Dimension limits on the columns.
     * @param options Further options for extraction.
     * 
     * @return A `DimensionAccess` object for random access to dense rows.
     */
    virtual std::unique_ptr<DenseFormat<Value_, Index_> > dense_row(DimensionLimit<Index_> row_limits, DimensionLimit<Index_> column_limits, ExtractOptions options) const = 0;

    /**
     * @param row_limits Dimension limits on the rows.
     * @param column_limits Dimension limits on the columns.
     * @param options Further options for extraction.
     * 
     * @return A `DimensionAccess` object for random access to dense columns.
     */
    virtual std::unique_ptr<DenseFormat<Value_, Index_> > dense_column(DimensionLimit<Index_> row_limits, DimensionLimit<Index_> column_limits, ExtractOptions options) const = 0;

    /***************************************
     **** Sparse access virtual methods ****
     ***************************************/
public:
    /**
     * @param row_limits Dimension limits on the rows.
     * @param column_limits Dimension limits on the columns.
     * @param options Further options for extraction.
     * 
     * @return A `DimensionAccess` object for random access to sparse rows.
     */
    virtual std::unique_ptr<SparseFormat<Value_, Index_> > sparse_row(DimensionLimit<Index_> row_limits, DimensionLimit<Index_> column_limits, ExtractOptions options) const = 0;

    /**
     * @param row_limits Dimension limits on the rows.
     * @param column_limits Dimension limits on the columns.
     * @param options Further options for extraction.
     * 
     * @return A `DimensionAccess` object for random access to sparse columns.
     */
    virtual std::unique_ptr<SparseFormat<Value_, Index_> > sparse_column(DimensionLimit<Index_> row_limits, DimensionLimit<Index_> column_limits, ExtractOptions options) const = 0;

    /***************************************
     **** Dense access options overload ****
     ***************************************/
public:
    /**
     * @param options Options for extraction.
     * @return A `DimensionAccess` object for random access to dense rows, without any dimension limits.
     */
    std::unique_ptr<DenseFormat<Value_, Index_> > dense_row(ExtractOptions options) const {
        return dense_row(DimensionLimit<Index_>(), DimensionLimit<Index_>(), std::move(options));
    }

    /**
     * @param options Options for extraction.
     * @return A `DimensionAccess` object for random access to dense columns, without any dimension limits.
     */
    std::unique_ptr<DenseFormat<Value_, Index_> > dense_column(ExtractOptions options) const {
        return dense_column(DimensionLimit<Index_>(), DimensionLimit<Index_>(), std::move(options));
    }

    /****************************************
     **** Sparse access options overload ****
     ****************************************/
public:
    /**
     * @param options Options for extraction.
     * @return A `DimensionAccess` object for random access to sparse rows, without any dimension limits.
     */
    std::unique_ptr<SparseFormat<Value_, Index_> > sparse_row(ExtractOptions options) const {
        return sparse_row(DimensionLimit<Index_>(), DimensionLimit<Index_>(), std::move(options));
    }

    /**
     * @param options Options for extraction.
     * @return A `DimensionAccess` object for random access to sparse columns, without any dimension limits.
     */
    std::unique_ptr<SparseFormat<Value_, Index_> > sparse_column(ExtractOptions options) const {
        return sparse_column(DimensionLimit<Index_>(), DimensionLimit<Index_>(), std::move(options));
    }

    /***************************************
     **** Dense access default overload ****
     ***************************************/
public:
    /**
     * @return A `DimensionAccess` object for random access to dense rows, without any dimension limits or non-default options.
     */
    std::unique_ptr<DenseFormat<Value_, Index_> > dense_row() const {
        return dense_row(ExtractOptions());
    }

    /**
     * @return A `DimensionAccess` object for random access to dense columns, without any dimension limits or non-default options.
     */
    std::unique_ptr<DenseFormat<Value_, Index_> > dense_column() const {
        return dense_column(ExtractOptions());
    }

    /****************************************
     **** Sparse access options overload ****
     ****************************************/
public:
    /**
     * @return A `DimensionAccess` object for random access to sparse rows, without any dimension limits or non-default options.
     */
    std::unique_ptr<SparseFormat<Value_, Index_> > sparse_row() const {
        return sparse_row(ExtractOptions());
    }

    /**
     * @return A `DimensionAccess` object for random access to sparse columns, without any dimension limits or non-default options.
     */
    std::unique_ptr<SparseFormat<Value_, Index_> > sparse_column() const {
        return sparse_column(ExtractOptions());
    }
};

/**
 * A convenient shorthand for the most common use case of double-precision matrices.
 */
using NumericMatrix = Matrix<double, int>;

}

#endif
