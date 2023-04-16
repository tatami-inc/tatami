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
     * @param iopt Iteration options.
     * @param eopt Extraction options.
     * 
     * @return A `DimensionAccess` object for random access to dense rows.
     */
    virtual std::unique_ptr<DenseExtractor<Value_, Index_> > dense_row(IterationOptions<Index_> iopt, ExtractionOptions<Index_> eopt) const = 0;

    /**
     * @param iopt Iteration options.
     * @param eopt Extraction options.
     * 
     * @return A `DimensionAccess` object for random access to dense columns.
     */
    virtual std::unique_ptr<DenseExtractor<Value_, Index_> > dense_column(IterationOptions<Index_> iopt, ExtractionOptions<Index_> eopt) const = 0;

    /***************************************
     **** Sparse access virtual methods ****
     ***************************************/
public:
    /**
     * @param iopt Iteration options.
     * @param eopt Extraction options.
     * 
     * @return A `DimensionAccess` object for random access to sparse rows.
     */
    virtual std::unique_ptr<SparseExtractor<Value_, Index_> > sparse_row(IterationOptions<Index_> iopt, ExtractionOptions<Index_> eopt) const = 0;

    /**
     * @param iopt Iteration options.
     * @param eopt Extraction options.
     * 
     * @return A `DimensionAccess` object for random access to sparse columns.
     */
    virtual std::unique_ptr<SparseExtractor<Value_, Index_> > sparse_column(IterationOptions<Index_> iopt, ExtractionOptions<Index_> eopt) const = 0;

    /***************************************
     **** Dense access default overload ****
     ***************************************/
public:
    /**
     * @return A `DimensionAccess` object for random access to dense rows, without any dimension limits or non-default options.
     */
    std::unique_ptr<DenseExtractor<Value_, Index_> > dense_row() const {
        return dense_row(IterationOptions<Index_>(), ExtractionOptions<Index_>());
    }

    /**
     * @return A `DimensionAccess` object for random access to dense columns, without any dimension limits or non-default options.
     */
    std::unique_ptr<DenseExtractor<Value_, Index_> > dense_column() const {
        return dense_column(IterationOptions<Index_>(), ExtractionOptions<Index_>());
    }

    /****************************************
     **** Sparse access options overload ****
     ****************************************/
public:
    /**
     * @return A `DimensionAccess` object for random access to sparse rows, without any dimension limits or non-default options.
     */
    std::unique_ptr<SparseExtractor<Value_, Index_> > sparse_row() const {
        return sparse_row(IterationOptions<Index_>(), ExtractionOptions<Index_>());
    }

    /**
     * @return A `DimensionAccess` object for random access to sparse columns, without any dimension limits or non-default options.
     */
    std::unique_ptr<SparseExtractor<Value_, Index_> > sparse_column() const {
        return sparse_column(IterationOptions<Index_>(), ExtractionOptions<Index_>());
    }
};

/**
 * A convenient shorthand for the most common use case of double-precision matrices.
 */
using NumericMatrix = Matrix<double, int>;

}

#endif
