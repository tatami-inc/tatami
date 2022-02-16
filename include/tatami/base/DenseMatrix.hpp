#ifndef TATAMI_DENSE_MATRIX_H
#define TATAMI_DENSE_MATRIX_H

#include "Matrix.hpp"
#include "has_data.hpp"

#include <vector>
#include <algorithm>

/**
 * @file DenseMatrix.hpp
 *
 * Dense matrix representation, with `typedef`s for the usual row- and column-major formats.
 */

namespace tatami {

/**
 * @brief Dense matrix representation.
 *
 * @tparam ROW Whether this is a row-major representation.
 * If `false`, a column-major representation is assumed instead.
 * @tparam T Type of the matrix values.
 * @tparam IDX Type of the row/column indices.
 * @tparam V Vector class used to store the matrix values internally.
 * This does not necessarily have to contain `T`, as long as the type is convertible to `T`.
 * Methods should be available for `size()`, `begin()`, `end()` and `[]`.
 * If a method is available for `data()` that returns a `const T*`, it will also be used.
 */
template<bool ROW, typename T, typename IDX = int, class V = std::vector<T> >
class DenseMatrix : public Matrix<T, IDX> {
public: 
    /**
     * @param nr Number of rows.
     * @param nc Number of columns.
     * @param source Vector of values, or length equal to the product of `nr` and `nc`.
     */
    DenseMatrix(size_t nr, size_t nc, const V& source) : nrows(nr), ncols(nc), values(source) {
        if (nrows * ncols != values.size()) {
            throw std::runtime_error("length of 'values' should be equal to product of 'nrows' and 'ncols'");
        }
        return;
    }

    /**
     * @param nr Number of rows.
     * @param nc Number of columns.
     * @param source Vector of values, or length equal to the product of `nr` and `nc`.
     */
    DenseMatrix(size_t nr, size_t nc, V&& source) : nrows(nr), ncols(nc), values(source) {
        if (nrows * ncols != values.size()) {
            throw std::runtime_error("length of 'values' should be equal to product of 'nrows' and 'ncols'");
        }
        return;
    }

public:
    size_t nrow() const { return nrows; }

    size_t ncol() const { return ncols; }

    /**
     * @return `true` if `ROW = true` (for row-major matrices), otherwise returns `false` (for column-major matrices).
     */
    bool prefer_rows() const { return ROW; }

public:
    const T* row(size_t r, T* buffer, size_t start, size_t end, Workspace* work=nullptr) const {
        if constexpr(ROW) {
            return primary(r, buffer, start, end, work, ncols);
        } else {
            secondary(r, buffer, start, end, work, nrows);
            return buffer;
        }
    }

    const T* column(size_t c, T* buffer, size_t start, size_t end, Workspace* work=nullptr) const {
        if constexpr(ROW) {
            secondary(c, buffer, start, end, work, ncols);
            return buffer;
        } else {
            return primary(c, buffer, start, end, work, nrows);
        }
    }

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::column;

private: 
    size_t nrows, ncols;
    V values;

    const T* primary(size_t c, T* buffer, size_t start, size_t end, Workspace* work, size_t dim_secondary) const {
        size_t shift = c * dim_secondary;
        if constexpr(has_data<T, V>::value) {
            return values.data() + shift + start;
        } else {
            end = std::min(end, dim_secondary);
            std::copy(values.begin() + shift + start, values.begin() + shift + end, buffer);
            return buffer;
        }
    }

    void secondary(size_t r, T* buffer, size_t start, size_t end, Workspace* work, size_t dim_secondary) const {
        auto it = values.begin() + r + start * dim_secondary;
        for (size_t i = start; i < end; ++i, ++buffer, it+=dim_secondary) {
            *buffer = *it; 
        }
        return;
    }
};

/**
 * Column-major matrix.
 * See `tatami::DenseMatrix` for details on the template parameters.
 */
template<typename T, typename IDX = int, class V = std::vector<T> >
using DenseColumnMatrix = DenseMatrix<false, T, IDX, V>;

/**
 * Row-major matrix.
 * See `tatami::DenseMatrix` for details on the template parameters.
 */
template<typename T, typename IDX = int, class V = std::vector<T> >
using DenseRowMatrix = DenseMatrix<true, T, IDX, V>;

}

#endif
