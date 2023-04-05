#ifndef TATAMI_DENSE_MATRIX_H
#define TATAMI_DENSE_MATRIX_H

#include "Matrix.hpp"
#include "has_data.hpp"

#include <vector>
#include <algorithm>
#include <stdexcept>

/**
 * @file DenseMatrix.hpp
 *
 * @brief Dense matrix representation.
 *
 * `typedef`s are provided for the usual row- and column-major formats.
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

private: 
    size_t nrows, ncols;
    V values;

public:
    size_t nrow() const { return nrows; }

    size_t ncol() const { return ncols; }

    /**
     * @return `true` if `ROW = true` (for row-major matrices), otherwise returns `false` (for column-major matrices).
     */
    bool prefer_rows() const { return ROW; }

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::column;

public:
    std::shared_ptr<RowWorkspace> new_row_workspace(bool = false) const { return nullptr; }

    std::shared_ptr<ColumnWorkspace> new_column_workspace(bool = false) const { return nullptr; }

    const T* row(size_t r, T* buffer, RowWorkspace* work) const {
        if constexpr(ROW) {
            return primary(r, buffer, 0, ncols, ncols);
        } else {
            secondary(r, buffer, 0, ncols, nrows);
            return buffer;
        }
    }

    const T* column(size_t c, T* buffer, ColumnWorkspace* work=nullptr) const {
        if constexpr(ROW) {
            secondary(c, buffer, 0, nrows, ncols);
            return buffer;
        } else {
            return primary(c, buffer, 0, nrows, nrows);
        }
    }

public:
    /**
     * @cond
     */
    template<bool WORKROW>
    struct DenseBlockWorkspace : public BlockWorkspace<WORKROW> {
        DenseBlockWorkspace(size_t s, size_t l) : details(s, l) {};
        std::pair<size_t, size_t> details;
        const std::pair<size_t, size_t>& block() const { return details; }
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowBlockWorkspace> new_row_workspace(size_t start, size_t len, bool = false) const {
        return std::shared_ptr<RowBlockWorkspace>(new DenseBlockWorkspace<true>(start, len));
    }

    std::shared_ptr<ColumnBlockWorkspace> new_column_workspace(size_t start, size_t len, bool = false) const {
        return std::shared_ptr<ColumnBlockWorkspace>(new DenseBlockWorkspace<false>(start, len));
    }

    const T* row(size_t r, T* buffer, RowBlockWorkspace* work) const {
        const auto& deets = work->block();
        size_t start = deets.first, end = start + deets.second;
        if constexpr(ROW) {
            return primary(r, buffer, start, end, ncols);
        } else {
            secondary(r, buffer, start, end, nrows);
            return buffer;
        }
    }

    const T* column(size_t c, T* buffer, ColumnBlockWorkspace* work) const {
        const auto& deets = work->block();
        size_t start = deets.first, end = start + deets.second;
        if constexpr(ROW) {
            secondary(c, buffer, start, end, ncols);
            return buffer;
        } else {
            return primary(c, buffer, start, end, nrows);
        }
    }

private:
    const T* primary(size_t c, T* buffer, size_t start, size_t end, size_t dim_secondary) const {
        size_t shift = c * dim_secondary;
        if constexpr(has_data<T, V>::value) {
            return values.data() + shift + start;
        } else {
            std::copy(values.begin() + shift + start, values.begin() + shift + end, buffer);
            return buffer;
        }
    }

    void secondary(size_t r, T* buffer, size_t start, size_t end, size_t dim_secondary) const {
        auto it = values.begin() + r + start * dim_secondary;
        for (size_t i = start; i < end; ++i, ++buffer, it += dim_secondary) {
            *buffer = *it; 
        }
        return;
    }

public:
    /**
     * @cond
     */
    template<bool WORKROW>
    struct DenseIndexWorkspace : public IndexWorkspace<IDX, WORKROW> {
        DenseIndexWorkspace(std::vector<IDX> i) : indices_(std::move(i)) {}
        std::vector<IDX> indices_;
        const std::vector<IDX>& indices() const { return indices_; }
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(std::vector<IDX> i, bool = false) const {
        return std::shared_ptr<RowIndexWorkspace<IDX> >(new DenseIndexWorkspace<true>(std::move(i)));
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX> > new_column_workspace(std::vector<IDX> i, bool = false) const {
        return std::shared_ptr<ColumnIndexWorkspace<IDX> >(new DenseIndexWorkspace<false>(std::move(i)));
    }

    const T* row(size_t r, T* buffer, RowIndexWorkspace<IDX>* work) const {
        if constexpr(ROW) {
            return primary_indexed(r, buffer, work->indices(), ncols);
        } else {
            secondary_indexed(r, buffer, work->indices(), nrows);
            return buffer;
        }
    }

    const T* column(size_t c, T* buffer, ColumnIndexWorkspace<IDX>* work) const {
        if constexpr(ROW) {
            secondary_indexed(c, buffer, work->indices(), ncols);
            return buffer;
        } else {
            return primary_indexed(c, buffer, work->indices(), nrows);
        }
    }

private:
    const T* primary_indexed(size_t c, T* buffer, const std::vector<IDX>& indices, size_t dim_secondary) const {
        auto offset = c * dim_secondary;
        for (size_t i = 0, end = indices.size(); i < end; ++i) {
            buffer[i] = values[indices[i] + offset];
        }
        return buffer;
    }

    void secondary_indexed(size_t r, T* buffer, const std::vector<IDX>& indices, size_t dim_secondary) const {
        for (size_t i = 0, end = indices.size(); i < end; ++i, ++buffer) {
            *buffer = values[indices[i] * dim_secondary + r]; 
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
