#ifndef DENSE_MATRIX_H
#define DENSE_MATRIX_H

#include "typed_matrix.hpp"
#include "../utils/is_contiguous.hpp"

#include <vector>
#include <algorithm>

namespace bioc {

template<bool ROW, typename T, class V = std::vector<T>, typename IDX = int>
class DenseMatrix : public typed_matrix<T, IDX> {
public: 
    DenseMatrix(size_t nr, size_t nc) : nrows(nr), ncols(nc), values(nr * nc) {}

    DenseMatrix(size_t nr, size_t nc, const V& source) : nrows(nr), ncols(nc), values(source) {
        if (nrows * ncols != values.size()) {
            throw std::runtime_error("length of 'values' should be equal to product of 'nrows' and 'ncols'");
        }
        return;
    }

    DenseMatrix(size_t nr, size_t nc, V&& source) : nrows(nr), ncols(nc), values(source) {
        if (nrows * ncols != values.size()) {
            throw std::runtime_error("length of 'values' should be equal to product of 'nrows' and 'ncols'");
        }
        return;
    }

    ~DenseMatrix() {}

public:
    size_t nrow() const { return nrows; }

    size_t ncol() const { return ncols; }

public:
    const T* get_row(size_t r, T* buffer, size_t start, size_t end, workspace * wrk=NULL) const {
        if constexpr(ROW) {
            return get_primary(r, buffer, start, end, wrk, ncols);
        } else {
            get_secondary(r, buffer, start, end, wrk, nrows);
            return buffer;
        }
    }

    const T* get_column(size_t c, T* buffer, size_t start, size_t end, workspace* wrk=NULL) const {
        if constexpr(ROW) {
            get_secondary(c, buffer, start, end, wrk, ncols);
            return buffer;
        } else {
            return get_primary(c, buffer, start, end, wrk, nrows);
        }
    }

    using typed_matrix<T, IDX>::get_row;

    using typed_matrix<T, IDX>::get_column;

private: 
    size_t nrows, ncols;
    V values;

    const T* get_primary(size_t c, T* buffer, size_t start, size_t end, workspace* wrk, size_t dim_secondary) const {
        size_t shift = c * dim_secondary;
        if constexpr(has_data<T, V>::value) {
            return values.data() + shift + start;
        } else {
            end = std::min(end, dim_secondary);
            std::copy(values.begin() + start, values.begin() + end, buffer);
            return buffer;
        }
    }

    void get_secondary(size_t r, T* buffer, size_t start, size_t end, workspace * wrk, size_t dim_secondary) const {
        auto it = values.begin() + r + start * dim_secondary;
        for (size_t i = start; i < end; ++i, ++buffer, it+=dim_secondary) {
            *buffer = *it; 
        }
        return;
    }
};

template<typename T, class V = std::vector<T> >
using DenseColumnMatrix = DenseMatrix<false, T, V>;

template<typename T, class V = std::vector<T> >
using DenseRowMatrix = DenseMatrix<true, T, V>;

}

#endif
