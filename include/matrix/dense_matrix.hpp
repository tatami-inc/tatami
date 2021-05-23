#ifndef DENSE_MATRIX_H
#define DENSE_MATRIX_H

#include "matrix.hpp"

#include <vector>
#include <algorithm>

namespace bioc {

template<typename T, class V = std::vector<T> >
class dense_matrix : public typed_matrix<T> {
public: 
    dense_matrix(size_t nr, size_t nc) : nrows(nr), ncols(nc), values(nr * nc) {}

    dense_matrix(size_t nr, size_t nc, V&& source) : nrows(nr), ncols(nc), values(source) {}

    ~dense_matrix() {}

    size_t nrow() const { return nrows; }

    size_t ncol() const { return ncols; }

    const T* get_row(size_t r, T* buffer, size_t start=0, size_t end=-1) const {
        auto it = values.begin() + r + start * nrows;
        end = std::min(end, ncols);
        for (size_t i = start; i < end; ++i, ++buffer, it += nrows) {
            *buffer = *it; 
        }
        return buffer;
    }

    const T* get_column(size_t c, T* buffer, size_t start=0, size_t end=-1) const {
        auto it = values.begin() + c * nrows;
        if constexpr(std::is_same<V, std::vector<T> >::value) {
            return it + start;
        } else {
            end = std::min(end, nrows);
            std::copy(it + start, it + end, buffer);
            return buffer;
        }
    }

    void set_row(size_t r, const T* buffer, size_t start=0, size_t end=-1) {
        auto it = values.begin() + r + start * nrows;
        end = std::min(end, ncols);
        for (size_t i = start; i < end; ++i, ++buffer, it += nrows) {
            *it = *buffer; 
        }
        return;
    }

    void set_column(size_t c, const T* buffer, size_t start=0, size_t end=-1) {
        auto it = values.begin() + c * nrows;
        end = std::min(end, nrows);
        std::copy(buffer, buffer + end - start, it + start);
        return;
    }
private: 
    size_t nrows, ncols;
    V values;
};

}

#endif
