#ifndef TYPED_MATRIX_H
#define TYPED_MATRIX_H

#include "matrix.hpp"
#include "../utils/sparse_range.hpp"
#include "../utils/workspace.hpp"

namespace bioc {

template <typename T, typename IDX = int>
class typed_matrix : public matrix {
public:
    ~typed_matrix() {}

    virtual const T* get_row(size_t, T*, size_t, size_t, workspace* work=NULL) const = 0;

    virtual const T* get_column(size_t, T*, size_t, size_t, workspace* work=NULL) const = 0;

    const T* get_row(size_t i, T* out, workspace* work=NULL) const {
        return get_row(i, out, 0, this->ncol(), work);
    }

    const T* get_column(size_t i, T* out, workspace* work=NULL) const {
        return get_column(i, out, 0, this->nrow(), work);
    }

    virtual sparse_range<T, IDX> get_sparse_row(size_t i, T* outv, IDX* outi, size_t start, size_t end, workspace* work=NULL) const {
        const T* val = get_row(i, outv, start, end, work);
        for (size_t i = start; i < end; ++i) {
            outi[i - start] = i;
        }
        return sparse_range(end - start, outv, outi); 
    }

    virtual sparse_range<T, IDX> get_sparse_column(size_t i, T* outv, IDX* outi, size_t start, size_t end, workspace* work=NULL) const {
        const T* val = get_column(i, outv, start, end, work);
        for (size_t i = start; i < end; ++i) {
            outi[i - start] = i;
        }
        return sparse_range(end - start, outv, outi); 
    }

    sparse_range<T, IDX> get_sparse_row(size_t i, T* out_values, IDX* out_indices, workspace* work=NULL) const {
        return get_sparse_row(i, out_values, out_indices, 0, this->ncol(), work);
    }

    sparse_range<T, IDX> get_sparse_column(size_t i, T* out_values, IDX* out_indices, workspace* work=NULL) const {
        return get_sparse_column(i, out_values, out_indices, 0, this->nrow(), work);
    }

    content_type type() const { return determine_content_type<T>(); }
};

using numeric_matrix = typed_matrix<double>;

}

#endif
