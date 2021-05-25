#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "typed_matrix.hpp"
#include "../utils/workspace.hpp"

namespace bioc {

template<typename T, typename IDX=int>
class sparse_matrix : public typed_matrix<T> {
public:
    ~sparse_matrix() {}

    struct sparse_range {
        sparse_range(size_t n, const T* v=NULL, const IDX* i=NULL) : number(n), value(v), index(i) {}
        sparse_range() {}
        size_t number;
        const T* value;
        const IDX* index;
    };
    
    virtual sparse_range get_sparse_row(size_t, T*, IDX*, size_t=0, size_t=-1, workspace* work=NULL) const = 0;

    virtual sparse_range get_sparse_column(size_t, T*, IDX*, size_t=0, size_t=-1, workspace* work=NULL) const = 0;

    sparse_range get_sparse_row(size_t i, T* out_values, IDX* out_indices, workspace* work) const {
        return get_sparse_row(i, out_values, out_indices, 0, this->ncol(), work);
    }

    sparse_range get_sparse_column(size_t i, T* out_values, IDX* out_indices, workspace* work) const {
        return get_sparse_column(i, out_values, out_indices, 0, this->nrow(), work);
    }

    bool is_sparse() const { return true; }
};

}

#endif
