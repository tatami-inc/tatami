#ifndef TYPED_MATRIX_H
#define TYPED_MATRIX_H

#include "matrix.hpp"
#include "../utils/workspace.hpp"

namespace bioc {

template <typename T>
class typed_matrix : public matrix {
public:
    ~typed_matrix() {}

    virtual const T* get_row(size_t, T*, size_t=0, size_t=-1, workspace* work=NULL) const = 0;

    virtual const T* get_column(size_t, T*, size_t=0, size_t=-1, workspace* work=NULL) const = 0;

    const T* get_row(size_t i, T* out, workspace* work) const {
        return get_row(i, out, 0, this->ncol(), work);
    }

    const T* get_column(size_t i, T* out, workspace* work) const {
        return get_column(i, out, 0, this->nrow(), work);
    }

    content_type type() const { return determine_content_type<T>(); }
};

using numeric_matrix = typed_matrix<double>;

}

#endif
