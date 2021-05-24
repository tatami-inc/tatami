#ifndef TYPED_MATRIX_H
#define TYPED_MATRIX_H

#include "matrix.hpp"

namespace bioc {

template <typename T>
class typed_matrix : public matrix {
public:
    ~typed_matrix() {}

    virtual const T* get_row(size_t, T*, size_t=0, size_t=-1, void* work=NULL) const = 0;

    virtual const T* get_column(size_t, T*, size_t=0, size_t=-1, void * work=NULL) const = 0;

    content_type type() const { return determine_content_type<T>(); }
};

using numeric_matrix = typed_matrix<double>;

}

#endif
