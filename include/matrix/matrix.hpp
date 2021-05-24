#ifndef MATRIX_H
#define MATRIX_H

#include "../utils/types.hpp"

namespace bioc {

class matrix {
public:
    virtual ~matrix() {}

    virtual size_t nrow() const = 0;

    virtual size_t ncol() const = 0;

    virtual void * create_workspace() const = 0;

    virtual content_type type() const { return _unknown; }

    virtual bool is_sparse() const { return false; }
};

}

#endif
