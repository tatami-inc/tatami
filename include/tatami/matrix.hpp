#ifndef TATAMI_MATRIX_H
#define TATAMI_MATRIX_H

#include "types.hpp"
#include "workspace.hpp"

namespace tatami {

class matrix {
public:
    virtual ~matrix() {}

    virtual size_t nrow() const = 0;

    virtual size_t ncol() const = 0;

    virtual workspace* create_workspace(bool) const { return NULL; }

    virtual content_type type() const { return _unknown; }

    virtual bool is_sparse() const { return false; }

    virtual int preferred_dimension() const { return 1; }
};

}

#endif
