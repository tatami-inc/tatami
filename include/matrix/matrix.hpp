#ifndef MATRIX_H
#define MATRIX_H

#include "../utils/types.hpp"
#include "../utils/workspace.hpp"

namespace bioc {

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
