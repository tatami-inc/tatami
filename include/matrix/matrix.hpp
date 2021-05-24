#ifndef MATRIX_H
#define MATRIX_H

namespace bioc {

class matrix {
public:
    virtual ~matrix() {}

    virtual size_t nrow() const = 0;

    virtual size_t ncol() const = 0;

    virtual void * create_workspace() const = 0;
};

}

#endif
