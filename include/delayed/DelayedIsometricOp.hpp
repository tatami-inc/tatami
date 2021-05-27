#ifndef DELAYED_ISOMETRIC_OP_H
#define DELAYED_ISOMETRIC_OP_H

#include "../utils/types.hpp"

namespace bioc {

template<typename T, class OP, typename X = T, typename IDX = int>
class DelayedIsometricOp : public typed_matrix<X, IDX> {
public:
    DelayedIsometricOp(std::shared_ptr<const typed_matrix<X, IDX> > p, const OP& op) : mat(p), operation(op) {}

    DelayedIsometricOp(std::shared_ptr<const typed_matrix<X, IDX> > p, OP&& op) : mat(p), operation(op) {}

    ~DelayedIsometricOp() {}

    const T* get_row(size_t r, T* buffer, size_t start=0, size_t end=-1, workspace* work=NULL) const {
        const X* raw = mat.get_row(r, buffer, start, end, work);
        for (size_t i = start; i < end; ++i, ++raw) {
            buffer[i] = operation.modify(r, i, *raw);
        }
        return buffer;
    }

    const T* get_column(size_t c, T* buffer, size_t start=0, size_t end=-1, workspace* work=NULL) const {
        const X* raw = mat.get_column(c, buffer, start, end, work);
        for (size_t i = start; i < end; ++i, ++raw) {
            buffer[i] = operation.modify(i, c, *raw);
        }
        return buffer;
    }

    sparse_range<T, IDX> get_sparse_row(size_t r, T* out_values, IDX* out_indices, size_t start=0, size_t end=-1, workspace* work=NULL) const {
        auto raw = mat.get_sparse_row(r, out_values, out_indices, start, end, work);
        for (size_t i = 0; i < raw.number; ++i) {
            out_values[i] = operation(r, i, raw.value[i]);
        }
        return sparse_range<T, IDX>(raw.number, out_values, raw.index);
    }

    sparse_range<T, IDX> get_sparse_column(size_t c, T* out_values, IDX* out_indices, size_t start=0, size_t end=-1, workspace* work=NULL) const {
        auto raw = mat.get_sparse_column(c, out_values, out_indices, start, end, work);
        for (size_t i = 0; i < raw.number; ++i) {
            out_values[i] = operation(i, c, raw.value[i]);
        }
        return sparse_range<T, IDX>(raw.number, out_values, raw.index);
    }

    size_t nrow() const {
        return mat->nrow();
    }
    
    size_t ncol() const {
        return mat->ncol();
    }

    workspace* create_workspace() const {
        return mat->create_workspace();
    }

    bool is_sparse() const {
        return mat->is_sparse() && OP.sparse();
    }
protected:
    std::shared_ptr<typed_matrix<X, IDX> > mat;
    OP operation;
};

}

#endif
