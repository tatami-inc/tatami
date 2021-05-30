#ifndef DELAYED_ISOMETRIC_OP_H
#define DELAYED_ISOMETRIC_OP_H

#include "../utils/types.hpp"

namespace bioc {

template<typename T, class OP, typename IDX = int>
class DelayedIsometricOp : public typed_matrix<T, IDX> {
public:
    DelayedIsometricOp(std::shared_ptr<const typed_matrix<T, IDX> > p, const OP& op) : mat(p), operation(op) {}

    DelayedIsometricOp(std::shared_ptr<const typed_matrix<T, IDX> > p, OP&& op) : mat(p), operation(op) {}

    ~DelayedIsometricOp() {}

public:
    const T* get_row(size_t r, T* buffer, size_t start, size_t end, workspace* work=NULL) const {
        const T* raw = mat->get_row(r, buffer, start, end, work);
        for (size_t i = start; i < end; ++i, ++raw) {
            buffer[i - start] = operation(r, i, *raw);
        }
        return buffer;
    }

    const T* get_column(size_t c, T* buffer, size_t start, size_t end, workspace* work=NULL) const {
        const T* raw = mat->get_column(c, buffer, start, end, work);
        for (size_t i = start; i < end; ++i, ++raw) {
            buffer[i - start] = operation(i, c, *raw);
        }
        return buffer;
    }

    using typed_matrix<T, IDX>::get_column;

    using typed_matrix<T, IDX>::get_row;

public:
    sparse_range<T, IDX> get_sparse_row(size_t r, T* out_values, IDX* out_indices, size_t start, size_t end, workspace* work=NULL) const {
        if (OP::sparse) {
            auto raw = mat->get_sparse_row(r, out_values, out_indices, start, end, work);
            for (size_t i = 0; i < raw.number; ++i) {
                out_values[i] = operation(r, raw.index[i], raw.value[i]);
            }
            return sparse_range<T, IDX>(raw.number, out_values, raw.index);
        } else {
            auto ptr = mat->get_row(r, out_values, start, end, work);
            auto original_values = out_values;
            auto original_indices = out_indices;
            for (size_t i = start; i < end; ++i, ++out_values, ++out_indices, ++ptr) {
                (*out_values) = operation(r, i, *ptr);
                (*out_indices) = i;
            }
            return sparse_range<T, IDX>(end - start, original_values, original_indices); 
        }
    }

    sparse_range<T, IDX> get_sparse_column(size_t c, T* out_values, IDX* out_indices, size_t start, size_t end, workspace* work=NULL) const {
        if (OP::sparse) {
            auto raw = mat->get_sparse_column(c, out_values, out_indices, start, end, work);
            for (size_t i = 0; i < raw.number; ++i) {
                out_values[i] = operation(raw.index[i], c, raw.value[i]);
            }
            return sparse_range<T, IDX>(raw.number, out_values, raw.index);
        } else {
            auto ptr = mat->get_column(c, out_values, start, end, work);
            auto original_values = out_values;
            auto original_indices = out_indices;
            for (size_t i = start; i < end; ++i, ++out_values, ++out_indices, ++ptr) {
                (*out_values) = operation(i, c, *ptr);
                (*out_indices) = i;
            }
            return sparse_range<T, IDX>(end - start, original_values, original_indices); 
        }
    }

    using typed_matrix<T, IDX>::get_sparse_column;

    using typed_matrix<T, IDX>::get_sparse_row;

public:
    size_t nrow() const {
        return mat->nrow();
    }
    
    size_t ncol() const {
        return mat->ncol();
    }

    workspace* create_workspace(bool row) const {
        return mat->create_workspace(row);
    }

    bool is_sparse() const {
        return mat->is_sparse() && OP::sparse;
    }
protected:
    std::shared_ptr<const typed_matrix<T, IDX> > mat;
    OP operation;
};

}

#include "arith_scalar_helpers.hpp"
#include "arith_vector_helpers.hpp"
#include "math_helpers.hpp"

#endif
