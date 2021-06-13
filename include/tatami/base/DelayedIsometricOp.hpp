#ifndef TATAMI_DELAYED_ISOMETRIC_OP_H
#define TATAMI_DELAYED_ISOMETRIC_OP_H

#include <memory>
#include "typed_matrix.hpp"

/**
 * @file DelayedIsometricOp.hpp
 *
 * Delayed isometric operations, equivalent to the class of the same name in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @brief Delayed isometric operations on a matrix.
 *
 * Implements any operation that preserves the shape of the matrix and operates on each matrix value independently.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam T Type of matrix value.
 * @tparam OP Functor class implementing the operation.
 * This should accept the row index, column index and value, and return the modified value after applying the operation. 
 * @tparam IDX Type of index value.
 */
template<typename T, typename IDX, class OP>
class DelayedIsometricOp : public typed_matrix<T, IDX> {
public:
    /**
     * @param p Pointer to the underlying matrix.
     * @param op Instance of the functor class.
     */
    DelayedIsometricOp(std::shared_ptr<const typed_matrix<T, IDX> > p, OP op) : mat(p), operation(std::move(op)) {}

    ~DelayedIsometricOp() {}

public:
    const T* row(size_t r, T* buffer, size_t start, size_t end, workspace* work=nullptr) const {
        const T* raw = mat->row(r, buffer, start, end, work);
        for (size_t i = start; i < end; ++i, ++raw) {
            buffer[i - start] = operation(r, i, *raw);
        }
        return buffer;
    }

    const T* column(size_t c, T* buffer, size_t start, size_t end, workspace* work=nullptr) const {
        const T* raw = mat->column(c, buffer, start, end, work);
        for (size_t i = start; i < end; ++i, ++raw) {
            buffer[i - start] = operation(i, c, *raw);
        }
        return buffer;
    }

    using typed_matrix<T, IDX>::column;

    using typed_matrix<T, IDX>::row;

public:
    sparse_range<T, IDX> sparse_row(size_t r, T* vbuffer, IDX* ibuffer, size_t start, size_t end, workspace* work=nullptr, bool sorted=true) const {
        if (OP::sparse) {
            auto raw = mat->sparse_row(r, vbuffer, ibuffer, start, end, work, sorted);
            for (size_t i = 0; i < raw.number; ++i) {
                vbuffer[i] = operation(r, raw.index[i], raw.value[i]);
            }
            return sparse_range<T, IDX>(raw.number, vbuffer, raw.index);
        } else {
            auto ptr = mat->row(r, vbuffer, start, end, work);
            auto original_values = vbuffer;
            auto original_indices = ibuffer;
            for (size_t i = start; i < end; ++i, ++vbuffer, ++ibuffer, ++ptr) {
                (*vbuffer) = operation(r, i, *ptr);
                (*ibuffer) = i;
            }
            return sparse_range<T, IDX>(end - start, original_values, original_indices); 
        }
    }

    sparse_range<T, IDX> sparse_column(size_t c, T* vbuffer, IDX* ibuffer, size_t start, size_t end, workspace* work=nullptr, bool sorted=true) const {
        if (OP::sparse) {
            auto raw = mat->sparse_column(c, vbuffer, ibuffer, start, end, work, sorted);
            for (size_t i = 0; i < raw.number; ++i) {
                vbuffer[i] = operation(raw.index[i], c, raw.value[i]);
            }
            return sparse_range<T, IDX>(raw.number, vbuffer, raw.index);
        } else {
            auto ptr = mat->column(c, vbuffer, start, end, work);
            auto original_values = vbuffer;
            auto original_indices = ibuffer;
            for (size_t i = start; i < end; ++i, ++vbuffer, ++ibuffer, ++ptr) {
                (*vbuffer) = operation(i, c, *ptr);
                (*ibuffer) = i;
            }
            return sparse_range<T, IDX>(end - start, original_values, original_indices); 
        }
    }

    using typed_matrix<T, IDX>::sparse_column;

    using typed_matrix<T, IDX>::sparse_row;

public:
    size_t nrow() const {
        return mat->nrow();
    }
    
    size_t ncol() const {
        return mat->ncol();
    }

    /**
     * @return A null pointer or a shared pointer to a `workspace` object, depending on the underlying (pre-operation) matrix.
     */
    std::shared_ptr<workspace> new_workspace(bool row) const {
        return mat->new_workspace(row);
    }

    /**
     * @return `true` if both the underlying (pre-operation) matrix is sparse and the operation preserves sparsity.
     * Otherwise returns `false`.
     */
    bool sparse() const {
        return mat->sparse() && OP::sparse;
    }

    /**
     * @return `true` if row-wise extraction is preferred by the underlying (pre-operation) matrix, otherwise returns `false`.
     */
    bool prefer_rows() const { return mat->prefer_rows(); }

private:
    std::shared_ptr<const typed_matrix<T, IDX> > mat;
    OP operation;
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam MAT A specialized `typed_matrix`, to be automatically deducted.
 * @tparam OP Helper class defining the operation.
 *
 * @param p Pointer to a `typed_matrix`.
 * @param op Instance of the operation helper class.
 */
template<class MAT, class OP>
std::shared_ptr<MAT> make_DelayedIsometricOp(std::shared_ptr<MAT> p, OP op) {
    return std::shared_ptr<MAT>(
        new DelayedIsometricOp<typename MAT::value, typename MAT::index, typename std::remove_reference<OP>::type>(
            p,
            std::move(op)
        )
    );
}

}

#include "arith_scalar_helpers.hpp"

#include "arith_vector_helpers.hpp"

#include "math_helpers.hpp"

#endif
