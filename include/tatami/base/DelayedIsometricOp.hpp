#ifndef TATAMI_DELAYED_ISOMETRIC_OP_H
#define TATAMI_DELAYED_ISOMETRIC_OP_H

#include <memory>
#include "Matrix.hpp"
#include "Workspace.hpp"

/**
 * @file DelayedIsometricOp.hpp
 *
 * @brief Delayed isometric operations.
 *
 * This is equivalent to the class of the same name in the **DelayedArray** package.
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
class DelayedIsometricOp : public Matrix<T, IDX> {
public:
    /**
     * @param p Pointer to the underlying matrix.
     * @param op Instance of the functor class.
     */
    DelayedIsometricOp(std::shared_ptr<const Matrix<T, IDX> > p, OP op) : mat(p), operation(std::move(op)) {}

private:
    std::shared_ptr<const Matrix<T, IDX> > mat;
    OP operation;
    static_assert(std::is_same<T, decltype(operation(0, 0, 0))>::value);

public:
    size_t nrow() const {
        return mat->nrow();
    }
    
    size_t ncol() const {
        return mat->ncol();
    }

    /**
     * @return `true` if both the underlying (pre-operation) matrix is sparse and the operation preserves sparsity.
     * Otherwise returns `false`.
     */
    bool sparse() const {
        if constexpr(OP::sparse) {
            return mat->sparse();
        } else {
            return false;
        }
    }

    /**
     * @return `true` if row-wise extraction is preferred by the underlying (pre-operation) matrix, otherwise returns `false`.
     */
    bool prefer_rows() const { 
        return mat->prefer_rows();
    }

    using Matrix<T, IDX>::column;

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::sparse_column;

    using Matrix<T, IDX>::sparse_row;

public:
    std::shared_ptr<RowWorkspace> new_row_workspace() const {
        return mat->new_row_workspace();
    }

    std::shared_ptr<ColumnWorkspace> new_column_workspace() const {
        return mat->new_column_workspace();
    }

    const T* row(size_t r, T* buffer, RowWorkspace* work) const {
        return operate_on_row(r, 0, mat->ncol(), buffer, work);
    }

    const T* column(size_t c, T* buffer, ColumnWorkspace* work) const {
        return operate_on_column(c, 0, mat->nrow(), buffer, work);
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* vbuffer, IDX* ibuffer, RowWorkspace* work, bool sorted=true) const {
        return operate_on_row(r, vbuffer, ibuffer, work, 0, mat->ncol(), sorted);
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* vbuffer, IDX* ibuffer, ColumnWorkspace* work, bool sorted=true) const {
        return operate_on_column(c, vbuffer, ibuffer, work, 0, mat->nrow(), sorted);
    }

public:
    std::shared_ptr<RowBlockWorkspace> new_row_workspace(size_t start, size_t length) const {
        return mat->new_row_workspace(start, length);
    }

    std::shared_ptr<ColumnBlockWorkspace> new_column_workspace(size_t start, size_t length) const {
        return mat->new_column_workspace(start, length);
    }

    const T* row(size_t r, T* buffer, RowBlockWorkspace* work) const {
        const auto& deets = work->block();
        size_t start = deets.first, end = start + deets.second;
        return operate_on_row(r, start, end, buffer, work);
    }

    const T* column(size_t c, T* buffer, ColumnBlockWorkspace* work) const {
        const auto& deets = work->block();
        size_t start = deets.first, end = start + deets.second;
        return operate_on_column(c, start, end, buffer, work);
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* vbuffer, IDX* ibuffer, RowBlockWorkspace* work, bool sorted=true) const {
        const auto& deets = work->block();
        size_t start = deets.first, end = start + deets.second;
        return operate_on_row(r, vbuffer, ibuffer, work, start, end, sorted);
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* vbuffer, IDX* ibuffer, ColumnBlockWorkspace* work, bool sorted=true) const {
        const auto& deets = work->block();
        size_t start = deets.first, end = start + deets.second;
        return operate_on_column(c, vbuffer, ibuffer, work, start, end, sorted);
    }

private:
    template<class SomeWorkspace>
    const T* operate_on_row(size_t r, size_t start, size_t end, T* buffer, SomeWorkspace* work) const {
        const T* raw = mat->row(r, buffer, work);
        for (size_t i = start; i < end; ++i, ++raw) {
            buffer[i - start] = operation(r, i, *raw);
        }
        return buffer;
    }

    template<class SomeWorkspace>
    const T* operate_on_column(size_t c, size_t start, size_t end, T* buffer, SomeWorkspace* work) const {
        const T* raw = mat->column(c, buffer, work);
        for (size_t i = start; i < end; ++i, ++raw) {
            buffer[i - start] = operation(i, c, *raw);
        }
        return buffer;
    }

    template<class SomeWorkspace> 
    SparseRange<T, IDX> operate_on_row(size_t r, T* vbuffer, IDX* ibuffer, SomeWorkspace* work, bool sorted) const {
        auto raw = mat->sparse_row(r, vbuffer, ibuffer, work, sorted);
        for (size_t i = 0; i < raw.number; ++i) {
            vbuffer[i] = operation(r, raw.index[i], raw.value[i]);
        }
        return SparseRange<T, IDX>(raw.number, vbuffer, raw.index);
    }

    template<class SomeWorkspace>
    SparseRange<T, IDX> operate_on_row(size_t r, T* vbuffer, IDX* ibuffer, SomeWorkspace* work, size_t start, size_t end, bool sorted) const {
        if constexpr(OP::sparse) {
            return operate_on_row(r, vbuffer, ibuffer, work, sorted);
        } else {
            auto ptr = mat->row(r, vbuffer, work);
            auto original_vbuffer = vbuffer;
            auto original_ibuffer = ibuffer;
            for (size_t i = start; i < end; ++i, ++vbuffer, ++ibuffer, ++ptr) {
                (*vbuffer) = operation(r, i, *ptr);
                (*ibuffer) = i;
            }
            return SparseRange<T, IDX>(end - start, original_vbuffer, original_ibuffer); 
        }
    }

    template<class SomeWorkspace>
    SparseRange<T, IDX> operate_on_column(size_t c, T* vbuffer, IDX* ibuffer, SomeWorkspace* work, bool sorted) const {
        auto raw = mat->sparse_column(c, vbuffer, ibuffer, work, sorted);
        for (size_t i = 0; i < raw.number; ++i) {
            vbuffer[i] = operation(raw.index[i], c, raw.value[i]);
        }
        return SparseRange<T, IDX>(raw.number, vbuffer, raw.index);
    }

    template<class SomeWorkspace>
    SparseRange<T, IDX> operate_on_column(size_t c, T* vbuffer, IDX* ibuffer, SomeWorkspace* work, size_t start, size_t end, bool sorted) const {
        if constexpr(OP::sparse) {
            return operate_on_column(c, vbuffer, ibuffer, work, sorted);
        } else {
            auto ptr = mat->column(c, vbuffer, work);
            auto original_vbuffer = vbuffer;
            auto original_ibuffer = ibuffer;
            for (size_t i = start; i < end; ++i, ++vbuffer, ++ibuffer, ++ptr) {
                (*vbuffer) = operation(i, c, *ptr);
                (*ibuffer) = i;
            }
            return SparseRange<T, IDX>(end - start, original_vbuffer, original_ibuffer); 
        }
    }

public:
    std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(std::vector<IDX> indices) const {
        return mat->new_row_workspace(std::move(indices));
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX> > new_column_workspace(std::vector<IDX> indices) const {
        return mat->new_column_workspace(std::move(indices));
    }

    const T* row(size_t r, T* buffer, RowIndexWorkspace<IDX>* work) const {
        const T* raw = mat->row(r, buffer, work);
        const auto& indices = work->indices();
        for (size_t i = 0, end = indices.size(); i < end; ++i, ++raw) {
            buffer[i] = operation(r, indices[i], *raw);
        }
        return buffer;
    }

    const T* column(size_t c, T* buffer, ColumnIndexWorkspace<IDX>* work) const {
        const T* raw = mat->column(c, buffer, work);
        const auto& indices = work->indices();
        for (size_t i = 0, end = indices.size(); i < end; ++i, ++raw) {
            buffer[i] = operation(indices[i], c, *raw);
        }
        return buffer;
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* vbuffer, IDX* ibuffer, RowIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(OP::sparse) {
            return operate_on_row(r, vbuffer, ibuffer, work, sorted);
        } else {
            auto ptr = mat->row(r, vbuffer, work);
            const auto& indices = work->indices();
            size_t len = indices.size();

            auto original_vbuffer = vbuffer;
            for (size_t i = 0; i < len; ++i, ++vbuffer, ++ptr) {
                (*vbuffer) = operation(r, indices[i], *ptr);
            }
            std::copy(indices.begin(), indices.end(), ibuffer); // avoid lifetime issues with RowIndexWorkspace's indices.

            return SparseRange<T, IDX>(len, original_vbuffer, ibuffer);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* vbuffer, IDX* ibuffer, ColumnIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(OP::sparse) {
            return operate_on_column(c, vbuffer, ibuffer, work, sorted);
        } else {
            auto ptr = mat->column(c, vbuffer, work);
            const auto& indices = work->indices();
            size_t len = indices.size();

            auto original_vbuffer = vbuffer;
            for (size_t i = 0; i < len; ++i, ++vbuffer, ++ptr) {
                (*vbuffer) = operation(indices[i], c, *ptr);
            }
            std::copy(indices.begin(), indices.end(), ibuffer); // avoid lifetime issues with RowIndexWorkspace's indices.

            return SparseRange<T, IDX>(len, original_vbuffer, ibuffer); 
        }
    }
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam MAT A specialized `Matrix`, to be automatically deducted.
 * @tparam OP Helper class defining the operation.
 *
 * @param p Pointer to a `Matrix`.
 * @param op Instance of the operation helper class.
 */
template<class MAT, class OP>
std::shared_ptr<MAT> make_DelayedIsometricOp(std::shared_ptr<MAT> p, OP op) {
    return std::shared_ptr<MAT>(
        new DelayedIsometricOp<typename MAT::data_type, typename MAT::index_type, typename std::remove_reference<OP>::type>(
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
