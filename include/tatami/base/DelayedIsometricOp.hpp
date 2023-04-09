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

    using Matrix<T, IDX>::dense_row_workspace;

    using Matrix<T, IDX>::dense_column_workspace;

    using Matrix<T, IDX>::sparse_row_workspace;

    using Matrix<T, IDX>::sparse_column_workspace;

public:
    std::shared_ptr<DenseRowWorkspace> dense_row_workspace(const WorkspaceOptions& opt) const {
        return mat->dense_row_workspace(opt);
    }

    std::shared_ptr<DenseColumnWorkspace> dense_column_workspace(const WorkspaceOptions& opt) const {
        return mat->dense_column_workspace(opt);
    }

    const T* row(size_t r, T* buffer, DenseRowWorkspace* work) const {
        return operate_on_row(r, 0, mat->ncol(), buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnWorkspace* work) const {
        return operate_on_column(c, 0, mat->nrow(), buffer, work);
    }

    std::shared_ptr<SparseRowWorkspace> dense_row_workspace(const WorkspaceOptions& opt) const {
        return new_workspace<true, false>(opt);
    }

    std::shared_ptr<SparseColumnWorkspace> dense_column_workspace(const WorkspaceOptions& opt) const {
        return new_workspace<false, false>(opt);
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowWorkspace* work) const {
        return operate_on_row(r, vbuffer, ibuffer, work, 0, mat->ncol(), sorted);
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnWorkspace* work) const {
        return operate_on_column(c, vbuffer, ibuffer, work, 0, mat->nrow(), sorted);
    }

private:
    template<bool WORKROW, bool SPARSE>
    struct IntermediateWorkspace : public Workspace<WORKROW, SPARSE> {
        std::shared_ptr<Workspace<WORKROW, SPARSE> > inner;
        std::vector<IDX> ibuffer;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<Workspace<WORKROW, SPARSE> > new_intermediate_workspace(const WorkspaceOptions& opt) const {
        if (sparse_extract_index(opt.mode)) {
            return new_workspace<WORKROW, SPARSE>(mat.get(), opt);
        }

        auto copy = opt;
        copy.mode = SparseExtractMode::BOTH; // OP::needs_index<true>;
        return new_workspace<WORKROW, SPARSE>(mat.get(), copy);
    }

public:
    std::shared_ptr<DenseRowBlockWorkspace> new_row_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return mat->new_row_workspace(start, length, opt);
    }

    std::shared_ptr<DenseColumnBlockWorkspace> new_column_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return mat->new_column_workspace(start, length, opt);
    }

    const T* row(size_t r, T* buffer, DenseRowBlockWorkspace* work) const {
        const auto& deets = work->block();
        size_t start = deets.first, end = start + deets.second;
        return operate_on_row(r, start, end, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnBlockWorkspace* work) const {
        const auto& deets = work->block();
        size_t start = deets.first, end = start + deets.second;
        return operate_on_column(c, start, end, buffer, work);
    }

    std::shared_ptr<SparseRowBlockWorkspace> new_row_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return mat->new_row_workspace(start, length, opt);
    }

    std::shared_ptr<SparseColumnBlockWorkspace> new_column_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return mat->new_column_workspace(start, length, opt);
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowBlockWorkspace* work) const {
        const auto& deets = work->block();
        size_t start = deets.first, end = start + deets.second;
        return operate_on_row(r, vbuffer, ibuffer, work, start, end, sorted);
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnBlockWorkspace* work) const {
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
    SparseRange<T, IDX> operate_on_row(size_t r, T* vbuffer, IDX* ibuffer, SomeWorkspace* work) const {
        auto raw = mat->row(r, vbuffer, ibuffer, work, sorted);
        if (raw.value) {
            for (size_t i = 0; i < raw.number; ++i) {
                vbuffer[i] = operation(r, raw.index[i], raw.value[i]);
            }
        }
        return SparseRange<T, IDX>(raw.number, vbuffer, raw.index);
    }

    template<class SomeWorkspace>
    SparseRange<T, IDX> operate_on_row(size_t r, T* vbuffer, IDX* ibuffer, SomeWorkspace* work, size_t start, size_t end, bool sorted) const {
        if constexpr(OP::sparse) {
            return operate_on_row(r, vbuffer, ibuffer, work, sorted);
        }

        if (vbuffer) {
            auto ptr = mat->row(r, vbuffer, work);
            auto vcopy = vbuffer;
            for (size_t i = start; i < end; ++i, ++vcopy, ++ptr) {
                (*vcopy) = operation(r, i, *ptr);
            }
        }

        if (ibuffer) {
            std::iota(ibuffer, ibuffer + (end - start), start);
        }

        return SparseRange<T, IDX>(end - start, vbuffer, ibuffer);
    }

    template<class SomeWorkspace>
    SparseRange<T, IDX> operate_on_column(size_t c, T* vbuffer, IDX* ibuffer, SomeWorkspace* work, bool sorted) const {
        auto raw = mat->sparse_column(c, vbuffer, ibuffer, work, sorted);
        if (vbuffer) {
            for (size_t i = 0; i < raw.number; ++i) {
                vbuffer[i] = operation(raw.index[i], c, raw.value[i]);
            }
        }
        return SparseRange<T, IDX>(raw.number, vbuffer, raw.index);
    }

    template<class SomeWorkspace>
    SparseRange<T, IDX> operate_on_column(size_t c, T* vbuffer, IDX* ibuffer, SomeWorkspace* work, size_t start, size_t end, bool sorted) const {
        if constexpr(OP::sparse) {
            return operate_on_column(c, vbuffer, ibuffer, work, sorted);
        }

        if (vbuffer) {
            auto ptr = mat->column(c, vbuffer, work);
            auto vcopy = vbuffer;
            for (size_t i = start; i < end; ++i, ++vcopy, ++ptr) {
                (*vcopy) = operation(i, c, *ptr);
            }
        }

        if (ibuffer) {
            std::iota(ibuffer, ibuffer + (end - start), start);
        }

        return SparseRange<T, IDX>(end - start, vbuffer, ibuffer);
    }

public:
    std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(std::vector<IDX> indices, const WorkspaceOptions& opt) const {
        return mat->new_row_workspace(std::move(indices), opt);
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX> > new_column_workspace(std::vector<IDX> indices, const WorkspaceOptions& opt) const {
        return mat->new_column_workspace(std::move(indices), opt);
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
        }

        const auto& indices = work->indices();
        size_t len = indices.size();

        if (vbuffer) {
            auto ptr = mat->row(r, vbuffer, work);
            auto vcopy = vbuffer;
            for (size_t i = 0; i < len; ++i, ++vcopy, ++ptr) {
                (*vcopy) = operation(r, indices[i], *ptr);
            }
        }

        if (ibuffer) {
            std::copy(indices.begin(), indices.end(), ibuffer); // avoid lifetime issues with RowIndexWorkspace's indices.
        }

        return SparseRange<T, IDX>(len, vbuffer, ibuffer);
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* vbuffer, IDX* ibuffer, ColumnIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(OP::sparse) {
            return operate_on_column(c, vbuffer, ibuffer, work, sorted);
        }

        const auto& indices = work->indices();
        size_t len = indices.size();

        if (vbuffer) {
            auto ptr = mat->column(c, vbuffer, work);
            auto vcopy = vbuffer;
            for (size_t i = 0; i < len; ++i, ++vcopy, ++ptr) {
                (*vcopy) = operation(indices[i], c, *ptr);
            }
        }

        if (ibuffer) {
            std::copy(indices.begin(), indices.end(), ibuffer); // avoid lifetime issues with RowIndexWorkspace's indices.
        }

        return SparseRange<T, IDX>(len, vbuffer, ibuffer);
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
