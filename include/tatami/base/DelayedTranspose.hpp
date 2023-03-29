#ifndef TATAMI_DELAYED_TRANSPOSE
#define TATAMI_DELAYED_TRANSPOSE

#include "Matrix.hpp"
#include <memory>

/**
 * @file DelayedTranspose.hpp
 *
 * @brief Delayed transposition.
 *
 * This is equivalent to the `DelayedAperm` class in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @brief Delayed transposition of a matrix.
 *
 * Implements delayed transposition of a matrix.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam T Type of matrix value.
 * @tparam IDX Type of index value.
 */
template<typename T, typename IDX>
class DelayedTranspose : public Matrix<T, IDX> {
public:
    /**
     * @param p Pointer to the underlying (pre-transpose) matrix.
     */
    DelayedTranspose(std::shared_ptr<const Matrix<T, IDX> > p) : mat(std::move(p)) {}

private:
    std::shared_ptr<const Matrix<T, IDX> > mat;

public:
    size_t nrow() const {
        return mat->ncol();
    }

    size_t ncol() const {
        return mat->nrow();
    }

    bool sparse() const {
        return mat->sparse();
    }

    bool prefer_rows() const {
        return !mat->prefer_rows();
    }

    std::pair<double, double> dimension_preference () const {
        auto temp = mat->dimension_preference();
        return std::pair<double, double>(temp.second, temp.first);
    }

    using Matrix<T, IDX>::column;

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::sparse_column;

    using Matrix<T, IDX>::sparse_row;

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct TransposedWorkspace : public Workspace<ROW> {
        TransposedWorkspace(std::shared_ptr<Workspace<!ROW> > w) : base(std::move(w)) {}
        std::shared_ptr<Workspace<!ROW> > base;
    };

    typedef TransposedWorkspace<true> TransposedRowWorkspace;
    typedef TransposedWorkspace<false> TransposedColumnWorkspace;
    /**
     * @endcond
     */

    std::shared_ptr<RowWorkspace> new_row_workspace() const {
        return std::shared_ptr<RowWorkspace>(new TransposedRowWorkspace(mat->new_column_workspace()));
    }

    std::shared_ptr<ColumnWorkspace> new_column_workspace() const {
        return std::shared_ptr<ColumnWorkspace>(new TransposedColumnWorkspace(mat->new_row_workspace()));
    }

    const T* row(size_t r, T* buffer, RowWorkspace* work) const {
        return mat->column(r, buffer, static_cast<TransposedRowWorkspace*>(work)->base.get());
    }

    const T* column(size_t c, T* buffer, ColumnWorkspace* work) const {
        return mat->row(c, buffer, static_cast<TransposedColumnWorkspace*>(work)->base.get());
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowWorkspace* work, bool sorted=true) const {
        return mat->sparse_column(r, out_values, out_indices, static_cast<TransposedRowWorkspace*>(work)->base.get(), sorted);
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnWorkspace* work, bool sorted=true) const {
        return mat->sparse_row(c, out_values, out_indices, static_cast<TransposedColumnWorkspace*>(work)->base.get(), sorted);
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct TransposedBlockWorkspace : public BlockWorkspace<ROW> {
        TransposedBlockWorkspace(std::shared_ptr<BlockWorkspace<!ROW> > w) : base(std::move(w)) {}
        std::shared_ptr<BlockWorkspace<!ROW> > base;
        const std::pair<size_t, size_t>& block() const { return base->block(); }
    };

    typedef TransposedBlockWorkspace<true> TransposedRowBlockWorkspace;
    typedef TransposedBlockWorkspace<false> TransposedColumnBlockWorkspace;
    /**
     * @endcond
     */

    std::shared_ptr<RowBlockWorkspace> new_row_workspace(size_t start, size_t length) const {
        return std::shared_ptr<RowBlockWorkspace>(new TransposedRowBlockWorkspace(mat->new_column_workspace(start, length)));
    }

    std::shared_ptr<ColumnBlockWorkspace> new_column_workspace(size_t start, size_t length) const {
        return std::shared_ptr<ColumnBlockWorkspace>(new TransposedColumnBlockWorkspace(mat->new_row_workspace(start, length)));
    }

    const T* row(size_t r, T* buffer, RowBlockWorkspace* work) const {
        return mat->column(r, buffer, static_cast<TransposedRowBlockWorkspace*>(work)->base.get());
    }

    const T* column(size_t c, T* buffer, ColumnBlockWorkspace* work) const {
        return mat->row(c, buffer, static_cast<TransposedColumnBlockWorkspace*>(work)->base.get());
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowBlockWorkspace* work, bool sorted=true) const {
        return mat->sparse_column(r, out_values, out_indices, static_cast<TransposedRowBlockWorkspace*>(work)->base.get(), sorted);
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnBlockWorkspace* work, bool sorted=true) const {
        return mat->sparse_row(c, out_values, out_indices, static_cast<TransposedColumnBlockWorkspace*>(work)->base.get(), sorted);
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct TransposedIndexWorkspace : public IndexWorkspace<IDX, ROW> {
        TransposedIndexWorkspace(std::shared_ptr<IndexWorkspace<IDX, !ROW> > w) : base(std::move(w)) {}
        std::shared_ptr<IndexWorkspace<IDX, !ROW> > base;
        const std::vector<IDX>& indices() const { return base->indices(); }
    };

    typedef TransposedIndexWorkspace<true> TransposedRowIndexWorkspace;
    typedef TransposedIndexWorkspace<false> TransposedColumnIndexWorkspace;
    /**
     * @endcond
     */

    std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(std::vector<IDX> subset) const {
        return std::shared_ptr<RowIndexWorkspace<IDX> >(new TransposedRowIndexWorkspace(mat->new_column_workspace(std::move(subset))));
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX> > new_column_workspace(std::vector<IDX> subset) const {
        return std::shared_ptr<ColumnIndexWorkspace<IDX> >(new TransposedColumnIndexWorkspace(mat->new_row_workspace(std::move(subset))));
    }

    const T* row(size_t r, T* buffer, RowIndexWorkspace<IDX>* work) const {
        return mat->column(r, buffer, static_cast<TransposedRowIndexWorkspace*>(work)->base.get());
    }

    const T* column(size_t c, T* buffer, ColumnIndexWorkspace<IDX>* work) const {
        return mat->row(c, buffer, static_cast<TransposedColumnIndexWorkspace*>(work)->base.get());
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowIndexWorkspace<IDX>* work, bool sorted=true) const {
        return mat->sparse_column(r, out_values, out_indices, static_cast<TransposedRowIndexWorkspace*>(work)->base.get(), sorted);
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnIndexWorkspace<IDX>* work, bool sorted=true) const {
        return mat->sparse_row(c, out_values, out_indices, static_cast<TransposedColumnIndexWorkspace*>(work)->base.get(), sorted);
    }
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam MAT A specialized `Matrix`, to be automatically deducted.
 *
 * @param p Pointer to a `Matrix`.
 *
 * @return A pointer to a `DelayedTranspose` instance.
 */
template<class MAT>
std::shared_ptr<MAT> make_DelayedTranspose(std::shared_ptr<MAT> p) {
    return std::shared_ptr<MAT>(
        new DelayedTranspose<typename MAT::data_type, typename MAT::index_type>(std::move(p))
    );
}

}

#endif
