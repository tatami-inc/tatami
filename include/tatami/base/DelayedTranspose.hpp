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

    using Matrix<T, IDX>::dense_row_workspace;

    using Matrix<T, IDX>::dense_column_workspace;

    using Matrix<T, IDX>::sparse_row_workspace;

    using Matrix<T, IDX>::sparse_column_workspace;

private:
    template<bool ROW, bool SPARSE>
    struct TransposedWorkspace : public Workspace<ROW, SPARSE> {
        TransposedWorkspace(std::shared_ptr<Workspace<!ROW, SPARSE> > w) : base(std::move(w)) {}
        std::shared_ptr<Workspace<!ROW, SPARSE> > base;
    };

    typedef TransposedWorkspace<true, false> DenseTransposedRowWorkspace;
    typedef TransposedWorkspace<false, false> DenseTransposedColumnWorkspace;
    typedef TransposedWorkspace<true, true> SparseTransposedRowWorkspace;
    typedef TransposedWorkspace<false, true> SparseTransposedColumnWorkspace;

public:
    std::shared_ptr<DenseRowWorkspace> dense_row_workspace(const WorkspaceOptions& opt) const {
        return std::shared_ptr<DenseRowWorkspace>(new DenseTransposedRowWorkspace(mat->dense_column_workspace(opt)));
    }

    std::shared_ptr<DenseColumnWorkspace> dense_column_workspace(const WorkspaceOptions& opt) const {
        return std::shared_ptr<DenseColumnWorkspace>(new DenseTransposedColumnWorkspace(mat->dense_row_workspace(opt)));
    }

    const T* row(size_t r, T* buffer, DenseRowWorkspace* work) const {
        return mat->column(r, buffer, static_cast<DenseTransposedRowWorkspace*>(work)->base.get());
    }

    const T* column(size_t c, T* buffer, DenseColumnWorkspace* work) const {
        return mat->row(c, buffer, static_cast<DenseTransposedColumnWorkspace*>(work)->base.get());
    }

    std::shared_ptr<SparseRowWorkspace> sparse_row_workspace(const WorkspaceOptions& opt) const {
        return std::shared_ptr<SparseRowWorkspace>(new SparseTransposedRowWorkspace(mat->sparse_column_workspace(opt)));
    }

    std::shared_ptr<SparseColumnWorkspace> sparse_column_workspace(const WorkspaceOptions& opt) const {
        return std::shared_ptr<SparseColumnWorkspace>(new SparseTransposedColumnWorkspace(mat->sparse_row_workspace(opt)));
    }

    SparseRange<T, IDX> row(size_t r, T* out_values, IDX* out_indices, SparseRowWorkspace* work) const {
        return mat->column(r, out_values, out_indices, static_cast<SparseTransposedRowWorkspace*>(work)->base.get());
    }

    SparseRange<T, IDX> column(size_t c, T* out_values, IDX* out_indices, SparseColumnWorkspace* work) const {
        return mat->row(c, out_values, out_indices, static_cast<SparseTransposedColumnWorkspace*>(work)->base.get());
    }

private:
    template<bool ROW, bool SPARSE>
    struct TransposedBlockWorkspace : public BlockWorkspace<ROW, SPARSE> {
        TransposedBlockWorkspace(std::shared_ptr<BlockWorkspace<!ROW, SPARSE> > w) : 
            BlockWorkspace<ROW, SPARSE>(w->start, w->length), base(std::move(w)) {}
        std::shared_ptr<BlockWorkspace<!ROW, SPARSE> > base;
    };

    typedef TransposedBlockWorkspace<true, false> DenseTransposedRowBlockWorkspace;
    typedef TransposedBlockWorkspace<false, false> DenseTransposedColumnBlockWorkspace;
    typedef TransposedBlockWorkspace<true, true> SparseTransposedRowBlockWorkspace;
    typedef TransposedBlockWorkspace<false, true> SparseTransposedColumnBlockWorkspace;

public:
    std::shared_ptr<DenseRowBlockWorkspace> dense_row_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return std::shared_ptr<DenseRowBlockWorkspace>(new DenseTransposedRowBlockWorkspace(mat->dense_column_workspace(start, length, opt)));
    }

    std::shared_ptr<DenseColumnBlockWorkspace> dense_column_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return std::shared_ptr<DenseColumnBlockWorkspace>(new DenseTransposedColumnBlockWorkspace(mat->dense_row_workspace(start, length, opt)));
    }

    const T* row(size_t r, T* buffer, DenseRowBlockWorkspace* work) const {
        return mat->column(r, buffer, static_cast<DenseTransposedRowBlockWorkspace*>(work)->base.get());
    }

    const T* column(size_t c, T* buffer, DenseColumnBlockWorkspace* work) const {
        return mat->row(c, buffer, static_cast<DenseTransposedColumnBlockWorkspace*>(work)->base.get());
    }

    std::shared_ptr<SparseRowBlockWorkspace> sparse_row_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return std::shared_ptr<SparseRowBlockWorkspace>(new SparseTransposedRowBlockWorkspace(mat->sparse_column_workspace(start, length, opt)));
    }

    std::shared_ptr<SparseColumnBlockWorkspace> sparse_column_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return std::shared_ptr<SparseColumnBlockWorkspace>(new SparseTransposedColumnBlockWorkspace(mat->sparse_row_workspace(start, length, opt)));
    }

    SparseRange<T, IDX> row(size_t r, T* out_values, IDX* out_indices, SparseRowBlockWorkspace* work) const {
        return mat->column(r, out_values, out_indices, static_cast<SparseTransposedRowBlockWorkspace*>(work)->base.get());
    }

    SparseRange<T, IDX> column(size_t c, T* out_values, IDX* out_indices, SparseColumnBlockWorkspace* work) const {
        return mat->row(c, out_values, out_indices, static_cast<SparseTransposedColumnBlockWorkspace*>(work)->base.get());
    }

private:
    template<bool ROW, bool SPARSE>
    struct TransposedIndexWorkspace : public IndexWorkspace<IDX, ROW, SPARSE> {
        TransposedIndexWorkspace(std::shared_ptr<IndexWorkspace<IDX, !ROW, SPARSE> > w) : 
            IndexWorkspace<IDX, ROW, SPARSE>(w->length), base(std::move(w)) {}
        std::shared_ptr<IndexWorkspace<IDX, !ROW, SPARSE> > base;
        const std::vector<IDX>& indices() const { return base->indices(); }
    };

    typedef TransposedIndexWorkspace<true, false> DenseTransposedRowIndexWorkspace;
    typedef TransposedIndexWorkspace<false, false> DenseTransposedColumnIndexWorkspace;
    typedef TransposedIndexWorkspace<true, true> SparseTransposedRowIndexWorkspace;
    typedef TransposedIndexWorkspace<false, true> SparseTransposedColumnIndexWorkspace;

public:
    std::shared_ptr<DenseRowIndexWorkspace<IDX> > dense_row_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return std::shared_ptr<DenseRowIndexWorkspace<IDX> >(new DenseTransposedRowIndexWorkspace(mat->dense_column_workspace(std::move(subset), opt)));
    }

    std::shared_ptr<DenseColumnIndexWorkspace<IDX> > dense_column_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return std::shared_ptr<DenseColumnIndexWorkspace<IDX> >(new DenseTransposedColumnIndexWorkspace(mat->dense_row_workspace(std::move(subset), opt)));
    }

    const T* row(size_t r, T* buffer, DenseRowIndexWorkspace<IDX>* work) const {
        return mat->column(r, buffer, static_cast<DenseTransposedRowIndexWorkspace*>(work)->base.get());
    }

    const T* column(size_t c, T* buffer, DenseColumnIndexWorkspace<IDX>* work) const {
        return mat->row(c, buffer, static_cast<DenseTransposedColumnIndexWorkspace*>(work)->base.get());
    }

    std::shared_ptr<SparseRowIndexWorkspace<IDX> > sparse_row_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return std::shared_ptr<SparseRowIndexWorkspace<IDX> >(new SparseTransposedRowIndexWorkspace(mat->sparse_column_workspace(std::move(subset), opt)));
    }

    std::shared_ptr<SparseColumnIndexWorkspace<IDX> > sparse_column_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return std::shared_ptr<SparseColumnIndexWorkspace<IDX> >(new SparseTransposedColumnIndexWorkspace(mat->sparse_row_workspace(std::move(subset), opt)));
    }

    SparseRange<T, IDX> row(size_t r, T* out_values, IDX* out_indices, SparseRowIndexWorkspace<IDX>* work) const {
        return mat->column(r, out_values, out_indices, static_cast<SparseTransposedRowIndexWorkspace*>(work)->base.get());
    }

    SparseRange<T, IDX> column(size_t c, T* out_values, IDX* out_indices, SparseColumnIndexWorkspace<IDX>* work) const {
        return mat->row(c, out_values, out_indices, static_cast<SparseTransposedColumnIndexWorkspace*>(work)->base.get());
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
