#ifndef TATAMI_DELAYED_TRANSPOSE
#define TATAMI_DELAYED_TRANSPOSE

#include "Matrix.hpp"
#include <memory>

/**
 * @file DelayedTranspose.hpp
 *
 * Delayed transpose, equivalent to the `DelayedAperm` class in the **DelayedArray** package.
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

public:
    const T* row(size_t r, T* buffer, size_t start, size_t end, Workspace* work=nullptr) const {
        return mat->column(r, buffer, start, end, work);
    }

    const T* column(size_t c, T* buffer, size_t start, size_t end, Workspace* work=nullptr) const {
        return mat->row(c, buffer, start, end, work);
    }

    using Matrix<T, IDX>::column;

    using Matrix<T, IDX>::row;

public:
    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, size_t start, size_t end, Workspace* work=nullptr, bool sorted=true) const {
        return mat->sparse_column(r, out_values, out_indices, start, end, work, sorted);
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, size_t start, size_t end, Workspace* work=nullptr, bool sorted=true) const {
        return mat->sparse_row(c, out_values, out_indices, start, end, work, sorted);
    }

    using Matrix<T, IDX>::sparse_column;

    using Matrix<T, IDX>::sparse_row;

public:
    /**
     * @return Number of rows after transposition.
     */
    size_t nrow() const {
        return mat->ncol();
    }
    
    /**
     * @return Number of columns after transposition.
     */
    size_t ncol() const {
        return mat->nrow();
    }

    /**
     * @return A null pointer or a shared pointer to a `Workspace` object, depending on the underlying (pre-subsetted) matrix.
     */
    std::shared_ptr<Workspace> new_workspace(bool row) const {
        return mat->new_workspace(!row);
    }

    /**
     * @return The sparsity status of the underlying (pre-subsetted) matrix.
     */
    bool sparse() const {
        return mat->sparse();
    }

    /**
     * @return Whether the underlying (pre-subsetted) matrix prefers row access.
     */
    bool prefer_rows() const {
        return !mat->prefer_rows();
    }
private:
    std::shared_ptr<const Matrix<T, IDX> > mat;
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
