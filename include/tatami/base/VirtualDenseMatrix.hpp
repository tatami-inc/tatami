#ifndef TATAMI_VIRTUAL_DENSE_MATRIX_H
#define TATAMI_VIRTUAL_DENSE_MATRIX_H

#include "SparseRange.hpp"
#include "Workspace.hpp"
#include "Matrix.hpp"
#include "utils.hpp"
#include <algorithm>
#include <numeric>

/**
 * @file VirtualDenseMatrix.hpp
 *
 * @brief Virtual class for a dense matrix of some numeric type.
 */

namespace tatami {

/**
 * @brief Virtual class for a dense matrix with a defined type.
 * 
 * @tparam T Type of the matrix data.
 * @tparam IDX Type of the row/column indices.
 *
 * This virtual class provides default methods for sparse extraction that just wrap the dense methods.
 * By inheriting from this class, implementers of dense matrices can skip the implementation of irrelevant methods for `Matrix::sparse_row()`, `Matrix::sparse_column_workspace()`, etc.,
 */
template <typename T, typename IDX = int>
class VirtualDenseMatrix : public Matrix<T, IDX> {
protected:
    VirtualDenseMatrix() = default;

public:
    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::column;

    using Matrix<T, IDX>::dense_row_workspace;

    using Matrix<T, IDX>::dense_column_workspace;

    using Matrix<T, IDX>::sparse_row_workspace;

    using Matrix<T, IDX>::sparse_column_workspace;

    bool sparse() const { return false; }

    /****************************************
     ***** Workspace-generating methods *****
     ****************************************/
public:
    /**
     * @cond
     */
    template<bool ROW>
    struct DefaultSparseWorkspace : public SparseWorkspace<ROW> {
        DefaultSparseWorkspace(const WorkspaceOptions& options, const Matrix<T, IDX>* self) : mode(options.mode) {
            if (sparse_extract_value(mode)) {
                dwork = new_workspace<ROW, false>(self);
            }
        }

        SparseExtractMode mode; 
        std::shared_ptr<DenseWorkspace<ROW> > dwork; // Need this for the dense extraction.
    };
    /**
     * @endcond
     */

    virtual std::shared_ptr<SparseRowWorkspace> sparse_row_workspace(const WorkspaceOptions& options) const {
        return std::shared_ptr<SparseRowWorkspace>(new DefaultSparseWorkspace<true>(options, this));
    }

    virtual std::shared_ptr<SparseColumnWorkspace> sparse_column_workspace(const WorkspaceOptions& options) const {
        return std::shared_ptr<SparseColumnWorkspace>(new DefaultSparseWorkspace<false>(options, this));
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct DefaultSparseBlockWorkspace : public SparseBlockWorkspace<ROW> {
        DefaultSparseBlockWorkspace(size_t start, size_t length, const WorkspaceOptions& options, const Matrix<T, IDX>* self) : SparseBlockWorkspace<ROW>(start, length), mode(options.mode) {
            if (sparse_extract_value(mode)) {
                dwork = new_workspace<ROW, false>(self, start, length);
            }
        }

        SparseExtractMode mode; 
        std::shared_ptr<DenseBlockWorkspace<ROW> > dwork;
    };
    /**
     * @endcond
     */

    virtual std::shared_ptr<SparseRowBlockWorkspace> sparse_row_workspace(size_t start, size_t length, const WorkspaceOptions& options) const {
        return std::shared_ptr<SparseRowBlockWorkspace>(new DefaultSparseBlockWorkspace<true>(start, length, options, this));
    }

    virtual std::shared_ptr<SparseColumnBlockWorkspace> sparse_column_workspace(size_t start, size_t length, const WorkspaceOptions& options) const {
        return std::shared_ptr<SparseColumnBlockWorkspace>(new DefaultSparseBlockWorkspace<false>(start, length, options, this));
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct DefaultSparseIndexWorkspace : public SparseIndexWorkspace<IDX, ROW> {
        DefaultSparseIndexWorkspace(std::vector<IDX> i, const WorkspaceOptions& options, const Matrix<T, IDX>* self) : SparseIndexWorkspace<IDX, ROW>(i.size()), mode(options.mode) {
            // Avoid holding an unnecessary copy of the indices.
            if (sparse_extract_value(mode)) {
                dwork = new_workspace<ROW, false>(self, std::move(i));
            } else {
                indices_ = std::move(i);
            }
        }

        SparseExtractMode mode;
        std::shared_ptr<DenseIndexWorkspace<IDX, ROW> > dwork;

        const std::vector<IDX>& indices() const { 
            if (dwork) {
                const auto& output = dwork->indices();
                return output;
            } else {
                return indices_; 
            }
        };

    private:
        std::vector<IDX> indices_;
    };
    /**
     * @endcond
     */

    virtual std::shared_ptr<SparseRowIndexWorkspace<IDX> > sparse_row_workspace(std::vector<IDX> indices, const WorkspaceOptions& options) const {
        return std::shared_ptr<SparseRowIndexWorkspace<IDX> >(new DefaultSparseIndexWorkspace<true>(std::move(indices), options, this));
    }

    virtual std::shared_ptr<SparseColumnIndexWorkspace<IDX> > sparse_column_workspace(std::vector<IDX> indices, const WorkspaceOptions& options) const {
        return std::shared_ptr<SparseColumnIndexWorkspace<IDX> >(new DefaultSparseIndexWorkspace<false>(std::move(indices), options, this));
    }

    /**********************************
     ***** Sparse virtual methods *****
     **********************************/
public:
    virtual SparseRange<T, IDX> row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowWorkspace* work) const {
        return default_dim<true>(r, vbuffer, ibuffer, work);
    }

    virtual SparseRange<T, IDX> column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnWorkspace* work) const {
        return default_dim<false>(c, vbuffer, ibuffer, work);
    }

    virtual SparseRange<T, IDX> row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowBlockWorkspace* work) const {
        return default_dim<true>(r, vbuffer, ibuffer, work);
    }

    virtual SparseRange<T, IDX> column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnBlockWorkspace* work) const {
        return default_dim<false>(c, vbuffer, ibuffer, work);
    }

    virtual SparseRange<T, IDX> row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowIndexWorkspace<IDX>* work) const {
        return default_dim<true>(r, vbuffer, ibuffer, work);
    }

    virtual SparseRange<T, IDX> column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnIndexWorkspace<IDX>* work) const {
        return default_dim<false>(c, vbuffer, ibuffer, work);
    }

private:
    template<bool ROW>
    SparseRange<T, IDX> default_dim(size_t i, T* vbuffer, IDX* ibuffer, Workspace<ROW, true>* work) const {
        auto wptr = static_cast<DefaultSparseWorkspace<ROW>*>(work);

        const T* val = NULL;
        if (sparse_extract_value(wptr->mode)) {
            val = extract_dense<ROW>(this, i, vbuffer, wptr->dwork.get());
        }

        size_t length = (ROW ? this->ncol() : this->nrow());
        if (sparse_extract_index(wptr->mode)) {
            std::iota(ibuffer, ibuffer + length, 0);
        } else {
            ibuffer = NULL;
        }

        return SparseRange(length, val, ibuffer); 
    }

    template<bool ROW>
    SparseRange<T, IDX> default_dim(size_t i, T* vbuffer, IDX* ibuffer, BlockWorkspace<ROW, true>* work) const {
        auto wptr = static_cast<DefaultSparseBlockWorkspace<ROW>*>(work);

        const T* val = NULL;
        if (sparse_extract_value(wptr->mode)) {
            val = extract_dense<ROW>(this, i, vbuffer, wptr->dwork.get());
        }

        if (sparse_extract_index(wptr->mode)) {
            std::iota(ibuffer, ibuffer + work->length, static_cast<IDX>(work->start));
        } else {
            ibuffer = NULL;
        }

        return SparseRange(work->length, val, ibuffer); 
    }

    template<bool ROW>
    SparseRange<T, IDX> default_dim(size_t i, T* vbuffer, IDX* ibuffer, IndexWorkspace<IDX, ROW, true>* work) const {
        auto wptr = static_cast<DefaultSparseIndexWorkspace<ROW>*>(work);

        const T* val = NULL;
        if (sparse_extract_value(wptr->mode)) {
            val = extract_dense<ROW>(this, i, vbuffer, wptr->dwork.get());
        }

        const auto& indices = wptr->indices();
        if (sparse_extract_index(wptr->mode)) {
            std::copy(indices.begin(), indices.end(), ibuffer); // avoid lifetime issues with the workspace.
        } else {
            ibuffer = NULL;
        }

        return SparseRange(indices.size(), val, ibuffer); 
    }
};

}

#endif
