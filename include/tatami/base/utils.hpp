#ifndef TATAMI_BASE_UTILS_HPP
#define TATAMI_BASE_UTILS_HPP

#include "Matrix.hpp"
#include <memory>

namespace tatami {

/**
 * Create a new workspace for dense extraction of full rows or columns.
 * This is a convenience wrapper for switching between rows/columns at compile time.
 *
 * @tparam ROW Whether to create a workspace for row extraction.
 * @tparam T Type of the matrix data.
 * @tparam IDX Type of the row/column indices.
 *
 * @param ptr Pointer to a `Matrix` instance.
 *
 * @return A workspace for extracting rows from `ptr` if `ROW = true`, or columns otherwise.
 */
template<bool ROW, typename T, typename IDX>
std::shared_ptr<DenseWorkspace<ROW> > dense_workspace(const Matrix<T, IDX>* ptr) {
    if constexpr(ROW) {
        return ptr->dense_row_workspace();
    } else {
        return ptr->dense_column_workspace();
    }
}

/**
 * Create a new workspace for extraction of full rows or columns.
 * This is a convenience wrapper for switching between rows/columns at compile time.
 *
 * @tparam ROW Whether to create a workspace for row extraction.
 * @tparam SPARSE Whether to create a workspace for sparse extraction.
 * @tparam T Type of the matrix data.
 * @tparam IDX Type of the row/column indices.
 *
 * @param ptr Pointer to a `Matrix` instance.
 * @param options Optional parameters for workspace construction.
 *
 * @return A workspace for extracting rows from `ptr` if `ROW = true`, or columns otherwise.
 */
template<bool ROW, bool SPARSE, typename T, typename IDX>
auto new_workspace(const Matrix<T, IDX>* ptr, const WorkspaceOptions& options) {
    if constexpr(ROW) {
        if constexpr(SPARSE) {
            return ptr->sparse_row_workspace(options);
        } else {
            return ptr->dense_row_workspace(options);
        }
    } else {
        if constexpr(SPARSE) {
            return ptr->sparse_column_workspace(options);
        } else {
            return ptr->dense_column_workspace(options);
        }
    }
}

/**
 * Create a new workspace for extraction of full rows or columns.
 * This is a convenience wrapper for switching between rows/columns at compile time.
 *
 * @tparam ROW Whether to create a workspace for row extraction.
 * @tparam SPARSE Whether to create a workspace for sparse extraction.
 * @tparam T Type of the matrix data.
 * @tparam IDX Type of the row/column indices.
 *
 * @param ptr Pointer to a `Matrix` instance.
 *
 * @return A workspace for extracting rows from `ptr` if `ROW = true`, or columns otherwise.
 */
template<bool ROW, bool SPARSE, typename T, typename IDX>
auto new_workspace(const Matrix<T, IDX>* ptr) {
    return new_workspace<ROW, SPARSE>(ptr, WorkspaceOptions());
}

/**
 * Create a new workspace for extraction of a contiguous block from each row or column. 
 * This is a convenience wrapper for switching between rows/columns at compile time.
 *
 * @tparam ROW Whether to create a workspace for row extraction.
 * @tparam SPARSE Whether to create a workspace for sparse extraction.
 * @tparam T Type of the matrix data.
 * @tparam IDX Type of the row/column indices.
 *
 * @param ptr Pointer to a `Matrix` instance.
 * @param start Index of the first column (if `ROW = true`) or column (otherwise) in the block.
 * @param length Number of columns (if `ROW = true`) or rows (otherwise) in the block.
 * @param options Optional parameters for workspace construction.
 *
 * @return A workspace for extracting rows from `ptr` if `ROW = true`, or columns otherwise.
 */
template<bool ROW, bool SPARSE, typename T, typename IDX>
auto new_workspace(const Matrix<T, IDX>* ptr, size_t start, size_t length, const WorkspaceOptions& options) {
    if constexpr(ROW) {
        if constexpr(SPARSE) {
            return ptr->sparse_row_workspace(start, length, options);
        } else {
            return ptr->dense_row_workspace(start, length, options);
        }
    } else {
        if constexpr(SPARSE) {
            return ptr->sparse_column_workspace(start, length, options);
        } else {
            return ptr->dense_column_workspace(start, length, options);
        }
    }
}

/**
 * Create a new workspace for extraction of a contiguous block from each row or column. 
 * This is a convenience wrapper for switching between rows/columns at compile time.
 *
 * @tparam ROW Whether to create a workspace for row extraction.
 * @tparam SPARSE Whether to create a workspace for sparse extraction.
 * @tparam T Type of the matrix data.
 * @tparam IDX Type of the row/column indices.
 *
 * @param ptr Pointer to a `Matrix` instance.
 * @param start Index of the first column (if `ROW = true`) or column (otherwise) in the block.
 * @param length Number of columns (if `ROW = true`) or rows (otherwise) in the block.
 *
 * @return A workspace for extracting rows from `ptr` if `ROW = true`, or columns otherwise.
 */
template<bool ROW, bool SPARSE, typename T, typename IDX>
auto new_workspace(const Matrix<T, IDX>* ptr, size_t start, size_t length) {
    return new_workspace<ROW, SPARSE>(ptr, start, length, WorkspaceOptions());
}

/**
 * Create a new workspace for dense extraction of a subset of entries from each row or column. 
 * This is a convenience wrapper for switching between rows/columns at compile time.
 *
 * @tparam ROW Whether to create a workspace for row extraction.
 * @tparam SPARSE Whether to create a workspace for sparse extraction.
 * @tparam T Type of the matrix data.
 * @tparam IDX Type of the row/column indices.
 *
 * @param ptr Pointer to a `Matrix` instance.
 * @param indices Vector of unique and sorted indices for columns (if `ROW = true`) or rows (otherwise) in the subset.
 * @param options Optional parameters for workspace construction.
 *
 * @return A workspace for extracting rows from `ptr` if `ROW = true`, or columns otherwise.
 */
template<bool ROW, bool SPARSE, typename T, typename IDX>
auto new_workspace(const Matrix<T, IDX>* ptr, std::vector<IDX> indices, const WorkspaceOptions& options) {
    if constexpr(ROW) {
        if constexpr(SPARSE) {
            return ptr->sparse_row_workspace(std::move(indices), options);
        } else {
            return ptr->dense_row_workspace(std::move(indices), options);
        }
    } else {
        if constexpr(SPARSE) {
            return ptr->sparse_column_workspace(std::move(indices), options);
        } else {
            return ptr->dense_column_workspace(std::move(indices), options);
        }
    }
}

/**
 * Create a new workspace for dense extraction of a subset of entries from each row or column. 
 * This is a convenience wrapper for switching between rows/columns at compile time.
 *
 * @tparam ROW Whether to create a workspace for row extraction.
 * @tparam SPARSE Whether to create a workspace for sparse extraction.
 * @tparam T Type of the matrix data.
 * @tparam IDX Type of the row/column indices.
 *
 * @param ptr Pointer to a `Matrix` instance.
 * @param indices Vector of unique and sorted indices for columns (if `ROW = true`) or rows (otherwise) in the subset.
 *
 * @return A workspace for extracting rows from `ptr` if `ROW = true`, or columns otherwise.
 */
template<bool ROW, bool SPARSE, typename T, typename IDX>
auto new_workspace(const Matrix<T, IDX>* ptr, std::vector<IDX> indices) {
    return new_workspace<ROW, SPARSE>(ptr, std::move(indices), WorkspaceOptions());    
}

template<bool ROW, typename T, typename IDX, class SomeWorkspace>
const T* extract_dense(const Matrix<T, IDX>* ptr, size_t i, T* buffer, SomeWorkspace* work) {
    if constexpr(ROW) {
        return ptr->row(i, buffer, work);
    } else {
        return ptr->column(i, buffer, work);
    }
}

template<bool ROW, typename T, typename IDX, class SomeWorkspace>
SparseRange<T, IDX> extract_sparse(const Matrix<T, IDX>* ptr, size_t i, T* vbuffer, IDX* ibuffer, SomeWorkspace* work) {
    if constexpr(ROW) {
        return ptr->row(i, vbuffer, ibuffer, work);
    } else {
        return ptr->column(i, vbuffer, ibuffer, work);
    }
}

}

#endif
