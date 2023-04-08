#ifndef TATAMI_MATRIX_H
#define TATAMI_MATRIX_H

#include "SparseRange.hpp"
#include "Workspace.hpp"
#include <algorithm>
#include <numeric>

/**
 * @file Matrix.hpp
 *
 * @brief Virtual class for a matrix of some numeric type.
 */

namespace tatami {

/**
 * @brief Virtual class for a matrix with a defined type.
 * 
 * @tparam T Type of the matrix data.
 * @tparam IDX Type of the row/column indices.
 */
template <typename T, typename IDX = int>
class Matrix {
public:
    virtual ~Matrix() = default;

    // Defining the other constructors for rule of 5. The move constructors
    // don't get auto-made when a destructor is declared, and if I'm defining
    // them, I might as well define the copy constructors.  Technically I
    // suppose I don't need them because this class is just an interface, but
    // who knows if I'll add some move-able stuff here later.

    /**
     * Default move constructor.
     */
    Matrix(Matrix&&) = default;

    /**
     * Default move assignment operator.
     */
    Matrix& operator=(Matrix&&) = default;

    /**
     * Default copy constructor.
     */
    Matrix(const Matrix&) = default;

    /**
     * Default copy assignment operator.
     */
    Matrix& operator=(const Matrix&) = default;

protected:
    Matrix() = default;

public:
    /** 
     * Type of data to be returned by getters.
     */
    typedef T data_type;

    /** 
     * Type of index to be returned by the sparse getters.
     */
    typedef IDX index_type;
public:
    /**
     * @return Number of rows.
     */
    virtual size_t nrow() const = 0;

    /**
     * @return Number of columns.
     */
    virtual size_t ncol() const = 0;

    /**
     * @return Is this matrix sparse?
     * Defaults to `false` if no specialized method is provided in derived classes.
     */
    virtual bool sparse() const = 0;

    /**
     * @return The preferred dimension for extracting values.
     * If `true`, row-wise extraction is preferred; if `false`, column-wise extraction is preferred.
     */
    virtual bool prefer_rows() const = 0;

    /**
     * @return A `pair` containing the number of matrix elements that prefer row-level access (`first`) or column-level access (`second`).
     *
     * This method is useful for determining the return value of `prefer_rows()` in combined matrices consisting of both row- and column-preferred submatrices.
     * In such cases, the net preference can be determined based on the combined size of the submatrices for each preference.
     *
     * For simpler matrices, the return value contains the total size of the matrix in one of the `double`s and zero in the other.
     */
    virtual std::pair<double, double> dimension_preference () const {
        double size = static_cast<double>(nrow()) * static_cast<double>(ncol());
        if (prefer_rows()) {
            return std::make_pair(size, 0.0);
        } else {
            return std::make_pair(0.0, size);
        }
    }

    /****************************************
     ***** Workspace-generating methods *****
     ****************************************/
public:
    /**
     * @param cache Whether to cache information from each call to `row()` with this workspace.
     * This enables faster iterations if the same row is extracted during multiple passes over the matrix.
     * @return A shared pointer to a `DenseRowWorkspace` for row-wise data extraction, or a null pointer if no workspace is required.
     */
    virtual std::shared_ptr<DenseRowWorkspace> dense_row_workspace(bool cache = false) const = 0;

    /**
     * @param cache Whether to cache information from each call to `column()` with this workspace.
     * This enables faster iterations if the same column is extracted during multiple passes over the matrix.
     * @return A shared pointer to a `DenseColumnWorkspace` for column-wise data extraction, or a null pointer if no workspace is required.
     */
    virtual std::shared_ptr<DenseColumnWorkspace> dense_column_workspace(bool cache = false) const = 0;

    /**
     * @param cache Whether to cache information from each call to `sparse_row()` with this workspace.
     * This enables faster iterations if the same row is extracted during multiple passes over the matrix.
     * @param mode Whether to extract the values, or the indices, or both.
     * This can be used to avoid unnecessary computation and copying.
     * @param sorted Whether the sparse values should be sorted by increasing indices.
     * Setting this to `false` can reduce computational work in situations where the order of non-zero elements does not matter.
     *
     * @return A shared pointer to a `SparseRowWorkspace` for row-wise data extraction, or a null pointer if no workspace is required.
     */
    virtual std::shared_ptr<SparseRowWorkspace> sparse_row_workspace(bool cache = false, SparseExtractMode mode = SparseExtractMode::BOTH, bool sorted = true) const = 0;

    /**
     * @param cache Whether to cache information from each call to `sparse_column()` with this workspace.
     * This enables faster iterations if the same column is extracted during multiple passes over the matrix.
     * @param mode Whether to extract the values, or the indices, or both.
     * This can be used to avoid unnecessary computation and copying.
     * @param sorted Whether the sparse values should be sorted by increasing indices.
     * Setting this to `false` can reduce computational work in situations where the order of non-zero elements does not matter.
     *
     * @return A shared pointer to a `SparseColumnWorkspace` for column-wise data extraction, or a null pointer if no workspace is required.
     */
    virtual std::shared_ptr<SparseColumnWorkspace> sparse_column_workspace(bool cache = false, SparseExtractMode mode = SparseExtractMode::BOTH, bool sorted = true) const = 0;

    /**
     * @param start Index of the first column in the block.
     * @param length Number of columns in the block.
     * @param cache Whether to cache information from each call to `row()` with this workspace.
     * This enables faster iterations if the same row is extracted during multiple passes over the matrix.
     *
     * @return A shared pointer to a `DenseRowBlockWorkspace` for row-wise extraction of data from a contiguous block of columns from `[start, start + length)`. 
     */
    virtual std::shared_ptr<DenseRowBlockWorkspace> dense_row_workspace(size_t start, size_t length, bool cache = false) const = 0;

    /**
     * @param start Index of the first row in the block.
     * @param length Number of rows in the block.
     * @param cache Whether to cache information from each call to `column()` with this workspace.
     * This enables faster iterations if the same column is extracted during multiple passes over the matrix.
     *
     * @return A shared pointer to a `DenseRowBlockWorkspace` for column-wise extraction of data from a contiguous block of rows from `[start, start + length)`.
     */
    virtual std::shared_ptr<DenseColumnBlockWorkspace> dense_column_workspace(size_t start, size_t length, bool cache = false) const = 0;

    /**
     * @param start Index of the first column in the block.
     * @param length Number of columns in the block.
     * @param cache Whether to cache information from each call to `sparse_row()` with this workspace.
     * This enables faster iterations if the same row is extracted during multiple passes over the matrix.
     * @param mode Whether to extract the values, or the indices, or both.
     * This can be used to avoid unnecessary computation and copying.
     * @param sorted Whether the sparse values should be sorted by increasing indices.
     * Setting this to `false` can reduce computational work in situations where the order of non-zero elements does not matter.
     *
     * @return A shared pointer to a `SparseRowBlockWorkspace` for row-wise extraction of data from a contiguous block of columns from `[start, start + length)`. 
     */
    virtual std::shared_ptr<SparseRowBlockWorkspace> sparse_row_workspace(size_t start, size_t length, 
        bool cache = false, SparseExtractMode mode = SparseExtractMode::BOTH, bool sorted = true) const = 0;

    /**
     * @param start Index of the first row in the block.
     * @param length Number of rows in the block.
     * @param cache Whether to cache information from each call to `sparse_column()` with this workspace.
     * This enables faster iterations if the same column is extracted during multiple passes over the matrix.
     * @param mode Whether to extract the values, or the indices, or both.
     * This can be used to avoid unnecessary computation and copying.
     * @param sorted Whether the sparse values should be sorted by increasing indices.
     * Setting this to `false` can reduce computational work in situations where the order of non-zero elements does not matter.
     *
     * @return A shared pointer to a `SparseRowBlockWorkspace` for column-wise extraction of data from a contiguous block of rows from `[start, start + length)`.
     */
    virtual std::shared_ptr<SparseColumnBlockWorkspace> sparse_column_workspace(size_t start, size_t length,
        bool cache = false, SparseExtractMode mode = SparseExtractMode::BOTH, bool sorted = true) const = 0;

    // Note to self: we use a vector rather than taking a pointer to the
    // indices. This is because workspaces of wrapper matrices may need to
    // construct their own vectors, in which case the validity of the pointers
    // is ambiguous, e.g., is the saved pointer to a vector still valid after
    // the vector is moved to a new location? See also discussion at:
    // https://stackoverflow.com/questions/11021764/does-moving-a-vector-invalidate-iterators

    /**
     * @param indices Vector containing sorted and unique column indices.
     * @param cache Whether to cache information from each call to `row()` with this workspace.
     * This enables faster iterations if the same row is extracted during multiple passes over the matrix.
     *
     * @return A shared pointer to a `DenseRowIndexWorkspace` for row-wise extraction of data from a subset of columns defined by `indices`.
     */
    virtual std::shared_ptr<DenseRowIndexWorkspace<IDX> > dense_row_workspace(std::vector<IDX> indices, bool cache = false) const = 0;

    /**
     * @param indices Vector containing sorted and unique row indices.
     * @param cache Whether to cache information from each call to `column()` with this workspace.
     * This enables faster iterations if the same row is extracted during multiple passes over the matrix.
     *
     * @return A shared pointer to a `DenseColumnIndexWorkspace` for column-wise extraction of data from a subset of rows defined by `indices`
     */
    virtual std::shared_ptr<DenseColumnIndexWorkspace<IDX> > dense_column_workspace(std::vector<IDX> indices, bool cache = false) const = 0;

    /**
     * @param indices Vector containing sorted and unique column indices.
     * @param cache Whether to cache information from each call to `sparse_row()` with this workspace.
     * This enables faster iterations if the same row is extracted during multiple passes over the matrix.
     * @param mode Whether to extract the values, or the indices, or both.
     * This can be used to avoid unnecessary computation and copying.
     * @param sorted Whether the sparse values should be sorted by increasing indices.
     * Setting this to `false` can reduce computational work in situations where the order of non-zero elements does not matter.
     *
     * @return A shared pointer to a `SparseRowIndexWorkspace` for row-wise extraction of data from a subset of columns defined by `indices`.
     */
    virtual std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(std::vector<IDX> indices, 
        bool cache = false, SparseExtractMode mode = SparseExtractMode::BOTH, bool sorted = true) const = 0;

    /**
     * @param indices Vector containing sorted and unique row indices.
     * @param cache Whether to cache information from each call to `sparse_column()` with this workspace.
     * This enables faster iterations if the same column is extracted during multiple passes over the matrix.
     * @param mode Whether to extract the values, or the indices, or both.
     * This can be used to avoid unnecessary computation and copying.
     * @param sorted Whether the sparse values should be sorted by increasing indices.
     * Setting this to `false` can reduce computational work in situations where the order of non-zero elements does not matter.
     *
     * @return A shared pointer to a `SparseColumnIndexWorkspace` for column-wise extraction of data from a subset of rows defined by `indices`.
     */
    virtual std::shared_ptr<ColumnIndexWorkspace<IDX> > new_column_workspace(std::vector<IDX> indices, 
        bool cache = false, SparseExtractMode mode = SparseExtractMode::BOTH, bool sorted = true) const = 0;

    /*********************************
     ***** Dense virtual methods *****
     *********************************/
public:
    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param r Index of the row.
     * @param buffer Pointer to an array with enough space for at least `ncol()` values.
     * @param work Pointer to the workspace created with `dense_row_workspace()`.
     *
     * @return Pointer to the values of row `r`, containing `ncol()` valid entries.
     */
    virtual const T* row(size_t r, T* buffer, DenseRowWorkspace* work) const = 0;

    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param c Index of the column.
     * @param buffer Pointer to an array with enough space for at least `nrow()` values.
     * @param work Pointer to the workspace created with `dense_column_workspace()`.
     *
     * @return Pointer to the values of column `c`, containing `ncol()` valid entries.
     */
    virtual const T* column(size_t c, T* buffer, DenseColumnWorkspace* work) const = 0;

    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param r Index of the row.
     * @param buffer Pointer to an array with enough space for at least `DenseRowBlockWorkspace::length` values.
     * @param work Pointer to the workspace created with `dense_row_workspace()`.
     *
     * @return Pointer to an array containing a block of values from row `r`, where the array is of length `RowBlockWorkspace::length`.
     */
    virtual const T* row(size_t r, T* buffer, DenseRowBlockWorkspace* work) const = 0;

    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param c Index of the column.
     * @param buffer Pointer to an array with enough space for at least `DenseColumnBlockWorkspace::length` values.
     * @param work Pointer to the workspace created with `dense_column_workspace()`.
     *
     * @return Pointer to an array containing a block of values from column `c`, where the array is of length `ColumnBlockWorkspace::length`.
     */
    virtual const T* column(size_t c, T* buffer, DenseColumnBlockWorkspace* work) const = 0;

    /**
     * @param r Index of the row.
     * @param buffer Pointer to an array with enough space for at least `DenseRowIndexWorkspace::length` values.
     * @param work Pointer to the workspace created with `dense_row_workspace()`.
     *
     * @return The array pointed to by `buffer` is filled with a subset of values from column `r`.
     * `buffer` itself is returned.
     */
    virtual const T* row(size_t r, T* buffer, DenseRowIndexWorkspace<IDX>* work) const = 0;

    /**
     * @param c Index of the column.
     * @param buffer Pointer to an array with enough space for at least `DenseColumnIndexWorkspace::length` values.
     * @param work Pointer to the workspace created with `dense_column_workspace()`.
     *
     * @return The array pointed to by `buffer` is filled with a subset of values from column `c`.
     * `buffer` itself is returned.
     */
    virtual const T* column(size_t c, T* buffer, DenseColumnIndexWorkspace<IDX>* work) const = 0;

    /*************************************
     ***** Dense non-virtual methods *****
     *************************************/
private:
    template<typename X>
    static void copy_over(const X* src, X* dest, size_t n) {
        if (src!=dest) {
            std::copy(src, src + n, dest);
        }
        return;
    }

public:
    /**
     * @param r Index of the row.
     * @param buffer Pointer to an array with enough space for at least `ncol()` values.
     * @param work Pointer to the workspace created with `dense_row_workspace()`.
     *
     * @return The array at `buffer` is filled with the values of row `r`.
     * `buffer` itself is returned.
     */
    const T* row_copy(size_t r, T* buffer, DenseRowWorkspace* work) const {
        auto ptr = row(r, buffer, work);
        copy_over(ptr, buffer, ncol());
        return buffer;
    }

    /**
     * @param c Index of the column.
     * @param buffer Pointer to an array with enough space for at least `nrow()` values.
     * @param work Pointer to the workspace created with `dense_column_workspace()`.
     *
     * @return The array at `buffer` is filled with the values of column `c`.
     * `buffer` itself is returned.
     */
    const T* column_copy(size_t c, T* buffer, DenseColumnWorkspace* work) const {
        auto ptr = column(c, buffer, work);
        copy_over(ptr, buffer, nrow());
        return buffer;
    }

    /**
     * @param r Index of the row.
     * @param buffer Pointer to an array with enough space for at least `DenseRowBlockWorkspace::length` values.
     * @param work Pointer to the workspace created with `dense_row_workspace()`.
     *
     * @return The array pointed to by `buffer` is filled with a block of values from row `r`.
     * `buffer` itself is returned.
     */
    const T* row_copy(size_t r, T* buffer, DenseRowBlockWorkspace* work) const {
        auto ptr = row(r, buffer, work);
        copy_over(ptr, buffer, work->length());
        return buffer;
    }

    /**
     * @param c Index of the column.
     * @param buffer Pointer to an array with enough space for at least `DenseColumnBlockWorkspace::length` values.
     * @param work Pointer to the workspace created with `dense_column_workspace()`.
     *
     * @return The array pointed to by `buffer` is filled with a block of values from column `c`.
     * `buffer` itself is returned.
     */
    const T* column_copy(size_t c, T* buffer, DenseColumnBlockWorkspace* work) const {
        auto ptr = column(c, buffer, work);
        copy_over(ptr, buffer, work->length());
        return buffer;
    }

    /**
     * @param r Index of the row.
     * @param buffer Pointer to an array with enough space for at least `DenseRowIndexWorkspace::length` values.
     * @param work Pointer to the workspace created with `dense_row_workspace()`.
     *
     * @return The array pointed to by `buffer` is filled with a block of values from row `r`.
     * `buffer` itself is returned.
     */
    const T* row_copy(size_t r, T* buffer, DenseRowIndexWorkspace<IDX>* work) const {
        // Indexed extraction should almost certainly copy, but we'll just make sure.
        auto ptr = row(r, buffer, work);
        copy_over(ptr, buffer, work->length());
        return buffer;
    }

    /**
     * @param c Index of the column.
     * @param buffer Pointer to an array with enough space for at least `ColumnIndexWorkspace<IDX>::length` values.
     * @param work Pointer to the workspace created with `dense_column_workspace()`.
     *
     * @return The array pointed to by `buffer` is filled with a block of values from column `c`.
     * `buffer` itself is returned.
     */
    const T* column_copy(size_t c, T* buffer, DenseColumnIndexWorkspace<IDX>* work) const {
        // Indexed extraction should almost certainly copy, but we'll just make sure.
        auto ptr = column(c, buffer, work);
        copy_over(ptr, buffer, work->length());
        return buffer;
    }

public:
    /**
     * A more convenient but less efficient version of the `row()` method.
     * Callers do not have to supply `buffer`; instead a new allocation is performed every time.
     *
     * @param r Index of the row.
     * @param work Pointer to the workspace created with `dense_row_workspace()`.
     *
     * @return A vector containing the values of row `r`.
     */
    std::vector<T> row(size_t r, DenseRowWorkspace* work) const {
        std::vector<T> output(ncol());
        row_copy(r, output.data(), work);
        return output;
    }

    /**
     * A more convenient but (slightly) less efficient version of the `column()` method.
     * Callers do not have to supply `buffer`; instead a new allocation is performed every time.
     *
     * @param c Index of the column.
     * @param work Pointer to the workspace created with `dense_column_workspace()`.
     *
     * @return A vector containing the values of column `c`.
     */
    std::vector<T> column(size_t c, DenseColumnWorkspace* work) const {
        std::vector<T> output(nrow());
        column_copy(c, output.data(), work);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `row()` method.
     * Callers do not have to supply `buffer`; instead a new allocation is performed every time.
     *
     * @param r Index of the row.
     * @param work Pointer to the workspace created with `dense_row_workspace()`.
     *
     * @return A vector containing all values of row `r`.
     */
    std::vector<T> row(size_t r, DenseRowBlockWorkspace* work) const {
        std::vector<T> output(work->length());
        row_copy(r, output.data(), work);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `column()` method.
     * Callers do not have to supply `buffer`; instead a new allocation is performed every time.
     *
     * @param c Index of the column.
     * @param work Pointer to the workspace created with `dense_column_workspace()`.
     *
     * @return A vector containing all values of column `c`.
     */
    std::vector<T> column(size_t c, DenseColumnBlockWorkspace* work) const {
        std::vector<T> output(work->length());
        column_copy(c, output.data(), work);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `row()` method.
     * Callers do not have to supply `buffer`; instead a new allocation is performed every time.
     *
     * @param r Index of the row.
     * @param work Pointer to the workspace created with `dense_row_workspace()`.
     *
     * @return A vector containing all values of row `r`.
     */
    std::vector<T> row(size_t r, DenseRowIndexWorkspace<IDX>* work) const {
        std::vector<T> output(work->length());
        row_copy(r, output.data(), work);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `column()` method.
     * Callers do not have to supply `buffer`; instead a new allocation is performed every time.
     *
     * @param c Index of the column.
     * @param work Pointer to the workspace created with `dense_column_workspace()`.
     *
     * @return A vector containing all values of column `c`.
     */
    std::vector<T> column(size_t c, DenseColumnIndexWorkspace<IDX>* work) const {
        std::vector<T> output(work->length());
        column_copy(c, output.data(), work);
        return output;
    }

    /**********************************
     ***** Sparse virtual methods *****
     **********************************/
public:
    /**
     * `vbuffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `vbuffer`; if they are the same, `vbuffer` has been filled. 
     * The same applies for `ibuffer` and the returned `SparseRange::index` pointer.
     *
     * Entries in the output `value` array are not guaranteed to be non-zero.
     * If zeroes are explicitly initialized in the underlying representation, they will be reported here.
     * However, one can safely assume that all values _not_ in `value` are zero.
     *
     * @param r Index of the row.
     * @param vbuffer Pointer to an array with enough space for at least `ncol()` values.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::INDEX` or `SparseExtractMode::NONE`.
     * @param ibuffer Pointer to an array with enough space for at least `ncol()` indices.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::VALUE` or `SparseExtractMode::NONE`.
     * @param work Pointer to a workspace created with `sparse_row_workspace()`.
     *
     * @return A `SparseRange` object describing the contents of row `r`.
     * Either or both of `value` or `index` may be `NULL`, depending on `work->mode`.
     */
    virtual SparseRange<T, IDX> sparse_row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowWorkspace* work) const {
        const T* val = (vbuffer ? row(r, vbuffer, work) : NULL);
        if (sparse_extract_index(work->mode)) {
            std::iota(ibuffer, ibuffer + ncol(), static_cast<IDX>(0));
        }
        return SparseRange(ncol(), val, ibuffer); 
    }

    /**
     * `vbuffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `vbuffer`; if they are the same, `vbuffer` has been filled. 
     * The same applies for `ibuffer` and the returned `SparseRange::index` pointer.
     *
     * Entries in the output `value` array are not guaranteed to be non-zero.
     * If zeroes are explicitly initialized in the underlying representation, they will be reported here.
     * However, one can safely assume that all values _not_ in `value` are zero.
     *
     * @param c Index of the column.
     * @param vbuffer Pointer to an array with enough space for at least `nrow()` values.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::INDEX` or `SparseExtractMode::NONE`.
     * @param ibuffer Pointer to an array with enough space for at least `nrow()` indices.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::VALUE` or `SparseExtractMode::NONE`.
     * @param work Pointer to a workspace.
     *
     * @return A `SparseRange` object describing the contents of column `c`.
     * Either or both of `value` or `index` may be `NULL`, depending on `work->mode`.
     */
    virtual SparseRange<T, IDX> sparse_column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnWorkspace* work) const {
        const T* val = (vbuffer ? column(c, vbuffer, work) : NULL);
        if (sparse_extract_index(work->mode)) {
            std::iota(ibuffer, ibuffer + nrow(), static_cast<IDX>(0));
        }
        return SparseRange(nrow(), val, ibuffer); 
    }

    /**
     * `vbuffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `vbuffer`; if they are the same, `vbuffer` has been filled. 
     * The same applies for `ibuffer` and the returned `SparseRange::index` pointer.
     *
     * Entries in the output `value` array are not guaranteed to be non-zero.
     * If zeroes are explicitly initialized in the underlying representation, they will be reported here.
     * However, one can safely assume that all values _not_ in `value` are zero.
     *
     * @param r Index of the row.
     * @param vbuffer Pointer to an array with enough space for at least `SparseRowBlockWorkspace::length` values.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::INDEX` or `SparseExtractMode::NONE`.
     * @param ibuffer Pointer to an array with enough space for at least `SparseRowBlockWorkspace::length` indices.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::VALUE` or `SparseExtractMode::NONE`.
     * @param work Pointer to a workspace created with `sparse_row_workspace()`.
     * @param sorted Should the non-zero elements be sorted by their indices?
     *
     * @return A `SparseRange` object describing the contents of a block of columns in row `r`.
     * Either or both of `value` or `index` may be `NULL`, depending on `work->mode`.
     */
    virtual SparseRange<T, IDX> sparse_row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowBlockWorkspace* work) const {
        const T* val = (vbuffer ? row(r, vbuffer, work) : NULL);
        if (sparse_extract_index(work->mode)) {
            std::iota(ibuffer, ibuffer + work->length, static_cast<IDX>(work->start));
        }
        return SparseRange(work->length, val, ibuffer); 
    }

    /**
     * `vbuffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `vbuffer`; if they are the same, `vbuffer` has been filled. 
     * The same applies for `ibuffer` and the returned `SparseRange::index` pointer.
     *
     * Entries in the output `value` array are not guaranteed to be non-zero.
     * If zeroes are explicitly initialized in the underlying representation, they will be reported here.
     * However, one can safely assume that all values _not_ in `value` are zero.
     *
     * @param c Index of the column.
     * @param vbuffer Pointer to an array with enough space for at least `SparseColumnBlockWorkspace::length` values.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::INDEX` or `SparseExtractMode::NONE`.
     * @param ibuffer Pointer to an array with enough space for at least `SparseColumnBlockWorkspace::length` indices.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::VALUE` or `SparseExtractMode::NONE`.
     * @param work Pointer to a workspace created with `sparse_column_workspace()`.
     * @param sorted Should the non-zero elements be sorted by their indices?
     *
     * @return A `SparseRange` object describing the contents of a block of rows in column `c`.
     * Either or both of `value` or `index` may be `NULL`, depending on `work->mode`.
     */
    virtual SparseRange<T, IDX> sparse_column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnBlockWorkspace* work) const {
        const T* val = (vbuffer ? column(c, vbuffer, work) : NULL);
        if (sparse_extract_index(work->mode)) {
            std::iota(ibuffer, ibuffer + work->length, static_cast<IDX>(work->start));
        }
        return SparseRange(work->length, val, ibuffer); 
    }

    /**
     * Entries in the output `value` array are not guaranteed to be non-zero.
     * If zeroes are explicitly initialized in the underlying representation, they will be reported here.
     * However, one can safely assume that all values _not_ in `value` are zero.
     *
     * @param r Index of the row.
     * @param vbuffer Pointer to an array with enough space for at least `SparseRowIndexWorkspace::length` values.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::INDEX` or `SparseExtractMode::NONE`.
     * @param ibuffer Pointer to an array with enough space for at least `SparseRowIndexWorkspace::length` indices.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::VALUE` or `SparseExtractMode::NONE`.
     * @param work Pointer to a workspace created with `sparse_row_workspace()`.
     * @param sorted Should the non-zero elements be sorted by their indices?
     *
     * @return A `SparseRange` object describing the contents of a subset of columns in row `r`.
     * Either or both of `value` or `index` may be `NULL`, depending on `work->mode`.
     */
    virtual SparseRange<T, IDX> sparse_row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowIndexWorkspace<IDX>* work) const {
        const T* val = (vbuffer ? row(r, vbuffer, work) : NULL);
        const auto& indices = work->indices();
        if (sparse_extract_index(work->mode)) {
            std::copy(indices.begin(), indices.end(), ibuffer); // avoid lifetime issues with the workspace.
        }
        return SparseRange(indices.size(), val, ibuffer); 
    }

    /**
     * Entries in the output `value` array are not guaranteed to be non-zero.
     * If zeroes are explicitly initialized in the underlying representation, they will be reported here.
     * However, one can safely assume that all values _not_ in `value` are zero.
     *
     * @param c Index of the column.
     * @param vbuffer Pointer to an array with enough space for at least `SparseColumnIndexWorkspace::length` values.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::INDEX` or `SparseExtractMode::NONE`.
     * @param ibuffer Pointer to an array with enough space for at least `SparseColumnIndexWorkspace::length` indices.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::VALUE` or `SparseExtractMode::NONE`.
     * @param work Pointer to a workspace created with `sparse_column_workspace()`.
     *
     * @return A `SparseRange` object describing the contents of a subset of rows in column `c`.
     * Either or both of `value` or `index` may be `NULL`, depending on `work->mode`.
     */
    virtual SparseRange<T, IDX> sparse_column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnIndexWorkspace<IDX>* work) const {
        const T* val = (vbuffer ? column(c, vbuffer, work) : NULL);
        const auto& indices = work->indices();
        if (sparse_extract_index(work->mode)) {
            std::copy(indices.begin(), indices.end(), ibuffer); // avoid lifetime issues with the workspace.
        }
        return SparseRange(indices.size(), val, ibuffer); 
    }

    /**************************************
     ***** Sparse non-virtual methods *****
     **************************************/
private:
    static void copy_over(SparseRange<T, IDX>& output, T* vbuffer, IDX* ibuffer, SparseExtractMode mode) {
        if (sparse_extract_index(mode) && output.index != ibuffer) {
            copy_over(output.index, ibuffer, output.number);
            output.index = ibuffer;
        }

        if (sparse_extract_value(mode) && output.value != vbuffer) {
            copy_over(output.value, vbuffer, output.number);
            output.value = vbuffer;
        }
    }
    
public:
    /**
     * @param r Index of the row.
     * @param vbuffer Pointer to an array with enough space for at least `ncol()` values.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::INDEX` or `SparseExtractMode::NONE`.
     * @param ibuffer Pointer to an array with enough space for at least `ncol()` indices.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::VALUE` or `SparseExtractMode::NONE`.
     * @param work Pointer to a workspace created with `sparse_row_workspace()`.
     *
     * @return A `SparseRange` object describing the contents of row `r`.
     * Either or both of `value` or `index` may be `NULL`, depending on `work->mode`.
     */
    SparseRange<T, IDX> sparse_row_copy(size_t r, T* vbuffer, IDX* ibuffer, SparseRowWorkspace* work) const {
        auto output = sparse_row(r, vbuffer, ibuffer, work);
        copy_over(output, vbuffer, ibuffer, work->mode);
        return output;
    }

    /**
     * @param c Index of the column.
     * @param vbuffer Pointer to an array with enough space for at least `nrow()` values.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::INDEX` or `SparseExtractMode::NONE`.
     * @param ibuffer Pointer to an array with enough space for at least `nrow()` indices.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::VALUE` or `SparseExtractMode::NONE`.
     * @param work Pointer to a workspace created with `sparse_column_workspace()`.
     *
     * @return A `SparseRange` object describing the contents of column `c`.
     * Either or both of `value` or `index` may be `NULL`, depending on `work->mode`.
     */
    SparseRange<T, IDX> sparse_column_copy(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnWorkspace* work) const {
        auto output = sparse_column(c, vbuffer, ibuffer, work);
        copy_over(output, vbuffer, ibuffer, work->mode);
        return output;
    }

    /**
     * @param r Index of the row.
     * @param vbuffer Pointer to an array with enough space for at least `SparseRowBlockWorkspace::length` values.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::INDEX` or `SparseExtractMode::NONE`.
     * @param ibuffer Pointer to an array with enough space for at least `SparseRowBlockWorkspace::length` indices.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::VALUE` or `SparseExtractMode::NONE`.
     * @param work Pointer to a workspace created with `sparse_row_workspace()`.
     *
     * @return A `SparseRange` object describing the contents of a block of columns in row `r`.
     * Either or both of `value` or `index` may be `NULL`, depending on `work->mode`.
     */
    SparseRange<T, IDX> sparse_row_copy(size_t r, T* vbuffer, IDX* ibuffer, SparseRowBlockWorkspace* work) const {
        auto output = sparse_row(r, vbuffer, ibuffer, work);
        copy_over(output, vbuffer, ibuffer, work->mode);
        return output;
    }

    /**
     * @param c Index of the column.
     * @param vbuffer Pointer to an array with enough space for at least `SparseColumnBlockWorkspace::length` values.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::INDEX` or `SparseExtractMode::NONE`.
     * @param ibuffer Pointer to an array with enough space for at least `SparseColumnBlockWorkspace::length` indices.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::VALUE` or `SparseExtractMode::NONE`.
     * @param work Pointer to a workspace created with `sparse_row_workspace()`.
     *
     * @return A `SparseRange` object describing the contents of a block of rows in column `c`.
     * Either or both of `value` or `index` may be `NULL`, depending on `work->mode`.
     */
    SparseRange<T, IDX> sparse_column_copy(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnBlockWorkspace* work) const {
        auto output = sparse_column(c, vbuffer, ibuffer, work);
        copy_over(output, vbuffer, ibuffer, work->mode);
        return output;
    }

    /**
     * @param r Index of the row.
     * @param vbuffer Pointer to an array with enough space for at least `SparseRowIndexWorkspace::length` values.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::INDEX` or `SparseExtractMode::NONE`.
     * @param ibuffer Pointer to an array with enough space for at least `SparseRowIndexWorkspace::length` indices.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::VALUE` or `SparseExtractMode::NONE`.
     * @param work Pointer to a workspace created with `sparse_row_workspace()`.
     *
     * @return A `SparseRange` object describing the contents of a subset of columns in row `r`.
     * Either or both of `value` or `index` may be `NULL`, depending on `work->mode`.
     */
    SparseRange<T, IDX> sparse_row_copy(size_t r, T* vbuffer, IDX* ibuffer, SparseRowIndexWorkspace<IDX>* work) const {
        auto output = sparse_row(r, vbuffer, ibuffer, work);
        copy_over(output, vbuffer, ibuffer, work->mode);
        return output;
    }

    /**
     * @param c Index of the column.
     * @param vbuffer Pointer to an array with enough space for at least `SparseColumnIndexWorkspace::length` values.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::INDEX` or `SparseExtractMode::NONE`.
     * @param ibuffer Pointer to an array with enough space for at least `SparseColumnIndexWorkspace::length` indices.
     * Alternatively a null pointer, if `work->mode` is `SparseExtractMode::VALUE` or `SparseExtractMode::NONE`.
     * @param work Pointer to a workspace created with `sparse_row_workspace()`.
     *
     * @return A `SparseRange` object describing the contents of a subset of rows in column `c`.
     * Either or both of `value` or `index` may be `NULL`, depending on `work->mode`.
     */
    SparseRange<T, IDX> sparse_column_copy(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnIndexWorkspace<IDX>* work) const {
        auto output = sparse_column(c, vbuffer, ibuffer, work);
        copy_over(output, vbuffer, ibuffer, work->mode);
        return output;
    }

public:
    /**
     * A more convenient but less efficient version of the `sparse_row()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param r Index of the row.
     * @param work Pointer to a workspace created with `sparse_row_workspace()`.
     *
     * @return A `SparseRangeCopy` object containing the contents of row `r`.
     */
    SparseRangeCopy<T, IDX> sparse_row(size_t r, SparseRowWorkspace* work) const {
        SparseRangeCopy<T, IDX> output(ncol());
        auto ret = sparse_row_copy(r, output.value.data(), output.index.data(), work);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `sparse_column()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param c Index of the column.
     * @param work Pointer to a workspace created with `sparse_column_workspace()`.
     *
     * @return A `SparseRangeCopy` object containing the contents of column `c`.
     */
    SparseRangeCopy<T, IDX> sparse_column(size_t c, SparseColumnWorkspace* work) const {
        SparseRangeCopy<T, IDX> output(nrow());
        auto ret = sparse_column_copy(c, output.value.data(), output.index.data(), work);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `sparse_row()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param r Index of the row.
     * @param work Pointer to a workspace created with `sparse_row_workspace()`.
     *
     * @return A `SparseRangeCopy` object containing the contents of a block of columns in row `r`.
     */
    SparseRangeCopy<T, IDX> sparse_row(size_t r, SparseRowBlockWorkspace* work) const {
        SparseRangeCopy<T, IDX> output(work->length());
        auto ret = sparse_row_copy(r, output.value.data(), output.index.data(), work);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `sparse_column()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param c Index of the column.
     * @param work Pointer to a workspace created with `sparse_column_workspace()`.
     *
     * @return A `SparseRangeCopy` object containing the contents of a block of rows in column `c`.
     */
    SparseRangeCopy<T, IDX> sparse_column(size_t c, SparseColumnBlockWorkspace* work) const {
        SparseRangeCopy<T, IDX> output(work->length());
        auto ret = sparse_column_copy(c, output.value.data(), output.index.data(), work);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `sparse_row()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param r Index of the row.
     * @param work Pointer to a workspace created with `sparse_row_workspace()`.
     *
     * @return A `SparseRangeCopy` object containing the contents of a subset of columns in row `r`.
     */
    SparseRangeCopy<T, IDX> sparse_row(size_t r, SparseRowIndexWorkspace<IDX>* work) const {
        SparseRangeCopy<T, IDX> output(work->length());
        auto ret = sparse_row_copy(r, output.value.data(), output.index.data(), work);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `sparse_column()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param c Index of the column.
     * @param work Pointer to a workspace created with `sparse_column_workspace()`.
     *
     * @return A `SparseRangeCopy` object containing the contents of a subset of rows in column `c`.
     */
    SparseRangeCopy<T, IDX> sparse_column(size_t c, SparseColumnIndexWorkspace<IDX>* work) const {
        SparseRangeCopy<T, IDX> output(work->length());
        auto ret = sparse_column_copy(c, output.value.data(), output.index.data(), work);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }
};

/**
 * A convenient shorthand for the most common use case of double-precision matrices.
 */
using NumericMatrix = Matrix<double, int>;

/**
 * Create a new workspace for extraction of full rows or columns.
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
std::shared_ptr<Workspace<ROW> > new_workspace(const Matrix<T, IDX>* ptr) {
    if constexpr(ROW) {
        return ptr->new_row_workspace();
    } else {
        return ptr->new_column_workspace();
    }
}

/**
 * Create a new workspace for extraction of a contiguous block from each row or column. 
 * This is a convenience wrapper for switching between rows/columns at compile time.
 *
 * @tparam ROW Whether to create a workspace for row extraction.
 * @tparam T Type of the matrix data.
 * @tparam IDX Type of the row/column indices.
 *
 * @param ptr Pointer to a `Matrix` instance.
 * @param start Index of the first column (if `ROW = true`) or column (otherwise) in the block.
 * @param length Number of columns (if `ROW = true`) or rows (otherwise) in the block.
 *
 * @return A workspace for extracting rows from `ptr` if `ROW = true`, or columns otherwise.
 */
template<bool ROW, typename T, typename IDX>
std::shared_ptr<BlockWorkspace<ROW> > new_workspace(const Matrix<T, IDX>* ptr, size_t start, size_t length) {
    if constexpr(ROW) {
        return ptr->new_row_workspace(start, length);
    } else {
        return ptr->new_column_workspace(start, length);
    }
}

/**
 * Create a new workspace for extraction of a subset of entries from each row or column. 
 * This is a convenience wrapper for switching between rows/columns at compile time.
 *
 * @tparam ROW Whether to create a workspace for row extraction.
 * @tparam T Type of the matrix data.
 * @tparam IDX Type of the row/column indices.
 *
 * @param ptr Pointer to a `Matrix` instance.
 * @param indices Vector of unique and sorted indices for columns (if `ROW = true`) or rows (otherwise) in the subset.
 *
 * @return A workspace for extracting rows from `ptr` if `ROW = true`, or columns otherwise.
 */
template<bool ROW, typename T, typename IDX>
std::shared_ptr<BlockWorkspace<ROW> > new_workspace(const Matrix<T, IDX>* ptr, std::vector<IDX> indices) {
    if constexpr(ROW) {
        return ptr->new_row_workspace(std::move(indices));
    } else {
        return ptr->new_column_workspace(std::move(indices));
    }
}

}

#endif
