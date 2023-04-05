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
    virtual bool sparse() const { return false; }

    /**
     * @return The preferred dimension for extracting values.
     * If `true`, row-wise extraction is preferred; if `false`, column-wise extraction is preferred.
     * Defaults to `false` if no specialized method is provided in derived classes.
     */
    virtual bool prefer_rows() const { return false; }

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
     * @param cache Whether to cache information from each call to `row()` with this workspace, for faster iterations if the same row is extracted during multiple passes over the matrix.
     * @return A shared pointer to a `RowWorkspace` for row-wise data extraction, or a null pointer if no workspace is required.
     */
    virtual std::shared_ptr<RowWorkspace> new_row_workspace(bool cache = false) const = 0;

    /**
     * @param cache Whether to cache information from each call to `column()` with this workspace, for faster iterations if the same column is extracted during multiple passes over the matrix.
     * @return A shared pointer to a `ColumnWorkspace` for column-wise data extraction, or a null pointer if no workspace is required.
     */
    virtual std::shared_ptr<ColumnWorkspace> new_column_workspace(bool cache = false) const = 0;

    /**
     * @param start Index of the first column in the block.
     * @param length Number of columns in the block.
     * @param cache Whether to cache information from each call to `row()` with this workspace, for faster iterations if the same row is extracted during multiple passes over the matrix.
     *
     * @return A shared pointer to a `RowBlockWorkspace` for row-wise extraction of data from a contiguous block of columns from `[start, start + length)`. 
     */
    virtual std::shared_ptr<RowBlockWorkspace> new_row_workspace(size_t start, size_t length, bool cache = false) const = 0;

    /**
     * @param start Index of the first row in the block.
     * @param length Number of rows in the block.
     * @param cache Whether to cache information from each call to `column()` with this workspace, for faster iterations if the same column is extracted during multiple passes over the matrix.
     *
     * @return A shared pointer to a `RowBlockWorkspace` for column-wise extraction of data from a contiguous block of rows from `[start, start + length)`.
     */
    virtual std::shared_ptr<ColumnBlockWorkspace> new_column_workspace(size_t start, size_t length, bool cache = false) const = 0;

    // Note to self: we use a vector rather than taking a pointer to the
    // indices. This is because workspaces of wrapper matrices may need to
    // construct their own vectors, in which case the validity of the pointers
    // is ambiguous, e.g., is the saved pointer to a vector still valid after
    // the vector is moved to a new location? See also discussion at:
    // https://stackoverflow.com/questions/11021764/does-moving-a-vector-invalidate-iterators

    /**
     * @param indices Vector containing sorted and unique column indices.
     * @param cache Whether to cache information from each call to `row()` with this workspace, for faster iterations if the same row is extracted during multiple passes over the matrix.
     *
     * @return A shared pointer to a `RowIndexWorkspace` for row-wise extraction of data from a subset of columns defined by `indices`, or a null pointer if no workspace is required.
     */
    virtual std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(std::vector<IDX> indices, bool cache = false) const = 0;

    /**
     * @param indices Vector containing sorted and unique row indices.
     * @param cache Whether to cache information from each call to `column()` with this workspace, for faster iterations if the same column is extracted during multiple passes over the matrix.
     *
     * @return A shared pointer to a `ColumnIndexWorkspace` for column-wise extraction of data from a subset of rows defined by `indices`, or a null pointer if no workspace is required.
     */
    virtual std::shared_ptr<ColumnIndexWorkspace<IDX> > new_column_workspace(std::vector<IDX> indices, bool cache = false) const = 0;

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
     * @param work Pointer to the workspace created with `new_row_workspace()`.
     *
     * @return Pointer to the values of row `r`, containing `ncol()` valid entries.
     */
    virtual const T* row(size_t r, T* buffer, RowWorkspace* work) const = 0;

    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param c Index of the column.
     * @param buffer Pointer to an array with enough space for at least `last - first` values.
     * @param work Pointer to the workspace created with `new_column_workspace()`.
     *
     * @return Pointer to the values of column `c`, containing `ncol()` valid entries.
     */
    virtual const T* column(size_t c, T* buffer, ColumnWorkspace* work) const = 0;

    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param r Index of the row.
     * @param buffer Pointer to an array with enough space for at least `RowBlockWorkspace::length` values.
     * @param work Pointer to the workspace created with `new_row_workspace()`.
     *
     * @return Pointer to an array containing a block of values from row `r`, where the array is of length `RowBlockWorkspace::length`.
     */
    virtual const T* row(size_t r, T* buffer, RowBlockWorkspace* work) const = 0;

    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param c Index of the column.
     * @param buffer Pointer to an array with enough space for at least `ColumnBlockWorkspace::length` values.
     * @param work Pointer to the workspace created with `new_column_workspace()`.
     *
     * @return Pointer to an array containing a block of values from column `c`, where the array is of length `ColumnBlockWorkspace::length`.
     */
    virtual const T* column(size_t c, T* buffer, ColumnBlockWorkspace* work) const = 0;

    /**
     * @param r Index of the row.
     * @param buffer Pointer to an array with enough space for at least `RowIndexWorkspace::length` values.
     * @param work Pointer to the workspace created with `new_row_workspace()`.
     *
     * @return The array pointed to by `buffer` is filled with a subset of values from column `r`.
     * `buffer` itself is returned.
     */
    virtual const T* row(size_t r, T* buffer, RowIndexWorkspace<IDX>* work) const = 0;

    /**
     * @param c Index of the column.
     * @param buffer Pointer to an array with enough space for at least `ColumnIndexWorkspace::length` values.
     * @param work Pointer to the workspace created with `new_column_workspace()`.
     *
     * @return The array pointed to by `buffer` is filled with a subset of values from column `c`.
     * `buffer` itself is returned.
     */
    virtual const T* column(size_t c, T* buffer, ColumnIndexWorkspace<IDX>* work) const = 0;

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
     * @param work Pointer to the workspace created with `new_row_workspace()`.
     *
     * @return The array at `buffer` is filled with the values of row `r`.
     * `buffer` itself is returned.
     */
    const T* row_copy(size_t r, T* buffer, RowWorkspace* work) const {
        auto ptr = row(r, buffer, work);
        copy_over(ptr, buffer, ncol());
        return buffer;
    }

    /**
     * @param c Index of the column.
     * @param buffer Pointer to an array with enough space for at least `nrow()` values.
     * @param work Pointer to the workspace created with `new_column_workspace()`.
     *
     * @return The array at `buffer` is filled with the values of column `c`.
     * `buffer` itself is returned.
     */
    const T* column_copy(size_t c, T* buffer, ColumnWorkspace* work) const {
        auto ptr = column(c, buffer, work);
        copy_over(ptr, buffer, nrow());
        return buffer;
    }

    /**
     * @param r Index of the row.
     * @param buffer Pointer to an array with enough space for at least `RowBlockWorkspace::length` values.
     * @param work Pointer to the workspace created with `new_row_workspace()`.
     *
     * @return The array pointed to by `buffer` is filled with a block of values from row `r`.
     * `buffer` itself is returned.
     */
    const T* row_copy(size_t r, T* buffer, RowBlockWorkspace* work) const {
        auto ptr = row(r, buffer, work);
        copy_over(ptr, buffer, work->length());
        return buffer;
    }

    /**
     * @param c Index of the column.
     * @param buffer Pointer to an array with enough space for at least `ColumnBlockWorkspace::length` values.
     * @param work Pointer to the workspace created with `new_column_workspace()`.
     *
     * @return The array pointed to by `buffer` is filled with a block of values from column `c`.
     * `buffer` itself is returned.
     */
    const T* column_copy(size_t c, T* buffer, ColumnBlockWorkspace* work) const {
        auto ptr = column(c, buffer, work);
        copy_over(ptr, buffer, work->length());
        return buffer;
    }

    /**
     * @param r Index of the row.
     * @param buffer Pointer to an array with enough space for at least `RowIndexWorkspace<IDX>::length` values.
     * @param work Pointer to the workspace created with `new_row_workspace()`.
     *
     * @return The array pointed to by `buffer` is filled with a block of values from row `r`.
     * `buffer` itself is returned.
     */
    const T* row_copy(size_t r, T* buffer, RowIndexWorkspace<IDX>* work) const {
        // Indexed extraction should almost certainly copy, but we'll just make sure.
        auto ptr = row(r, buffer, work);
        copy_over(ptr, buffer, work->length());
        return buffer;
    }

    /**
     * @param c Index of the column.
     * @param buffer Pointer to an array with enough space for at least `ColumnIndexWorkspace<IDX>::length` values.
     * @param work Pointer to the workspace created with `new_column_workspace()`.
     *
     * @return The array pointed to by `buffer` is filled with a block of values from column `c`.
     * `buffer` itself is returned.
     */
    const T* column_copy(size_t c, T* buffer, ColumnIndexWorkspace<IDX>* work) const {
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
     * @param work Pointer to the workspace created with `new_row_workspace()`.
     *
     * @return A vector containing the values of row `r`.
     */
    std::vector<T> row(size_t r, RowWorkspace* work) const {
        std::vector<T> output(ncol());
        row_copy(r, output.data(), work);
        return output;
    }

    /**
     * A more convenient but (slightly) less efficient version of the `column()` method.
     * Callers do not have to supply `buffer`; instead a new allocation is performed every time.
     *
     * @param c Index of the column.
     * @param work Pointer to the workspace created with `new_column_workspace()`.
     *
     * @return A vector containing the values of column `c`.
     */
    std::vector<T> column(size_t c, ColumnWorkspace* work) const {
        std::vector<T> output(nrow());
        column_copy(c, output.data(), work);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `row()` method.
     * Callers do not have to supply `buffer`; instead a new allocation is performed every time.
     *
     * @param r Index of the row.
     * @param work Pointer to the workspace created with `new_row_workspace()`.
     *
     * @return A vector containing all values of row `r`.
     */
    std::vector<T> row(size_t r, RowBlockWorkspace* work) const {
        std::vector<T> output(work->length());
        row_copy(r, output.data(), work);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `column()` method.
     * Callers do not have to supply `buffer`; instead a new allocation is performed every time.
     *
     * @param c Index of the column.
     * @param work Pointer to the workspace created with `new_column_workspace()`.
     *
     * @return A vector containing all values of column `c`.
     */
    std::vector<T> column(size_t c, ColumnBlockWorkspace* work) const {
        std::vector<T> output(work->length());
        column_copy(c, output.data(), work);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `row()` method.
     * Callers do not have to supply `buffer`; instead a new allocation is performed every time.
     *
     * @param r Index of the row.
     * @param work Pointer to the workspace created with `new_row_workspace()`.
     *
     * @return A vector containing all values of row `r`.
     */
    std::vector<T> row(size_t r, RowIndexWorkspace<IDX>* work) const {
        std::vector<T> output(work->length());
        row_copy(r, output.data(), work);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `column()` method.
     * Callers do not have to supply `buffer`; instead a new allocation is performed every time.
     *
     * @param c Index of the column.
     * @param work Pointer to the workspace created with `new_column_workspace()`.
     *
     * @return A vector containing all values of column `c`.
     */
    std::vector<T> column(size_t c, ColumnIndexWorkspace<IDX>* work) const {
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
     * Values in `vbuffer` are not guaranteed to be non-zero.
     * If zeroes are explicitly initialized in the underlying representation, they will be reported here.
     * However, one can safely assume that all values _not_ in `vbuffer` are zero.
     *
     * Setting `sorted = false` can reduce computational work in situations where the order of non-zero elements does not matter.
     *
     * @param r Index of the row.
     * @param vbuffer Pointer to an array with enough space for at least `ncol()` values.
     * @param ibuffer Pointer to an array with enough space for at least `ncol()` indices.
     * @param work Pointer to a workspace created with `new_row_workspace()`.
     * @param sorted Should the non-zero elements be sorted by their indices?
     *
     * @return A `SparseRange` object describing the contents of row `r`.
     */
    virtual SparseRange<T, IDX> sparse_row(size_t r, T* vbuffer, IDX* ibuffer, RowWorkspace* work, bool sorted = true) const {
        const T* val = row(r, vbuffer, work);
        std::iota(ibuffer, ibuffer + ncol(), static_cast<IDX>(0));
        return SparseRange(ncol(), val, ibuffer); 
    }

    /**
     * `vbuffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `vbuffer`; if they are the same, `vbuffer` has been filled. 
     * The same applies for `ibuffer` and the returned `SparseRange::index` pointer.
     *
     * Values in `vbuffer` are not guaranteed to be non-zero.
     * If zeroes are explicitly initialized in the underlying representation, they will be reported here.
     * However, one can safely assume that all values _not_ in `vbuffer` are zero.
     *
     * Setting `sorted = false` can reduce computational work in situations where the order of non-zero elements does not matter.
     *
     * @param c Index of the column.
     * @param vbuffer Pointer to an array with enough space for at least `nrow()` values.
     * @param ibuffer Pointer to an array with enough space for at least `nrow()` indices.
     * @param work Pointer to a workspace.
     * @param sorted Should the non-zero elements be sorted by their indices?
     *
     * @return A `SparseRange` object describing the contents of column `c`.
     */
    virtual SparseRange<T, IDX> sparse_column(size_t c, T* vbuffer, IDX* ibuffer, ColumnWorkspace* work, bool sorted = true) const {
        const T* val = column(c, vbuffer, work);
        std::iota(ibuffer, ibuffer + nrow(), static_cast<IDX>(0));
        return SparseRange(nrow(), val, ibuffer); 
    }

    /**
     * `vbuffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `vbuffer`; if they are the same, `vbuffer` has been filled. 
     * The same applies for `ibuffer` and the returned `SparseRange::index` pointer.
     *
     * Values in `vbuffer` are not guaranteed to be non-zero.
     * If zeroes are explicitly initialized in the underlying representation, they will be reported here.
     * However, one can safely assume that all values _not_ in `vbuffer` are zero.
     *
     * Setting `sorted = false` can reduce computational work in situations where the order of non-zero elements does not matter.
     *
     * @param r Index of the row.
     * @param vbuffer Pointer to an array with enough space for at least `RowBlockWorkspace::length` values.
     * @param ibuffer Pointer to an array with enough space for at least `RowBlockWorkspace::length` indices.
     * @param work Pointer to a workspace created with `new_row_workspace()`.
     * @param sorted Should the non-zero elements be sorted by their indices?
     *
     * @return A `SparseRange` object describing the contents of a block of columns in row `r`.
     */
    virtual SparseRange<T, IDX> sparse_row(size_t r, T* vbuffer, IDX* ibuffer, RowBlockWorkspace* work, bool sorted = true) const {
        const T* val = row(r, vbuffer, work);
        const auto& deets = work->block();
        size_t start = deets.first, len = deets.second;
        std::iota(ibuffer, ibuffer + len, static_cast<IDX>(start));
        return SparseRange(len, val, ibuffer); 
    }

    /**
     * `vbuffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `vbuffer`; if they are the same, `vbuffer` has been filled. 
     * The same applies for `ibuffer` and the returned `SparseRange::index` pointer.
     *
     * Values in `vbuffer` are not guaranteed to be non-zero.
     * If zeroes are explicitly initialized in the underlying representation, they will be reported here.
     * However, one can safely assume that all values _not_ in `vbuffer` are zero.
     *
     * Setting `sorted = false` can reduce computational work in situations where the order of non-zero elements does not matter.
     *
     * @param c Index of the column.
     * @param vbuffer Pointer to an array with enough space for at least `ColumnBlockWorkspace::length` values.
     * @param ibuffer Pointer to an array with enough space for at least `ColumnBlockWorkspace::length` indices.
     * @param work Pointer to a workspace created with `new_column_workspace()`.
     * @param sorted Should the non-zero elements be sorted by their indices?
     *
     * @return A `SparseRange` object describing the contents of a block of rows in column `c`.
     */
    virtual SparseRange<T, IDX> sparse_column(size_t c, T* vbuffer, IDX* ibuffer, ColumnBlockWorkspace* work, bool sorted = true) const {
        const T* val = column(c, vbuffer, work);
        const auto& deets = work->block();
        size_t start = deets.first, len = deets.second;
        std::iota(ibuffer, ibuffer + len, static_cast<IDX>(start));
        return SparseRange(len, val, ibuffer); 
    }

    /**
     * Values in `vbuffer` are not guaranteed to be non-zero.
     * If zeroes are explicitly initialized in the underlying representation, they will be reported here.
     * However, one can safely assume that all values _not_ in `vbuffer` are zero.
     *
     * Setting `sorted = false` can reduce computational work in situations where the order of non-zero elements does not matter.
     *
     * @param r Index of the row.
     * @param vbuffer Pointer to an array with enough space for at least `RowIndexWorkspace::length` values.
     * @param ibuffer Pointer to an array with enough space for at least `RowIndexWorkspace::length` indices.
     * @param work Pointer to a workspace created with `new_row_workspace()`.
     * @param sorted Should the non-zero elements be sorted by their indices?
     *
     * @return A `SparseRange` object describing the contents of a subset of columns in row `r`.
     * `vbuffer` is set as `SparseRange::value`, while `ibuffer` is set as `SparseRange::index`.
     */
    virtual SparseRange<T, IDX> sparse_row(size_t r, T* vbuffer, IDX* ibuffer, RowIndexWorkspace<IDX>* work, bool sorted = true) const {
        const T* val = row(r, vbuffer, work);
        const auto& indices = work->indices();
        std::copy(indices.begin(), indices.end(), ibuffer); // avoid lifetime issues with the workspace.
        return SparseRange(indices.size(), val, ibuffer); 
    }

    /**
     * Values in `vbuffer` are not guaranteed to be non-zero.
     * If zeroes are explicitly initialized in the underlying representation, they will be reported here.
     * However, one can safely assume that all values _not_ in `vbuffer` are zero.
     *
     * Setting `sorted = false` can reduce computational work in situations where the order of non-zero elements does not matter.
     *
     * @param c Index of the column.
     * @param vbuffer Pointer to an array with enough space for at least `ColumnIndexWorkspace::length` values.
     * @param ibuffer Pointer to an array with enough space for at least `ColumnIndexWorkspace::length` indices.
     * @param work Pointer to a workspace created with `new_column_workspace()`.
     * @param sorted Should the non-zero elements be sorted by their indices?
     *
     * @return A `SparseRange` object describing the contents of a subset of rows in column `c`.
     * `vbuffer` is set as `SparseRange::value`, while `ibuffer` is set as `SparseRange::index`.
     */
    virtual SparseRange<T, IDX> sparse_column(size_t c, T* vbuffer, IDX* ibuffer, ColumnIndexWorkspace<IDX>* work, bool sorted = true) const {
        const T* val = column(c, vbuffer, work);
        const auto& indices = work->indices();
        std::copy(indices.begin(), indices.end(), ibuffer); // avoid lifetime issues with the workspace.
        return SparseRange(indices.size(), val, ibuffer); 
    }

    /**************************************
     ***** Sparse non-virtual methods *****
     **************************************/
private:
    static void copy_over(SparseRange<T, IDX>& output, T* vbuffer, IDX* ibuffer, SparseCopyMode copy) {
        if ((copy == SPARSE_COPY_BOTH || copy == SPARSE_COPY_INDEX) && output.index != ibuffer) {
            copy_over(output.index, ibuffer, output.number);
            output.index = ibuffer;
        }

        if ((copy == SPARSE_COPY_BOTH || copy == SPARSE_COPY_VALUE) && output.value != vbuffer) {
            copy_over(output.value, vbuffer, output.number);
            output.value = vbuffer;
        }
    }
    
public:
    /**
     * @param r Index of the row.
     * @param vbuffer Pointer to an array with enough space for at least `ncol()` values.
     * @param ibuffer Pointer to an array with enough space for at least `ncol()` indices.
     * @param work Pointer to a workspace created with `new_row_workspace()`.
     * @param copy Whether the non-zero values and/or indices should be copied into `vbuffer` and `ibuffer`, respectively.
     * @param sorted Should the non-zero elements be sorted by their indices? See `sparse_row()` for details.
     *
     * @return A `SparseRange` object describing the contents of row `r`.
     * Depending on `copy`, values and indices will be copied into `vbuffer` and/or `ibuffer`.
     */
    SparseRange<T, IDX> sparse_row_copy(size_t r, T* vbuffer, IDX* ibuffer, RowWorkspace* work, SparseCopyMode copy = SPARSE_COPY_BOTH, bool sorted = true) const {
        auto output = sparse_row(r, vbuffer, ibuffer, work, sorted);
        copy_over(output, vbuffer, ibuffer, copy);
        return output;
    }

    /**
     * @param c Index of the column.
     * @param vbuffer Pointer to an array with enough space for at least `nrow()` values.
     * @param ibuffer Pointer to an array with enough space for at least `nrow()` indices.
     * @param work Pointer to a workspace created with `new_column_workspace()`.
     * @param copy Whether the non-zero values and/or indices should be copied into `vbuffer` and `ibuffer`, respectively.
     * @param sorted Should the non-zero elements be sorted by their indices? See `sparse_column()` for details.
     *
     * @return A `SparseRange` object describing the contents of column `c`.
     * Depending on `copy`, values and indices will be copied into `vbuffer` and/or `ibuffer`.
     */
    SparseRange<T, IDX> sparse_column_copy(size_t c, T* vbuffer, IDX* ibuffer, ColumnWorkspace* work, SparseCopyMode copy = SPARSE_COPY_BOTH, bool sorted = true) const {
        auto output = sparse_column(c, vbuffer, ibuffer, work, sorted);
        copy_over(output, vbuffer, ibuffer, copy);
        return output;
    }

    /**
     * @param r Index of the row.
     * @param vbuffer Pointer to an array with enough space for at least `RowBlockWorkspace::length` values.
     * @param ibuffer Pointer to an array with enough space for at least `RowBlockWorkspace::length` indices.
     * @param work Pointer to a workspace created with `new_row_workspace()`.
     * @param copy Whether the non-zero values and/or indices should be copied into `vbuffer` and `ibuffer`, respectively.
     * @param sorted Should the non-zero elements be sorted by their indices? See `sparse_row()` for details.
     *
     * @return A `SparseRange` object describing the contents of a block of columns in row `r`.
     * Depending on `copy`, values and indices will be copied into `vbuffer` and/or `ibuffer`.
     */
    SparseRange<T, IDX> sparse_row_copy(size_t r, T* vbuffer, IDX* ibuffer, RowBlockWorkspace* work, SparseCopyMode copy = SPARSE_COPY_BOTH, bool sorted = true) const {
        auto output = sparse_row(r, vbuffer, ibuffer, work, sorted);
        copy_over(output, vbuffer, ibuffer, copy);
        return output;
    }

    /**
     * @param c Index of the column.
     * @param vbuffer Pointer to an array with enough space for at least `ColumnBlockWorkspace::length` values.
     * @param ibuffer Pointer to an array with enough space for at least `ColumnBlockWorkspace::length` indices.
     * @param work Pointer to a workspace created with `new_row_workspace()`.
     * @param copy Whether the non-zero values and/or indices should be copied into `vbuffer` and `ibuffer`, respectively.
     * @param sorted Should the non-zero elements be sorted by their indices? See `sparse_column()` for details.
     *
     * @return A `SparseRange` object describing the contents of a block of rows in column `c`.
     * Depending on `copy`, values and indices will be copied into `vbuffer` and/or `ibuffer`.
     */
    SparseRange<T, IDX> sparse_column_copy(size_t c, T* vbuffer, IDX* ibuffer, ColumnBlockWorkspace* work, SparseCopyMode copy = SPARSE_COPY_BOTH, bool sorted = true) const {
        auto output = sparse_column(c, vbuffer, ibuffer, work, sorted);
        copy_over(output, vbuffer, ibuffer, copy);
        return output;
    }

    /**
     * @param r Index of the row.
     * @param vbuffer Pointer to an array with enough space for at least `RowIndexWorkspace::length` values.
     * @param ibuffer Pointer to an array with enough space for at least `RowIndexWorkspace::length` indices.
     * @param work Pointer to a workspace created with `new_row_workspace()`.
     * @param copy Whether the non-zero values and/or indices should be copied into `vbuffer` and `ibuffer`, respectively.
     * @param sorted Should the non-zero elements be sorted by their indices? See `sparse_row()` for details.
     *
     * @return A `SparseRange` object describing the contents of a subset of columns in row `r`.
     * Depending on `copy`, values and indices will be copied into `vbuffer` and/or `ibuffer`.
     */
    SparseRange<T, IDX> sparse_row_copy(size_t r, T* vbuffer, IDX* ibuffer, RowIndexWorkspace<IDX>* work, SparseCopyMode copy = SPARSE_COPY_BOTH, bool sorted = true) const {
        auto output = sparse_row(r, vbuffer, ibuffer, work, sorted);
        copy_over(output, vbuffer, ibuffer, copy);
        return output;
    }

    /**
     * @param c Index of the column.
     * @param vbuffer Pointer to an array with enough space for at least `ColumnIndexWorkspace::length` values.
     * @param ibuffer Pointer to an array with enough space for at least `ColumnIndexWorkspace::length` indices.
     * @param work Pointer to a workspace created with `new_row_workspace()`.
     * @param copy Whether the non-zero values and/or indices should be copied into `vbuffer` and `ibuffer`, respectively.
     * @param sorted Should the non-zero elements be sorted by their indices? See `sparse_column()` for details.
     *
     * @return A `SparseRange` object describing the contents of a subset of rows in column `c`.
     * Depending on `copy`, values and indices will be copied into `vbuffer` and/or `ibuffer`.
     */
    SparseRange<T, IDX> sparse_column_copy(size_t c, T* vbuffer, IDX* ibuffer, ColumnIndexWorkspace<IDX>* work, SparseCopyMode copy = SPARSE_COPY_BOTH, bool sorted = true) const {
        auto output = sparse_column(c, vbuffer, ibuffer, work, sorted);
        copy_over(output, vbuffer, ibuffer, copy);
        return output;
    }

public:
    /**
     * A more convenient but less efficient version of the `sparse_row()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param r Index of the row.
     * @param work Pointer to a workspace created with `new_row_workspace()`.
     * @param sorted Should the non-zero elements be sorted by their indices? See `sparse_row()` for details.
     *
     * @return A `SparseRangeCopy` object containing the contents of row `r`.
     */
    SparseRangeCopy<T, IDX> sparse_row(size_t r, RowWorkspace* work, bool sorted = true) const {
        SparseRangeCopy<T, IDX> output(ncol());
        auto ret = sparse_row_copy(r, output.value.data(), output.index.data(), work, SPARSE_COPY_BOTH, sorted);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `sparse_column()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param c Index of the column.
     * @param work Pointer to a workspace created with `new_column_workspace()`.
     * @param sorted Should the non-zero elements be sorted by their indices? See `sparse_column()` for details.
     *
     * @return A `SparseRangeCopy` object containing the contents of column `c`.
     */
    SparseRangeCopy<T, IDX> sparse_column(size_t c, ColumnWorkspace* work, bool sorted = true) const {
        SparseRangeCopy<T, IDX> output(nrow());
        auto ret = sparse_column_copy(c, output.value.data(), output.index.data(), work, SPARSE_COPY_BOTH, sorted);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `sparse_row()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param r Index of the row.
     * @param work Pointer to a workspace created with `new_row_workspace()`.
     * @param sorted Should the non-zero elements be sorted by their indices? See `sparse_row()` for details.
     *
     * @return A `SparseRangeCopy` object containing the contents of a block of columns in row `r`.
     */
    SparseRangeCopy<T, IDX> sparse_row(size_t r, RowBlockWorkspace* work, bool sorted = true) const {
        SparseRangeCopy<T, IDX> output(work->length());
        auto ret = sparse_row_copy(r, output.value.data(), output.index.data(), work, SPARSE_COPY_BOTH, sorted);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `sparse_column()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param c Index of the column.
     * @param work Pointer to a workspace created with `new_column_workspace()`.
     * @param sorted Should the non-zero elements be sorted by their indices? See `sparse_column()` for details.
     *
     * @return A `SparseRangeCopy` object containing the contents of a block of rows in column `c`.
     */
    SparseRangeCopy<T, IDX> sparse_column(size_t c, ColumnBlockWorkspace* work, bool sorted = true) const {
        SparseRangeCopy<T, IDX> output(work->length());
        auto ret = sparse_column_copy(c, output.value.data(), output.index.data(), work, SPARSE_COPY_BOTH, sorted);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `sparse_row()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param r Index of the row.
     * @param work Pointer to a workspace created with `new_row_workspace()`.
     * @param sorted Should the non-zero elements be sorted by their indices? See `sparse_row()` for details.
     *
     * @return A `SparseRangeCopy` object containing the contents of a subset of columns in row `r`.
     */
    SparseRangeCopy<T, IDX> sparse_row(size_t r, RowIndexWorkspace<IDX>* work, bool sorted = true) const {
        SparseRangeCopy<T, IDX> output(work->length());
        auto ret = sparse_row_copy(r, output.value.data(), output.index.data(), work, SPARSE_COPY_BOTH, sorted);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `sparse_column()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param c Index of the column.
     * @param work Pointer to a workspace created with `new_column_workspace()`.
     * @param sorted Should the non-zero elements be sorted by their indices? See `sparse_column()` for details.
     *
     * @return A `SparseRangeCopy` object containing the contents of a subset of rows in column `c`.
     */
    SparseRangeCopy<T, IDX> sparse_column(size_t c, ColumnIndexWorkspace<IDX>* work, bool sorted = true) const {
        SparseRangeCopy<T, IDX> output(work->length());
        auto ret = sparse_column_copy(c, output.value.data(), output.index.data(), work, SPARSE_COPY_BOTH, sorted);
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
