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
 * @cond
 */
// Forward declarations.
template<typename T, typename IDX>
class Matrix;

template<bool ROW, typename T, typename IDX>
std::shared_ptr<DenseWorkspace<ROW> > dense_workspace(const Matrix<T, IDX>* ptr);

template<bool ROW, typename T, typename IDX>
std::shared_ptr<DenseBlockWorkspace<ROW> > dense_workspace(const Matrix<T, IDX>*, size_t, size_t);

template<bool ROW, typename T, typename IDX>
std::shared_ptr<DenseIndexWorkspace<IDX, ROW> > dense_workspace(const Matrix<T, IDX>*, std::vector<IDX>);
/**
 * @endcond
 */

/**
 * @brief Options for workspace construction.
 * 
 * This allows users to pass in more options to methods like `Matrix::dense_row_workspace()`, 
 * in order to fine-tune how the extraction of data out of the `Matrix` is performed.
 */
struct WorkspaceOptions {
    /** 
     * Whether to extract the sparse values, indices, both or neither.
     * This can be used to avoid unnecessary computation and copying.
     * Only used in the sparse workspace methods.
     */
    SparseExtractMode mode = SparseExtractMode::BOTH;

    /**
     * Whether the sparse values should be sorted by increasing indices.
     * Setting this to `false` can reduce computational work in situations where the order of non-zero elements does not matter.
     * Only used in the sparse workspace methods.
     */
    bool sorted = true;

    /** 
     * Whether to cache information from each call to `Matrix::row()`/`Matrix::column()` with this workspace.
     * This enables faster iterations if the same row/column is extracted during multiple passes over the matrix.
     */
    bool cache = false;
};

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
     * @param options Optional parameters for workspace construction.
     * @return A shared pointer to a `DenseRowWorkspace` for row-wise data extraction.
     */
    virtual std::shared_ptr<DenseRowWorkspace> dense_row_workspace(const WorkspaceOptions& options) const = 0;

    /**
     * @param options Optional parameters for workspace construction.
     * @return A shared pointer to a `DenseColumnWorkspace` for column-wise data extraction.
     */
    virtual std::shared_ptr<DenseColumnWorkspace> dense_column_workspace(const WorkspaceOptions& options) const = 0;

    /**
     * @cond
     */
    template<bool ROW>
    struct DefaultSparseWorkspace : public SparseWorkspace<ROW> {
        DefaultSparseWorkspace(const WorkspaceOptions& options, const Matrix* self) : SparseWorkspace<ROW>(options.mode, options.sorted) {
            if (sparse_extract_value(options.mode)) {
                dwork = dense_workspace<ROW>(self);
            }
        }

        // Need this for the dense extraction.
        std::shared_ptr<DenseWorkspace<ROW> > dwork;
    };
    /**
     * @endcond
     */

    /**
     * @param options Optional parameters for workspace construction.
     * @return A shared pointer to a `SparseRowWorkspace` for row-wise data extraction.
     */
    virtual std::shared_ptr<SparseRowWorkspace> sparse_row_workspace(const WorkspaceOptions& options) const {
        return std::shared_ptr<SparseRowWorkspace>(new DefaultSparseWorkspace<true>(options, this));
    }

    /**
     * @param options Optional parameters for workspace construction.
     * @return A shared pointer to a `SparseColumnWorkspace` for column-wise data extraction.
     */
    virtual std::shared_ptr<SparseColumnWorkspace> sparse_column_workspace(const WorkspaceOptions& options) const {
        return std::shared_ptr<SparseColumnWorkspace>(new DefaultSparseWorkspace<false>(options, this));
    }

public:
    /**
     * @param options Optional parameters for workspace construction.
     * @return A shared pointer to a `DenseRowBlockWorkspace` for row-wise extraction of data from a contiguous block of columns from `[start, start + length)`. 
     */
    virtual std::shared_ptr<DenseRowBlockWorkspace> dense_row_workspace(size_t start, size_t length, const WorkspaceOptions& options) const = 0;

    /**
     * @param start Index of the first row in the block.
     * @param length Number of rows in the block.
     * @param options Optional parameters for workspace construction.
     *
     * @return A shared pointer to a `DenseRowBlockWorkspace` for column-wise extraction of data from a contiguous block of rows from `[start, start + length)`.
     */
    virtual std::shared_ptr<DenseColumnBlockWorkspace> dense_column_workspace(size_t start, size_t length, const WorkspaceOptions& options) const = 0;

    /**
     * @cond
     */
    template<bool ROW>
    struct DefaultSparseBlockWorkspace : public SparseBlockWorkspace<ROW> {
        DefaultSparseBlockWorkspace(size_t start, size_t length, const WorkspaceOptions& options, const Matrix* self) :
            SparseBlockWorkspace<ROW>(start, length, options.mode, options.sorted) 
        {
            if (sparse_extract_value(options.mode)) {
                dwork = dense_workspace<ROW>(self, start, length);
            }
        }

        // Need this for the dense extraction.
        std::shared_ptr<DenseBlockWorkspace<ROW> > dwork;
    };
    /**
     * @endcond
     */

    /**
     * @param start Index of the first column in the block.
     * @param length Number of columns in the block.
     * @param options Optional parameters for workspace construction.
     *
     * @return A shared pointer to a `SparseRowBlockWorkspace` for row-wise extraction of data from a contiguous block of columns from `[start, start + length)`. 
     */
    virtual std::shared_ptr<SparseRowBlockWorkspace> sparse_row_workspace(size_t start, size_t length, const WorkspaceOptions& options) const {
        return std::shared_ptr<SparseRowBlockWorkspace>(new DefaultSparseBlockWorkspace<true>(start, length, options, this));
    }

    /**
     * @param start Index of the first row in the block.
     * @param length Number of rows in the block.
     * @param options Optional parameters for workspace construction.
     *
     * @return A shared pointer to a `SparseColumnBlockWorkspace` for column-wise extraction of data from a contiguous block of rows from `[start, start + length)`.
     */
    virtual std::shared_ptr<SparseColumnBlockWorkspace> sparse_column_workspace(size_t start, size_t length, const WorkspaceOptions& options) const {
        return std::shared_ptr<SparseColumnBlockWorkspace>(new DefaultSparseBlockWorkspace<false>(start, length, options, this));
    }

public:
    // Note to self: we use a vector rather than taking a pointer to the
    // indices. This is because workspaces of wrapper matrices may need to
    // construct their own vectors, in which case the validity of the pointers
    // is ambiguous, e.g., is the saved pointer to a vector still valid after
    // the vector is moved to a new location? See also discussion at:
    // https://stackoverflow.com/questions/11021764/does-moving-a-vector-invalidate-iterators

    /**
     * @param indices Vector containing sorted and unique column indices.
     * @param options Optional parameters for workspace construction.
     *
     * @return A shared pointer to a `DenseRowIndexWorkspace` for row-wise extraction of data from a subset of columns defined by `indices`.
     */
    virtual std::shared_ptr<DenseRowIndexWorkspace<IDX> > dense_row_workspace(std::vector<IDX> indices, const WorkspaceOptions& options) const = 0;

    /**
     * @param indices Vector containing sorted and unique row indices.
     * @param options Optional parameters for workspace construction.
     *
     * @return A shared pointer to a `DenseColumnIndexWorkspace` for column-wise extraction of data from a subset of rows defined by `indices`
     */
    virtual std::shared_ptr<DenseColumnIndexWorkspace<IDX> > dense_column_workspace(std::vector<IDX> indices, const WorkspaceOptions& options) const = 0;

    /**
     * @cond
     */
    template<bool ROW>
    struct DefaultSparseIndexWorkspace : public SparseIndexWorkspace<IDX, ROW> {
        DefaultSparseIndexWorkspace(std::vector<IDX> i, const WorkspaceOptions& options, const Matrix* self) : 
            SparseIndexWorkspace<IDX, ROW>(i.size(), options.mode, options.sorted)
        {
            if (sparse_extract_value(options.mode)) {
                dwork = dense_workspace<ROW>(self, std::move(i));
            } else {
                indices_ = std::move(i);
            }
        }

        std::vector<IDX> indices_;

        // Need this for the dense extraction.
        std::shared_ptr<DenseIndexWorkspace<IDX, ROW> > dwork;

        const std::vector<IDX>& indices() const { 
            if (indices_.empty()) {
                return dwork->indices();
            } else {
                return indices_; 
            }
        };
    };
    /**
     * @endcond
     */

    /**
     * @param indices Vector containing sorted and unique column indices.
     * @param options Optional parameters for workspace construction.
     *
     * @return A shared pointer to a `SparseRowIndexWorkspace` for row-wise extraction of data from a subset of columns defined by `indices`.
     */
    virtual std::shared_ptr<SparseRowIndexWorkspace<IDX> > sparse_row_workspace(std::vector<IDX> indices, const WorkspaceOptions& options) const {
        return std::shared_ptr<SparseRowIndexWorkspace<IDX> >(new DefaultSparseIndexWorkspace<true>(std::move(indices), options, this));
    }

    /**
     * @param indices Vector containing sorted and unique row indices.
     * @param options Optional parameters for workspace construction.
     *
     * @return A shared pointer to a `SparseColumnIndexWorkspace` for column-wise extraction of data from a subset of rows defined by `indices`.
     */
    virtual std::shared_ptr<SparseColumnIndexWorkspace<IDX> > sparse_column_workspace(std::vector<IDX> indices, const WorkspaceOptions& options) const {
        return std::shared_ptr<SparseColumnIndexWorkspace<IDX> >(new DefaultSparseIndexWorkspace<false>(std::move(indices), options, this));
    }

    /*****************************************
     ***** Non-virtual workspace methods *****
     *****************************************/
public:
    /**
     * @return A shared pointer to a `DenseRowWorkspace` for row-wise data extraction, or a null pointer if no workspace is required.
     */
    std::shared_ptr<DenseRowWorkspace> dense_row_workspace() const {
        return dense_row_workspace(WorkspaceOptions());
    }

    /**
     * @return A shared pointer to a `DenseColumnWorkspace` for row-wise data extraction, or a null pointer if no workspace is required.
     */
    std::shared_ptr<DenseColumnWorkspace> dense_column_workspace() const {
        return dense_column_workspace(WorkspaceOptions()); 
    }

    /**
     * @return A shared pointer to a `SparseRowWorkspace` for row-wise data extraction, or a null pointer if no workspace is required.
     */
    std::shared_ptr<SparseRowWorkspace> sparse_row_workspace() const {
        return sparse_row_workspace(WorkspaceOptions()); 
    }

    /**
     * @return A shared pointer to a `SparseColumnWorkspace` for row-wise data extraction, or a null pointer if no workspace is required.
     */
    std::shared_ptr<SparseColumnWorkspace> sparse_column_workspace() const {
        return sparse_column_workspace(WorkspaceOptions());
    }

    /**
     * @param start Index of the first row in the block.
     * @param length Number of rows in the block.
     *
     * @return A shared pointer to a `DenseRowBlockWorkspace` for row-wise extraction of data from a contiguous block of columns from `[start, start + length)`. 
     */
    std::shared_ptr<DenseRowBlockWorkspace> dense_row_workspace(size_t start, size_t length) const {
        return dense_row_workspace(start, length, WorkspaceOptions());
    }

    /**
     * @param start Index of the first row in the block.
     * @param length Number of rows in the block.
     *
     * @return A shared pointer to a `DenseRowBlockWorkspace` for column-wise extraction of data from a contiguous block of rows from `[start, start + length)`.
     */
    std::shared_ptr<DenseColumnBlockWorkspace> dense_column_workspace(size_t start, size_t length) const {
        return dense_column_workspace(start, length, WorkspaceOptions()); 
    }

    /**
     * @param start Index of the first column in the block.
     * @param length Number of columns in the block.
     *
     * @return A shared pointer to a `SparseRowBlockWorkspace` for row-wise extraction of data from a contiguous block of columns from `[start, start + length)`. 
     */
    std::shared_ptr<SparseRowBlockWorkspace> sparse_row_workspace(size_t start, size_t length) const {
        return sparse_row_workspace(start, length, WorkspaceOptions());
    }

    /**
     * @param start Index of the first row in the block.
     * @param length Number of rows in the block.
     *
     * @return A shared pointer to a `SparseRowBlockWorkspace` for column-wise extraction of data from a contiguous block of rows from `[start, start + length)`.
     */
    std::shared_ptr<SparseColumnBlockWorkspace> sparse_column_workspace(size_t start, size_t length) const {
        return sparse_column_workspace(start, length, WorkspaceOptions());
    }

    /**
     * @param indices Vector containing sorted and unique column indices.
     *
     * @return A shared pointer to a `DenseRowIndexWorkspace` for row-wise extraction of data from a subset of columns defined by `indices`.
     */
    std::shared_ptr<DenseRowIndexWorkspace<IDX> > dense_row_workspace(std::vector<IDX> indices) const {
        return dense_row_workspace(std::move(indices), WorkspaceOptions());
    }

    /**
     * @param indices Vector containing sorted and unique row indices.
     *
     * @return A shared pointer to a `DenseColumnIndexWorkspace` for column-wise extraction of data from a subset of rows defined by `indices`
     */
    std::shared_ptr<DenseColumnIndexWorkspace<IDX> > dense_column_workspace(std::vector<IDX> indices) const {
        return dense_column_workspace(std::move(indices), WorkspaceOptions());
    }

    /**
     * @param indices Vector containing sorted and unique column indices.
     * @param options Optional parameters for workspace construction.
     *
     * @return A shared pointer to a `SparseRowIndexWorkspace` for row-wise extraction of data from a subset of columns defined by `indices`.
     */
    std::shared_ptr<SparseRowIndexWorkspace<IDX> > sparse_row_workspace(std::vector<IDX> indices) const {
        return sparse_row_workspace(std::move(indices), WorkspaceOptions());
    }

    /**
     * @param indices Vector containing sorted and unique row indices.
     * @param options Optional parameters for workspace construction.
     *
     * @return A shared pointer to a `SparseColumnIndexWorkspace` for column-wise extraction of data from a subset of rows defined by `indices`.
     */
    virtual std::shared_ptr<SparseColumnIndexWorkspace<IDX> > sparse_column_workspace(std::vector<IDX> indices) const {
        return sparse_column_workspace(std::move(indices), WorkspaceOptions());
    }

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
        copy_over(ptr, buffer, work->length);
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
        copy_over(ptr, buffer, work->length);
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
        copy_over(ptr, buffer, work->length);
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
        copy_over(ptr, buffer, work->length);
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
        std::vector<T> output(work->length);
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
        std::vector<T> output(work->length);
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
        std::vector<T> output(work->length);
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
        std::vector<T> output(work->length);
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
     * @param work Pointer to a workspace created with `row_workspace()`.
     *
     * @return A `SparseRange` object describing the contents of row `r`.
     * Either or both of `value` or `index` may be `NULL`, depending on `work->mode`.
     */
    virtual SparseRange<T, IDX> row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowWorkspace* work) const {
        const T* val = NULL;
        if (sparse_extract_value(work->mode)) {
            val = row(r, vbuffer, static_cast<DefaultSparseWorkspace<true>*>(work)->dwork.get());
        }
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
    virtual SparseRange<T, IDX> column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnWorkspace* work) const {
        const T* val = NULL;
        if (sparse_extract_value(work->mode)) {
            val = column(c, vbuffer, static_cast<DefaultSparseWorkspace<false>*>(work)->dwork.get());
        }
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
    virtual SparseRange<T, IDX> row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowBlockWorkspace* work) const {
        const T* val = NULL;
        if (sparse_extract_value(work->mode)) {
            val = row(r, vbuffer, static_cast<DefaultSparseBlockWorkspace<true>*>(work)->dwork.get());
        }
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
    virtual SparseRange<T, IDX> column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnBlockWorkspace* work) const {
        const T* val = NULL;
        if (sparse_extract_value(work->mode)) {
            val = column(c, vbuffer, static_cast<DefaultSparseBlockWorkspace<false>*>(work)->dwork.get());
        }
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
    virtual SparseRange<T, IDX> row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowIndexWorkspace<IDX>* work) const {
        const T* val = NULL;
        if (sparse_extract_value(work->mode)) {
            val = row(r, vbuffer, static_cast<DefaultSparseIndexWorkspace<true>*>(work)->dwork.get());
        }
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
    virtual SparseRange<T, IDX> column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnIndexWorkspace<IDX>* work) const {
        const T* val = NULL;
        if (sparse_extract_value(work->mode)) {
            val = column(c, vbuffer, static_cast<DefaultSparseIndexWorkspace<false>*>(work)->dwork.get());
        }
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
     * If non-`NULL`, they will be set to `vbuffer` and `ibuffer`, which will be filled with values and indices respectively.
     */
    SparseRange<T, IDX> row_copy(size_t r, T* vbuffer, IDX* ibuffer, SparseRowWorkspace* work) const {
        auto output = row(r, vbuffer, ibuffer, work);
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
     * If non-`NULL`, they will be set to `vbuffer` and `ibuffer`, which will be filled with values and indices respectively.
     */
    SparseRange<T, IDX> column_copy(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnWorkspace* work) const {
        auto output = column(c, vbuffer, ibuffer, work);
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
     * If non-`NULL`, they will be set to `vbuffer` and `ibuffer`, which will be filled with values and indices respectively.
     */
    SparseRange<T, IDX> row_copy(size_t r, T* vbuffer, IDX* ibuffer, SparseRowBlockWorkspace* work) const {
        auto output = row(r, vbuffer, ibuffer, work);
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
     * If non-`NULL`, they will be set to `vbuffer` and `ibuffer`, which will be filled with values and indices respectively.
     */
    SparseRange<T, IDX> column_copy(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnBlockWorkspace* work) const {
        auto output = column(c, vbuffer, ibuffer, work);
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
     * If non-`NULL`, they will be set to `vbuffer` and `ibuffer`, which will be filled with values and indices respectively.
     */
    SparseRange<T, IDX> row_copy(size_t r, T* vbuffer, IDX* ibuffer, SparseRowIndexWorkspace<IDX>* work) const {
        auto output = row(r, vbuffer, ibuffer, work);
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
     * If non-`NULL`, they will be set to `vbuffer` and `ibuffer`, which will be filled with values and indices respectively.
     */
    SparseRange<T, IDX> column_copy(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnIndexWorkspace<IDX>* work) const {
        auto output = column(c, vbuffer, ibuffer, work);
        copy_over(output, vbuffer, ibuffer, work->mode);
        return output;
    }

public:
    /**
     * A more convenient but less efficient version of the `row()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param r Index of the row.
     * @param work Pointer to a workspace created with `sparse_row_workspace()`.
     *
     * @return A `SparseRangeCopy` object containing the contents of row `r`.
     */
    SparseRangeCopy<T, IDX> row(size_t r, SparseRowWorkspace* work) const {
        SparseRangeCopy<T, IDX> output(ncol());
        auto ret = row_copy(r, output.value.data(), output.index.data(), work);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `column()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param c Index of the column.
     * @param work Pointer to a workspace created with `sparse_column_workspace()`.
     *
     * @return A `SparseRangeCopy` object containing the contents of column `c`.
     */
    SparseRangeCopy<T, IDX> column(size_t c, SparseColumnWorkspace* work) const {
        SparseRangeCopy<T, IDX> output(nrow());
        auto ret = column_copy(c, output.value.data(), output.index.data(), work);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `row()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param r Index of the row.
     * @param work Pointer to a workspace created with `sparse_row_workspace()`.
     *
     * @return A `SparseRangeCopy` object containing the contents of a block of columns in row `r`.
     */
    SparseRangeCopy<T, IDX> row(size_t r, SparseRowBlockWorkspace* work) const {
        SparseRangeCopy<T, IDX> output(work->length);
        auto ret = row_copy(r, output.value.data(), output.index.data(), work);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `column()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param c Index of the column.
     * @param work Pointer to a workspace created with `sparse_column_workspace()`.
     *
     * @return A `SparseRangeCopy` object containing the contents of a block of rows in column `c`.
     */
    SparseRangeCopy<T, IDX> column(size_t c, SparseColumnBlockWorkspace* work) const {
        SparseRangeCopy<T, IDX> output(work->length);
        auto ret = column_copy(c, output.value.data(), output.index.data(), work);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `row()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param r Index of the row.
     * @param work Pointer to a workspace created with `sparse_row_workspace()`.
     *
     * @return A `SparseRangeCopy` object containing the contents of a subset of columns in row `r`.
     */
    SparseRangeCopy<T, IDX> row(size_t r, SparseRowIndexWorkspace<IDX>* work) const {
        SparseRangeCopy<T, IDX> output(work->length);
        auto ret = row_copy(r, output.value.data(), output.index.data(), work);
        output.index.resize(ret.number);
        output.value.resize(ret.number);
        return output;
    }

    /**
     * A more convenient but less efficient version of the `column()` method.
     * Callers do not have to supply the buffers; instead a new allocation is performed every time.
     *
     * @param c Index of the column.
     * @param work Pointer to a workspace created with `sparse_column_workspace()`.
     *
     * @return A `SparseRangeCopy` object containing the contents of a subset of rows in column `c`.
     */
    SparseRangeCopy<T, IDX> column(size_t c, SparseColumnIndexWorkspace<IDX>* work) const {
        SparseRangeCopy<T, IDX> output(work->length);
        auto ret = column_copy(c, output.value.data(), output.index.data(), work);
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
 * Create a new workspace for sparse extraction of full rows or columns.
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
std::shared_ptr<SparseWorkspace<ROW> > sparse_workspace(const Matrix<T, IDX>* ptr) {
    if constexpr(ROW) {
        return ptr->sparse_row_workspace();
    } else {
        return ptr->sparse_column_workspace();
    }
}

/**
 * Create a new workspace for dense extraction of a contiguous block from each row or column. 
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
std::shared_ptr<DenseBlockWorkspace<ROW> > dense_workspace(const Matrix<T, IDX>* ptr, size_t start, size_t length) {
    if constexpr(ROW) {
        return ptr->dense_row_workspace(start, length);
    } else {
        return ptr->dense_column_workspace(start, length);
    }
}

/**
 * Create a new workspace for sparse extraction of a contiguous block from each row or column. 
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
std::shared_ptr<SparseBlockWorkspace<ROW> > sparse_workspace(const Matrix<T, IDX>* ptr, size_t start, size_t length) {
    if constexpr(ROW) {
        return ptr->sparse_row_workspace(start, length);
    } else {
        return ptr->sparse_column_workspace(start, length);
    }
}

/**
 * Create a new workspace for dense extraction of a subset of entries from each row or column. 
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
std::shared_ptr<DenseIndexWorkspace<IDX, ROW> > dense_workspace(const Matrix<T, IDX>* ptr, std::vector<IDX> indices) {
    if constexpr(ROW) {
        return ptr->dense_row_workspace(std::move(indices));
    } else {
        return ptr->dense_column_workspace(std::move(indices));
    }
}

/**
 * Create a new workspace for sparse extraction of a subset of entries from each row or column. 
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
std::shared_ptr<SparseIndexWorkspace<IDX, ROW> > sparse_workspace(const Matrix<T, IDX>* ptr, std::vector<IDX> indices) {
    if constexpr(ROW) {
        return ptr->sparse_row_workspace(std::move(indices));
    } else {
        return ptr->sparse_column_workspace(std::move(indices));
    }
}

}

#endif
