#ifndef TATAMI_WORKSPACE_H
#define TATAMI_WORKSPACE_H

#include <memory>

/**
 * @file Workspace.hpp
 *
 * @brief Base workspace classes.
 *
 * Workspaces allow `Matrix` instances to persist information across certain method calls within the same thread.
 * `Matrix` subclass implementations can then cache bits and pieces of information for repeated accesses.
 */

namespace tatami {

/**
 * @brief Base row workspace class.
 *
 * This provides a workspace for row-wise extraction of data from `Matrix` instances.
 * Instances of this class cannot be constructed; it is only provided as a base class for `Matrix::new_row_workspace()` methods to return appropriate objects.
 */
class RowWorkspace {
protected:
    /**
     * @cond
     */
    RowWorkspace() = default;
    virtual ~RowWorkspace() = default;

    // Defining the other constructors for rule of 5. The move constructors
    // don't get auto-made when a destructor is declared, and if I'm defining
    // them, I might as well define the copy constructors.  Technically I
    // suppose I don't need them because this class is just an interface, but
    // who knows if I'll add some move-able stuff here later.
    RowWorkspace(RowWorkspace&&) = default;
    RowWorkspace& operator=(RowWorkspace&&) = default;
    RowWorkspace(const RowWorkspace&) = default;
    RowWorkspace& operator=(const RowWorkspace&) = default;
    /**
     * @endcond
     */
};

/**
 * @brief Base column workspace class.
 *
 * This provides a workspace for column-wise extraction of data from `Matrix` instances.
 * Instances of this class cannot be constructed; it is only provided as a base class for `Matrix::new_column_workspace()` methods to return appropriate objects.
 */
class ColumnWorkspace {
protected:
    /**
     * @cond
     */
    ColumnWorkspace() = default;
    virtual ~Columnkspace() = default;
    ColumnWorkspace(ColumnWorkspace&&) = default;
    ColumnWorkspace& operator=(ColumnWorkspace&&) = default;
    ColumnWorkspace(const ColumnWorkspace&) = default;
    ColumnWorkspace& operator=(const ColumnWorkspace&) = default;
    /**
     * @endcond
     */
};

/**
 * @brief Base row-block workspace class.
 *
 * This provides a workspace for row-wise extraction of a block of data from `Matrix` instances,
 * where the block is defined as a constant and contiguous set of columns from `[start, end)`.
 * Instances of this class cannot be constructed; it is only provided as a base class for `Matrix::new_row_workspace()` methods to return appropriate objects.
 */
class RowBlockWorkspace {
protected:
    /**
     * @cond
     */
    RowBlockWorkspace() = default;
    virtual ~RowBlockWorkspace() = default;
    RowBlockWorkspace(RowBlockWorkspace&&) = default;
    RowBlockWorkspace& operator=(RowBlockWorkspace&&) = default;
    RowBlockWorkspace(const RowBlockWorkspace&) = default;
    RowBlockWorkspace& operator=(const RowBlockWorkspace&) = default;
    /**
     * @endcond
     */

public:
    /**
     * Index of the first column in the block.
     */
    size_t start;

    /**
     * Number of columns in the block.
     */
    size_t length;
};

/**
 * @brief Base column-block workspace class.
 *
 * This provides a workspace for column-wise extraction of a block of data from `Matrix` instances,
 * where the block is defined as a constant and contiguous set of rows.
 * Instances of this class cannot be constructed; it is only provided as a base class for `Matrix::new_column_workspace()` methods to return appropriate objects.
 */
class ColumnBlockWorkspace {
protected:
    /**
     * @cond
     */
    ColumnBlockWorkspace() = default;
    virtual ~ColumnBlockWorkspace() = default;
    ColumnBlockWorkspace(ColumnBlockWorkspace&&) = default;
    ColumnBlockWorkspace& operator=(ColumnBlockWorkspace&&) = default;
    ColumnBlockWorkspace(const ColumnBlockWorkspace&) = default;
    ColumnBlockWorkspace& operator=(const ColumnBlockWorkspace&) = default;
    /**
     * @endcond
     */

public:
    /**
     * Index of the first row in the block.
     */
    size_t start;

    /**
     * Number of rows in the block.
     */
    size_t length;
};

/**
 * @tparam IDX Type of the column indices.
 * @brief Base row-index workspace class.
 *
 * This provides a workspace for row-wise extraction of a subset of data from `Matrix` instances,
 * where the subset is defined from sorted and unique column indices.
 * Instances of this class cannot be constructed; it is only provided as a base class for `Matrix::new_row_workspace()` methods to return appropriate objects.
 */
template<typename IDX>
class RowIndexWorkspace {
protected:
    /**
     * @cond
     */
    RowIndexWorkspace() = default;
    virtual ~RowIndexWorkspace() = default;
    RowIndexWorkspace(RowIndexWorkspace&&) = default;
    RowIndexWorkspace& operator=(RowIndexWorkspace&&) = default;
    RowIndexWorkspace(const RowIndexWorkspace&) = default;
    RowIndexWorkspace& operator=(const RowIndexWorkspace&) = default;
    /**
     * @endcond
     */

public:
    /**
     * Number of columns in the subset.
     */
    size_t length;

    /**
     * Pointer to an array of length `length`, containing sorted and unique indices for columns in the subset.
     */
    const IDX* indices;
};

/**
 * @tparam IDX Type of the row indices.
 * @brief Base column-index workspace class.
 *
 * This provides a workspace for column-wise extraction of a subset of data from `Matrix` instances,
 * where the subset is defined from sorted and unique row indices.
 * Instances of this class cannot be constructed; it is only provided as a base class for `Matrix::new_column_workspace()` methods to return appropriate objects.
 */
template<typename IDX>
class ColumnIndexWorkspace {
protected:
    /**
     * @cond
     */
    ColumnIndexWorkspace() = default;
    virtual ~ColumnIndexWorkspace() = default;
    ColumnIndexWorkspace(ColumnIndexWorkspace&&) = default;
    ColumnIndexWorkspace& operator=(ColumnIndexWorkspace&&) = default;
    ColumnIndexWorkspace(const ColumnIndexWorkspace&) = default;
    ColumnIndexWorkspace& operator=(const ColumnIndexWorkspace&) = default;
    /**
     * @endcond
     */

public:
    /**
     * Number of rows in the subset.
     */
    size_t length;

    /**
     * Pointer to an array of length `length`, containing sorted and unique indices for rows in the subset.
     */
    const IDX* indices;
};



}

#endif
