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
 * @brief Base workspace class.
 * @tparam ROW Whether to perform row-wise extraction.
 *
 * This provides a workspace for extraction of entire vectors from `Matrix` instances, either for each row (if `ROW = true`) or column (otherwise).
 * Instances of this class cannot be constructed; it is only provided as a base class for `Matrix::new_row_workspace()` and `Matrix::new_column_workspace()` methods to return appropriate objects.
 */
template<bool ROW>
class Workspace {
protected:
    /**
     * @cond
     */
    virtual ~Workspace() = default;

    // Defining the other constructors for rule of 5. The move constructors
    // don't get auto-made when a destructor is declared, and if I'm defining
    // them, I might as well define the copy constructors.  Technically I
    // suppose I don't need them because this class is just an interface, but
    // who knows if I'll add some move-able stuff here later.
    Workspace(Workspace&&) = default;
    Workspace& operator=(Workspace&&) = default;
    Workspace(const Workspace&) = default;
    Workspace& operator=(const Workspace&) = default;
    /**
     * @endcond
     */

public:
    /**
     * Default constructor.
     */
    Workspace() = default;
};

/**
 * Base workspace class for extracting entire rows.
 */
typedef Workspace<true> RowWorkspace;

/**
 * Base workspace class for extracting entire columns.
 */
typedef Workspace<false> ColumnWorkspace;

/**
 * @brief Base workspace class for block extraction.
 * @tparam ROW Whether to perform row-wise extraction.
 *
 * This provides a workspace for extraction of a block of data from `Matrix` instances.
 * If `ROW = true`, each row-wise vector is extracted and the block is defined as a constant and contiguous set of columns from `[start, start + length)`.
 * Otherwise, each column-wise vector is extracted and the block is defined as a constant and contiguous set of rows from `[start, start + length)`.
 * Instances of this class cannot be constructed; it is only provided as a base class for `Matrix::new_row_workspace()` and `Matrix::new_column_workspace()` methods to return appropriate objects.
 */
template<bool ROW>
class BlockWorkspace {
protected:
    /**
     * @cond
     */
    virtual ~BlockWorkspace() = default;
    BlockWorkspace(BlockWorkspace&&) = default;
    BlockWorkspace& operator=(BlockWorkspace&&) = default;
    BlockWorkspace(const BlockWorkspace&) = default;
    BlockWorkspace& operator=(const BlockWorkspace&) = default;
    /**
     * @endcond
     */

public:
    /**
     * Default constructor.
     */
    BlockWorkspace() = default;

    /**
     * @param s Index of the first row/column in the block.
     * @param l Number of rows/columns in the block.
     */
    BlockWorkspace(size_t s, size_t l) : start(s), length(l) {}

    /**
     * Index of the first row/column in the block.
     */
    size_t start;

    /**
     * Number of rows/columns in the block.
     */
    size_t length;
};

/**
 * Base workspace class for extracting blocks of data from each row.
 */
typedef BlockWorkspace<true> RowBlockWorkspace;

/**
 * Base workspace class for extracting blocks of data from each column.
 */
typedef BlockWorkspace<false> ColumnBlockWorkspace;

/**
 * @brief Base workspace class for indexed extraction.
 *
 * @tparam IDX Type of the column indices.
 * @tparam ROW Whether to perform row-wise extraction.
 *
 * This provides a workspace for extraction of a subset of data from `Matrix` instances.
 * If `ROW = true`, each row-wise vector is extracted and the subset is defined from sorted and unique column indices.
 * Otherwise, each column-wise vector is extracted and the subset is defined from sorted and unique row indices.
 * Instances of this class cannot be constructed; it is only provided as a base class for `Matrix::new_row_workspace()` and `Matrix::new_column_workspace()` methods to return appropriate objects.
 */
template<typename IDX, bool ROW>
class IndexWorkspace {
protected:
    /**
     * @cond
     */
    virtual ~IndexWorkspace() = default;
    IndexWorkspace(IndexWorkspace&&) = default;
    IndexWorkspace& operator=(IndexWorkspace&&) = default;
    IndexWorkspace(const IndexWorkspace&) = default;
    IndexWorkspace& operator=(const IndexWorkspace&) = default;
    /**
     * @endcond
     */

public:
    /**
     * Default constructor.
     */
    IndexWorkspace() = default;

    /**
     * @param l Number of rows/columns in the subset.
     * @param i Pointer to an array of sorted and unique row/column indices.
     */
    IndexWorkspace(size_t l, const IDX* i) : length(l), indices(i) {}

    /**
     * Number of row/columns in the subset.
     */
    size_t length;

    /**
     * Pointer to an array of length `length`, containing sorted and unique indices for rows/columns in the subset.
     * It is assumed that the lifetime of the array exceeds that of this object.
     */
    const IDX* indices;
};

/**
 * Base workspace class for extracting subsets of data from each row.
 */
template<typename IDX>
using RowIndexWorkspace = IndexWorkspace<IDX, true>;

/**
 * Base workspace class for extracting subsets of data from each column.
 */
template<typename IDX>
using ColumnIndexWorkspace = IndexWorkspace<IDX, false>;


}

#endif
