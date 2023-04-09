#ifndef TATAMI_WORKSPACE_H
#define TATAMI_WORKSPACE_H

#include <memory>
#include "SparseRange.hpp"

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
 *
 * @tparam ROW Whether to perform row-wise extraction.
 *
 * This provides a workspace for dense extraction of entire vectors from `Matrix` instances, either for each row (if `ROW = true`) or column (otherwise).
 * Instances of this class cannot be directly constructed, and should instead be created by calling the appropriate `Matrix` methods.
 * Developers should inherit from this class when creating workspaces for their own `Matrix` implementation.
 */
template<bool ROW, bool SPARSE>
class Workspace {
protected:
   /**
     * @cond
     */
    Workspace() = default;

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
     * Is this class used for sparse extraction?
     */
    static constexpr bool sparse = SPARSE;

    /**
     * Is this class used for row-wise extraction?
     */
    static constexpr bool row = ROW;
};

/**
 * @tparam ROW Whether to perform row-wise extraction.
 * Workspace for extracting rows/columns in dense form.
 */
template<bool ROW>
using DenseWorkspace = Workspace<ROW, false>;

/**
 * Workspace for extracting entire rows in dense form.
 */
typedef DenseWorkspace<true> DenseRowWorkspace;

/**
 * Workspace for extracting entire columns in dense form.
 */
typedef DenseWorkspace<false> DenseColumnWorkspace;

/**
 * @tparam ROW Whether to perform row-wise extraction.
 * Workspace for extracting rows/columns in sparse form.
 */
template<bool ROW>
using SparseWorkspace = Workspace<ROW, true>;

/**
 * Workspace for extracting entire rows in sparse form.
 */
typedef SparseWorkspace<true> SparseRowWorkspace;

/**
 * Workspace for extracting entire columns in sparse form.
 */
typedef SparseWorkspace<false> SparseColumnWorkspace;

/**
 * @brief Base workspace class for block extraction.
 * @tparam ROW Whether to perform row-wise extraction.
 * @tparam SPARSE Whether to perform sparse extraction.
 *
 * This provides a workspace for dense extraction of a block of data from `Matrix` instances.
 * If `ROW = true`, each row-wise vector is extracted and the block is defined as a constant and contiguous set of columns from `[start, start + length)`.
 * Otherwise, each column-wise vector is extracted and the block is defined as a constant and contiguous set of rows from `[start, start + length)`.
 *
 * Instances of this class cannot be directly constructed, and should instead be created by calling the appropriate `Matrix` methods.
 * Developers should inherit from this class when creating workspaces for their own `Matrix` implementation.
 */
template<bool ROW, bool SPARSE>
class BlockWorkspace {
protected:
    /**
     * @cond
     */
    BlockWorkspace(size_t s, size_t l) : start(s), length(l) {}

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
     * Is this class used for sparse extraction?
     */
    static constexpr bool sparse = SPARSE;

    /**
     * Is this class used for row-wise extraction?
     */
    static constexpr bool row = ROW;

public:
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
 * @tparam ROW Whether to perform row-wise extraction.
 * Workspace for extracting blocks from rows/columns in dense form.
 */
template<bool ROW>
using DenseBlockWorkspace = BlockWorkspace<ROW, false>;

/**
 * Workspace for extracting entire rows in dense form.
 */
typedef DenseBlockWorkspace<true> DenseRowBlockWorkspace;

/**
 * Workspace for extracting entire columns in dense form.
 */
typedef DenseBlockWorkspace<false> DenseColumnBlockWorkspace;

/**
 * @tparam ROW Whether to perform row-wise extraction.
 * Workspace for extracting blocks from rows/columns in sparse form.
 */
template<bool ROW>
using SparseBlockWorkspace = BlockWorkspace<ROW, true>;

/**
 * Workspace for extracting entire rows in sparse form.
 */
typedef SparseBlockWorkspace<true> SparseRowBlockWorkspace;

/**
 * Workspace for extracting entire columns in sparse form.
 */
typedef SparseBlockWorkspace<false> SparseColumnBlockWorkspace;

/**
 * @brief Base workspace class for indexed extraction.
 *
 * @tparam IDX Type of the row/column indices.
 * @tparam ROW Whether to perform row-wise extraction.
 * @tparam SPARSE Whether to perform sparse extraction.
 *
 * This provides a workspace for extraction of a subset of data from `Matrix` instances.
 * If `ROW = true`, each row-wise vector is extracted and the subset is defined from sorted and unique column indices.
 * Otherwise, each column-wise vector is extracted and the subset is defined from sorted and unique row indices.
 *
 * Developers of `Matrix` implementations should define workspaces that inherit from this class's derived subclasses for dense/sparse extraction.
 */
template<typename IDX, bool ROW, bool SPARSE>
class IndexWorkspace {
protected:
    /**
     * @cond
     */
    IndexWorkspace(size_t l) : length(l) {}

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
     * Is this class used for sparse extraction?
     */
    static constexpr bool sparse = SPARSE;

    /**
     * Is this class used for row-wise extraction?
     */
    static constexpr bool row = ROW;

public:
    /**
     * @return Vector containing sorted and unique indices for rows/columns in the subset.
     */
    virtual const std::vector<IDX>& indices() const = 0;

    /**
     * @param Number of rows/columns in the subset.
     */
    size_t length;
};

/**
 * @tparam IDX Type of the row/column indices.
 * @tparam ROW Whether to perform row-wise extraction.
 * Workspace for extracting indexed subsets from rows/columns in dense form.
 */
template<typename IDX, bool ROW>
using DenseIndexWorkspace = IndexWorkspace<IDX, ROW, false>;

/**
 * @tparam IDX Type of the row/column indices.
 * Workspace for extracting subsets of data from each row in dense form.
 */
template<typename IDX>
using DenseRowIndexWorkspace = DenseIndexWorkspace<IDX, true>;

/**
 * @tparam IDX Type of the row/column indices.
 * Workspace for extracting subsets of data from each column in dense form.
 */
template<typename IDX>
using DenseColumnIndexWorkspace = DenseIndexWorkspace<IDX, false>;

/**
 * @tparam IDX Type of the row/column indices.
 * @tparam ROW Whether to perform row-wise extraction.
 * Workspace for extracting indexed subsets from rows/columns in sparse form.
 */
template<typename IDX, bool ROW>
using SparseIndexWorkspace = IndexWorkspace<IDX, ROW, true>;

/**
 * @tparam IDX Type of the row/column indices.
 * Workspace for extracting subsets of data from each row in sparse form.
 */
template<typename IDX>
using SparseRowIndexWorkspace = SparseIndexWorkspace<IDX, true>;

/**
 * @tparam IDX Type of the row/column indices.
 * Workspace for extracting subsets of data from each column in sparse form.
 */
template<typename IDX>
using SparseColumnIndexWorkspace = SparseIndexWorkspace<IDX, false>;

}

#endif
