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
 * This only exists to define the various constructors for a virtual base class,
 * given that we need to define a virtual destructor.
 * Developers should inherit from the derived classes.
 */
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
};

/**
 * @brief Dense workspace class.
 * @tparam ROW Whether to perform row-wise extraction.
 *
 * This provides a workspace for dense extraction of entire vectors from `Matrix` instances, either for each row (if `ROW = true`) or column (otherwise).
 * Instances of this class cannot be directly constructed, and should instead be created by calling the appropriate `Matrix` methods.
 */
template<bool ROW>
class DenseWorkspace : public Workspace {
protected:
    /**
     * Default constructor.
     */
    DenseWorkspace() = default;
};

/**
 * Workspace for extracting entire rows in dense form.
 */
typedef DenseWorkspace<true> DenseRowWorkspace;

/**
 * Workspace for extracting entire columns in dense form.
 */
typedef DenseWorkspace<false> DenseColumnWorkspace;

}

namespace tatami {

/**
 * @brief Sparse workspace class.
 * @tparam ROW Whether to perform row-wise extraction.
 *
 * This provides a workspace for sparse extraction of entire vectors from `Matrix` instances, either for each row (if `ROW = true`) or column (otherwise).
 * Instances of this class cannot be directly constructed, and should instead be created by calling the appropriate `Matrix` methods.
 */
template<bool ROW>
class SparseWorkspace : public Workspace {
protected:
    /**
     * @param m Extraction mode.
     * @param s Whether to sort by indices.
     */
    SparseWorkspace(SparseExtractMode m, bool so) : mode(m), sorted(so) {}

public:
    /**
     * Extraction mode - whether to extract just the indices, values, or both.
     */
    SparseExtractMode mode;

    /**
     * Whether to sort sparse values by index in the output.
     */
    bool sorted;
};

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
 *
 * This provides a workspace for dense extraction of a block of data from `Matrix` instances.
 * If `ROW = true`, each row-wise vector is extracted and the block is defined as a constant and contiguous set of columns from `[start, start + length)`.
 * Otherwise, each column-wise vector is extracted and the block is defined as a constant and contiguous set of rows from `[start, start + length)`.
 *
 * Developers of `Matrix` implementations should define workspaces that inherit from this class's derived subclasses for dense/sparse extraction.
 */
template<bool ROW>
class BlockWorkspace : public Workspace {
protected:
    /**
     * @cond
     */
    BlockWorkspace(size_t s, size_t l) : start(s), length(l) {}
    /**
     * @endcond
     */

public:
    /**
     * Index of the orst row/column in the block.
     */
    size_t start;

    /**
     * Number of rows/columns in the block.
     */
    size_t length;
};

/**
 * @brief Dense workspace class for block extraction.
 * @tparam ROW Whether to perform row-wise extraction.
 *
 * This provides a workspace for dense extraction of a block of data from `Matrix` instances.
 * Instances of this class cannot be directly constructed, and should instead be created by calling the appropriate `Matrix` methods.
 */
template<bool ROW>
class DenseBlockWorkspace : public BlockWorkspace<ROW> {
protected:
    /**
     * @param s Start of the block.
     * @param l Length of the block.
     */
    DenseBlockWorkspace(size_t s, size_t l) : BlockWorkspace<ROW>(s, l) {}
};

/**
 * Workspace for extracting blocks of data from each row in dense form.
 */
typedef DenseBlockWorkspace<true> DenseRowBlockWorkspace;

/**
 * Workspace for extracting blocks of data from each column in dense form.
 */
typedef DenseBlockWorkspace<false> DenseColumnBlockWorkspace;

/**
 * @brief Sparse workspace class for block extraction.
 * @tparam ROW Whether to perform row-wise extraction.
 *
 * This provides a workspace for sparse extraction of a block of data from `Matrix` instances.
 * Instances of this class cannot be directly constructed, and should instead be created by calling the appropriate `Matrix` methods.
 */
template<bool ROW>
class SparseBlockWorkspace : public BlockWorkspace<ROW> {
protected:
    /**
     * @param s Start of the block.
     * @param l Length of the block.
     * @param m Extraction mode.
     * @param so Whether to sort by indices.
     */
    SparseBlockWorkspace(size_t s, size_t l, SparseExtractMode m, bool so) : BlockWorkspace<ROW>(s, l), mode(m), sorted(so) {}

public:
    /**
     * Extraction mode - whether to extract just the indices, values, or both.
     */
    SparseExtractMode mode;

    /**
     * Whether to sort sparse values by index in the output.
     */
    bool sorted;
};

/**
 * Workspace for extracting blocks of data from each row in sparse form.
 */
typedef SparseBlockWorkspace<true> SparseRowBlockWorkspace;

/**
 * Workspace for extracting blocks of data from each column in sparse form.
 */
typedef SparseBlockWorkspace<false> SparseColumnBlockWorkspace;

/**
 * @brief Base workspace class for indexed extraction.
 *
 * @tparam IDX Type of the column indices.
 * @tparam ROW Whether to perform row-wise extraction.
 *
 * This provides a workspace for extraction of a subset of data from `Matrix` instances.
 * If `ROW = true`, each row-wise vector is extracted and the subset is defined from sorted and unique column indices.
 * Otherwise, each column-wise vector is extracted and the subset is defined from sorted and unique row indices.
 *
 * Developers of `Matrix` implementations should define workspaces that inherit from this class's derived subclasses for dense/sparse extraction.
 */
template<typename IDX, bool ROW>
class IndexWorkspace : public Workspace<ROW> {
protected:
    /**
     * @cond
     */
    IndexWorkspace(size_t l) : length(l) {}
    /**
     * @endcond
     */

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
 * @brief Dense workspace class for indexed extraction.
 *
 * @tparam IDX Type of the column indices.
 * @tparam ROW Whether to perform row-wise extraction.
 *
 * This provides a workspace for dense extraction of a subset of data from `Matrix` instances.
 * Instances of this class cannot be directly constructed, and should instead be created by calling the appropriate `Matrix` methods.
 */
template<typename IDX, bool ROW>
class DenseIndexWorkspace : public IndexWorkspace<IDX, ROW> {
protected:
    /**
     * @param l Length of the subset.
     */
    DenseIndexWorkspace(size_t l) : IndexWorkspace<IDX, ROW>(l) {}
};

/**
 * Workspace for extracting subsets of data from each row in dense form.
 */
template<typename IDX>
using DenseRowIndexWorkspace = DenseIndexWorkspace<IDX, true>;

/**
 * Workspace for extracting subsets of data from each column in dense form.
 */
template<typename IDX>
using DenseColumnIndexWorkspace = DenseIndexWorkspace<IDX, false>;

/**
 * @brief Sparse workspace class for indexed extraction.
 *
 * @tparam IDX Type of the column indices.
 * @tparam ROW Whether to perform row-wise extraction.
 *
 * This provides a workspace for sparse extraction of a subset of data from `Matrix` instances.
 * Instances of this class cannot be directly constructed, and should instead be created by calling the appropriate `Matrix` methods.
 */
template<typename IDX, bool ROW>
class SparseIndexWorkspace : public IndexWorkspace<IDX, ROW> {
protected:
    /**
     * @param l Length of the subset.
     * @param m Extraction mode.
     * @param so Whether to sort by indices.
     */
    SparseIndexWorkspace(size_t l, SparseExtractMode m, bool so) : IndexWorkspace<IDX, ROW>(l), mode(m), sorted(so) {}

public:
    /**
     * Extraction mode - whether to extract just the indices, values, or both.
     */
    SparseExtractMode mode;

    /**
     * Whether to sort sparse values by index in the output.
     */
    bool sorted;
};

/**
 * Workspace for extracting subsets of data from each row in sparse form.
 */
template<typename IDX>
using SparseRowIndexWorkspace = SparseIndexWorkspace<IDX, true>;

/**
 * Workspace for extracting subsets of data from each column in sparse form.
 */
template<typename IDX>
using SparseColumnIndexWorkspace = SparseIndexWorkspace<IDX, false>;

}

#endif
