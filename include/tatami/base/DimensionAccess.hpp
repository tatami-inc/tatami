#ifndef TATAMI_GRID_HPP
#define TATAMI_GRID_HPP

#include <type_traits>
#include "SparseRange.hpp"

/**
 * @file DimensionAccess.hpp
 *
 * @brief Virtual class for acessing matrix dimensions.
 */

namespace tatami {

/**
 * @tparam Index Row/column index type, should be integer.
 * @brief Virtual base class for consecutive iteration.
 */
template<typename Value, typename Index>
class SortedOrder {
protected:
    /**
     * @cond
     */
    virtual ~SortedOrder() {}
    /**
     * @endcond
     */

public:
    /**
     * Advance to the next element in the dimension.
     * For the last element, this should behave like `reset()`.
     */
    virtual void next() = 0;

    /**
     * Reset to the first element in the dimension.
     */
    virtual void reset() = 0; 
};

/**
 * @tparam Index Row/column index type, should be integer.
 * @brief Virtual base class for random access.
 */
template<typename Index>
class RandomOrder {
protected:
    /**
     * @cond
     */
    virtual ~RandomOrder() {}
    /**
     * @endcond
     */

public:
    /**
     * @param i Index of the dimension.
     * The position of this instance is set to `i`.
     */
    virtual void set(Index i) = 0;
};

/**
 * Access pattern for the dimension elements.
 * This can be sorted (i.e., non-decreasing) or random.
 */
enum class AccessOrder char { SORTED, RANDOM };

/**
 * @tparam order Access order.
 * @tparam Index Row/column index type, should be integer.
 *
 * Conditional class definition for different access patterns.
 */
template<AccessOrder order, typename Index>
using ConditionalOrder = typename std::conditional<order == AccessOrder::RANDOM, RandomOrder<Index>, SortedOrder<Index> >::type;

/**
 * @tparam Value Data value type, should be numeric.
 * @tparam Index Row/column index type, should be integer.
 * @brief Virtual base class for dense extraction.
 */
template<typename Value, typename Index>
class DenseFormat {
protected:
    /**
     * @cond
     */
    virtual ~DenseFormat() {}
    /**
     * @endcond
     */

public:
    /**
     * @return Minimum length of the buffer to use in `fetch()`.
     */
    virtual Index buffer_length() const = 0;

    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param buffer Pointer to an array of length specified by `buffer_length()`.
     *
     * @return Pointer to the values of the current element of the dimension, containing `length()` valid entries.
     */
    virtual const Value* fetch(Value* buffer) const = 0;
};

/**
 * @tparam Value Data value type, should be numeric.
 * @tparam Index Row/column index type, should be integer.
 * @brief Virtual base class for dense extraction.
 */
class SparseFormat {
protected:
    /**
     * @cond
     */
    virtual ~SparseFormat() {}
    /**
     * @endcond
     */

public:
    /**
     * @return Minimum length of the buffers to use in `fetch()`.
     */
    virtual Index buffer_length() const = 0;

    /**
     * `vbuffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `vbuffer`; if they are the same, `vbuffer` has been filled. 
     * The same applies for `ibuffer` and the returned `SparseRange::index` pointer.
     *
     * Entries in the output `value` array are not guaranteed to be non-zero.
     * If zeroes are explicitly initialized in the underlying representation, they will be reported here.
     * However, one can safely assume that all values _not_ in `value` are zero.
     *
     * @param vbuffer Pointer to an array with enough space for at least `buffer_length()` values.
     * Ignored if `WorkspaceOptions::sparse_extract_value` was set to `false`.
     * @param ibuffer Pointer to an array with enough space for at least `buffer_length()` indices.
     * Ignored if `WorkspaceOptions::sparse_extract_index` was set to `false`.
     *
     * @return A `SparseRange` object describing the contents of the current dimension element.
     * Either or both of `value` or `index` is set to `NULL` if extraction of that field is skipped, based on the setting of `WorkspaceOptions::sparse_extract_mode` used to construct this object.
     */
    virtual SparseRange<T, IDX> fetch(Value* vbuffer, Index* ibuffer) const = 0;
};

/**
 * @tparam sparse Whether to perform sparse retrieval.
 * @tparam Value Data value type, should be numeric.
 * @tparam Index Row/column index type, should be integer.
 *
 * Conditional class definition for dense/sparse formats.
 */
template<bool sparse, typename Value, typename Index>
using ConditionalFormat = typename std::conditional<sparse, SparseFormat<Value, Index>, DenseFormat<Value, Index> >::type;

/**
 * @tparam order Access order.
 * @tparam sparse Whether to perform sparse retrieval.
 * @tparam Value Data value type, should be numeric.
 * @tparam Index Row/column index type, should be integer.
 * 
 * @brief Virtual base class for accessing dimension elements.
 */
template<AccessOrder order, bool sparse, typename Value, typename Index>
struct DimensionAccess: public ConditionalOrder<order, Index>, public ConditionalFormat<sparse, Value, Index> {};

/**
 * @tparam Value Data value type, should be numeric.
 * @tparam Index Row/column index type, should be integer.
 *
 * Access random elements in dense format.
 */
template<typename Value, typename Index>
typedef DimensionAccess<AccessOrder::RANDOM, false, Value, Index> DenseRandomDimensionAccess;

/**
 * @tparam Value Data value type, should be numeric.
 * @tparam Index Row/column index type, should be integer.
 *
 * Access random elements in sparse format.
 */
template<typename Value, typename Index>
typedef DimensionAccess<AccessOrder::RANDOM, true, Value, Index> SparseRandomDimensionAccess;

/**
 * @tparam Value Data value type, should be numeric.
 * @tparam Index Row/column index type, should be integer.
 *
 * Access ordered elements in dense format.
 */
template<typename Value, typename Index>
typedef DimensionAccess<AccessOrder::ORDERED, false, Value, Index> DenseOrderedDimensionAccess;

/**
 * @tparam Value Data value type, should be numeric.
 * @tparam Index Row/column index type, should be integer.
 *
 * Access ordered elements in sparse format.
 */
template<typename Value, typename Index>
typedef DimensionAccess<AccessOrder::ORDERED, true, Value, Index> SparseOrderedDimensionAccess;

}

#endif
