#ifndef TATAMI_SPARSE_RANGE_H
#define TATAMI_SPARSE_RANGE_H

#include <cstddef>
#include <vector>

/**
 * @file SparseRange.hpp
 *
 * @brief Store information about extracted sparse elements.
 */

namespace tatami {

/**
 * @brief A range of a sparse vector.
 *
 * This class defines a range along a sparse vector.
 * It records the number of structural "non-zero" elements within this range, and contains pointers to their values and row/column indices.
 *
 * @note
 * Note that the elements in `value` are not guaranteed to be non-zero.
 * If zeroes are explicitly initialized in the underlying structure, they will be reported here.
 * However, one can safely assume that all indices _not_ reported in `index` have values of zero.
 *
 * @tparam Value Data value type, should be numeric.
 * @tparam Index Row/column index type, should be integer.
 */
template <typename Value, typename Index>
struct SparseRange {
    /**
     * @param n Number of structural non-zero values.
     * @param v Pointer to the values. This should have at least `n` addressible elements.
     * @param i Pointer to the indices. This should have at least `n` addressible elements.
     */ 
    SparseRange(Index n, const Value* v=NULL, const Index* i=NULL) : number(n), value(v), index(i) {}

    /**
     * Default constructor.
     */
    SparseRange() {}

    /**
     * Number of structural non-zero elements.
     */
    Index number = 0;

    /**
     * Pointer to the values of the structural non-zeros.
     * Has at least `number` addressible entries. 
     */
    const Value* value = NULL;

    /**
     * Pointer to the (row/column) indices of the structural non-zeros.
     * Has at least `number` addressible entries. 
     */
    const Index* index = NULL;
};

/**
 * @brief A range of a sparse vector with copying 
 *
 * This class defines a range along a sparse vector, where the values and indices of the structural non-zero elements are copied into internal `vector`s.
 * It provides more memory safety than the `SparseRange` class, at the cost of extra allocation and copying.
 *
 * @tparam Value Data value type, should be numeric.
 * @tparam Index Row/column index type, should be integer.
 */
template<typename Value, typename Index>
struct SparseRangeCopy {
    /**
     * Number of structural non-zero elements.
     */
    Index number = 0;

    /** 
     * Values of the structural non-zeros.
     *
     * Note that this may not be filled if extraction was performed with a `SparseExtractMode` that did not extract the values.
     */
    std::vector<Value> value;

    /** 
     * Indices of the structural non-zeros.
     * This should be of the same length as `value`.
     * 
     * Note that this may not be filled if extraction was performed with a `SparseExtractMode` that did not extract the indices.
     */
    std::vector<Index> index;

    /**
     * @param n Number of structural non-zeros.
     */
    SparseRangeCopy(Index n) : number(n), index(n), value(n) {}
};

}

#endif
