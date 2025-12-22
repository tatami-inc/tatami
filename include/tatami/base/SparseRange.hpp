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
 * It records the number of structural "non-zero" elements within this range, and contains pointers to their values and indices.
 *
 * When a `SparseRange` is created by methods like `MyopicSparseExtractor::fetch()`, each of its indices refers to a position on the non-target dimension of the original `Matrix`.
 * For example, when iterating over the rows, the indices would correspond to the columns.
 *
 * Note that the elements in `value` are not guaranteed to be non-zero.
 * If zeroes are explicitly initialized in the underlying structure, they will be reported here.
 * However, one can safely assume that all indices _not_ reported in `index` have values of zero.
 *
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 */
template <typename Value_, typename Index_>
struct SparseRange {
    /**
     * @param number Number of structural non-zero values.
     * @param value Pointer to the values.
     * This should have at least `number` addressible elements.
     * @param index Pointer to the indices.
     * This should have at least `number` addressible elements.
     */ 
    SparseRange(Index_ number, const Value_* value = NULL, const Index_* index = NULL) : number(number), value(value), index(index) {}

    /**
     * Default constructor.
     */
    SparseRange() = default;

    /**
     * Number of structural non-zero elements.
     */
    Index_ number = 0;

    /**
     * Pointer to an array containing the values of the structural non-zeros.
     * If non-`NULL`, this has at least `number` addressible entries. 
     */
    const Value_* value = NULL;

    /**
     * Pointer to an array containing the indices of the structural non-zeros.
     * If non-`NULL`, this has at least `number` addressible entries. 
     */
    const Index_* index = NULL;
};

}

#endif
