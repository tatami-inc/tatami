#ifndef TATAMI_SPARSE_RANGE_H
#define TATAMI_SPARSE_RANGE_H

#include <cstddef>
#include <vector>

/**
 * @file SparseRange.hpp
 *
 * Defines the `SparseRange` class to hold information about extracted sparse elements.
 */

using std::size_t;

namespace tatami {

/**
 * @brief A range of a sparse vector.
 *
 * This class defines a range along a sparse vector.
 * It records the number of "non-zero" elements within this range, and contains pointers to their values and row/column indices.
 * This is most commonly returned by `Matrix::sparse_row()` and `Matrix::sparse_column()` methods to obtain the contents of a row or column of a sparse matrix.
 *
 * @note
 * Note that the elements in `value` are not guaranteed to be non-zero.
 * Zero values are usually not explicitly filtered out.
 *
 * @tparam T Type of value.
 * @tparam IDX Type of index.
 */
template <typename T, typename IDX>
struct SparseRange {
    /**
     * @param n Number of non-zero values.
     * @param v Pointer to the values. This should have at least `n` addressible elements.
     * @param i Pointer to the indices. This should have at least `n` addressible elements.
     */ 
    SparseRange(size_t n, const T* v=NULL, const IDX* i=NULL) : number(n), value(v), index(i) {}

    /**
     * Default constructor.
     */
    SparseRange() {}

    /**
     * Number of non-zero elements.
     */
    size_t number = 0;

    /**
     * Pointer to the values of the non-zero elements.
     * Has at least `number` addressible entries. 
     */
    const T* value = NULL;

    /**
     * Pointer to the (row/column) indices of the non-zero elements.
     * Has at least `number` addressible entries. 
     */
    const IDX* index = NULL;
};

/**
 * @brief A range of a sparse vector with copying 
 *
 * This class defines a range along a sparse vector, where the values and indices of the non-zero elements are copied into internal `vector`s.
 * It provides more memory safety than the `SparseRange` class, at the cost of extra allocation and copying.
 *
 * @tparam T Type of value.
 * @tparam IDX Type of index.
 */
template<typename T, typename IDX>
struct SparseRangeCopy {
    /** 
     * Values of the non-zero elements.
     */
    std::vector<T> value;

    /** 
     * Indices of the non-zero elements.
     * This should be of the same length as `value`.
     */
    std::vector<IDX> index;

    /**
     * @param n Number of non-zero elements.
     */
    SparseRangeCopy(size_t n) : index(n), value(n) {}
};

/**
 * What components of the non-zero elements should be copied?
 * Just the indices (`INDEX`), the values (`VALUE`) or both (`BOTH`).
 */
enum SparseCopyMode { SPARSE_COPY_INDEX, SPARSE_COPY_VALUE, SPARSE_COPY_BOTH };

}

#endif
