#ifndef TATAMI_SPARSE_RANGE_H
#define TATAMI_SPARSE_RANGE_H

/**
 * @file workspace.hpp
 *
 * Defines the `sparse_range` class to hold information about extracted sparse values.
 */

namespace tatami {

/**
 * @brief A range of sparse values.
 *
 * More specifically, this class defines a range of sparse (i.e., "non-zero") values, e.g., along a slice of a row or column.
 * The aim is to hold (pointers to) the values and row/column indices, as well as the number of non-zero values within this range.
 * This is most commonly returned by `typed_matrix::get_sparse_row()` and `typed_matrix::get_sparse_column()` methods.
 *
 * @tparam T Type of value.
 * @tparam IDX Type of index.
 */
template <typename T, typename IDX>
struct sparse_range {
    /**
     * @param n Number of non-zero values.
     * @param v Pointer to the values. This should have at least `n` addressible elements.
     * @param i Pointer to the indices. This should have at least `n` addressible elements.
     */ 
    sparse_range(size_t n, const T* v=NULL, const IDX* i=NULL) : number(n), value(v), index(i) {}

    /**
     * Default constructor.
     */
    sparse_range() {}

    size_t number = 0;
    const T* value = NULL;
    const IDX* index = NULL;
};

}

#endif
