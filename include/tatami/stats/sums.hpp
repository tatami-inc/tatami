#ifndef TATAMI_STATS_SUMS_HPP
#define TATAMI_STATS_SUMS_HPP

#include "../base/typed_matrix.hpp"
#include "apply.hpp"
#include <vector>
#include <algorithm>

/**
 * @file sums.hpp
 *
 * Compute row and column sums from a `tatami::typed_matrix`.
 */

namespace tatami {

template<typename T, bool SPARSE, bool RUNNABLE> 
struct StatsSumHelper {
    StatsSumHelper(size_t n, size_t d) : store(n), dim(d) {}

    static const bool sparse = SPARSE;
    static const bool runnable = RUNNABLE;
    typedef std::vector<T> value;

    void direct(size_t i, const T* ptr) {
        static_assert(!SPARSE && !RUNNABLE);
        store[i] = std::accumulate(ptr, ptr + dim, static_cast<T>(0));
        return;
    }

    template<typename IDX>
    void direct(size_t i, const sparse_range<T, IDX>& range) {
        static_assert(SPARSE && !RUNNABLE);
        store[i] = std::accumulate(range.value, range.value + range.number, static_cast<T>(0));
        return;
    }

    void running(const T* ptr) {
        static_assert(!SPARSE && RUNNABLE);
        for (auto sIt = store.begin(); sIt != store.end(); ++sIt, ++ptr) {
            *sIt += *ptr;
        }
        return;
    }

    template<typename IDX>
    void running(sparse_range<T, IDX> range) {
        static_assert(SPARSE && RUNNABLE);
        for (size_t j = 0; j < range.number; ++j, ++range.index, ++range.value) {
            store[*range.index] += *range.value;
        }
        return;
    }

    value yield() {
        return store;
    }
private:
    value store;
    size_t dim;
};

/**
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::typed_matrix`.
 *
 * @return A vector of length equal to the number of columns, containing the column sums.
 */
template<typename T, typename IDX>
inline std::vector<T> column_sums(const typed_matrix<T, IDX>* p) {
    return apply<1, T, IDX, StatsSumHelper>(p);
}

/**
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::typed_matrix`.
 *
 * @return A vector of length equal to the number of rows, containing the row sums.
 */
template<typename T, typename IDX>
inline std::vector<T> row_sums(const typed_matrix<T, IDX>* p) {
    return apply<0, T, IDX, StatsSumHelper>(p);
}

}

#endif
