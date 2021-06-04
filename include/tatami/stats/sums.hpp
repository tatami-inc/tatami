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

template<typename T> 
struct StatsSumHelper {
    StatsSumHelper(size_t n) : store(n) {}

    static const bool sparse = true;
    static const bool runnable = true;

    void direct(size_t i, const T* ptr, size_t dim) {
        store[i] = std::accumulate(ptr, ptr + dim, static_cast<T>(0));
        return;
    }

    template<typename IDX>
    void direct(size_t i, const sparse_range<T, IDX>& range, size_t dim) {
        store[i] = std::accumulate(range.value, range.value + range.number, static_cast<T>(0));
        return;
    }

    void running(const T* ptr) {
        for (auto sIt = store.begin(); sIt != store.end(); ++sIt, ++ptr) {
            *sIt += *ptr;
        }
        return;
    }

    template<typename IDX>
    void running(const sparse_range<T, IDX>& range) {
        auto vptr = range.value;
        auto iptr = range.index;
        for (size_t j = 0; j < range.number; ++j, ++iptr, ++vptr) {
            store[*iptr] += *vptr;
        }
        return;
    }

    std::vector<T> yield() {
        return store;
    }
private:
    std::vector<T> store;
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
    StatsSumHelper<T> stats(p->ncol());
    apply<1, T, IDX>(p, stats);
    return stats.yield();
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
    StatsSumHelper<T> stats(p->nrow());
    apply<0, T, IDX>(p, stats);
    return stats.yield();
}

}

#endif
