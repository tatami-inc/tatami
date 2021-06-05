#ifndef TATAMI_STATS_MEDIANS_HPP
#define TATAMI_STATS_MEDIANS_HPP

#include "../base/typed_matrix.hpp"
#include "apply.hpp"
#include <vector>
#include <algorithm>
#include <limits>

/**
 * @file medians.hpp
 *
 * Compute row and column medians from a `tatami::typed_matrix`.
 */

namespace tatami {

template<typename T, bool SPARSE, bool RUNNABLE> 
struct StatsMedianHelper {
    StatsMedianHelper(size_t n, size_t d) : store(n), dim(d), halfway(dim/2), is_even(dim%2==0) {}

    static const bool sparse = SPARSE;
    static const bool runnable = false;
    typedef std::vector<double> value;

    void direct(size_t i, const T* ptr, T* buffer) {
        if (dim == 0) {
            store[i] = std::numeric_limits<double>::quiet_NaN();
            return;
        }
        if (ptr != buffer) {
            std::copy(ptr, ptr + dim, buffer);
        }

        // At some point, I found two nth_element calls to be faster than partial_sort.
        std::nth_element(buffer, buffer + halfway, buffer + dim);
        double medtmp = *(buffer + halfway);
        if (is_even) {
            std::nth_element(buffer, buffer + halfway - 1, buffer + dim);
            store[i] = (medtmp + *(buffer + halfway - 1))/2;
        } else {
            store[i] = medtmp;
        }
        return;
    }

    template<typename IDX>
    void direct(size_t i, const sparse_range<T, IDX>& range, T* vbuffer, IDX* ibuffer) {
        if (range.number == dim) {
            direct(i, range.value, vbuffer);
        } else if (range.number * 2 < dim) {
            // zero is the median.
        } else {
            if (range.value != vbuffer) {
                std::copy(range.value, range.value + range.number, vbuffer);
            }

            auto vend = vbuffer + range.number;
            std::sort(vbuffer, vend);
            size_t zeropos = std::lower_bound(vbuffer, vend, 0) - vbuffer;
            size_t nzero = dim - range.number;

            if (!is_even) {
                if (zeropos > halfway) {
                    store[i] = vbuffer[halfway];
                } else if (halfway >= zeropos + nzero) {
                    store[i] = vbuffer[halfway - nzero];
                } else {
                    ; // zero is the median.
                }
            } else {
                double& tmp = store[i];
                if (zeropos > halfway) {
                    tmp = vbuffer[halfway] + vbuffer[halfway - 1];
                } else if (zeropos == halfway) {
                    // guaranteed to be at least 1 zero.
                    tmp += vbuffer[halfway - 1];
                } else if (zeropos < halfway && zeropos + nzero > halfway) {
                    ; // zero is the median.
                } else if (zeropos + nzero == halfway) {
                    // guaranteed to be at least 1 zero.
                    tmp += vbuffer[halfway - nzero];
                } else {
                    tmp = vbuffer[halfway - nzero] + vbuffer[halfway - nzero - 1];
                }
                tmp /= 2;
            }
        }
        return;
    }

    value yield() {
        return store;
    }
private:
    value store;
    size_t dim, halfway;
    bool is_even;
};

/**
 * @tparam T Type of the matrix value.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Shared pointer to a `tatami::typed_matrix`.
 *
 * @return A vector of length equal to the number of columns, containing the column medians.
 */
template<typename T, typename IDX>
inline std::vector<T> column_medians(const typed_matrix<T, IDX>* p) {
    return apply<1, T, IDX, StatsMedianHelper>(p);
}

/**
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Shared pointer to a `tatami::typed_matrix`.
 *
 * @return A vector of length equal to the number of rows, containing the row medians.
 */
template<typename T, typename IDX>
inline std::vector<T> row_medians(const typed_matrix<T, IDX>* p) {
    return apply<0, T, IDX, StatsMedianHelper>(p);
}

}

#endif
