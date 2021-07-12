#ifndef TATAMI_STATS_MEDIANS_HPP
#define TATAMI_STATS_MEDIANS_HPP

#include "../base/Matrix.hpp"
#include "apply.hpp"

#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

/**
 * @file medians.hpp
 *
 * Compute row and column medians from a `tatami::Matrix`.
 */

namespace tatami {

namespace stats {

/**
 * @brief Helper to compute the median along each dimension.
 */
struct MedianHelper {
    /**
     * This statistic can be computed from sparse inputs.
     */
    static const bool supports_sparse = true;

    /**
     * This statistic cannot be computed in a running manner.
     */
    static const bool supports_running = false;

public:
    /**
     * Compute the median along a vector.
     *
     * @tparam T Type of the input data.
     *
     * @param ptr Pointer to an array of values of length `n`.
     * @param n Size of the array.
     * @param buffer Pointer to an array with `n` addressible elements, to be used as a workspace.
     *
     * @return The median of values in `[ptr, ptr + n)`.
     */
    template<typename T = double>
    static double compute(const T* ptr, size_t n, T* buffer) {
        if (n == 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (ptr != buffer) {
            std::copy(ptr, ptr + n, buffer);
        }

        size_t halfway = n / 2;
        bool is_even = (n % 2 == 0);

        // At some point, I found two nth_element calls to be faster than partial_sort.
        std::nth_element(buffer, buffer + halfway, buffer + n);
        double medtmp = *(buffer + halfway);
        if (is_even) {
            std::nth_element(buffer, buffer + halfway - 1, buffer + n);
            return (medtmp + *(buffer + halfway - 1))/2;
        } else {
            return medtmp;
        }
    }

    /**
     * Compute the median along a sparse vector.
     * This achieves faster processing by only sorting the non-zero elements.
     *
     * @tparam T Type of the input data.
     * @tparam IDX Type of the indices.
     *
     * @param range A `SparseRange` object specifying the number and values of all non-zero indices.
     * @param n Total length of the vector, including zero values.
     * @param vbuffer,ibuffer Pointer to arrays with `range.number` addressible elements, to be used as a workspace.
     *
     * @return The median of values in the vector.
     */
    template<typename T = double, typename IDX = int>
    static double compute(const SparseRange<T, IDX>& range, size_t n, T* vbuffer, IDX* ibuffer) {
        if (range.number == n) {
            return compute(range.value, n, vbuffer);
        } else if (range.number * 2 < n) {
            return 0; // zero is the median if there are too many zeroes.
        } else {
            if (range.value != vbuffer) {
                std::copy(range.value, range.value + range.number, vbuffer);
            }

            size_t halfway = n / 2;
            bool is_even = (n % 2 == 0);

            auto vend = vbuffer + range.number;
            std::sort(vbuffer, vend);
            size_t zeropos = std::lower_bound(vbuffer, vend, 0) - vbuffer;
            size_t nzero = n - range.number;

            if (!is_even) {
                if (zeropos > halfway) {
                    return vbuffer[halfway];
                } else if (halfway >= zeropos + nzero) {
                    return vbuffer[halfway - nzero];
                } else {
                    return 0; // zero is the median.
                }
            } else {
                double tmp = 0;
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
                return tmp / 2;
            }
        }
    }
};

}

/**
 * @tparam T Type of the matrix value.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Shared pointer to a `tatami::Matrix`.
 *
 * @return A vector of length equal to the number of columns, containing the column medians.
 */
template<typename T, typename IDX>
inline std::vector<T> column_medians(const Matrix<T, IDX>* p) {
    return apply<1, T, IDX, stats::MedianHelper>(p);
}

/**
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Shared pointer to a `tatami::Matrix`.
 *
 * @return A vector of length equal to the number of rows, containing the row medians.
 */
template<typename T, typename IDX>
inline std::vector<T> row_medians(const Matrix<T, IDX>* p) {
    return apply<0, T, IDX, stats::MedianHelper>(p);
}

}

#endif
