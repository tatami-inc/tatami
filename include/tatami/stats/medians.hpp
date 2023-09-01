#ifndef TATAMI_STATS_MEDIANS_HPP
#define TATAMI_STATS_MEDIANS_HPP

#include "../base/Matrix.hpp"
#include "./utils.hpp"

#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>

/**
 * @file medians.hpp
 *
 * @brief Compute row and column medians from a `tatami::Matrix`.
 */

namespace tatami {

namespace stats {

/**
 * @cond
 */
template<typename Output_ = double, typename Value_>
Output_ compute_median(Value_* buffer, size_t n) {
    if (n == 0) {
        return std::numeric_limits<Output_>::quiet_NaN();
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

template<bool row_, typename Value_, typename Index_, typename Output_>
void dimension_medians(const Matrix<Value_, Index_>* p, Output_* output, int threads) {
    auto dim = (row_ ? p->nrow() : p->ncol());
    auto otherdim = (row_ ? p->ncol() : p->nrow());

    if (p->sparse()) {
        Options opt;
        opt.sparse_extract_index = false;
        opt.sparse_ordered_index = false; // we'll be sorting by value anyway.

        parallelize([&](int, Index_ s, Index_ l) -> void {
            auto ext = consecutive_extractor<row_, true>(p, s, l, opt);

            std::vector<Value_> buffer(otherdim);
            auto vbuffer = buffer.data();
            for (Index_ i = s, e = s + l; i < e; ++i) {
                auto range = ext->fetch_copy(i, vbuffer, NULL);
                auto n = range.number;

                if (n == otherdim) {
                    output[i] = compute_median<Output_>(vbuffer, otherdim);
                } else if (n * 2 < otherdim) {
                    output[i] = 0; // zero is the median if there are too many zeroes.
                } else {
                    size_t halfway = otherdim / 2;
                    bool is_even = (otherdim % 2 == 0);

                    auto vend = vbuffer  + n;
                    std::sort(vbuffer, vend);
                    size_t zeropos = std::lower_bound(vbuffer, vend, 0) - vbuffer;
                    size_t nzero = otherdim - n;

                    if (!is_even) {
                        if (zeropos > halfway) {
                            output[i] = vbuffer[halfway];
                        } else if (halfway >= zeropos + nzero) {
                            output[i] = vbuffer[halfway - nzero];
                        } else {
                            output[i] = 0; // zero is the median.
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
                        output[i] = tmp / 2;
                    }
                }
            }
        }, dim, threads);

    } else {
        parallelize([&](int, Index_ s, Index_ l) -> void {
            std::vector<Value_> buffer(otherdim);
            auto ext = consecutive_extractor<row_, false>(p, s, l);
            for (Index_ i = s, e = s + l; i < e; ++i) {
                ext->fetch_copy(i, buffer.data());
                output[i] = compute_median<Output_>(buffer.data(), otherdim);
            }
        }, dim, threads);
    }
}
/**
 * @endcond
 */

}

/**
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output_ Type of the output.
 *
 * @param p Shared pointer to a `tatami::Matrix`.
 * @param[out] output Pointer to an array of length equal to the number of columns.
 * On output, this will contain the column medians.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Output_>
void column_medians(const Matrix<Value_, Index_>* p, Output_* output, int threads = 1) {
    stats::dimension_medians<false>(p, output, threads);
}

/**
 * @tparam Output_ Type of the output.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Shared pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A vector of length equal to the number of columns, containing the column medians.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> column_medians(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> output(p->ncol());
    column_medians(p, output.data(), threads);
    return output;
}

/**
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output_ Type of the output.
 *
 * @param p Shared pointer to a `tatami::Matrix`.
 * @param[out] output Pointer to an array of length equal to the number of rows.
 * On output, this will contain the row medians.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Output_>
void row_medians(const Matrix<Value_, Index_>* p, Output_* output, int threads = 1) {
    stats::dimension_medians<true>(p, output, threads);
}

/**
 * @tparam Output_ Type of the output.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Shared pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A vector of length equal to the number of rows, containing the row medians.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> row_medians(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> output(p->nrow());
    row_medians(p, output.data(), threads);
    return output;
}

}

#endif
