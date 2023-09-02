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
 * Compute the median from a dense vector.
 *
 * @param[in] ptr Pointer to an array of values.
 * This may be modified on output.
 * @param n Length of the array.
 *
 * @tparam Output_ Type of the output value.
 * This should be a floating-point value for possible averaging.
 * @tparam Value_ Type of the input values.
 *
 * @return The median of values in `[ptr, ptr + n)`.
 */
template<typename Output_ = double, typename Value_>
Output_ compute_median(Value_* ptr, size_t n) {
    if (n == 0) {
        return std::numeric_limits<Output_>::quiet_NaN();
    }

    size_t halfway = n / 2;
    bool is_even = (n % 2 == 0);

    // At some point, I found two nth_element calls to be faster than partial_sort.
    std::nth_element(ptr, ptr + halfway, ptr + n);
    double medtmp = *(ptr + halfway);
    if (is_even) {
        std::nth_element(ptr, ptr + halfway - 1, ptr + n);
        return (medtmp + *(ptr + halfway - 1))/2;
    } else {
        return medtmp;
    }
}

/**
 * Compute the median from a sparse vector.
 *
 * @param[in] ptr Pointer to an array of non-zero values.
 * This may be modified on output.
 * @param nnz Number of non-zero elements, i.e., the length of the array referenced by `ptr`.
 * @param nall Total number of elements in the set,
 * i.e., `nall - nnz` is the number of zeros.
 *
 * @tparam Output_ Type of the output value.
 * This should be a floating-point value for possible averaging.
 * @tparam Value_ Type of the input values.
 *
 * @return The median of values in the sparse vector.
 */
template<typename Output_ = double, typename Value_>
Output_ compute_median(Value_* ptr, size_t nnz, size_t nall) {
    if (nnz == nall) {
        return compute_median<Output_>(ptr, nall);
    }

    if (nnz * 2 < nall) {
        return 0; // zero is the median if there are too many zeroes.
    } 
    
    size_t halfway = nall / 2;
    bool is_even = (nall % 2 == 0);

    auto vend = ptr + nnz;
    std::sort(ptr, vend);
    size_t zeropos = std::lower_bound(ptr, vend, 0) - ptr;
    size_t nzero = nall - nnz;

    if (!is_even) {
        if (zeropos > halfway) {
            return ptr[halfway];
        } else if (halfway >= zeropos + nzero) {
            return ptr[halfway - nzero];
        } else {
            return 0; // zero is the median.
        }
    }

    double tmp = 0;
    if (zeropos > halfway) {
        tmp = ptr[halfway] + ptr[halfway - 1];
    } else if (zeropos == halfway) {
        // guaranteed to be at least 1 zero.
        tmp += ptr[halfway - 1];
    } else if (zeropos < halfway && zeropos + nzero > halfway) {
        ; // zero is the median.
    } else if (zeropos + nzero == halfway) {
        // guaranteed to be at least 1 zero.
        tmp += ptr[halfway - nzero];
    } else {
        tmp = ptr[halfway - nzero] + ptr[halfway - nzero - 1];
    }
    return tmp / 2;
}

/**
 * @cond
 */
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
                output[i] = compute_median<Output_>(vbuffer, range.number, otherdim);
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
