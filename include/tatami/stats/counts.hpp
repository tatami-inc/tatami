#ifndef TATAMI_STATS_COUNTS_HPP
#define TATAMI_STATS_COUNTS_HPP

#include "../base/Matrix.hpp"
#include "utils.hpp"
#include <vector>
#include <algorithm>
#include <cmath>

/**
 * @file counts.hpp
 *
 * @brief Compute row and column counts from a `tatami::Matrix`.
 */

namespace tatami {

namespace stats {

/**
 * @cond
 */
template<bool row_, typename Output_, typename Value_, typename Index_, class Function_>
void dimension_counts(const Matrix<Value_, Index_>* p, int threads, Output_* output, Function_ fun) {
    auto dim = (row_ ? p->nrow() : p->ncol());
    auto otherdim = (row_ ? p->ncol() : p->nrow());
    std::fill(output, output + dim, 0);

    if (p->prefer_rows() == row_) {
        if (p->sparse()) {
            Options opt;
            opt.sparse_ordered_index = false;
            Output_ zerocount = fun(0);

            tatami::parallelize([&](int, Index_ start, Index_ len) -> void {
                std::vector<Value_> xbuffer(otherdim);
                std::vector<Index_> ibuffer(otherdim);
                auto ext = tatami::consecutive_extractor<row_, true>(p, start, len, opt);

                for (Index_ i = start, end = start + len; i < end; ++i) {
                    auto range = ext->fetch(i, xbuffer.data(), ibuffer.data());
                    Output_ target = 0;
                    for (Index_ j = 0; j < range.number; ++j) {
                        target += fun(range.value[j]);
                    }
                    output[i] = target + zerocount * (otherdim - range.number);
                }
            }, dim, threads);

        } else {
            tatami::parallelize([&](int, Index_ start, Index_ len) -> void {
                std::vector<Value_> xbuffer(otherdim);
                auto ext = tatami::consecutive_extractor<row_, false>(p, start, len);

                for (Index_ i = start, end = start + len; i < end; ++i) {
                    auto ptr = ext->fetch(i, xbuffer.data());
                    Output_ target = 0;
                    for (Index_ j = 0; j < otherdim; ++j) {
                        target += fun(ptr[j]);
                    }
                    output[i] = target;
                }
            }, dim, threads);
        }

    } else {
        std::vector<Output_*> threaded_output_ptrs(threads, output);
        std::vector<std::vector<Output_> > threaded_output(threads - 1);
        for (int t = 1; t < threads; ++t) {
            auto& curout = threaded_output[t - 1];
            curout.resize(dim);
            threaded_output_ptrs[t] = curout.data();
        }

        if (p->sparse()) {
            Options opt;
            opt.sparse_ordered_index = false;
            Output_ zerocount = fun(0);

            tatami::parallelize([&](int t, Index_ start, Index_ len) -> void {
                std::vector<Value_> xbuffer(dim);
                std::vector<Index_> ibuffer(dim);
                auto ext = tatami::consecutive_extractor<!row_, true>(p, start, len, opt);

                auto curoutput = threaded_output_ptrs[t];
                std::vector<Index_> nonzeros(dim);

                for (Index_ i = start, end = start + len; i < end; ++i) {
                    auto range = ext->fetch(i, xbuffer.data(), ibuffer.data());
                    for (Index_ j = 0; j < range.number; ++j) {
                        curoutput[range.index[j]] += fun(range.value[j]);
                        ++(nonzeros[range.index[j]]);
                    }
                }

                for (int d = 0; d < dim; ++d) {
                    curoutput[d] += zerocount * (len - nonzeros[d]);
                }
            }, otherdim, threads);

        } else {
            tatami::parallelize([&](int t, Index_ start, Index_ len) -> void {
                std::vector<Value_> xbuffer(dim);
                auto ext = tatami::consecutive_extractor<!row_, false>(p, start, len);
                auto curoutput = threaded_output_ptrs[t];

                for (Index_ i = start, end = start + len; i < end; ++i) {
                    auto ptr = ext->fetch(i, xbuffer.data());
                    for (Index_ j = 0; j < dim; ++j) {
                        curoutput[j] += fun(ptr[j]);
                    }
                }
            }, otherdim, threads);
        }

        for (int t = 1; t < threads; ++t) {
            auto curoutput = threaded_output_ptrs[t];
            for (Index_ d = 0; d < dim; ++d) {
                output[d] += curoutput[d];
            }
        }
    }
}
/**
 * @endcond
 */

}

/**
 * @tparam Value_ Type of the matrix value, should be summable.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output_ Type of the output value.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[out] output Pointer to an array of length equal to the number of rows.
 * On output, this will store the number of NaNs in each row.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Output_>
void row_nan_counts(const Matrix<Value_, Index_>* p, Output_* output, int threads = 1) {
    stats::dimension_counts<true>(p, threads, output, [](Value_ x) -> bool { return std::isnan(x); });
}

/**
 * @tparam Output_ Type of the output value.
 * @tparam Value_ Type of the matrix value, should be summable.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A vector of length equal to the number of rows, containing the number of NaNs in each row.
 */
template<typename Output_ = int, typename Value_, typename Index_>
std::vector<Output_> row_nan_counts(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> output(p->nrow());
    row_nan_counts(p, output.data(), threads);
    return output;
}

/**
 * @tparam Value_ Type of the matrix value, should be summable.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output_ Type of the output value.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[out] output Pointer to an array of length equal to the number of columns.
 * On output, this will store the number of NaNs in each column.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Output_>
void column_nan_counts(const Matrix<Value_, Index_>* p, Output_* output, int threads = 1) {
    stats::dimension_counts<false>(p, threads, output, [](Value_ x) -> bool { return std::isnan(x); });
}

/**
 * @tparam Output_ Type of the output value.
 * @tparam Value_ Type of the matrix value, should be summable.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A vector of length equal to the number of columns, containing the number of NaNs in each column.
 */
template<typename Output_ = int, typename Value_, typename Index_>
std::vector<Output_> column_nan_counts(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> output(p->ncol());
    column_nan_counts(p, output.data(), threads);
    return output;
}

/**
 * @tparam Value_ Type of the matrix value, should be summable.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output_ Type of the output value.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[out] output Pointer to an array of length equal to the number of rows.
 * On output, this will store the number of zeros in each row.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Output_>
void row_zero_counts(const Matrix<Value_, Index_>* p, Output_* output, int threads = 1) {
    stats::dimension_counts<true>(p, threads, output, [](Value_ x) -> bool { return x == 0; });
}

/**
 * @tparam Output_ Type of the output value.
 * @tparam Value_ Type of the matrix value, should be summable.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A vector of length equal to the number of rows, containing the number of zeros in each row.
 */
template<typename Output_ = int, typename Value_, typename Index_>
std::vector<Output_> row_zero_counts(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> output(p->nrow());
    row_zero_counts(p, output.data(), threads);
    return output;
}

/**
 * @tparam Value_ Type of the matrix value, should be summable.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output_ Type of the output value.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[out] output Pointer to an array of length equal to the number of columns.
 * On output, this will store the number of zeros in each column.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Output_>
void column_zero_counts(const Matrix<Value_, Index_>* p, Output_* output, int threads = 1) {
    stats::dimension_counts<false>(p, threads, output, [](Value_ x) -> bool { return x == 0; });
}

/**
 * @tparam Output_ Type of the output value.
 * @tparam Value_ Type of the matrix value, should be summable.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A vector of length equal to the number of columns, containing the number of zeros in each column.
 */
template<typename Output_ = int, typename Value_, typename Index_>
std::vector<Output_> column_zero_counts(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> output(p->ncol());
    column_zero_counts(p, output.data(), threads);
    return output;
}

}

#endif
