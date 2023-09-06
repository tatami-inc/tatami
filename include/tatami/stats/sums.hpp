#ifndef TATAMI_STATS_SUMS_HPP
#define TATAMI_STATS_SUMS_HPP

#include "../base/Matrix.hpp"
#include "utils.hpp"
#include <vector>
#include <numeric>
#include <algorithm>

/**
 * @file sums.hpp
 *
 * @brief Compute row and column sums from a `tatami::Matrix`.
 */

namespace tatami {

namespace stats {

/**
 * @cond
 */
template<bool row_, typename Value_, typename Index_, typename Output_>
void dimension_sums(const Matrix<Value_, Index_>* p, Output_* output, int threads) {
    auto dim = (row_ ? p->nrow() : p->ncol());
    auto otherdim = (row_ ? p->ncol() : p->nrow());
    const bool direct = p->prefer_rows() == row_;

    if (p->sparse()) {
        if (direct) {
            Options opt;
            opt.sparse_extract_index = false;
            parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = consecutive_extractor<row_, true>(p, s, l, opt);
                std::vector<Value_> vbuffer(otherdim);
                for (Index_ i = s, e = s + l; i < e; ++i) {
                    auto out = ext->fetch(i, vbuffer.data(), NULL);
                    output[i] = std::accumulate(out.value, out.value + out.number, static_cast<Output_>(0));
                }
            }, dim, threads);

        } else {
            std::fill(output, output + dim, static_cast<Output_>(0));

            parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = consecutive_extractor<!row_, true>(p, 0, otherdim, s, l);
                auto len = ext->block_length;
                std::vector<Value_> vbuffer(len);
                std::vector<Index_> ibuffer(len);

                for (Index_ i = 0; i < otherdim; ++i) {
                    auto out = ext->fetch(i, vbuffer.data(), ibuffer.data());
                    for (Index_ j = 0; j < out.number; ++j) {
                        output[out.index[j]] += out.value[j];
                    }
                }
            }, dim, threads);
        }

    } else {
        if (direct) {
            parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = consecutive_extractor<row_, false>(p, s, l);
                std::vector<Value_> buffer(otherdim);
                for (Index_ i = s, e = s + l; i < e; ++i) {
                    auto out = ext->fetch(i, buffer.data());
                    output[i] = std::accumulate(out, out + otherdim, static_cast<Output_>(0));
                }
            }, dim, threads);

        } else {
            std::fill(output, output + dim, static_cast<Output_>(0));

            parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = consecutive_extractor<!row_, false>(p, 0, otherdim, s, l);
                std::vector<Value_> buffer(ext->block_length);
                auto len = ext->block_length;
                for (Index_ i = 0; i < otherdim; ++i) {
                    auto out = ext->fetch(i, buffer.data());
                    for (Index_ j = 0; j < len; ++j) {
                        output[s + j] += out[j];
                    }
                }
            }, dim, threads);
        }
    }

    return;
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
 * @param[out] output Pointer to an array of length equal to the number of columns.
 * On output, this will store the sum of values for each column.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Output_>
void column_sums(const Matrix<Value_, Index_>* p, Output_* output, int threads = 1) {
    stats::dimension_sums<false>(p, output, threads);
    return;
}

/**
 * @tparam Output_ Type of the output value.
 * @tparam Value_ Type of the matrix value, should be summable.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A vector of length equal to the number of columns, containing the column sums.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> column_sums(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> output(p->ncol());
    column_sums(p, output.data(), threads);
    return output;
}

/**
 * @tparam Output_ Type of the output value.
 * @tparam Value_ Type of the matrix value, should be summable.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[out] output Pointer to an array of length equal to the number of rows.
 * On output, this will contain the row sums.
 * @param threads Number of threads to use.
 */
template<typename Output_ = double, typename Value_, typename Index_>
void row_sums(const Matrix<Value_, Index_>* p, Output_* output, int threads = 1) {
    stats::dimension_sums<true>(p, output, threads);
    return;
}

/**
 * @tparam Output_ Type of the output value.
 * @tparam Value_ Type of the matrix value, should be summable.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A vector of length equal to the number of rows, containing the row sums.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> row_sums(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> output(p->nrow());
    row_sums(p, output.data(), threads);
    return output;
}

}

#endif
