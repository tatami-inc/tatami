#ifndef TATAMI_STATS_SUMS_HPP
#define TATAMI_STATS_SUMS_HPP

#include "../base/Matrix.hpp"
#include "utils.hpp"
#include <vector>
#include <numeric>

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
template<typename Output_, bool row_, typename Value_, typename Index_>
std::vector<Output_> dimension_sums(const Matrix<Value_, Index_>* p, int threads) {
    Index_ dim = (row_ ? p->nrow() : p->ncol());
    Index_ otherdim = (row_ ? p->ncol() : p->nrow());
    std::vector<Output_> output(dim);

    const bool prefrow = p->prefer_rows();
    const bool direct = prefrow == row_;
    const bool oracled = p->uses_oracle(prefrow);
    Options opt;

    if (p->sparse()) {
        if (direct) {
            opt.sparse_extract_index = false;
            parallelize([&](size_t, Index_ s, Index_ e) {
                auto ext = direct_extractor<row_, true>(p, oracled, s, e, opt);
                std::vector<Value_> vbuffer(otherdim);
                for (Index_ i = s; i < e; ++i) {
                    auto out = ext->fetch(i, vbuffer.data(), NULL);
                    output[i] = std::accumulate(out.value, out.value + out.number, static_cast<Output_>(0));
                }
            }, dim, threads);

        } else {
            parallelize([&](size_t, Index_ s, Index_ e) {
                auto len = e - s;
                auto ext = running_extractor<row_, true>(p, s, e, oracled, otherdim, opt);
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
            parallelize([&](size_t, Index_ s, Index_ e) {
                auto ext = direct_extractor<row_, false>(p, oracled, s, e, opt);
                std::vector<Value_> buffer(otherdim);
                for (Index_ i = s; i < e; ++i) {
                    auto out = ext->fetch(i, buffer.data());
                    output[i] = std::accumulate(out, out + otherdim, static_cast<Output_>(0));
                }
            }, dim, threads);

        } else {
            parallelize([&](size_t, Index_ s, Index_ e) {
                auto len = e - s;
                auto ext = running_extractor<row_, false>(p, s, e, oracled, otherdim, opt);
                std::vector<Value_> buffer(len);
                for (Index_ i = 0; i < otherdim; ++i) {
                    auto out = ext->fetch(i, buffer.data());
                    for (Index_ j = 0; j < len; ++j) {
                        output[s + j] += out[j];
                    }
                }
            }, dim, threads);
        }
    }

    return output;
}
/**
 * @endcond
 */

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
    return stats::dimension_sums<Output_, false>(p, threads);
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
    return stats::dimension_sums<Output_, true>(p, threads);
}

}

#endif
