#ifndef TATAMI_GROUPED_SUMS_HPP
#define TATAMI_GROUPED_SUMS_HPP

#include "../base/Matrix.hpp"
#include "sums.hpp"
#include "utils.hpp"
#include <vector>
#include <algorithm>

/**
 * @file grouped_sums.hpp
 *
 * @brief Compute group-wise sums from a `tatami::Matrix`.
 */

namespace tatami {

namespace stats {

/**
 * @cond
 */
template<bool row_, typename Value_, typename Index_, typename Group_, typename Output_>
void grouped_sums(const tatami::Matrix<Value_, Index_>* p, const Group_* groups, const std::vector<Index_>& group_sizes, Output_* output, int threads) {
    Index_ dim = (row_ ? p->nrow() : p->ncol());
    Index_ otherdim = (row_ ? p->ncol() : p->nrow());

    tatami::parallelize([&](int, Index_ start, Index_ len) -> void {
        // Always convert to size_t when doing any pointer arithmetic.
        size_t num_groups = group_sizes.size();
        auto curoutput = output + static_cast<size_t>(start) * num_groups;
        std::fill(curoutput, curoutput + num_groups * static_cast<size_t>(len), 0.0);

        if (p->sparse()) {
            Options opt;
            opt.sparse_ordered_index = false;

            if (p->prefer_rows() == row_) {
                auto ext = tatami::consecutive_extractor<row_, true>(p, start, len, opt);
                std::vector<Value_> xbuffer(otherdim);
                std::vector<Index_> ibuffer(otherdim);

                for (Index_ i = start, end = start + len; i < end; ++i) {
                    auto range = ext->fetch(i, xbuffer.data(), ibuffer.data());
                    for (int j = 0; j < range.number; ++j) {
                        curoutput[groups[range.index[j]]] += range.value[j];
                    }
                    curoutput += num_groups;
                }

            } else {
                auto ext = tatami::consecutive_extractor<!row_, true>(p, 0, otherdim, start, len, opt);
                std::vector<Value_> xbuffer(len);
                std::vector<Index_> ibuffer(len);

                for (int i = 0; i < otherdim; ++i) {
                    auto range = ext->fetch(i, xbuffer.data(), ibuffer.data());
                    auto outcopy = curoutput + groups[i];
                    for (int j = 0; j < range.number; ++j) {
                        outcopy[static_cast<size_t>(range.index[j] - start) * num_groups] += range.value[j];
                    }
                }
            }

        } else {
            if (p->prefer_rows() == row_) {
                std::vector<Value_> xbuffer(otherdim);
                auto ext = tatami::consecutive_extractor<row_, false>(p, start, len);

                for (Index_ i = start, end = start + len; i < end; ++i) {
                    auto ptr = ext->fetch(i, xbuffer.data());
                    for (Index_ j = 0; j < otherdim; ++j) {
                        curoutput[groups[j]] += ptr[j];
                    }
                    curoutput += num_groups;
                }

            } else {
                std::vector<double> xbuffer(len);
                auto ext = tatami::consecutive_extractor<!row_, false>(p, 0, otherdim, start, len);

                for (int i = 0; i < otherdim; ++i) {
                    auto ptr = ext->fetch(i, xbuffer.data());
                    auto outcopy = curoutput + groups[i];
                    for (int j = 0; j < len; ++j) {
                        *outcopy += ptr[j];
                        outcopy += num_groups;
                    }
                }
            }
        }
    }, dim, threads);
}
/**
 * @endcond
 */

}

/**
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Group_ Type of the group assignments for each column.
 * @tparam Output_ Type of the output.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of columns.
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in `[0, N)` where `N` is the number of unique groups.
 * @param group_sizes Vector of length `N`, specifying the number of columns assigned to each group.
 * @param[out] output Pointer to an array of length equal to the product of the number of rows and `N`.
 * On output, `output[g + N * r]` will contain the sum for row `r` and group `g`.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Group_, typename Output_>
void row_sums_by_group(const tatami::Matrix<Value_, Index_>* p, const Group_* group, const std::vector<Index_>& group_sizes, Output_* output, int threads = 1) {
    stats::grouped_sums<true>(p, group, group_sizes, output, threads);
}

/**
 * @tparam Output_ Type of the output.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Group_ Type of the group assignments for each column.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of columns.
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in `[0, N)` where `N` is the number of unique groups.
 * @param threads Number of threads to use.
 *
 * @return Vector of length equal to the product of the number of rows of `p` and `N`.
 * The entry at `g + N * r` will contain the sum for row `r` and group `g`.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<Output_> row_sums_by_group(const tatami::Matrix<Value_, Index_>* p, const Group_* group, int threads = 1) {
    auto group_sizes = tabulate_groups(group, p->ncol());
    std::vector<Output_> output(group_sizes.size() * static_cast<size_t>(p->nrow()));
    row_sums_by_group(p, group, group_sizes, output.data(), threads);
    return output;
}

/**
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Group_ Type of the group assignments for each column.
 * @tparam Output_ Type of the output.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of rows.
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in `[0, N)` where `N` is the number of unique groups.
 * @param group_sizes Vector of length `N`, specifying the number of columns assigned to each group.
 * @param[out] output Pointer to an array of length equal to the product of the number of columns and `N`.
 * On output, `output[g + N * c]` will contain the sum for column `c` and group `g`.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Group_, typename Output_>
void column_sums_by_group(const tatami::Matrix<Value_, Index_>* p, const Group_* group, const std::vector<Index_>& group_sizes, Output_* output, int threads = 1) {
    stats::grouped_sums<false>(p, group, group_sizes, output, threads);
}

/**
 * @tparam Output_ Type of the output.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the column/column indices.
 * @tparam Group_ Type of the group assignments for each column.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of columns.
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in `[0, N)` where `N` is the number of unique groups.
 * @param threads Number of threads to use.
 *
 * @return Vector of length equal to the product of the number of columns of `p` and `N`.
 * The entry at `g + N * c` will contain the sum for column `c` and group `g`.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<Output_> column_sums_by_group(const tatami::Matrix<Value_, Index_>* p, const Group_* group, int threads = 1) {
    auto group_sizes = tabulate_groups(group, p->nrow());
    std::vector<Output_> output(group_sizes.size() * static_cast<size_t>(p->ncol()));
    column_sums_by_group(p, group, group_sizes, output.data(), threads);
    return output;
}

}

#endif
