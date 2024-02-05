#ifndef TATAMI_GROUPED_MEDIANS_HPP
#define TATAMI_GROUPED_MEDIANS_HPP

#include "../base/Matrix.hpp"
#include "medians.hpp"
#include "utils.hpp"
#include <vector>
#include <algorithm>

/**
 * @file grouped_medians.hpp
 *
 * @brief Compute group-wise medians from a `tatami::Matrix`.
 */

namespace tatami {

namespace stats {

/**
 * @cond
 */
template<bool row_, typename Value_, typename Index_, typename Group_, class GroupSizes_, typename Output_>
void grouped_medians(const tatami::Matrix<Value_, Index_>* p, const Group_* groups, const GroupSizes_& group_sizes, Output_* output, int threads) {
    Index_ dim = (row_ ? p->nrow() : p->ncol());
    Index_ otherdim = (row_ ? p->ncol() : p->nrow());

    parallelize([&](int, Index_ start, Index_ len) -> void {
        std::vector<Value_> xbuffer(otherdim);

        size_t ngroups = group_sizes.size();
        std::vector<std::vector<double> > workspace(ngroups);
        for (size_t g = 0; g < ngroups; ++g) {
            workspace[g].reserve(group_sizes[g]);
        }

        // Always convert to size_t when doing any pointer arithmetic.
        auto curoutput = output + static_cast<size_t>(start) * ngroups;

        if (p->sparse()) {
            tatami::Options opt;
            opt.sparse_ordered_index = false;

            auto ext = tatami::consecutive_extractor<row_, true>(p, start, len, opt);
            std::vector<Index_> ibuffer(otherdim);

            for (Index_ i = start, end = start + len; i < end; ++i) {
                auto range = ext->fetch(i, xbuffer.data(), ibuffer.data());
                for (Index_ j = 0; j < range.number; ++j) {
                    workspace[groups[range.index[j]]].push_back(range.value[j]);
                }

                for (size_t g = 0; g < ngroups; ++g, ++curoutput) {
                    auto& w = workspace[g];
                    *curoutput = compute_median(w.data(), w.size(), group_sizes[g]);
                    w.clear();
                }
            }

        } else {
            auto ext = tatami::consecutive_extractor<row_, false>(p, start, len);
            for (Index_ i = start, end = start + len; i < end; ++i) {
                auto ptr = ext->fetch(i, xbuffer.data());
                for (Index_ j = 0; j < otherdim; ++j) {
                    workspace[groups[j]].push_back(ptr[j]);
                }

                for (auto& w : workspace) {
                    *curoutput = compute_median(w.data(), w.size());
                    ++curoutput;
                    w.clear();
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
 * @tparam GroupSizes_ Vector-like class that has `size()` and `[` methods.
 * @tparam Output_ Type of the output.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of columns.
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in `[0, N)` where `N` is the number of unique groups.
 * @param group_sizes Vector-like object of length `N`, specifying the number of columns assigned to each group.
 * This can be created by calling `stats::tabulate_groups()` on `group`.
 * @param[out] output Pointer to an array of length equal to the product of the number of rows and `N`.
 * On output, `output[g + N * r]` will contain the median for row `r` and group `g`.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Group_, typename GroupSizes_, typename Output_>
void row_medians_by_group(const tatami::Matrix<Value_, Index_>* p, const Group_* group, const GroupSizes_& group_sizes, Output_* output, int threads = 1) {
    stats::grouped_medians<true>(p, group, group_sizes, output, threads);
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
 * The entry at `g + N * r` will contain the median for row `r` and group `g`.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<Output_> row_medians_by_group(const tatami::Matrix<Value_, Index_>* p, const Group_* group, int threads = 1) {
    auto group_sizes = stats::tabulate_groups(group, p->ncol());
    std::vector<Output_> output(group_sizes.size() * static_cast<size_t>(p->nrow()));
    row_medians_by_group(p, group, group_sizes, output.data(), threads);
    return output;
}

/**
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Group_ Type of the group assignments for each column.
 * @tparam GroupSizes_ Vector-like class that has `size()` and `[` methods.
 * @tparam Output_ Type of the output.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of rows.
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in `[0, N)` where `N` is the number of unique groups.
 * @param group_sizes Vector-like object of length `N`, specifying the number of columns assigned to each group.
 * This can be created by calling `stats::tabulate_groups()` on `group`.
 * @param[out] output Pointer to an array of length equal to the product of the number of columns and `N`.
 * On output, `output[g + N * c]` will contain the median for column `c` and group `g`.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Group_, typename GroupSizes_, typename Output_>
void column_medians_by_group(const tatami::Matrix<Value_, Index_>* p, const Group_* group, const GroupSizes_& group_sizes, Output_* output, int threads = 1) {
    stats::grouped_medians<false>(p, group, group_sizes, output, threads);
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
 * The entry at `g + N * c` will contain the median for column `c` and group `g`.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<Output_> column_medians_by_group(const tatami::Matrix<Value_, Index_>* p, const Group_* group, int threads = 1) {
    auto group_sizes = stats::tabulate_groups(group, p->nrow());
    std::vector<Output_> output(group_sizes.size() * static_cast<size_t>(p->ncol()));
    column_medians_by_group(p, group, group_sizes, output.data(), threads);
    return output;
}

}

#endif
