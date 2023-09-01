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
template<bool row_, typename Value_, typename Index_, typename Group_, typename Output_>
void grouped_medians(const tatami::Matrix<Value_, Index_>* p, const Group_* groups, const std::vector<Index_>& group_sizes, Output_* output, int threads) {
    int dim = (row_ ? p->nrow() : p->ncol());
    int otherdim = (row_ ? p->ncol() : p->nrow());

    tatami::parallelize([&](int, int start, int len) -> void {
        std::vector<double> xbuffer(otherdim);

        std::vector<std::vector<double> > workspace;
        for (auto l : group_sizes) {
            workspace.emplace(l);
        }

        // Always convert to size_t when doing any pointer arithmetic.
        auto curoutput = output + static_cast<size_t>(start) * group_sizes.size();

        if (reference->sparse()) {
            tatami::Options opt;
            opt.sparse_ordered_index = false;

            auto ext = tatami::consecutive_extractor<true, true>(reference, start, len, opt);
            std::vector<int> ibuffer(NC);

            for (int r = start, end = start + len; r < end; ++r) {
                auto range = ext->fetch(r, xbuffer.data(), ibuffer.data());
                for (int i = 0; i < range.number; ++i) {
                    workspace[groups[range.index[i]]].push_back(range.value[i]);
                }

                auto lIt = group_sizes.begin();
                for (auto& w : workspace) {
                    *curoutput = compute_median(w.data(), w.size(), *lIt);
                    ++curoutput;
                    w.clear();
                }
            }

        } else {
            auto ext = tatami::consecutive_extractor<true, false>(reference, start, len);
            for (int r = start, end = start + len; r < end; ++r) {
                auto ptr = ext->fetch(r, xbuffer.data());
                for (int c = 0; c < NC; ++c) {
                    workspace[groups[c]].push_back(ptr[c]);
                }

                for (auto& w : workspace) {
                    *curoutput = compute_median(w.data(), w.size());
                    ++curoutput;
                    w.clear();
                }
            }
        }
    }, dim, num_threads);
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
 * On output, `output[g + N * r]` will contain the median for row `r` and group `g`.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Group_, typename Output_>
void row_medians_by_group(const tatami::Matrix<Value_, Index_>* p, const Group_* group, const std::vector<Index_>& group_sizes, Output_* output, int threads) {
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
std::vector<Output_> row_medians_by_group(const tatami::Matrix<Value_, Index_>* p, const Group_* group, int threads) const {
    Index_ NC = p->ncol();
    Index_ ngroups = (NC ? static_cast<Index_>(*std::max_element(group, group + NC)) + 1 : 0);

    std::vector<Index_> group_sizes(ngroups);
    for (Index_ c = 0; c < NC; ++c) {
        ++(group_sizes[group[c]]);
    }

    std::vector<Output_> output(static_cast<size_t>(ngroups) * static_cast<size_t>(p->nrow()));
    row_medians_by_group(p, group, group_sizes, output.data(), threads);
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
 * On output, `output[g + N * c]` will contain the median for column `c` and group `g`.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Group_, typename Output_>
void column_medians_by_group(const tatami::Matrix<Value_, Index_>* p, const Group_* group, const std::vector<Index_>& group_sizes, Output_* output, int threads) {
    stats::grouped_medians<true>(p, group, group_sizes, output, threads);
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
std::vector<Output_> column_medians_by_group(const tatami::Matrix<Value_, Index_>* p, const Group_* group, int threads) const {
    Index_ NC = p->ncol();
    Index_ ngroups = (NC ? static_cast<Index_>(*std::max_element(group, group + NC)) + 1 : 0);

    std::vector<Index_> group_sizes(ngroups);
    for (Index_ c = 0; c < NC; ++c) {
        ++(group_sizes[group[c]]);
    }

    std::vector<Output_> output(static_cast<size_t>(ngroups) * static_cast<size_t>(p->ncolumn()));
    column_medians_by_group(p, group, group_sizes, output.data(), threads);
    return output;
}


}

#endif
