#ifndef TATAMI_STATS_UTILS_HPP
#define TATAMI_STATS_UTILS_HPP

#include <vector>
#include <algorithm>

/**
 * @file utils.hpp
 *
 * @brief Utilities for computing matrix statistics.
 */

namespace tatami {

namespace stats {

/**
 * Count the total number of groups, typically for per-group memory allocations.
 *
 * @tparam Group_ Integer type for the group assignments.
 * @tparam Index_ Integer type for the number of observations.
 *
 * @param[in] group Pointer to an array of group assignments per observation.
 * Each assignment should be an integer in `[0, G)` where `G` is the total number of groups.
 * @param n Number of observations, i.e., the length of the array referenced by `group`.
 *
 * @return Total number of groups, i.e., `G`.
 * Note that not all groups may actually have non-zero occurrences in `group`.
 */
template<typename Group_, typename Size_>
size_t total_groups(const Group_* group, Size_ n) {
    if (n) {
        return static_cast<size_t>(*std::max_element(group, group + n)) + 1;
    } else {
        return 0;
    }
}

/**
 * Count the occurrences of each group.
 *
 * @tparam Group_ Integer type for the group assignments.
 * @tparam Index_ Integer type for the number of observations.
 *
 * @param[in] group Pointer to an array of group assignments per observation.
 * Each assignment should be an integer in `[0, G)` where `G` is the total number of groups.
 * @param n Number of observations, i.e., the length of the array referenced by `group`.
 *
 * @return Vector of length equal to `G`, containing the number of occurrences of each group.
 */
template<typename Group_, typename Size_>
std::vector<Size_> tabulate_groups(const Group_* group, Size_ n) {
    auto ngroups = total_groups(group, n);
    std::vector<Size_> group_sizes(ngroups);
    for (Size_ r = 0; r < n; ++r) {
        ++(group_sizes[group[r]]);
    }
    return group_sizes;
}

}

}

#endif
