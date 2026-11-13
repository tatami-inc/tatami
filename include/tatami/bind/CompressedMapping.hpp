#ifndef TATAMI_DELAYED_BIND_COMPRESSED_MAPPING_HPP
#define TATAMI_DELAYED_BIND_COMPRESSED_MAPPING_HPP

#include <vector>
#include <limits>

#include "sanisizer/sanisizer.hpp"

#include "../utils/copy.hpp"

namespace tatami {

namespace DelayedBind_internal {

/*
 * Consider the problem of mapping an index 'i' of the combined dimension of the DelayedBind object to the component matrix.
 * The simplest approach is to do a binary search on 'cumulative' but that has O(log(n)) scaling for 'n' component matrices.
 * Another approach is to create a look-up table that can be indexed directly by 'i', which has O(1) runtime but consumes a lot of memory.
 *
 * The CompressedMapping class is a compromise that breaks up the combined dimensions into bins of regular size no greater than the extent of the smallest component matrix.
 * Each bin stores the identity the component matrix in which it starts, plus the length of that component matrix covered by that bin.
 * Given 'i', we can divide by the bin size to find its matching bin; we then compare the remainder to the covered length to assign it to the starting or next matrix.
 * This has O(1) runtime while being much more memory efficient for a minimum bin size >> 1.
 *
 * In practice, the bin size is set to the largest power of 2 that is no greater than the small matrix extent.
 * This allows for (i) some very fast division via bit shifting, and (ii) bitpacking to store both the matrix identity and the covered length in a single integer.
 */
template<typename Index_>
class CompressedMapping {
public:
    CompressedMapping() = default;

    CompressedMapping(Index_ min_size, const std::vector<Index_>& cumulative) {
        assert(min_size > 0);
        assert(cumulative.size() > 1);
        assert(sanisizer::is_less_than_or_equal(min_size, compute_min_diff(cumulative))); 

        // Define bin size as the smallest power of 2 that is no greater than 'min_size'.
        // This means that 'my_shift' must be less than the number of binary digits of 'Index_'.
        while (min_size > 1) {
            ++my_shift;
            min_size /= 2;
        }
        assert(sanisizer::is_less_than(my_shift, std::numeric_limits<Index_>::digits));
        const auto bin_size = static_cast<Index_>(1) << my_shift;
        my_mask = (bin_size - 1);

        const auto total_dim = cumulative.back();
        assert(min_size <= total_dim);
        const auto total_bins = (total_dim >> my_shift) + ((total_dim & my_mask) > 0);
        my_bins.reserve(total_bins);

        Index_ pos = 1;
        for (I<decltype(total_bins)> b = 0; b < total_bins; ++b) {
            // Find the first cumulative boundary that exceeds the start of the bin.
            // 'pos - 1' is the index of the matrix covering the start of the bin.
            Index_ lower = b << my_shift;
            while (cumulative[pos] <= lower) {
                ++pos;
            }

            // We compute the gap between the start of the bin and the cumulative boundary immediately to its right.
            // This is the number of elements in the bin that are assigned to the 'pos - 1'-th matrix.
            const auto until_next = sanisizer::min(cumulative[pos] - lower, bin_size);

            // Now we pack both of these values into a single integer to save some space.
            // This is possible based on the following logic:
            //
            // - The total cumulative dimension is less than '2^X' where 'X' is the number of digits of 'Index_'.
            // - The minimum size of each matrix is greater than or equal to '2^Y', where 'Y' is 'my_shift'.
            //   We note that 'X > Y' as the minimum size cannot be greater than the total cumulative dimension.
            // - Thus, the maximum possible value of the number of matrices must be less than '2^(X - Y)'.
            //   This implies that the matrix index must fit into a '(X - Y)'-bit unsigned integer.
            // - We know that 'until_next' must be less than or equal to '2^Y'.
            //   We also know it must be positive as the cumulative boundary exceeds the bin start.
            //   This implies that 'until_next - 1' must fit into a Y-bit unsigned integer.
            //
            // Thus, we can shift the matrix index by 'Y' and use the cleared bits to store 'until_next - 1'.
            const Index_ encoding = (static_cast<Index_>(pos - 1) << my_shift) | static_cast<Index_>(until_next - 1);
            my_bins.push_back(encoding);
        }
    }

    Index_ operator[](const Index_ i) const {
        const auto encoded = my_bins[i >> my_shift];
        // 'encoded >> my_shift' recovers the matrix index.
        // 'encoded & my_mask' recovers 'until_next - 1'.
        // We want to know if 'i & my_mask >= until_next', which is the same as 'i & my_mask > until_next - 1', saving us an addition.
        return sanisizer::sum_unsafe<Index_>(encoded >> my_shift, sanisizer::is_greater_than(i & my_mask, encoded & my_mask)); 
    }

private:
    std::vector<Index_> my_bins;
    Index_ my_shift = 0; 
    Index_ my_mask;

    static Index_ compute_min_diff(const std::vector<Index_>& cumulative) {
        auto output = cumulative[1] - cumulative[0];
        const auto nmats = cumulative.size() - 1;
        for (I<decltype(nmats)> i = 1; i < nmats; ++i) {
            output = sanisizer::min(output, cumulative[i + 1] - cumulative[i]);
        }
        return output;
    }
};

}

}

#endif
