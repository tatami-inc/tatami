#ifndef TATAMI_COMPRESS_SPARSE_TRIPLETS_H
#define TATAMI_COMPRESS_SPARSE_TRIPLETS_H

#include <vector>
#include <algorithm>
#include <numeric>
#include <utility>
#include <cstddef>

#include "sanisizer/sanisizer.hpp"

#include "../utils/Index_to_container.hpp"
#include "../utils/copy.hpp"

/**
 * @file compress_sparse_triplets.hpp
 *
 * @brief Convert sparse data in triplet format to a compressed sparse row/column format.
 */

namespace tatami {

/**
 * @cond
 */
namespace compress_triplets {

template<class Primary_, class Secondary_>
int is_ordered(const Primary_& primary, const Secondary_& secondary) {
    if (!std::is_sorted(primary.begin(), primary.end())) {
        return 2;
    }

    const auto nprimary = primary.size();
    I<decltype(nprimary)> start = 0;
    while (start < nprimary) {
        I<decltype(nprimary)> end = start + 1;
        while (end < nprimary && primary[end] == primary[start]) {
            if (secondary[end] < secondary[end - 1]) {
                // Quit on first failure; we've seen enough.
                return 1;
            }
            ++end;
        }
        start = end;
    }

    return 0;
}

template<typename Size_, class Primary_, class Secondary_>
void order(const int status, std::vector<Size_>& indices, const Primary_& primary, const Secondary_& secondary) {
    if (status == 1) {
        const auto nprimary = primary.size();
        I<decltype(nprimary)> start = 0;
        while (start < nprimary) {
            I<decltype(nprimary)> end = start + 1;
            while (end < nprimary && primary[end] == primary[start]) {
                ++end;
            }

            // Checking if this particular run can be skipped.
            if (!std::is_sorted(secondary.begin() + start, secondary.begin() + end)) {
                std::sort(
                    indices.begin() + start,
                    indices.begin() + end,
                    [&](const Size_ left, const Size_ right) -> bool {
                        return secondary[left] < secondary[right];
                    }
                );
            }
            start = end;
        }

    } else if (status == 2) {
        std::sort(
            indices.begin(),
            indices.end(),
            [&](const Size_ left, const Size_ right) -> bool {
                if (primary[left] == primary[right]) {
                    return (secondary[left] < secondary[right]);
                }
                return (primary[left] < primary[right]);
            }
        );
    }
}

}
/**
 * @endcond
 */

/**
 * For compressed sparse matrices, we consider the "primary" dimension to be that along which the structural non-zero elements are grouped in memory,
 * e.g., the columns and rows for compressed sparse row and column matrices, respectively.
 * Thus, the matrix only needs to retain the indices of structural non-zeros along the other (i.e., "secondary") dimension.
 *
 * @tparam Values_ Random-access container for the values.
 * This should have a `size()` method and a `[` access operator.
 * @tparam Pointer_ Integer type of the index pointers in the output.
 * @tparam PrimaryIndices_ Random access container for the primary indices.
 * This should have a `size()` method and a `[` access operator.
 * @tparam SecondaryIndices_ Random access container for the secondary indices.
 * This should have a `size()` method and a `[` access operator.
 *
 * @param num_primary Extent of the primary dimension.
 * @param values Values of the structural non-zeros.
 * @param primary_indices Indices of the structural non-zeros along the primary dimension.
 * Values must be non-negative integers less than `num_primary`.
 * The length of the vector should be equal to `values`.
 * @param secondary_indices Indices of the structural non-zeros along the secondary dimension.
 * Values must be non-negative integers.
 * The length of the vector should be equal to `values`.
 *
 * @return `secondary_indices` and `values` are sorted in-place, as if all the structural non-zeros were sorted by increasing `primary_indices` and then `secondary_indices`.
 * A vector of index pointers is returned with length `num_primary + 1`, specifying the subarray of `values` and `secondary_indices` corresponding to each primary dimension element.
 *
 * Note that `primary_indices` is not modified inside this function for efficiency.
 * If needed, a sorted `primary_indices` can be easily created by simply filling `primary_indices` based on the differences between consecutive index pointers.
 */
template<class Values_, typename Pointer_ = I<decltype(std::declval<Values_>().size())>, class PrimaryIndices_, class SecondaryIndices_>
std::vector<Pointer_> compress_sparse_triplets(std::size_t num_primary, Values_& values, const PrimaryIndices_& primary_indices, SecondaryIndices_& secondary_indices) {
    const auto N = values.size();
    if (!safe_non_negative_equal(N, primary_indices.size()) || !safe_non_negative_equal(N, secondary_indices.size())) {
        throw std::runtime_error("'primary_indices', 'secondary_indices' and 'values' should have the same length");
    }

    const int order_status = compress_triplets::is_ordered(primary_indices, secondary_indices);
    if (order_status != 0) {
        auto indices = sanisizer::create<std::vector<I<decltype(N)> > >(N);
        std::iota(indices.begin(), indices.end(), static_cast<I<decltype(N)> >(0));
        compress_triplets::order(order_status, indices, primary_indices, secondary_indices);

        // Reordering values in place. This saves memory and allows us to work with Values_, RowIndices_, etc. that may not have well-defined copy constructors.
        auto used = sanisizer::create<std::vector<unsigned char> >(N);
        for (I<decltype(N)> i = 0; i < N; ++i) {
            if (used[i]) {
                continue;
            }
            auto current = i, replacement = indices[i];
            used[i] = 1;

            while (replacement != i) {
                std::swap(secondary_indices[current], secondary_indices[replacement]);
                std::swap(values[current], values[replacement]);
                current = replacement;
                used[current] = 1;
                replacement = indices[replacement];
            }
        }
    }

    typedef std::vector<Pointer_> Output;
    typedef typename Output::size_type OutputSize;
    Output output(sanisizer::sum<OutputSize>(num_primary, 1));
    for (const auto t : primary_indices) {
        ++(output[sanisizer::sum_unsafe<OutputSize>(t, 1)]);
    }
    std::partial_sum(output.begin(), output.end(), output.begin());

    return output;
}

/**
 * @cond
 */
// Back-compatibility.
template<class Values_, class RowIndices_, class ColumnIndices_>
std::vector<decltype(std::declval<Values_>().size())> compress_sparse_triplets(std::size_t nrow, std::size_t ncol, Values_& values, RowIndices_& row_indices, ColumnIndices_& column_indices, bool csr) {
    if (csr) {
        auto output = compress_sparse_triplets(nrow, values, row_indices, column_indices);
        for (decltype(nrow) r = 0; r < nrow; ++r) {
            std::fill_n(row_indices.begin() + output[r], output[r + 1] - output[r], r);
        }
        return output;
    } else {
        auto output = compress_sparse_triplets(ncol, values, column_indices, row_indices);
        for (decltype(ncol) c = 0; c < ncol; ++c) {
            std::fill_n(column_indices.begin() + output[c], output[c + 1] - output[c], c);
        }
        return output;
    }
}

template<bool row_, class Values_, class RowIndices_, class ColumnIndices_>
auto compress_sparse_triplets(std::size_t nrow, std::size_t ncol, Values_& values, RowIndices_& row_indices, ColumnIndices_& column_indices) {
    return compress_sparse_triplets(nrow, ncol, values, row_indices, column_indices, row_);
}
/**
 * @endcond
 */

}

#endif
