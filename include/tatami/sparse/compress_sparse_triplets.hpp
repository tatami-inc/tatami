#ifndef TATAMI_COMPRESS_SPARSE_TRIPLETS_H
#define TATAMI_COMPRESS_SPARSE_TRIPLETS_H

#include <vector>
#include <algorithm>
#include <numeric>
#include <utility>
#include <cstddef>

#include "sanisizer/sanisizer.hpp"

#include "../utils/copy.hpp"

/**
 * @file compress_sparse_triplets.hpp
 *
 * @brief Convert sparse data in triplet format to a compressed sparse row/column format.
 */

namespace tatami {

/**
 * For compressed sparse matrices, we consider the "primary" dimension to be that along which the structural non-zero elements are grouped in memory,
 * e.g., the columns and rows for compressed sparse row and column matrices, respectively.
 * Thus, the matrix only needs to retain the indices of structural non-zeros along the other (i.e., "secondary") dimension.
 *
 * @tparam Pointer_ Integer type of the output index pointers. 
 * @tparam Extent_ Integer type of the dimension extent.
 * @tparam Values_ Random-access container for the values.
 * This should have a `size()` method and a `[` access operator.
 * @tparam PrimaryIndices_ Random access container for the primary indices.
 * This should have a `size()` method, a `begin()` and `end()` method, and a `[` access operator.
 * @tparam SecondaryIndices_ Random access container for the secondary indices.
 * This should have a `size()` method, a `begin()` and `end()` method, and a `[` access operator.
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
template<typename Pointer_ = std::size_t, typename Extent_, class Values_, class PrimaryIndices_, class SecondaryIndices_>
std::vector<Pointer_> compress_sparse_triplets(const Extent_ num_primary, Values_& values, const PrimaryIndices_& primary_indices, SecondaryIndices_& secondary_indices) {
    const auto num_triplets = values.size();
    if (!sanisizer::is_equal(num_triplets, primary_indices.size()) || !sanisizer::is_equal(num_triplets, secondary_indices.size())) {
        throw std::runtime_error("'primary_indices', 'secondary_indices' and 'values' should have the same length");
    }
    sanisizer::cast<Pointer_>(num_triplets); // check that additions and cumulative sums will not overflow the Pointer_ type.

    typedef std::vector<Pointer_> Output;
    typedef typename Output::size_type OutputSize;
    Output ptrs(sanisizer::sum<OutputSize>(num_primary, 1));
    for (const auto p : primary_indices) {
        ++ptrs[sanisizer::sum_unsafe<OutputSize>(p, 1)];
    }
    std::partial_sum(ptrs.begin(), ptrs.end(), ptrs.begin());

    if (!std::is_sorted(primary_indices.begin(), primary_indices.end())) {
        std::vector<Pointer_> copy(ptrs.begin(), ptrs.begin() + num_primary); // don't need the last element.
        auto triplet_indices = sanisizer::create<std::vector<I<decltype(num_triplets)> > >(num_triplets);
        for (I<decltype(num_triplets)> t = 0; t < num_triplets; ++t) {
            auto& offset = copy[primary_indices[t]];
            triplet_indices[offset] = t;
            ++offset;
        }

        // Reordering in-place so that triplets are sorted by primary index.
        // Ordering of secondary_indices is not yet guaranteed and will be handled later.
        for (I<decltype(num_triplets)> i = 0; i < num_triplets; ++i) {
            if (triplet_indices[i] == num_triplets) {
                continue;
            }

            const auto cur_second = secondary_indices[i];
            const auto cur_value = values[i];
            auto current = i, replacement = triplet_indices[current];

            while (replacement != i) {
                secondary_indices[current] = secondary_indices[replacement];
                values[current] = values[replacement];
                triplet_indices[current] = num_triplets;
                current = replacement;
                replacement = triplet_indices[replacement];
            }

            secondary_indices[current] = cur_second;
            values[current] = cur_value;
            triplet_indices[current] = num_triplets;
        }
    }

    std::vector<std::pair<I<decltype(secondary_indices[0])>, I<decltype(values[0])> > > sortspace;
    for (Extent_ p = 0; p < num_primary; ++p) {
        const auto start = ptrs[p];
        const auto end = ptrs[p + 1]; // + 1 is safe as p < num_primary, which fits in an Extent_.
        if (std::is_sorted(secondary_indices.begin() + start, secondary_indices.begin() + end)) {
            continue; 
        }

        sortspace.clear();
        sortspace.reserve(end - start);
        for (Pointer_ x = start; x < end; ++x) {
            sortspace.emplace_back(secondary_indices[x], values[x]);
        }
        std::sort(sortspace.begin(), sortspace.end());

        auto ssIt = sortspace.begin();
        for (Pointer_ x = start; x < end; ++x, ++ssIt) {
            secondary_indices[x] = ssIt->first;
            values[x] = ssIt->second;
        }
    }

    return ptrs;
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
