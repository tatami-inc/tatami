#ifndef TATAMI_COMPRESS_SPARSE_TRIPLETS_H
#define TATAMI_COMPRESS_SPARSE_TRIPLETS_H

#include <vector>
#include <algorithm>
#include <numeric>

/**
 * @file compress_sparse_triplets.hpp
 *
 * @brief Convert sparse data in triplet format to a compressed sparse row/column format.
 */

namespace tatami {

namespace compress_triplets {

/**
 * @cond
 */
template<class Primary, class Secondary>
int is_ordered(const Primary& primary, const Secondary& secondary) {
    if (!std::is_sorted(primary.begin(), primary.end())) {
        return 2;
    }

    size_t start = 0;
    while (start < primary.size()) {
        size_t end = start + 1;
        while (end < primary.size() && primary[end] == primary[start]) {
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

template<class Primary, class Secondary>
void order(int status, std::vector<size_t>& indices, const Primary& primary, const Secondary& secondary) {
    if (status == 1) {
        size_t start = 0;
        while (start < primary.size()) {
            size_t end = start + 1;
            while (end < primary.size() && primary[end] == primary[start]) {
                ++end;
            }

            // Checking if this particular run can be skipped.
            if (!std::is_sorted(secondary.begin() + start, secondary.begin() + end)) {
                std::sort(indices.begin() + start, indices.begin() + end, [&](size_t left, size_t right) -> bool { 
                    return secondary[left] < secondary[right];
                });
            }
            start = end;
        }

    } else if (status == 2) {
        std::sort(indices.begin(), indices.end(), [&](size_t left, size_t right) -> bool {
            if (primary[left] == primary[right]) {
                return (secondary[left] < secondary[right]);
            }
            return (primary[left] < primary[right]);
        });
    }
}
/**
 * @endcond
 */

}

/**
 * @tparam Values_ Random-access container for the values.
 * @tparam RowIndices_ Random access container for the row indices.
 * @tparam ColumnIndices_ Random access container for the column indices.
 *
 * @param nrow Number of rows.
 * @param ncol Number of columns.
 * @param row_indices Row indices of the structural non-zeros.
 * Values must be non-negative integers less than `nrow`.
 * @param column_indices Column indices of the structural non-zeros.
 * This must be of the same length as `row_indices`, where corresponding entries contain data for a single structural non-zero.
 * Values must be non-negative integers less than `ncol`.
 * @param values Values of the structural non-zeros.
 * This must be of the same length as `row_indices` and `column_indices`, where corresponding entries contain data for a single structural non-zero.
 * @param csr Whether to create a compressed sparse row format.
 * If `false`, the compressed sparse column format is used instead.
 * 
 * @return `row_indices`, `column_indices` and `values` are sorted in-place by the row and column indices (if `csr = true`) or by the column and row indices (if `csr = false`).
 * A vector of index pointers is returned with length `nrow + 1` (if `csr = true`) or `ncol + 1` (if `csr = false`).
 */
template<class Values_, class RowIndices_, class ColumnIndices_>
std::vector<size_t> compress_sparse_triplets(size_t nrow, size_t ncol, Values_& values, RowIndices_& row_indices, ColumnIndices_& column_indices, bool csr) {
    const size_t N = row_indices.size();
    if (N != column_indices.size() || values.size() != N) { 
        throw std::runtime_error("'row_indices', 'column_indices' and 'values' should have the same length");
    }

    int order_status = 0;
    if (csr) {
        order_status = compress_triplets::is_ordered(row_indices, column_indices);
    } else {
        order_status = compress_triplets::is_ordered(column_indices, row_indices);
    }

    if (order_status != 0) {
        std::vector<size_t> indices(N);
        for (size_t i = 0; i < N; ++i) {
            indices[i] = i;
        }

        // Sorting without duplicating the data.
        if (csr) {
            compress_triplets::order(order_status, indices, row_indices, column_indices);
        } else {
            compress_triplets::order(order_status, indices, column_indices, row_indices);
        }

        // Reordering values in place. This (i) saves memory, and (ii) allows
        // us to work with Values_, RowIndices_, etc. that may not have well-defined copy
        // constructors (e.g., if they refer to external memory).
        for (size_t i = 0; i < indices.size(); ++i) {
            if (indices[i] == static_cast<size_t>(-1)) {
                continue;
            }

            size_t current = i, replacement = indices[i];
            indices[i] = -1;

            while (replacement != i) {
                std::swap(row_indices[current], row_indices[replacement]);
                std::swap(column_indices[current], column_indices[replacement]);
                std::swap(values[current], values[replacement]);

                current = replacement;
                auto next_replacement = indices[replacement]; 
                indices[replacement] = -1;
                replacement = next_replacement;
            } 
        }
    }

    // Collating the indices.
    std::vector<size_t> output(csr ? nrow + 1 : ncol + 1);
    if (csr) {
        for (auto t : row_indices) {
            ++(output[t+1]);
        } 
    } else {
        for (auto t : column_indices) {
            ++(output[t+1]);
        }
    }
    std::partial_sum(output.begin(), output.end(), output.begin());

    return output;
}

/**
 * @cond
 */
// Back-compatibility.
template<bool row_, class Values_, class RowIndices_, class ColumnIndices_>
std::vector<size_t> compress_sparse_triplets(size_t nrow, size_t ncol, Values_& values, RowIndices_& row_indices, ColumnIndices_& column_indices) {
    return compress_sparse_triplets(nrow, ncol, values, row_indices, column_indices, row_);
}
/**
 * @endcond
 */

}

#endif
