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
 * @tparam ROW Whether to create a compressed sparse row format.
 * If `false`, the compressed sparse column format is used instead.
 * @tparam U Random-access container for the values.
 * @tparam V Random access container for the row indices.
 * @tparam W Random access container for the column indices.
 *
 * @param nr Number of rows.
 * @param nc Number of columns.
 * @param rows Row indices. Values must be non-negative integers less than `nr`.
 * @param cols Column indices. Values must be non-negative integers less than `nc`.
 * @param values Non-zero values.
 * 
 * `rows`, `cols` and `values` must be of the same length.
 * Corresponding entries across these vectors are assumed to contain data for a single non-zero element.
 *
 * @return `rows`, `cols` and `values` are sorted in-place by the row and column indices (if `ROW = true`) or by the column and row indices (if `ROW = false`).
 * A vector of index pointers is returned with length `nr + 1` (if `ROW = true`) or `nc + 1` (if `ROW = false`).
 */
template <bool ROW, class U, class V, class W>
std::vector<size_t> compress_sparse_triplets(size_t nr, size_t nc, U& values, V& rows, W& cols) {
    const size_t N = rows.size();
    if (N != cols.size() || values.size() != N) { 
        throw std::runtime_error("'rows', 'cols' and 'values' should have the same length");
    }

    int order_status = 0;
    if constexpr(ROW) {
        order_status = compress_triplets::is_ordered(rows, cols);
    } else {
        order_status = compress_triplets::is_ordered(cols, rows);
    }

    if (order_status != 0) {
        std::vector<size_t> indices(N);
        for (size_t i = 0; i < N; ++i) {
            indices[i] = i;
        }

        // Sorting without duplicating the data.
        if constexpr(ROW) {
            compress_triplets::order(order_status, indices, rows, cols);
        } else {
            compress_triplets::order(order_status, indices, cols, rows);
        }

        // Reordering values in place. This (i) saves memory, and (ii) allows
        // us to work with classes U and V that may not have well-defined copy
        // constructors (e.g., if they refer to external memory).
        for (size_t i = 0; i < indices.size(); ++i) {
            if (indices[i] == -1) {
                continue;
            }

            size_t current = i, replacement = indices[i];
            indices[i] = -1;

            while (replacement != i) {
                std::swap(rows[current], rows[replacement]);
                std::swap(cols[current], cols[replacement]);
                std::swap(values[current], values[replacement]);

                current = replacement;
                auto next_replacement = indices[replacement]; 
                indices[replacement] = -1;
                replacement = next_replacement;
            } 
        }
    }

    // Collating the indices.
    std::vector<size_t> output(ROW ? nr + 1 : nc + 1);
    if constexpr(ROW) {
        for (auto t : rows) {
            ++(output[t+1]);
        } 
    } else {
        for (auto t : cols) {
            ++(output[t+1]);
        }
    }
    std::partial_sum(output.begin(), output.end(), output.begin());

    return output;
};

}

#endif
