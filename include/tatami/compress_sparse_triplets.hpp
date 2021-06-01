#ifndef TATAMI_COMPRESS_SPARSE_TRIPLETS_H
#define TATAMI_COMPRESS_SPARSE_TRIPLETS_H

#include <vector>
#include <algorithm>
#include <numeric>

/**
 * @file compress_sparse_triplets.hpp
 *
 * Convert sparse data in triplet format to a compressed sparse row/column format.
 */

namespace tatami {

/**
 * @tparam U Random-access container for the values.
 * @tparam V Random access container for the indices.
 * @tparam ROW Whether to compress into a row-based format, e.g., for `tatami::CompressedSparseRowMatrix` construction.
 * If `false`, compression to a column-based format is performed instead.
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
template <class U, class V, bool ROW = false>
std::vector<size_t> compress_sparse_triplets(size_t nr, size_t nc, U& values, V& rows, V& cols) {
    const size_t N = rows.size();
    if (N != cols.size() || values.size() != N) { 
        throw std::runtime_error("'rows', 'cols' and 'values' should have the same length");
    }

    std::vector<size_t> indices(N);
    for (size_t i = 0; i < N; ++i) {
        indices[i] = i;
    }

    // Sorting without duplicating the data.
    if constexpr(ROW) {
        std::sort(indices.begin(), indices.end(), [&](size_t left, size_t right) -> bool {
            if (rows[left] == rows[right]) {
                return (cols[left] < cols[right]);
            }
            return (rows[left] < rows[right]);
        });
    } else {
        std::sort(indices.begin(), indices.end(), [&](size_t left, size_t right) -> bool {
            if (cols[left] == cols[right]) {
                return (rows[left] < rows[right]);
            }
            return (cols[left] < cols[right]);
        });
    }

    // Reordering values in place.
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

    // Collating the indices.
    std::vector<size_t> output(ROW ? nr + 1 : nc + 1);
    const auto& target = (ROW ? rows : cols);
    for (auto t : target) {
        ++(output[t+1]);
    }
    std::partial_sum(output.begin(), output.end(), output.begin());

    return output;
};

}

#endif
