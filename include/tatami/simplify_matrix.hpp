#ifndef TATAMI_SIMPLIFY_MATRIX_H
#define TATAMI_SIMPLIFY_MATRIX_H

#include "DenseMatrix.hpp"
#include "CompressedSparseMatrix.hpp"

#include <memory>
#include <vector>
#include <deque>

/**
 * @file simplify_matrix.hpp
 *
 * Simplify a (possibly delayed) matrix into a dense or compressed sparse format.
 */

namespace tatami {

/**
 * @tparam T Type of the values in the matrix.
 * @tparam IDX Type of index values.
 *
 * @param incoming Pointer to a `tatami::typed_matrix`, possibly containing delayed operations.
 * @param row Whether the output matrix should be in a row-based format (i.e., row-major dense or compressed sparse row).
 * @param sparse Whether the output matrix should be sparse.
 *
 * @return A pointer to a new `tatami::DenseMatrix` or `tatami::CompressedSparseMatrix`, with the same dimensions and type as the matrix referenced by `incoming`.
 */
template <typename T, int IDX>
inline std::shared_ptr<typed_matrix<T, IDX> > simplify_matrix(std::shared_ptr<const typed_matrix<T, IDX> > incoming, bool row, bool sparse) {
    size_t NR = incoming->nrow();
    size_t NC = incoming->ncol();

    if (row) {
        if (sparse) {
            std::deque<T> buffer_v;
            std::deque<IDX> buffer_i;
            std::vector<size_t> indptrs(NC + 1);

            for (size_t c = 0; c < NR; ++c, optr+=NC) {
                auto range = incoming->sparse_row(c, buffer_v.data(), buffer_i.data());
                for (size_t i = 0; i < range.number; ++i) {
                    buffer_v.push_back(range.value[i]);
                    buffer_i.push_back(range.index[i]);
                }
                indptrs[c+1] = indptrs[c] + range.number;
            }
            return std::shared_ptr<typed_matrix<T, IDX>*>(new CompressedSparseRowMatrix(NR, NC, std::move(buffer_v), std::move(buffer_i), std::move(indptrs)));

        } else {
            std::vector<T> output(NR * NC);
            auto optr = output.data();
            for (size_t r = 0; r < NR; ++r, optr+=NC) {
                auto ptr = incoming->row(r, optr);
                if (ptr != optr) {
                    std::copy(ptr, ptr + NC, optr);
                }
            }
            return std::shared_ptr<typed_matrix<T, IDX>*>(new DenseRowMatrix(NR, NC, std::move(output)));
        }

    } else {
        if (sparse) {
            std::deque<T> buffer_v;
            std::deque<IDX> buffer_i;
            std::vector<size_t> indptrs(NC + 1);

            for (size_t c = 0; c < NC; ++c, optr+=NR) {
                auto range = incoming->sparse_column(c, buffer_v.data(), buffer_i.data());
                for (size_t i = 0; i < range.number; ++i) {
                    buffer_v.push_back(range.value[i]);
                    buffer_i.push_back(range.index[i]);
                }
                indptrs[c+1] = indptrs[c] + range.number;
            }
            return std::shared_ptr<typed_matrix<T, IDX>*>(new CompressedSparseColumnMatrix(NR, NC, std::move(buffer_v), std::move(buffer_i), std::move(indptrs)));

        } else {
            std::vector<T> output(NR * NC);
            auto optr = output.data();
            for (size_t c = 0; c < NC; ++c, optr+=NR) {
                auto ptr = incoming->column(c, optr);
                if (ptr != optr) {
                    std::copy(ptr, ptr + NR, optr);
                }
            }
            return std::shared_ptr<typed_matrix<T, IDX>*>(new DenseColumnMatrix(NR, NC, std::move(output)));
        }
    }
}

/**
 * This overload automatically sets the output row-format and sparsity (`sparse`) based on `tatami::matrix::prefer_rows()` and `tatami::matrix::sparse()`, respectively.
 *
 * @tparam T Type of the values in the matrix.
 * @tparam IDX Type of index values.
 *
 * @param incoming Pointer to a `tatami::typed_matrix`, possibly containing delayed operations.
 *
 * @return A pointer to a new `tatami::DenseMatrix` or `tatami::CompressedSparseMatrix`, with the same dimensions and type as the matrix referenced by `incoming`.
 */
template <typename T, int IDX>
inline std::shared_ptr<typed_matrix<T, IDX>*> simplify_matrix(std::shared_ptr<const typed_matrix<T, IDX>*> incoming) {
    return simplify_matrix(incoming, incoming->prefer_rows(), incoming->sparse());
}

}

#endif
