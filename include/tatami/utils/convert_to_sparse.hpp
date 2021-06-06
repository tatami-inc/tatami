#ifndef TATAMI_CONVERT_TO_SPARSE_H
#define TATAMI_CONVERT_TO_SPARSE_H

#include "../base/CompressedSparseMatrix.hpp"

#include <memory>
#include <vector>
#include <deque>

/**
 * @file convert_to_sparse.hpp
 *
 * Convert a matrix into a compressed sparse format.
 */

namespace tatami {

/**
 * @tparam T Type of the values in the matrix.
 * @tparam IDX Type of index values.
 *
 * @param incoming Pointer to a `tatami::typed_matrix`, possibly containing delayed operations.
 * @param row Whether the output matrix should be in compressed sparse row format.
 *
 * @return A pointer to a new `tatami::CompressedSparseMatrix`, with the same dimensions and type as the matrix referenced by `incoming`.
 */
template <typename T, typename IDX>
inline std::shared_ptr<typed_matrix<T, IDX> > convert_to_sparse(const typed_matrix<T, IDX>* incoming, bool row, double reserve = 0.2) {
    size_t NR = incoming->nrow();
    size_t NC = incoming->ncol();
    size_t reservation = static_cast<double>(NR * NC) * reserve;

    std::vector<T> output_v;
    output_v.reserve(reservation);
    std::vector<IDX> output_i;
    output_i.reserve(reservation);

    if (row) {
        auto wrk = incoming->new_workspace(true);
        std::vector<size_t> indptrs(NR + 1);
        std::vector<T> buffer_v(NC);

        if (incoming->sparse()) {
            std::vector<IDX> buffer_i(NC);
            for (size_t r = 0; r < NR; ++r) {
                auto range = incoming->sparse_row(r, buffer_v.data(), buffer_i.data(), wrk.get());
                for (size_t i = 0; i < range.number; ++i, ++range.value, ++range.index) {
                    if (*range.value) {
                        output_v.push_back(*range.value);
                        output_i.push_back(*range.index);
                    }
                }
                indptrs[r+1] = output_v.size();
            }
        } else {
            for (size_t r = 0; r < NR; ++r) {
                auto ptr = incoming->row(r, buffer_v.data(), wrk.get());
                for (size_t c = 0; c < NC; ++c, ++ptr) {
                    if (*ptr) {
                        output_v.push_back(*ptr);
                        output_i.push_back(c);
                    }
                }
                indptrs[r+1] = output_v.size();
            }
        }

        return std::shared_ptr<typed_matrix<T, IDX> >(new CompressedSparseRowMatrix<T, IDX>(NR, NC, std::move(output_v), std::move(output_i), std::move(indptrs)));

    } else {
        auto wrk = incoming->new_workspace(false);
        std::vector<size_t> indptrs(NC + 1);
        std::vector<T> buffer_v(NR);

        if (incoming->sparse()) {
            std::vector<IDX> buffer_i(NR);
            for (size_t c = 0; c < NC; ++c) {
                auto range = incoming->sparse_column(c, buffer_v.data(), buffer_i.data(), wrk.get());
                for (size_t i = 0; i < range.number; ++i, ++range.value, ++range.index) {
                    if (*range.value) {
                        output_v.push_back(*range.value);
                        output_i.push_back(*range.index);
                    }
                }
                indptrs[c+1] = output_v.size();
            }

        } else {
            for (size_t c = 0; c < NC; ++c) {
                auto ptr = incoming->column(c, buffer_v.data(), wrk.get());
                for (size_t r = 0; r < NR; ++r, ++ptr) {
                    if (*ptr) {
                        output_v.push_back(*ptr);
                        output_i.push_back(r);
                    }
                }
                indptrs[c+1] = output_v.size();
            }
        }

        return std::shared_ptr<typed_matrix<T, IDX> >(new CompressedSparseColumnMatrix<T, IDX>(NR, NC, std::move(output_v), std::move(output_i), std::move(indptrs)));
    }
}

/**
 * This overload automatically sets the output row format based on `tatami::matrix::prefer_rows()`.
 *
 * @tparam T Type of the values in the matrix.
 * @tparam IDX Type of index values.
 *
 * @param incoming Pointer to a `tatami::typed_matrix`, possibly containing delayed operations.
 *
 * @return A pointer to a new `tatami::CompressedSparseMatrix`, with the same dimensions and type as the matrix referenced by `incoming`.
 */
template <typename T, typename IDX>
inline std::shared_ptr<typed_matrix<T, IDX> > convert_to_sparse(const typed_matrix<T, IDX>* incoming) {
    return convert_to_sparse(incoming, incoming->prefer_rows(), incoming->sparse());
}

}

#endif
