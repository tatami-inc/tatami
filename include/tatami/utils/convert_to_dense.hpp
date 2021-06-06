#ifndef TATAMI_CONVERT_TO_DENSE_H
#define TATAMI_CONVERT_TO_DENSE_H

#include "../base/DenseMatrix.hpp"

#include <memory>
#include <vector>
#include <deque>

/**
 * @file convert_to_dense.hpp
 *
 * Convert a matrix into a dense format.
 */

namespace tatami {

/**
 * @tparam T Type of the values in the matrix.
 * @tparam IDX Type of index values.
 *
 * @param incoming Pointer to a `tatami::typed_matrix`.
 * @param row Whether the output matrix should be row-major.
 *
 * @return A pointer to a new `tatami::DenseMatrix` with the same dimensions and type as the matrix referenced by `incoming`.
 */
template <typename T, typename IDX>
inline std::shared_ptr<typed_matrix<T, IDX> > convert_to_dense(const typed_matrix<T, IDX>* incoming, bool row) {
    size_t NR = incoming->nrow();
    size_t NC = incoming->ncol();
    std::vector<T> output(NR * NC);
    auto optr = output.data();

    if (row) {
        auto wrk = incoming->new_workspace(true);
        for (size_t r = 0; r < NR; ++r, optr+=NC) {
            auto ptr = incoming->row(r, optr, wrk.get());
            if (ptr != optr) {
                std::copy(ptr, ptr + NC, optr);
            }
        }
        return std::shared_ptr<typed_matrix<T, IDX> >(new DenseRowMatrix<T, decltype(output), IDX>(NR, NC, std::move(output)));

    } else {
        auto wrk = incoming->new_workspace(false);
        for (size_t c = 0; c < NC; ++c, optr+=NR) {
            auto ptr = incoming->column(c, optr, wrk.get());
            if (ptr != optr) {
                std::copy(ptr, ptr + NR, optr);
            }
        }
        return std::shared_ptr<typed_matrix<T, IDX> >(new DenseColumnMatrix<T, decltype(output), IDX>(NR, NC, std::move(output)));
    }
}

/**
 * This overload automatically sets the output row-format based on `tatami::matrix::prefer_rows()`. 
 *
 * @tparam T Type of the values in the matrix.
 * @tparam IDX Type of index values.
 *
 * @param incoming Pointer to a `tatami::typed_matrix`.
 *
 * @return A pointer to a new `tatami::DenseMatrix` with the same dimensions and type as the matrix referenced by `incoming`.
 */
template <typename T, typename IDX>
inline std::shared_ptr<typed_matrix<T, IDX> > convert_to_dense(const typed_matrix<T, IDX>* incoming) {
    return convert_to_dense(incoming, incoming->prefer_rows());
}

}

#endif
