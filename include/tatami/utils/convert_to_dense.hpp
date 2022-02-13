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
 * @tparam row_ Whether to return a row-major matrix.
 * @tparam MatrixIn Input matrix class, most typically a `tatami::Matrix`.
 * @tparam Data Type of data values in the output interface.
 * @tparam Index Integer type for the indices in the output interface.
 *
 * @param incoming Pointer to a `tatami::Matrix`.
 *
 * @return A pointer to a new `tatami::DenseMatrix` with the same dimensions and type as the matrix referenced by `incoming`.
 * If `row = true`, the matrix is row-major, otherwise it is column-major.
 */
template <bool row_, class MatrixIn, typename DataOut = typename MatrixIn::data_type, typename IndexOut = typename MatrixIn::index_type>
inline std::shared_ptr<Matrix<DataOut, IndexOut> > convert_to_dense(const MatrixIn* incoming) {
    size_t NR = incoming->nrow();
    size_t NC = incoming->ncol();
    size_t primary = (row_ ? NR : NC);
    size_t secondary = (row_ ? NC : NR);

    typedef typename MatrixIn::data_type DataIn;
    std::vector<DataOut> buffer(NR * NC);
    std::vector<DataIn> temp(secondary);

    if (row_ == incoming->prefer_rows()) {
        auto bptr = buffer.data();
        auto wrk = incoming->new_workspace(row_);
        
        for (size_t p = 0; p < primary; ++p, bptr += secondary) {
            if constexpr(std::is_same<DataIn, DataOut>::value) {
                if constexpr(row_) {
                    incoming->row_copy(p, bptr, wrk.get());
                } else {
                    incoming->column_copy(p, bptr, wrk.get());
                }
            } else {
                const DataIn* ptr;
                if constexpr(row_) {
                    ptr = incoming->row(p, temp.data(), wrk.get());
                } else {
                    ptr = incoming->column(p, temp.data(), wrk.get());
                }
                std::copy(ptr, ptr + secondary, bptr);
            }
        }

    } else {
        // We iterate on the incoming matrix's preferred dimension,
        // under the assumption that it may be arbitrarily costly to
        // extract in the non-preferred dim; it is thus cheaper to
        // do cache-unfriendly inserts into the output buffer.
        auto wrk = incoming->new_workspace(!row_);
        for (size_t s = 0; s < secondary; ++s) {
            const DataIn* ptr;
            if constexpr(row_) {
                ptr = incoming->column(s, temp.data(), wrk.get());
            } else {
                ptr = incoming->row(s, temp.data(), wrk.get());
            }

            auto bptr = buffer.begin() + s;
            for (size_t p = 0; p < primary; ++p, bptr += secondary) {
                *bptr = ptr[p]; 
            }
        }
    }

    return std::shared_ptr<Matrix<DataOut, IndexOut> >(new DenseMatrix<row_, DataOut, IndexOut>(NR, NC, std::move(buffer)));
}

/**
 * This overload makes it easier to control the desired output order when it is not known at compile time.
 *
 * @tparam MatrixIn Input matrix class, most typically a `tatami::Matrix`.
 * @tparam Data Type of data values in the output interface.
 * @tparam Index Integer type for the indices in the output interface.
 *
 * @param incoming Pointer to a `tatami::Matrix`.
 * @param order Ordering of values in the output dense matrix - row-major (0) or column-major (1).
 * If set to -1, the ordering is chosen based on `tatami::Matrix::prefer_rows()`. 
 *
 * @return A pointer to a new `tatami::DenseMatrix` with the same dimensions and type as the matrix referenced by `incoming`.
 */
template <class MatrixIn, typename DataOut = typename MatrixIn::data_type, typename IndexOut = typename MatrixIn::index_type>
std::shared_ptr<Matrix<DataOut, IndexOut> > convert_to_dense(const MatrixIn* incoming, int order) {
    if (order < 0) {
        order = static_cast<int>(!incoming->prefer_rows()); 
    }
    if (order == 0) {
        return convert_to_dense<true, MatrixIn, DataOut, IndexOut>(incoming);
    } else {
        return convert_to_dense<false, MatrixIn, DataOut, IndexOut>(incoming);
    }
}

}

#endif
