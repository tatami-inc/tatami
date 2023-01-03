#ifndef TATAMI_CONVERT_TO_DENSE_H
#define TATAMI_CONVERT_TO_DENSE_H

#include "../base/DenseMatrix.hpp"

#include <memory>
#include <vector>
#include <deque>

/**
 * @file convert_to_dense.hpp
 *
 * @brief Convert a matrix into a dense format.
 */

namespace tatami {

/**
 * @tparam row Whether to return a row-major matrix.
 * @tparam DataInterface Type of data values in the output interface.
 * @tparam Index Integer type for the indices in the output interface.
 * @tparam DataStore Type of data values to be stored in the output.
 * @tparam MatrixIn Input matrix class, most typically a `tatami::Matrix`.
 *
 * @param incoming Pointer to a `tatami::Matrix`.
 *
 * @return A pointer to a new `tatami::DenseMatrix` with the same dimensions and type as the matrix referenced by `incoming`.
 * If `row = true`, the matrix is row-major, otherwise it is column-major.
 */
template <bool row, 
    typename DataInterface = double,
    typename Index = int,
    typename DataStore = DataInterface,
    class MatrixIn
>
inline std::shared_ptr<Matrix<DataInterface, Index> > convert_to_dense(const MatrixIn* incoming) {
    size_t NR = incoming->nrow();
    size_t NC = incoming->ncol();
    size_t primary = (row ? NR : NC);
    size_t secondary = (row ? NC : NR);

    typedef typename MatrixIn::data_type DataIn;
    std::vector<DataStore> buffer(NR * NC);
    std::vector<DataIn> temp(secondary);

    if (row == incoming->prefer_rows()) {
        auto bptr = buffer.data();
        auto wrk = incoming->new_workspace(row);
        
        for (size_t p = 0; p < primary; ++p, bptr += secondary) {
            if constexpr(std::is_same<DataIn, DataStore>::value) {
                if constexpr(row) {
                    incoming->row_copy(p, bptr, wrk.get());
                } else {
                    incoming->column_copy(p, bptr, wrk.get());
                }
            } else {
                const DataIn* ptr;
                if constexpr(row) {
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
        auto wrk = incoming->new_workspace(!row);
        for (size_t s = 0; s < secondary; ++s) {
            const DataIn* ptr;
            if constexpr(row) {
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

    return std::shared_ptr<Matrix<DataInterface, Index> >(new DenseMatrix<row, DataInterface, Index, decltype(buffer)>(NR, NC, std::move(buffer)));
}

/**
 * This overload makes it easier to control the desired output order when it is not known at compile time.
 *
 * @tparam DataInterface Type of data values in the output interface.
 * @tparam DataStore Type of data values to be stored in the output.
 * @tparam Index Integer type for the indices in the output interface.
 * @tparam MatrixIn Input matrix class, most typically a `tatami::Matrix`.
 *
 * @param incoming Pointer to a `tatami::Matrix`.
 * @param order Ordering of values in the output dense matrix - row-major (0) or column-major (1).
 * If set to -1, the ordering is chosen based on `tatami::Matrix::prefer_rows()`. 
 *
 * @return A pointer to a new `tatami::DenseMatrix` with the same dimensions and type as the matrix referenced by `incoming`.
 */
template <
    typename DataInterface = double,
    typename Index = int,
    typename DataStore = DataInterface,
    class MatrixIn
>
std::shared_ptr<Matrix<DataInterface, Index> > convert_to_dense(const MatrixIn* incoming, int order) {
    if (order < 0) {
        order = static_cast<int>(!incoming->prefer_rows()); 
    }
    if (order == 0) {
        return convert_to_dense<true, DataInterface, Index, DataStore, MatrixIn>(incoming);
    } else {
        return convert_to_dense<false, DataInterface, Index, DataStore, MatrixIn>(incoming);
    }
}

}

#endif
