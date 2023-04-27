#ifndef TATAMI_CONVERT_TO_DENSE_H
#define TATAMI_CONVERT_TO_DENSE_H

#include "../base/utils.hpp"
#include "../base/dense/DenseMatrix.hpp"

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
 * @tparam row_ Whether to return a row-major matrix.
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam Matrix_ Input matrix class, most typically a `tatami::Matrix`.
 *
 * @param incoming Pointer to a `tatami::Matrix`.
 * @param[out] store Pointer to an array of length equal to the product of the dimensions of `incoming`.
 * On output, this is filled with values from `incoming` in row- or column-major format depending on `row_`.
 */
template <bool row_, typename StoredValue_, class Matrix_>
void convert_to_dense(const Matrix_* incoming, StoredValue_* store) {
    typedef typename Matrix_::index_type Index_;
    typedef typename Matrix_::value_type Value_;

    size_t NR = incoming->nrow();
    size_t NC = incoming->ncol();
    size_t primary = (row_ ? NR : NC);
    size_t secondary = (row_ ? NC : NR);

    if (row_ == incoming->prefer_rows()) {
        auto wrk = new_extractor<row_, false, Value_, Index_>(incoming);
        constexpr bool same_type = std::is_same<Value_, StoredValue_>::value;
        std::vector<Value_> temp(same_type ? 0 : secondary);

        for (size_t p = 0; p < primary; ++p, store += secondary) {
            if constexpr(same_type) {
                wrk->fetch_copy(p, store);
            } else {
                auto ptr = wrk->fetch(p, temp.data());
                std::copy(ptr, ptr + secondary, store);
            }
        }

    } else {
        // We iterate on the incoming matrix's preferred dimension,
        // under the assumption that it may be arbitrarily costly to
        // extract in the non-preferred dim; it is thus cheaper to
        // do cache-unfriendly inserts into the output buffer.
        auto wrk = new_extractor<!row_, false, Value_, Index_>(incoming);
        std::vector<Value_> temp(primary);

        for (size_t s = 0; s < secondary; ++s) {
            auto ptr = wrk->fetch(s, temp.data());
            auto bptr = store + s;
            for (size_t p = 0; p < primary; ++p, bptr += secondary) {
                *bptr = ptr[p]; 
            }
        }
    }

    return;
}

/**
 * @tparam row Whether to return a row-major matrix.
 * @tparam Value_ Type of data values in the output interface.
 * @tparam Index Integer type for the indices in the output interface.
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam Matrix_ Input matrix class, most typically a `tatami::Matrix`.
 *
 * @param incoming Pointer to a `tatami::Matrix`.
 *
 * @return A pointer to a new `tatami::DenseMatrix` with the same dimensions and type as the matrix referenced by `incoming`.
 * If `row = true`, the matrix is row-major, otherwise it is column-major.
 */
template <
    bool row_, 
    typename Value_ = double, 
    typename Index = int, 
    typename StoredValue_ = Value_, 
    class Matrix_
>
inline std::shared_ptr<Matrix<Value_, Index> > convert_to_dense(const Matrix_* incoming) {
    size_t NR = incoming->nrow();
    size_t NC = incoming->ncol();
    std::vector<StoredValue_> buffer(NR * NC);
    convert_to_dense<row_>(incoming, buffer.data());
    return std::shared_ptr<Matrix<Value_, Index> >(new DenseMatrix<row_, Value_, Index, decltype(buffer)>(NR, NC, std::move(buffer)));
}

/**
 * This overload makes it easier to control the desired output order when it is not known at compile time.
 *
 * @tparam Value_ Type of data values in the output interface.
 * @tparam Index Integer type for the indices in the output interface.
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam Matrix_ Input matrix class, most typically a `tatami::Matrix`.
 *
 * @param incoming Pointer to a `tatami::Matrix`.
 * @param order Ordering of values in the output dense matrix - row-major (0) or column-major (1).
 * If set to -1, the ordering is chosen based on `tatami::Matrix::prefer_rows()`. 
 *
 * @return A pointer to a new `tatami::DenseMatrix` with the same dimensions and type as the matrix referenced by `incoming`.
 */
template <
    typename Value_ = double,
    typename Index_ = int,
    typename StoredValue_ = Value_,
    class Matrix_
>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_dense(const Matrix_* incoming, int order) {
    if (order < 0) {
        order = static_cast<int>(!incoming->prefer_rows()); 
    }
    if (order == 0) {
        return convert_to_dense<true, Value_, Index_, StoredValue_, Matrix_>(incoming);
    } else {
        return convert_to_dense<false, Value_, Index_, StoredValue_, Matrix_>(incoming);
    }
}

}

#endif
