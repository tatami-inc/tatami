#ifndef TATAMI_DIMSUM_HPP
#define TATAMI_DIMSUM_HPP

#include "../base/typed_matrix.hpp"
#include <vector>
#include <algorithm>

/**
 * @file sums.hpp
 *
 * Compute row and column sums from a `tatami::typed_matrix`.
 */

namespace tatami {

/**
 * @tparam ROW Whether to compute the row sums.
 * If false, computes the column sums instead.
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::typed_matrix`.
 *
 * @return A vector of length equal to the number of rows (if `ROW = true`) or columns (if `ROW = false`), containing the row or column sums respectively.
 */
template<bool ROW, typename T, typename IDX>
inline std::vector<T> sums(const typed_matrix<T, IDX>* p) {
    size_t NR = p->nrow(), NC = p->ncol();
    std::vector<T> obuffer(ROW ? NR : NC);
    auto output = obuffer.begin();

    // Deciding whether or not to perform row-wise extraction.
    if (p->prefer_rows()) {
        std::vector<T> buffer(NC);
        auto wrk = p->new_workspace(true);

        if (p->sparse()) {
            std::vector<IDX> ibuffer(NC);
            for (size_t i = 0; i < NR; ++i) {
                auto range = p->sparse_row(i, buffer.data(), ibuffer.data(), wrk.get());
                if constexpr(ROW) {
                    *output = std::accumulate(range.value, range.value + range.number, static_cast<T>(0));
                    ++output;
                } else {
                    for (size_t j = 0; j < range.number; ++j, ++range.index, ++range.value) {
                        *(output + *range.index) += *range.value;
                    }
                }
            }
        } else {
            for (size_t i = 0; i < NR; ++i) {
                auto ptr = p->row(i, buffer.data(), wrk.get());
                if constexpr(ROW) {
                    *output = std::accumulate(ptr, ptr + NC, static_cast<T>(0));
                    ++output;
                } else {
                    auto copy = output;
                    for (size_t j = 0; j < NC; ++j, ++copy, ++ptr) {
                        *copy += *ptr;
                    }
                }
            }
        }
    } else {
        std::vector<T> buffer(NR);
        auto wrk = p->new_workspace(false);

        if (p->sparse()) {
            std::vector<IDX> ibuffer(NR);
            for (size_t i = 0; i < NC; ++i) {
                auto range = p->sparse_column(i, buffer.data(), ibuffer.data(), wrk.get());
                if constexpr(ROW) {
                    for (size_t j = 0; j < range.number; ++j, ++range.index, ++range.value) {
                        *(output + *range.index) += *range.value;
                    }
                } else {
                    *output = std::accumulate(range.value, range.value + range.number, static_cast<T>(0));
                    ++output;
                }
            }
        } else {
            for (size_t i = 0; i < NC; ++i) {
                auto ptr = p->column(i, buffer.data(), wrk.get());
                if constexpr(ROW) {
                    auto copy = output;
                    for (size_t j = 0; j < NR; ++j, ++copy, ++ptr) {
                        *copy += *ptr;
                    }
                } else {
                    *output = std::accumulate(ptr, ptr + NR, static_cast<T>(0));
                    ++output;
                }
            }
        }
    }

    return obuffer;
}

/**
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::typed_matrix`.
 *
 * @return A vector of length equal to the number of columns, containing the column sums.
 */
template<typename T, typename IDX>
inline std::vector<T> column_sums(const typed_matrix<T, IDX>* p) {
    return sums<false, T, IDX>(p);
}

/**
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::typed_matrix`.
 *
 * @return A vector of length equal to the number of rows, containing the row sums.
 */
template<typename T, typename IDX>
inline std::vector<T> row_sums(const typed_matrix<T, IDX>* p) {
    return sums<true, T, IDX>(p);
}

}

#endif
