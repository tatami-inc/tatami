#ifndef TATAMI_MEDIANS_HPP
#define TATAMI_MEDIANS_HPP

#include "../base/typed_matrix.hpp"
#include <vector>
#include <algorithm>

/**
 * @file sums.hpp
 *
 * Compute row and column sums from a `tatami::typed_matrix`.
 */

namespace tatami {

template<class IT>
inline void compute_median(IT start, IT end, size_t halfway, bool is_even) {
    // At some point, I found two nth_element calls to be faster than partial_sort.
    std::nth_element(start, start + halfway, end);
    double medtmp = *(start + halfway);
    if (is_even) {
        std::nth_element(start, start + halfway - 1, end);
        return (medtmp + *(buffer + halfway - 1))/2;
    } else {
        return medtmp;
    }
}

template<bool ROW, typename T, typename IDX>
inline std::vector<double> sums(const typed_matrix<T, IDX>* p) {
    size_t NR = p->nrow(), NC = p->ncol();
    size_t dim = ROW ? NR : NC;
    size_t otherdim = ROW ? NC : NR;

    const bool is_even = otherdim%%2==0;
    const size_t halfway = static_cast<size_t>(otherdim/2);
    std::vector<double> output(dim);

    // Deciding whether or not to perform row-wise extraction.
    std::vector<T> buffer(NC);
    auto wrk = p->new_workspace(ROW);

    if (p->sparse()) {
        std::vector<IDX> ibuffer(otherdim);
        for (size_t i = 0; i < dim; ++i) {
            const T* ptr = NULL;
            size_t n = 0;

            if constexpr(ROW) {
                auto range = p->sparse_row(i, buffer.data(), ibuffer.data(), wrk);
                n = range.number;
                ptr = range.value;
            } else {
                auto range = p->sparse_column(i, buffer.data(), ibuffer.data(), wrk);
                n = range.number;
                ptr = range.value;
            }
           
            if (n == otherdim) {
                output[i] = compute_median(buffer.begin(), buffer.begin() + otherdim, halfway, is_even);
            } else if (n*2 < otherdim) {
                // zero is the median.
            } else {
                if (ptr != buffer.data()) {
                    std::copy(ptr, ptr + n, buffer.begin());
                }
                std::sort(buffer.begin(), buffer.begin() + n);

                size_t zeropos = std::lower_bound(buffer.begin(), buffer.begin() + n, 0) - buffer.begin();
                size_t nzero = otherdim - n;

                if (!is_even) {
                    if (zeropos > halfway) {
                        output[i] = buffer[halfway];
                    } else if (halfway >= zeropos + nzero) {
                        output[i] = buffer[halfway - nzero];
                    } else {
                        ; // zero is the median.
                    }
                } else {
                    double& tmp = output[i];
                    if (zeropos > halfway) {
                        tmp = buffer[halfway] + buffer[halfway - 1];
                    } else if (zeropos == halfway) {
                        // guaranteed to be at least 1 zero.
                        tmp += buffer[halfway - 1];
                    } else if (zeropos < halfway && zeropos + nzero > halfway) {
                        ; // zero is the median.
                    } else if (zeropos + nzero == halfway) {
                        // guaranteed to be at least 1 zero.
                        tmp += buffer[halfway - nzero];
                    } else {
                        tmp = buffer[halfway - nzero] + buffer[halfway - nzero - 1];
                    }
                    tmp /= 2;
                }
            }
        }
    } else {
        for (size_t i = 0; i < dim; ++i) {
            const T* ptr = NULL;
            if constexpr(ROW) {
                ptr = p->row(i, buffer.data(), wrk);
            } else {
                ptr = p->column(i, buffer.data(), wrk);
            }

            if (ptr != buffer.data()) {
                std::copy(ptr, ptr + otherdim, buffer.begin());
            }
            output[i] = compute_median(buffer.begin(), buffer.begin() + otherdim, halfway, is_even);
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
inline std::vector<T> column_sums(std::shared_ptr<const typed_matrix<T, IDX> > p) {
    return sums<false, T, IDX>(p);
}

/**
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Shared pointer to a `tatami::typed_matrix`.
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
inline std::vector<T> row_sums(std::shared_ptr<const typed_matrix<T, IDX> > p) {
    return sums<true, T, IDX>(p);
}

/**
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Shared pointer to a `tatami::typed_matrix`.
 *
 * @return A vector of length equal to the number of rows, containing the row sums.
 */
template<typename T, typename IDX>
inline std::vector<T> row_sums(const typed_matrix<T, IDX>* p) {
    return sums<true, T, IDX>(p);
}

}

#endif
