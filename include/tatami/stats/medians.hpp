#ifndef TATAMI_STATS_MEDIANS_HPP
#define TATAMI_STATS_MEDIANS_HPP

#include "../base/Matrix.hpp"
#include "apply.hpp"

#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

/**
 * @file medians.hpp
 *
 * Compute row and column medians from a `tatami::Matrix`.
 */

namespace tatami {

namespace stats {

template<typename O = double, typename T>
O compute_median(const T* src, T* buffer, size_t n) {
    if (n == 0) {
        return std::numeric_limits<O>::quiet_NaN();
    }

    if (src != buffer) {
        std::copy(src, src + n, buffer);
    }

    size_t halfway = n / 2;
    bool is_even = (n % 2 == 0);

    // At some point, I found two nth_element calls to be faster than partial_sort.
    std::nth_element(buffer, buffer + halfway, buffer + n);
    double medtmp = *(buffer + halfway);
    if (is_even) {
        std::nth_element(buffer, buffer + halfway - 1, buffer + n);
        return (medtmp + *(buffer + halfway - 1))/2;
    } else {
        return medtmp;
    }
}

template<typename O>
struct MedianFactory {
    MedianFactory(O* o, size_t d2) : output(o), otherdim(d2) {}

    static const bool supports_sparse = true;

    static const bool supports_running = false;
private:
    O* output;
    size_t otherdim;

public:
    struct Dense {
        Dense(O* o, size_t d2) : output(o), otherdim(d2) {}

        template<typename T>
        void compute(size_t i, const T* ptr, T* buffer) {
            output[i] = compute_median<O>(ptr, buffer, otherdim);
        }
    private:
        O* output;
        size_t otherdim;
    };

    Dense dense_direct() {
        return Dense(output, otherdim);
    }

public:
    struct Sparse {
        Sparse(O* o, size_t d2) : output(o), otherdim(d2) {}

        template<typename T, typename IDX>
        void compute(size_t i, const SparseRange<T, IDX>& range, T* vbuffer, IDX* ibuffer) {
            if (range.number == otherdim) {
                output[i] = compute_median<O>(range.value, vbuffer, otherdim);
            } else if (range.number * 2 < otherdim) {
                output[i] = 0; // zero is the median if there are too many zeroes.
            } else {
                if (range.value != vbuffer) {
                    std::copy(range.value, range.value + range.number, vbuffer);
                }

                size_t halfway = otherdim / 2;
                bool is_even = (otherdim % 2 == 0);

                auto vend = vbuffer + range.number;
                std::sort(vbuffer, vend);
                size_t zeropos = std::lower_bound(vbuffer, vend, 0) - vbuffer;
                size_t nzero = otherdim - range.number;

                if (!is_even) {
                    if (zeropos > halfway) {
                        output[i] = vbuffer[halfway];
                    } else if (halfway >= zeropos + nzero) {
                        output[i] = vbuffer[halfway - nzero];
                    } else {
                        output[i] = 0; // zero is the median.
                    }
                } else {
                    double tmp = 0;
                    if (zeropos > halfway) {
                        tmp = vbuffer[halfway] + vbuffer[halfway - 1];
                    } else if (zeropos == halfway) {
                        // guaranteed to be at least 1 zero.
                        tmp += vbuffer[halfway - 1];
                    } else if (zeropos < halfway && zeropos + nzero > halfway) {
                        ; // zero is the median.
                    } else if (zeropos + nzero == halfway) {
                        // guaranteed to be at least 1 zero.
                        tmp += vbuffer[halfway - nzero];
                    } else {
                        tmp = vbuffer[halfway - nzero] + vbuffer[halfway - nzero - 1];
                    }
                    output[i] = tmp / 2;
                }
            }
        }
    private:
        O* output;
        size_t otherdim;
    };

    Sparse sparse_direct() {
        return Sparse(output, otherdim);
    }
};

}

/**
 * @tparam Output Type of the output.
 * @tparam T Type of the matrix value.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Shared pointer to a `tatami::Matrix`.
 *
 * @return A vector of length equal to the number of columns, containing the column medians.
 */
template<typename Output = double, typename T, typename IDX>
inline std::vector<Output> column_medians(const Matrix<T, IDX>* p) {
    std::vector<Output> output(p->ncol());
    stats::MedianFactory factory(output.data(), p->nrow());
    apply<1>(p, factory);
    return output;
}

/**
 * @tparam Output Type of the output.
 * @tparam T Type of the matrix value.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Shared pointer to a `tatami::Matrix`.
 *
 * @return A vector of length equal to the number of rows, containing the row medians.
 */
template<typename Output = double, typename T, typename IDX>
inline std::vector<Output> row_medians(const Matrix<T, IDX>* p) {
    std::vector<Output> output(p->nrow());
    stats::MedianFactory factory(output.data(), p->ncol());
    apply<0>(p, factory);
    return output;
}

}

#endif
