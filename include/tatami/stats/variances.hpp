#ifndef TATAMI_VARS_HPP
#define TATAMI_VARS_HPP

#include "../base/typed_matrix.hpp"
#include <vector>
#include <algorithm>
#include <limits>

/**
 * @file variances.hpp
 *
 * Compute row and column variances from a `tatami::typed_matrix`.
 */

namespace tatami {

template<typename T, bool SPARSE, bool RUNNABLE> 
struct StatsVarianceHelper {
    StatsVarianceHelper(size_t n, size_t d) : store(n), running_mean(RUNNABLE ? n : 0), running_nnzero(RUNNABLE && SPARSE ? n : 1), dim(d) {}

    static const bool sparse = SPARSE;
    static const bool runnable = RUNNABLE;
    typedef std::vector<T> value;

    void direct(size_t i, const T* ptr, T* buffer) {
        static_assert(!SPARSE && !RUNNABLE);
        const double mean = std::accumulate(ptr, ptr + dim, 0.0)/dim;
        double& out = store[i];
        for (size_t j = 0; j < dim; ++j, ++ptr) {
            out += (*ptr - mean) * (*ptr - mean);
        }
        return;
    }

    template<typename IDX>
    void direct(size_t i, const sparse_range<T, IDX>& range, T* vbuffer, IDX* ibuffer) {
        static_assert(SPARSE && !RUNNABLE);
        const double mean = std::accumulate(range.value, range.value + range.number, 0.0)/dim;
        double& out = store[i];
        for (size_t j = 0; j < range.number; ++j) {
            out += (range.value[j] - mean) * (range.value[j] - mean);
        }
        out += mean * mean * (dim - range.number);
        return;
    }

    /* Using Welford's algorithm to compute the variance in the 'other' dimension.
     * which should be much more cache-friendly for large matrices. We run through
     * only non-zero counts when circumstances allow for it.
     */
    void running(const T* ptr, T* buffer) {
        static_assert(!SPARSE && RUNNABLE);
        auto mIt = running_mean.begin();
        auto& counter = running_nnzero[0];
        ++counter;

        for (auto sIt = store.begin(); sIt < store.end(); ++sIt, ++ptr, ++mIt) {
            const double delta=*ptr - *mIt;
            *mIt += delta/counter;
            *sIt += delta*(*ptr - *mIt);
        }
        return;
    }

    template<typename IDX>
    void running(sparse_range<T, IDX> range, T* vbuffer, IDX* ibuffer) {
        static_assert(SPARSE && RUNNABLE);
        auto sIt = store.begin();
        auto mIt = running_mean.begin();
        auto nzIt = running_nnzero.begin();

        for (size_t j = 0; j < range.number; ++j, ++range.index, ++range.value) {
            auto ri = *range.index;
            auto& curM = *(mIt + ri);
            auto& curS = *(sIt + ri);
            auto& curNZ = *(nzIt + ri);
            ++curNZ;

            const auto& curval = *range.value;
            const double delta = curval - curM;
            curM += delta / curNZ;
            curS += delta * (curval - curM);
        }

        return;
    }

    value yield() {
        if constexpr(SPARSE && RUNNABLE) {
            auto mIt = running_mean.begin();
            auto nzIt = running_nnzero.begin();
            for (auto sIt = store.begin(); sIt != store.end(); ++mIt, ++sIt, ++nzIt) {
                const double curNZ = *nzIt;
                const double ratio = curNZ / dim;
                auto& curM = *mIt;
                *sIt += curM * curM * ratio * (dim - curNZ);
                curM *= ratio;
            }
        }

        if (dim > 1) {
            for (auto& s : store) {
                s /= dim - 1;
            }
        } else {
            std::fill(store.begin(), store.end(), std::numeric_limits<double>::quiet_NaN());
        }
        return store;
    }
private:
    value store;
    std::vector<double> running_mean;
    std::vector<int> running_nnzero;
    size_t dim;
};

/**
 * This uses the usual algorithm for matrices where `tatami::matrix::prefer_rows()` is false, otherwise it uses Welford's algorithm.
 * As a result, the computed variances will be slightly different (within numerical precision) for row- and column-major matrices.
 *
 * @tparam T Type of the matrix value, should be numeric.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::typed_matrix`.
 *
 * @return A vector of length equal to the number of columns, containing the column variances.
 */
template<typename T, typename IDX>
inline std::vector<T> column_variances(const typed_matrix<T, IDX>* p) {
    return apply<1, T, IDX, StatsVarianceHelper>(p);
}

/**
 * This uses the usual algorithm for matrices where `tatami::matrix::prefer_rows()` is true, otherwise it uses Welford's algorithm.
 * As a result, the computed variances will be slightly different (within numerical precision) for row- and column-major matrices.
 *
 * @tparam T Type of the matrix value, should be numeric.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::typed_matrix`.
 *
 * @return A vector of length equal to the number of rows, containing the row variances.
 */
template<typename T, typename IDX>
inline std::vector<T> row_variances(const typed_matrix<T, IDX>* p) {
    return apply<0, T, IDX, StatsVarianceHelper>(p);
}

}

#endif
