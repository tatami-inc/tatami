#ifndef TATAMI_VARS_HPP
#define TATAMI_VARS_HPP

#include "../base/typed_matrix.hpp"
#include <vector>
#include <algorithm>

/**
 * @file variances.hpp
 *
 * Compute row and column variances from a `tatami::typed_matrix`.
 */

namespace tatami {

template<bool ROW, typename T, typename IDX>
inline std::vector<double> direct_variances(const typed_matrix<T, IDX>* p) {
    size_t NR = p->nrow(), NC = p->ncol();
    size_t dim = ROW ? NR : NC;
    size_t otherdim = ROW ? NC : NR;

    std::vector<double> output(dim);
    std::vector<T> buffer(otherdim);
    auto wrk = p->new_workspace(ROW);

    if (p->sparse()) {
        std::vector<IDX> ibuffer(otherdim);

        for (size_t i = 0; i < dim; ++i) {
            tatami::sparse_range<T, IDX> range;
            if (ROW) {
                range = p->sparse_row(i, buffer.data(), ibuffer.data(), wrk);
            } else {
                range = p->sparse_column(i, buffer.data(), ibuffer.data(), wrk);
            }

            const double mean = std::accumulate(range.value, range.value + range.number, 0.0)/otherdim;
            double& out = output[i];
            for (size_t j = 0; j < range.number; ++j) {
                out += (range.value[j] - mean) * (range.value[j] - mean);
            }
            out += mean * mean * (otherdim - range.number);
        }
    } else {
        for (size_t i = 0; i < dim; ++i) {
            const T* ptr;
            if (ROW) {
                ptr = p->row(i, buffer.data(), wrk);
            } else {
                ptr = p->column(i, buffer.data(), wrk);
            }

            const double mean = std::accumulate(ptr, ptr + otherdim, 0.0)/otherdim;
            double& out = output[i];
            for (size_t j = 0; j < otherdim; ++j, ++ptr) {
                out += (*ptr - mean) * (*ptr - mean);
            }
        }
    }

    return output;
}

// Using Welford's algorithm to compute the variance in the 'other' dimension.
// which should be much more cache-friendly for large matrices. We run through
// only non-zero counts when circumstances allow for it.
template<bool ROW, typename T, typename IDX>
inline std::vector<double> running_variances(const typed_matrix<T, IDX>* p) {
    size_t NR = p->nrow(), NC = p->ncol();
    size_t dim = ROW ? NR : NC;
    size_t otherdim = ROW ? NC : NR;

    std::vector<double> output(dim);
    std::vector<double> running_mean(dim);
    std::vector<T> buffer(dim);
    auto wrk = p->new_workspace(!ROW);

    if (p->sparse()) {
        std::vector<IDX> ibuffer(dim);
        std::vector<double> running_nzero(dim);

        for (size_t i = 0; i < otherdim; ++i) {
            tatami::sparse_range<T, IDX> range;
            if (!ROW) {
                range = p->sparse_row(i, buffer.data(), ibuffer.data(), wrk);
            } else {
                range = p->sparse_column(i, buffer.data(), ibuffer.data(), wrk);
            }

            auto sIt = output.begin();
            auto mIt = running_mean.begin();
            auto nzIt = running_nzero.begin();

            for (size_t j = 0; j < range.number; ++j) {
                auto ri = range.index[j];
                auto& curM = *(mIt + ri);
                auto& curS = *(sIt + ri);
                auto& curNZ = *(nzIt + ri);
                ++curNZ;
                const auto& curval = range.value[j];
                
                const double delta = curval - curM;
                curM += delta / curNZ;
                curS += delta * (curval - curM);
            }        
        }

        auto sIt = output.begin();
        auto mIt = running_mean.begin();
        auto nzIt = running_nzero.begin();
        for (size_t i = 0; i < dim; ++i, ++mIt, ++sIt, ++nzIt) {
            const double curNZ = *nzIt;
            const double ratio = curNZ / otherdim;
            auto& curM = *mIt;
            *sIt += curM * curM * ratio * (otherdim - curNZ);
            curM *= ratio;
        }
    } else {
        for (size_t i = 0; i < otherdim; ++i) {
            const T* ptr;
            if (!ROW) {
                ptr = p->row(i, buffer.data(), wrk);
            } else {
                ptr = p->column(i, buffer.data(), wrk);
            }

            auto oIt = output.begin();
            auto mIt = running_mean.begin();
            for (size_t j = 0; j < dim; ++j, ++ptr, ++oIt, ++mIt) {
                const double delta=*ptr - *mIt;
                *mIt += delta/(i+1);
                *oIt += delta*(*ptr - *mIt);
            }
        }
    }

    return output;
}


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
    if (p->prefer_rows()) {
        return running_variances<false, T, IDX>(p);
    } else {
        return direct_variances<false, T, IDX>(p);
    }
}

/**
 * @copydoc column_variances()
 */
template<typename T, typename IDX>
inline std::vector<T> column_variances(std::shared_ptr<const typed_matrix<T, IDX> > p) {
    return column_variances<T, IDX>(p.get());
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
    if (p->prefer_rows()) {
        return direct_variances<true, T, IDX>(p);
    } else {
        return running_variances<true, T, IDX>(p);
    }
}

/**
 * @copydoc row_variances()
 */
template<typename T, typename IDX>
inline std::vector<T> row_variances(std::shared_ptr<const typed_matrix<T, IDX> > p) {
    return row_variances<T, IDX>(p.get());
}

}

#endif
