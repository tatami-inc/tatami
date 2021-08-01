#ifndef TATAMI_STATS_APPLY_H
#define TATAMI_STATS_APPLY_H

#include "../base/Matrix.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace tatami {

/**
 * @tparam MARGIN The dimension over which to apply the calculation of statistics, i.e., rows (0) or columns (1).
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 * @tparam STAT Class implementing a statistics store.
 * This should have the `runnable` boolean member, indicating if it is capable of computing a running statistic;
 * the `sparse` boolean member, indicating if the statistic can be computed sparsely;
 * and the `direct()` method, to compute the statistics from a dense vector.
 * If `runnable = true`, there should be `running()` method to compute a running statistic;
 * and if `sparse = true`, there should be `*_sparse()` equivalents of the other methods.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param stat Instance of a statistics store.
 *
 * @return A vector of row- or column-wise statistics.
 */
template<int MARGIN, typename T, typename IDX, class Factory>
void apply(const Matrix<T, IDX>* p, Factory& factory) {
    size_t NR = p->nrow(), NC = p->ncol();

    /* One might question why we use MARGIN in the template if we just convert
     * it to ROW here. This is because 'MARGIN' is used when we're doing
     * something to the matrix; 'ROW' is used for the representation of the
     * matrix itself. I'm keeping the distinction clear in the interface.
     */
    constexpr bool ROW = (MARGIN == 0);

    const size_t dim = (ROW ? NR : NC);
    const size_t otherdim = (ROW ? NC : NR);

    /* If we support running calculations AND the preference 
     * is not consistent with the margin, we give it a shot.
     */
    if constexpr(Factory::supports_running) {
        if (p->prefer_rows() != ROW){

            if constexpr(Factory::supports_sparse) {
                if (p->sparse()) {
#ifdef _OPENMP
                    #pragma omp parallel
                    {
                        int nworkers = omp_get_num_threads();
                        size_t worker_size = std::ceil(static_cast<double>(dim) / nworkers);
                        size_t start = worker_size * omp_get_thread_num(), end = std::min(dim, start + worker_size);

                        std::vector<T> obuffer(end - start);
                        std::vector<IDX> ibuffer(obuffer.size());
                        auto wrk = p->new_workspace(!ROW);
                        auto stat = factory.sparse_running(start, end);

                        for (size_t i = 0; i < otherdim; ++i) {
                            if constexpr(ROW) { // flipped around; remember, we're trying to get the preferred dimension.
                                auto range = p->sparse_column(i, obuffer.data(), ibuffer.data(), start, end, wrk.get());
                                stat.add(range, obuffer.data(), ibuffer.data());
                            } else {
                                auto range = p->sparse_row(i, obuffer.data(), ibuffer.data(), start, end, wrk.get());
                                stat.add(range, obuffer.data(), ibuffer.data());
                            }
                        }
                        stat.finish();
                    }
#else
                    auto stat = factory.sparse_running();
                    std::vector<T> obuffer(dim);
                    std::vector<IDX> ibuffer(dim);
                    auto wrk = p->new_workspace(!ROW);

                    for (size_t i = 0; i < otherdim; ++i) {
                        if constexpr(ROW) { // flipped around; remember, we're trying to get the preferred dimension.
                            auto range = p->sparse_column(i, obuffer.data(), ibuffer.data(), wrk.get());
                            stat.add(range, obuffer.data(), ibuffer.data());
                        } else {
                            auto range = p->sparse_row(i, obuffer.data(), ibuffer.data(), wrk.get());
                            stat.add(range, obuffer.data(), ibuffer.data());
                        }
                    }
                    stat.finish();
#endif
                    return;
                }
            }

#ifdef _OPENMP
            #pragma omp parallel
            {
                int nworkers = omp_get_num_threads();
                size_t worker_size = std::ceil(static_cast<double>(dim) / nworkers);
                size_t start = worker_size * omp_get_thread_num(), end = std::min(dim, start + worker_size);

                auto stat = factory.dense_running(start, end);
                std::vector<T> obuffer(end - start);
                auto wrk = p->new_workspace(!ROW);

                for (size_t i = 0; i < otherdim; ++i) {
                    if constexpr(ROW) { // flipped around, see above.
                        auto ptr = p->column(i, obuffer.data(), start, end, wrk.get());
                        stat.add(ptr, obuffer.data());
                    } else {
                        auto ptr = p->row(i, obuffer.data(), start, end, wrk.get());
                        stat.add(ptr, obuffer.data());
                    }
                }
                stat.finish();
            }
#else
            auto stat = factory.dense_running();
            std::vector<T> obuffer(dim);
            auto wrk = p->new_workspace(!ROW);

            for (size_t i = 0; i < otherdim; ++i) {
                if constexpr(ROW) { // flipped around, see above.
                    auto ptr = p->column(i, obuffer.data(), wrk.get());
                    stat.add(ptr, obuffer.data());
                } else {
                    auto ptr = p->row(i, obuffer.data(), wrk.get());
                    stat.add(ptr, obuffer.data());
                }
            }
            stat.finish();
#endif
            return;
        }
    }

    if constexpr(Factory::supports_sparse) {
        if (p->sparse()) {
            #pragma omp parallel
            {
                std::vector<T> obuffer(otherdim);
                auto wrk = p->new_workspace(ROW);
                std::vector<IDX> ibuffer(otherdim);
                auto stat = factory.sparse_direct();

                #pragma omp for schedule(static) 
                for (size_t i = 0; i < dim; ++i) {
                    if constexpr(ROW) {
                        auto range = p->sparse_row(i, obuffer.data(), ibuffer.data(), wrk.get());
                        stat.compute(i, range, obuffer.data(), ibuffer.data());
                    } else {
                        auto range = p->sparse_column(i, obuffer.data(), ibuffer.data(), wrk.get());
                        stat.compute(i, range, obuffer.data(), ibuffer.data());
                    }
                }
            }
            return;
        }
    }

    #pragma omp parallel
    {
        std::vector<T> obuffer(otherdim);
        auto wrk = p->new_workspace(ROW);
        auto stat = factory.dense_direct();

        #pragma omp for schedule(static)
        for (size_t i = 0; i < dim; ++i) {
            if constexpr(ROW) {
                auto ptr = p->row(i, obuffer.data(), wrk.get());
                stat.compute(i, ptr, obuffer.data());
            } else {
                auto ptr = p->column(i, obuffer.data(), wrk.get());
                stat.compute(i, ptr, obuffer.data());
            }
        }
    }

    return; 
}

}

#endif
