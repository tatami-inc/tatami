#ifndef TATAMI_STATS_APPLY_H
#define TATAMI_STATS_APPLY_H

#include "../base/typed_matrix.hpp"

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
 * @param p Pointer to a `tatami::typed_matrix`.
 * @param stat Instance of a statistics store.
 *
 * @return The `direct()` or `running()` methods in `stat` are invoked to fill up the store with the relevant row- or column-wise statistics.
 */
template<int MARGIN, typename T, typename IDX, template<typename,bool,bool> class STAT>
inline typename STAT<T, false, false>::value apply(const typed_matrix<T, IDX>* p) {
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
    if constexpr(STAT<T, false, true>::runnable || STAT<T, true, true>::runnable) {
        if (p->prefer_rows() != ROW){
            std::vector<T> obuffer(dim);
            auto wrk = p->new_workspace(!ROW);

            if constexpr(STAT<T, true, true>::sparse) {
                if (p->sparse()) {
                    STAT<T, true, true> stat(dim, otherdim);
                    std::vector<IDX> ibuffer(dim);
                    for (size_t i = 0; i < otherdim; ++i) {
                        if constexpr(ROW) { // flipped around; remember, we're trying to get the preferred dimension.
                            auto range = p->sparse_column(i, obuffer.data(), ibuffer.data(), wrk.get());
                            stat.running(range);
                        } else {
                            auto range = p->sparse_row(i, obuffer.data(), ibuffer.data(), wrk.get());
                            stat.running(range);
                        }
                    }
                    return stat.yield();
                }
            }

            if constexpr(STAT<T, false, true>::runnable) {
                STAT<T, false, true> stat(dim, otherdim);
                for (size_t i = 0; i < otherdim; ++i) {
                    if constexpr(ROW) { // flipped around, see above.
                        auto ptr = p->column(i, obuffer.data(), wrk.get());
                        stat.running(ptr);
                    } else {
                        auto ptr = p->row(i, obuffer.data(), wrk.get());
                        stat.running(ptr);
                    }
                }
                return stat.yield();
            }
        }
    }

    std::vector<T> obuffer(otherdim);
    auto wrk = p->new_workspace(ROW);

    if constexpr(STAT<T, true, false>::sparse) {
        if (p->sparse()) {
            STAT<T, true, false> stat(dim, otherdim);
            std::vector<IDX> ibuffer(otherdim);
            for (size_t i = 0; i < dim; ++i) {
                if constexpr(ROW) {
                    auto range = p->sparse_row(i, obuffer.data(), ibuffer.data(), wrk.get());
                    stat.direct(i, range);
                } else {
                    auto range = p->sparse_column(i, obuffer.data(), ibuffer.data(), wrk.get());
                    stat.direct(i, range);
                }
            }
            return stat.yield();
        }
    }

    STAT<T, false, false> stat(dim, otherdim);
    for (size_t i = 0; i < dim; ++i) {
        if constexpr(ROW) {
            auto ptr = p->row(i, obuffer.data(), wrk.get());
            stat.direct(i, ptr);
        } else {
            auto ptr = p->column(i, obuffer.data(), wrk.get());
            stat.direct(i, ptr);
        }
    }
    return stat.yield();
}

}

#endif
