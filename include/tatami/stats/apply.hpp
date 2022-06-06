#ifndef TATAMI_STATS_APPLY_H
#define TATAMI_STATS_APPLY_H

#include "../base/Matrix.hpp"
#include "config.hpp"

#include <cmath>
#include <vector>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * @file apply.hpp
 *
 * @brief Apply arbitrary calculations along rows or columns.
 */

namespace tatami {

/**
 * @tparam MARGIN The dimension over which to apply the calculation of statistics, i.e., rows (0) or columns (1).
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 * @tparam Factory Factory class to create the statistic-calculating classes.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param factory Instance of a factory.
 *
 * @section apply_overview Overview
 * In this function, we consider the matrix to be a collection of "target vectors".
 * When `MARGIN = 0`, each row is a target vector, whereas when `MARGIN = 1`, each column is a target vector.
 * The goal is to compute some statistic along each of these target vectors, e.g., the per-row sums (`row_sums()`) or the per-column variances (`row_variances()`).
 * For brevity, we will refer to the vectors along the other dimension as "running vectors".
 *
 * @section apply_factory Factory requirements 
 * Arbitrary computations are supported via the supplied instance of the factory class.
 * The factory class should have a `dense_direct()` method that accepts no arguments and returns an arbitrary `struct`.
 * The returned `struct` should have a `compute()` method that accepts the index of the target vector and a `const` pointer to the target vector. 
 * This will be used in a "direct" calculation, i.e., the statistic for a target vector is directly computed from the target vector's values in memory.
 *
 * Optionally, the factory class may have one or more of:
 *
 * - A `sparse_direct()` method that accepts no arguments and returns an arbitrary `struct`.
 *   The returned `struct` should have a `compute()` method that accepts the index of the target vector and a `SparseRange` object containing the non-zero elements of the target vector.
 * - A `dense_running()` method that accepts no arguments and returns an arbitrary `struct`.
 *   The returned `struct` should have an `add()` method that accepts a pointer to a running vector.
 *   It should also have a `finish()` method, to be called to finalize any calculations after all running vectors are supplied.
 * - A `sparse_running()` method that accepts no arguments and returns an arbitrary `struct`.
 *   The returned `struct` should have an `add()` method that accepts a `SparseRange` object specifying the non-zero elements at a running vector.
 *   It should also have a `finish()` method, to be called to finalize any calculations after all running vectors are supplied.
 *
 * In the running calculations, the statistic for each target vector is computed incrementally as new values become available in successive running vectors.
 * The idea is that `apply()` will automatically choose the most appropriate calculation based on whether the matrix is sparse, whether it prefers row/column access,
 * whether we want row/column statistics, and whether the calculation supports sparse and/or running calculations (based on the availability of the methods above).
 * For example, if running calculations are feasible for a particular statistic, `apply()` can choose a more favorable access pattern.
 *
 * Each of these `compute()` and `add()` methods should modify the contents of `factory` by reference, 
 * which is usually achieved by passing pointers from `factory` to each returned `struct` and writing to those pointers in `compute()` and `add()`.
 * Computed statistics can then be extracted from `factory` once `apply()` has finished running.
 * Results should be agnostic to the choice of calculation, notwithstanding minor differences due to numerical precision.
 *
 * See the `VarianceFactory` in `variances.hpp` for an example of a valid factory class. 
 *
 * @section apply_mutable Mutating copies
 * In some applications, it may be necessary to mutate the buffers containing the contents of each row/column (e.g., sorting for quantile calculations).
 * This is accommodated by replacing the `compute()` method with a `compute_copy()` method for the `*_direct()` outputs.
 *
 * - In the dense case, the `compute_copy()` method should accept the index of the target vector and a non-`const` pointer to a buffer with the values of the row/column.
 * - In the sparse case, the `compute_copy()` should accept the index of the target vector, the number of non-zero elements, 
 *   a non-`const` pointer to a buffer with the non-zero values, and a non-`const` pointer to the indices of the non-zero elements.
 *   The `struct` may also contain a static `copy_mode` member specifying whether one or both the values/indices should be copied into the buffers, see the options for a `SparseCopyMode`.
 *
 * If a `compute_copy()` method is available, `apply()` will create a copy upon row/column extraction and the developer is free to modify the buffers inside the `compute_copy()` method.
 * Note that we do not provide a copyable option for the running calculations as this is rarely necessary.
 *
 * See the `MedianFactory` in `medians.hpp` for an example of a factory that mutates the buffer. 
 *
 * @section apply_parallel2 Caller parallelization
 * `apply()` supports parallelization via OpenMP by default, so callers can simply compile with `-fopenmp` to parallelize their code.
 * `apply()` will automatically distribute the calculations for each target vector across available threads.
 * The maximum number of threads can be controlled with `omp_set_max_threads()`.
 *
 * Advanced callers of `apply()` may specify their own parallelization scheme by defining a `TATAMI_CUSTOM_PARALLEL` macro.
 * This should be a function-like macro with two arguments: 
 *
 * - `n`, an integer specifying the total number of jobs.
 * - `f`, a lambda specifying the function to run on a range of jobs in `n`.
 *   This should take two `size_t` arguments specifying the start and one-past-the-end of the range.
 * 
 * The `TATAMI_CUSTOM_PARALLEL` function is expected to split `n` into non-overlapping ranges where each range is assigned to a worker.
 * It should then call `f` within each worker on the corresponding range of jobs.
 * For example:
 *
 * ```cpp
 * template<class Function>
 * void parallelize(size_t n, Function f) {
 *     size_t jobs_per_worker = std::ceil(static_cast<double>(n) / 3);
 *     size_t start = 0;
 *     std::vector<std::thread> jobs;
 *     
 *     for (size_t w = 0; w < 3; ++w) {
 *         size_t end = std::min(n, start + jobs_per_worker);
 *         if (start >= end) {
 *             break;
 *         }
 *         jobs.emplace_back(f, start, end);
 *         start += jobs_per_worker;
 *     }
 * 
 *     for (auto& job : jobs) {
 *         job.join();
 *     }
 * }
 * 
 * #define TATAMI_CUSTOM_PARALLEL parallelize
 * ```
 *
 * @section apply_parallel Factory parallelization
 * Thread safety is _not_ required in each call to `compute()` and `add()` within the same instance of a returned `struct`.
 * However, there should be thread safety across instances, which is most relevant when results are being written back to the shared memory in `factory`.
 * This is usually easy to achieve as long as the memory spaces for the results of two target vectors do not overlap.
 *
 * The direct calculations are trivially parallelizable - it is assumed that they can be called independently for different target indices.
 * It is also assumed that each `compute()` call will only write to the shared memory at its supplied target index.
 *
 * To support parallelization in the running calculations, we require some overloads to the `dense_running()` and `sparse_running()` methods:
 *
 * - `dense_running()` now accepts two arguments: namely, the indices of the first and one-past-the-last target vectors to be processed.
 *   This should return a `struct` with an `add()` method that accepts a pointer to the subinterval of the running vector, corresponding to the target vectors to be processed; 
 *   plus a buffer of the same length as that subinterval.
 *   It should also have a `finish()` method, to be called to finalize any calculations after all running vectors are supplied.
 * - `sparse_running()` now accepts two arguments: namely, the indices of the first and one-past-the-last target vectors to be processed.
 *   The returned `struct` should have an `add()` method that accepts a `SparseRange` object, specifying the non-zero elements in the subinterval of the running vector;
 *   plus two buffers of the same length as that subinterval (one for the non-zero values, another for their positional indices).
 *   It should also have a `finish()` method, to be called to finalize any calculations after all running vectors are supplied.
 *
 * These overloads are optional and the function will fall back to serial processing if they are not supplied (and the function decides perform a running calculation).
 */
template<int MARGIN, typename T, typename IDX, class Factory>
void apply(const Matrix<T, IDX>* p, Factory& factory) {
    /**
     * @cond
     */
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
    if constexpr(stats::has_sparse_running<Factory>::value || stats::has_dense_running<Factory>::value) {
        if (p->prefer_rows() != ROW){

            if constexpr(stats::has_sparse_running<Factory>::value) {
                if (p->sparse()) {
#if defined(_OPENMP) || defined(TATAMI_CUSTOM_PARALLEL)
                    if constexpr(stats::has_sparse_running_parallel<Factory>::value) {
#ifndef TATAMI_CUSTOM_PARALLEL
                        #pragma omp parallel
                        {
                            size_t worker_size = std::ceil(static_cast<double>(dim) / omp_get_num_threads());
                            size_t start = worker_size * omp_get_thread_num(), end = std::min(dim, start + worker_size);
#else
                        TATAMI_CUSTOM_PARALLEL(dim, [&](size_t start, size_t end) -> void {
#endif
                            if (start < end) {
                                std::vector<T> obuffer(end - start);
                                std::vector<IDX> ibuffer(obuffer.size());
                                auto wrk = p->new_workspace(!ROW);
                                auto stat = factory.sparse_running(start, end);

                                for (size_t i = 0; i < otherdim; ++i) {
                                    if constexpr(ROW) { // flipped around; remember, we're trying to get the preferred dimension.
                                        auto range = p->sparse_column(i, obuffer.data(), ibuffer.data(), start, end, wrk.get());
                                        stat.add(range);
                                    } else {
                                        auto range = p->sparse_row(i, obuffer.data(), ibuffer.data(), start, end, wrk.get());
                                        stat.add(range);
                                    }
                                }
                                stat.finish();
                            }
#ifndef TATAMI_CUSTOM_PARALLEL
                        }
#else
                        });
#endif
                        return;
                    }
#endif
                    auto stat = factory.sparse_running();
                    std::vector<T> obuffer(dim);
                    std::vector<IDX> ibuffer(dim);
                    auto wrk = p->new_workspace(!ROW);

                    for (size_t i = 0; i < otherdim; ++i) {
                        if constexpr(ROW) { // flipped around; remember, we're trying to get the preferred dimension.
                            auto range = p->sparse_column(i, obuffer.data(), ibuffer.data(), wrk.get());
                            stat.add(range);
                        } else {
                            auto range = p->sparse_row(i, obuffer.data(), ibuffer.data(), wrk.get());
                            stat.add(range);
                        }
                    }
                    stat.finish();
                    return;
                }
            }

#if defined(_OPENMP) || defined(TATAMI_CUSTOM_PARALLEL)
            if constexpr(stats::has_dense_running_parallel<Factory>::value) {
#ifndef TATAMI_CUSTOM_PARALLEL
                #pragma omp parallel
                {
                    size_t worker_size = std::ceil(static_cast<double>(dim) / omp_get_num_threads());
                    size_t start = worker_size * omp_get_thread_num(), end = std::min(dim, start + worker_size);
#else
                TATAMI_CUSTOM_PARALLEL(dim, [&](size_t start, size_t end) -> void {
#endif                
                    if (start < end) {
                        auto stat = factory.dense_running(start, end);
                        std::vector<T> obuffer(end - start);
                        auto wrk = p->new_workspace(!ROW);

                        for (size_t i = 0; i < otherdim; ++i) {
                            if constexpr(ROW) { // flipped around, see above.
                                auto ptr = p->column(i, obuffer.data(), start, end, wrk.get());
                                stat.add(ptr);
                            } else {
                                auto ptr = p->row(i, obuffer.data(), start, end, wrk.get());
                                stat.add(ptr);
                            }
                        }
                        stat.finish();
                    }
#ifndef TATAMI_CUSTOM_PARALLEL
                }
#else
                });
#endif
                return;
            }
#endif
            auto stat = factory.dense_running();
            std::vector<T> obuffer(dim);
            auto wrk = p->new_workspace(!ROW);

            for (size_t i = 0; i < otherdim; ++i) {
                if constexpr(ROW) { // flipped around, see above.
                    auto ptr = p->column(i, obuffer.data(), wrk.get());
                    stat.add(ptr);
                } else {
                    auto ptr = p->row(i, obuffer.data(), wrk.get());
                    stat.add(ptr);
                }
            }
            stat.finish();
            return;
        }
    }

    if constexpr(stats::has_sparse_direct<Factory>::value) {
        if (p->sparse()) {
#ifndef TATAMI_CUSTOM_PARALLEL
            #pragma omp parallel
            {
#else
            TATAMI_CUSTOM_PARALLEL(dim, [&](size_t start, size_t end) -> void {
#endif
                std::vector<T> obuffer(otherdim);
                auto wrk = p->new_workspace(ROW);
                std::vector<IDX> ibuffer(otherdim);
                auto stat = factory.sparse_direct();

                constexpr bool do_copy = stats::has_nonconst_sparse_compute<decltype(stat), T, IDX>::value;
                constexpr SparseCopyMode copy_mode = stats::nonconst_sparse_compute_copy_mode<decltype(stat)>::value;

#ifndef TATAMI_CUSTOM_PARALLEL
                #pragma omp for schedule(static) 
                for (size_t i = 0; i < dim; ++i) {
#else
                for (size_t i = start; i < end; ++i) {
#endif
                    if constexpr(ROW) {
                        if constexpr(do_copy) {
                            auto range = p->sparse_row_copy(i, obuffer.data(), ibuffer.data(), copy_mode, wrk.get());
                            stat.compute_copy(i, range.number, obuffer.data(), ibuffer.data());
                        } else {
                            auto range = p->sparse_row(i, obuffer.data(), ibuffer.data(), wrk.get());
                            stat.compute(i, range);
                        }
                    } else {
                        if constexpr(do_copy) {
                            auto range = p->sparse_column_copy(i, obuffer.data(), ibuffer.data(), copy_mode, wrk.get());
                            stat.compute_copy(i, range.number, obuffer.data(), ibuffer.data());
                        } else {
                            auto range = p->sparse_column(i, obuffer.data(), ibuffer.data(), wrk.get());
                            stat.compute(i, range);
                        }
                    }
                }
#ifndef TATAMI_CUSTOM_PARALLEL
            }
#else
            });
#endif
            return;
        }
    }

#ifndef TATAMI_CUSTOM_PARALLEL
    #pragma omp parallel
    {
#else
    TATAMI_CUSTOM_PARALLEL(dim, [&](size_t start, size_t end) -> void {
#endif
        std::vector<T> obuffer(otherdim);
        auto wrk = p->new_workspace(ROW);
        auto stat = factory.dense_direct();
        constexpr bool do_copy = stats::has_nonconst_dense_compute<decltype(stat), T>::value;

#ifndef TATAMI_CUSTOM_PARALLEL
        #pragma omp for schedule(static)
        for (size_t i = 0; i < dim; ++i) {
#else            
        for (size_t i = start; i < end; ++i) {
#endif
            if constexpr(ROW) {
                if constexpr(do_copy) {
                    auto ptr = p->row_copy(i, obuffer.data(), wrk.get());
                    stat.compute_copy(i, obuffer.data());
                } else {
                    auto ptr = p->row(i, obuffer.data(), wrk.get());
                    stat.compute(i, ptr);
                }
            } else {
                if constexpr(do_copy) {
                    auto ptr = p->column_copy(i, obuffer.data(), wrk.get());
                    stat.compute_copy(i, obuffer.data());
                } else {
                    auto ptr = p->column(i, obuffer.data(), wrk.get());
                    stat.compute(i, ptr);
                }
            }
        }
#ifndef TATAMI_CUSTOM_PARALLEL            
    }
#else
    });
#endif

    /**
     * @endcond
     */
    return; 
}

}

#endif
