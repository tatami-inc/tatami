#ifndef TATAMI_STATS_APPLY_H
#define TATAMI_STATS_APPLY_H

#include "../base/Matrix.hpp"
#include "config.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cmath>
#include <vector>
#include <algorithm>

/**
 * @file apply.hpp
 *
 * @brief Apply arbitrary calculations along rows or columns.
 */

namespace tatami {

namespace stats {

/**
 * @cond
 */
template<bool ROW, typename T, typename IDX, class Factory>
void apply_sparse_running(size_t dim, size_t otherdim, const Matrix<T, IDX>* p, Factory& factory, int threads) {
    if constexpr(stats::has_prepare_sparse_running<Factory>::value) {
        factory.prepare_sparse_running();
    }

#if defined(_OPENMP) || defined(TATAMI_CUSTOM_PARALLEL)
    if constexpr(stats::has_sparse_running_parallel<Factory>::value) {
#ifndef TATAMI_CUSTOM_PARALLEL
        #pragma omp parallel num_threads(threads)
        {
            // We use the omp methods here to get the total number of threads, just in case something
            // causes us to have fewer threads than the number we requested.
            size_t worker_size = std::ceil(static_cast<double>(dim) / omp_get_num_threads());
            size_t start = worker_size * omp_get_thread_num(), end = std::min(dim, start + worker_size);
#else
        TATAMI_CUSTOM_PARALLEL(dim, [&](size_t start, size_t end) -> void {
#endif
            if (start < end) {
                size_t len = end - start;
                std::vector<T> obuffer(len);
                std::vector<IDX> ibuffer(obuffer.size());
                auto stat = factory.sparse_running(start, end);
                auto wrk = new_workspace<!ROW>(p, start, len);

                for (size_t i = 0; i < otherdim; ++i) {
                    if constexpr(ROW) { // flipped around; remember, we're trying to get the preferred dimension.
                        auto range = p->sparse_column(i, obuffer.data(), ibuffer.data(), wrk.get());
                        stat.add(range);
                    } else {
                        auto range = p->sparse_row(i, obuffer.data(), ibuffer.data(), wrk.get());
                        stat.add(range);
                    }
                }

                if constexpr(stats::has_finish<decltype(stat)>::value) {
                    stat.finish();
                }
            }
#ifndef TATAMI_CUSTOM_PARALLEL
        }
#else
        }, threads);
#endif
        return;
    }
#endif

    auto stat = factory.sparse_running();
    std::vector<T> obuffer(dim);
    std::vector<IDX> ibuffer(dim);
    auto wrk = new_workspace<!ROW>(p);

    for (size_t i = 0; i < otherdim; ++i) {
        if constexpr(ROW) { // flipped around; remember, we're trying to get the preferred dimension.
            auto range = p->sparse_column(i, obuffer.data(), ibuffer.data(), wrk.get());
            stat.add(range);
        } else {
            auto range = p->sparse_row(i, obuffer.data(), ibuffer.data(), wrk.get());
            stat.add(range);
        }
    }

    if constexpr(stats::has_finish<decltype(stat)>::value) {
        stat.finish();
    }
    return;
}

template<bool ROW, typename T, typename IDX, class Factory>
void apply_dense_running(size_t dim, size_t otherdim, const Matrix<T, IDX>* p, Factory& factory, int threads) {
    if constexpr(stats::has_prepare_dense_running<Factory>::value) {
        factory.prepare_dense_running();
    }

#if defined(_OPENMP) || defined(TATAMI_CUSTOM_PARALLEL)
    if constexpr(stats::has_dense_running_parallel<Factory>::value) {
#ifndef TATAMI_CUSTOM_PARALLEL
        #pragma omp parallel num_threads(threads)
        {
            size_t worker_size = std::ceil(static_cast<double>(dim) / omp_get_num_threads());
            size_t start = worker_size * omp_get_thread_num(), end = std::min(dim, start + worker_size);
#else
        TATAMI_CUSTOM_PARALLEL(dim, [&](size_t start, size_t end) -> void {
#endif
            if (start < end) {
                size_t len = end - start;
                std::vector<T> obuffer(len);
                auto stat = factory.dense_running(start, end);
                auto wrk = new_workspace<!ROW>(p, start, len);

                for (size_t i = 0; i < otherdim; ++i) {
                    if constexpr(ROW) { // flipped around, see above.
                        auto ptr = p->column(i, obuffer.data(), wrk.get());
                        stat.add(ptr);
                    } else {
                        auto ptr = p->row(i, obuffer.data(), wrk.get());
                        stat.add(ptr);
                    }
                }

                if constexpr(stats::has_finish<decltype(stat)>::value) {
                    stat.finish();
                }
            }
#ifndef TATAMI_CUSTOM_PARALLEL
        }
#else
        }, threads);
#endif
        return;
    }
#endif
    auto stat = factory.dense_running();
    std::vector<T> obuffer(dim);
    auto wrk = new_workspace<!ROW>(p);

    for (size_t i = 0; i < otherdim; ++i) {
        if constexpr(ROW) { // flipped around, see above.
            auto ptr = p->column(i, obuffer.data(), wrk.get());
            stat.add(ptr);
        } else {
            auto ptr = p->row(i, obuffer.data(), wrk.get());
            stat.add(ptr);
        }
    }

    if constexpr(stats::has_finish<decltype(stat)>::value) {
        stat.finish();
    }
    return;
}

template<bool ROW, typename T, typename IDX, class Factory>
void apply_sparse_direct(size_t dim, size_t otherdim, const Matrix<T, IDX>* p, Factory& factory, int threads) {
    if constexpr(stats::has_prepare_sparse_direct<Factory>::value) {
        factory.prepare_sparse_direct();
    }

#ifndef TATAMI_CUSTOM_PARALLEL
    #pragma omp parallel num_threads(threads)
    {
#else
    TATAMI_CUSTOM_PARALLEL(dim, [&](size_t start, size_t end) -> void {
#endif
        std::vector<T> obuffer(otherdim);
        std::vector<IDX> ibuffer(otherdim);
        auto stat = factory.sparse_direct();
        auto wrk = new_workspace<ROW>(p);

        constexpr bool do_copy = stats::has_nonconst_sparse_compute<decltype(stat), T, IDX>::value;
        constexpr SparseCopyMode copy_mode = stats::nonconst_sparse_compute_copy_mode<decltype(stat)>::value;

#ifndef TATAMI_CUSTOM_PARALLEL
        // Without chunk size specification, static scheduling should distribute
        // one chunk to each job, so everything is as contiguous as possible.
        #pragma omp for schedule(static)
        for (size_t i = 0; i < dim; ++i) {
#else
        for (size_t i = start; i < end; ++i) {
#endif
            if constexpr(ROW) {
                if constexpr(do_copy) {
                    auto range = p->sparse_row_copy(i, obuffer.data(), ibuffer.data(), wrk.get(), copy_mode);
                    stat.compute_copy(i, range.number, obuffer.data(), ibuffer.data());
                } else {
                    auto range = p->sparse_row(i, obuffer.data(), ibuffer.data(), wrk.get());
                    stat.compute(i, range);
                }
            } else {
                if constexpr(do_copy) {
                    auto range = p->sparse_column_copy(i, obuffer.data(), ibuffer.data(), wrk.get(), copy_mode);
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
    }, threads);
#endif
}

template<bool ROW, typename T, typename IDX, class Factory>
void apply_dense_direct(size_t dim, size_t otherdim, const Matrix<T, IDX>* p, Factory& factory, int threads) {
    if constexpr(stats::has_prepare_dense_direct<Factory>::value) {
        factory.prepare_dense_direct();
    }

#ifndef TATAMI_CUSTOM_PARALLEL
    #pragma omp parallel num_threads(threads)
    {
#else
    TATAMI_CUSTOM_PARALLEL(dim, [&](size_t start, size_t end) -> void {
#endif
        std::vector<T> obuffer(otherdim);
        auto wrk = new_workspace<ROW>(p);
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
    }, threads);
#endif
    return;
}
/**
 * @endcond
 */

}

/**
 * @tparam MARGIN The dimension over which to apply the calculation of statistics, i.e., rows (0) or columns (1).
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 * @tparam Factory Factory class to create the statistic-calculating classes.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param factory Instance of a factory.
 * @param threads Number of threads to use for matrix traversal.
 *
 * @section apply_overview Overview
 * In this function, we consider the matrix to be a collection of "target vectors".
 * When `MARGIN = 0`, each row is a target vector, whereas when `MARGIN = 1`, each column is a target vector.
 * The goal is to compute some statistic along each of these target vectors, e.g., the per-row sums (`row_sums()`) or the per-column variances (`row_variances()`).
 * (For brevity, we will refer to the vectors along the other dimension as "non-target vectors".)
 *
 * `apply()` will automatically choose the most efficient iteration strategy for the supplied matrix.
 * This depends on whether the matrix is sparse or dense;
 * whether the matrix prefers iteration across rows or columns;
 * and whether the desired statistic must be computed directly (i.e., given one target vector at a time)
 * or can be computed in a "running" fashion (i.e., given one non-target vector at a time).
 *
 * @section apply_factory Factory requirements 
 * Arbitrary computations are supported via the supplied `factory` instance of the factory class.
 * The factory class should have a `dense_direct()` method that accepts no arguments and returns an arbitrary object (referred to here as a "worker").
 * The worker should have a `compute()` method that accepts the index of the target vector and a `const` pointer to the target vector. 
 * This will be used in the "direct" iteration strategy, i.e., the statistic for a target vector is directly computed from the target vector's values in memory.
 *
 * Optionally, the factory class may have one or more of:
 *
 * - A `sparse_direct()` method that accepts no arguments and returns an arbitrary worker (not necessarily of the same class as `dense_direct()`).
 *   The worker should have a `compute()` method that accepts the index of the target vector and a `SparseRange` object containing the non-zero elements of the target vector.
 *   This will be favored for direct iteration over a sparse matrix.
 * - A `dense_running()` method that accepts no arguments and returns an arbitrary worker.
 *   The worker should have an `add()` method that accepts a pointer to a non-target vector.
 *   It may also have a `finish()` method, to be called to finalize any calculations after all non-target vectors are supplied.
 *   This will be favored for running iteration over a dense matrix.
 * - A `sparse_running()` method that accepts no arguments and returns an arbitrary worker.
 *   The worker should have an `add()` method that accepts a `SparseRange` object specifying the non-zero elements of a non-target vector.
 *   It may also have a `finish()` method, to be called to finalize any calculations after all non-target vectors are supplied.
 *   This will be favored for running iteration over a sparse matrix.
 *
 * `apply()` will automatically choose the most appropriate iteration strategy based on the various inputs.
 * If `MARGIN` is consistent with the matrix's preferred dimension for iteration, a direct strategy is used; otherwise a running strategy will be chosen.
 * If the matrix is sparse, the `sparse_*()` versions of the methods above will be favored.
 * If the factory does not implement the method for the best iteration strategy, `apply()` will attempt other strategies, finally falling back to the `dense_direct()` method.
 * We assume that iteration in the non-preferred dimension is costly and will attempt to avoid this wherever possible.
 *
 * Each of these `compute()` and `add()` methods should modify the contents of `factory` by reference.
 * This is achieved by passing pointers from `factory` to each worker and writing to those pointers in `compute()` and `add()`.
 * Computed statistics can then be extracted from `factory` once `apply()` has finished running.
 * Results should be agnostic to the choice of calculation, notwithstanding minor differences due to numerical precision.
 *
 * In the running calculations, the statistic for each target vector should be computed incrementally as new values become available in successive non-target vectors.
 * For example, if the matrix prefers row access but we want to compute column sums, `apply()` can use a running strategy to exploit the more efficient access pattern.
 * This could be achieved in `add()` by iterating across the supplied non-target vector and adding each entry to the running statistic for the corresponding target vector.
 *
 * See the `VarianceFactory` in `variances.hpp` for an example of a valid factory class. 
 *
 * @section apply_mutable Mutating copies
 * In some applications, it may be necessary to mutate the buffers containing the contents of each target vector (e.g., sorting for quantile calculations).
 * This is accommodated by replacing the `compute()` method with a `compute_copy()` method for the `*_direct()` workers.
 *
 * - In the dense case, the `compute_copy()` method should accept the index of the target vector and a non-`const` pointer to a buffer with the values of the row/column.
 * - In the sparse case, the `compute_copy()` should accept the index of the target vector, the number of non-zero elements, 
 *   a non-`const` pointer to a buffer with the non-zero values, and a non-`const` pointer to the indices of the non-zero elements.
 *   The worker may also contain a static `copy_mode` member specifying whether one or both the values/indices should be copied into the buffers, see the options for a `SparseCopyMode`.
 *
 * If a `compute_copy()` method is available, `apply()` will create a copy upon row/column extraction and the developer is free to modify the buffers inside the `compute_copy()` method.
 * Note that we do not provide a copyable option for the `*_running()` workers, as there is nothing to mutate when the target vector is never instantiated.
 *
 * See the `MedianFactory` in `medians.hpp` for an example of a factory that mutates the buffer. 
 *
 * @section method_prep Method-specific preparation
 * The factory class may optionally implement any number of the following methods:
 *
 * - `prepare_dense_direct()`, which will be called once before any invocation of `dense_direct()`.
 * - `prepare_dense_running()`, which will be called once before any invocation of `dense_running()`.
 * - `prepare_sparse_direct()`, which will be called once before any invocation of `sparse_direct()`.
 * - `prepare_sparse_running()`, which will be called once before any invocation of `sparse_running()`.
 *
 * This can be used by developers to perform any necessary preparation before iteration over the input matrix.
 * Importantly, the manner of preparation can vary according to the chosen iteration strategy.
 * For example, running calculations will typically require more intermediate structures than their direct counterparts;
 * the `prepare_*_running()` methods can be used to set up those intermediates as needed,
 * without committing to the setup cost if `apply()` chooses to use a direct strategy.
 * 
 * @section apply_parallel2 Caller parallelization
 * `apply()` supports parallelization via OpenMP by default, so callers can simply compile with `-fopenmp` to parallelize their code.
 * `apply()` will automatically distribute the calculations for each target vector across available threads.
 * The maximum number of threads is defined by the `threads` argument.
 *
 * Advanced callers of `apply()` may specify their own parallelization scheme by defining a `TATAMI_CUSTOM_PARALLEL` macro.
 * This should be a function-like macro with three arguments: 
 *
 * - `n`, an integer specifying the total number of jobs.
 * - `f`, a lambda specifying the function to run on a range of jobs in `n`.
 *   This should take two `size_t` arguments specifying the start and one-past-the-end of the range.
 * - `t`, an integer specifying the number of available threads. 
 *   This is set to `threads`.
 * 
 * The `TATAMI_CUSTOM_PARALLEL` function is expected to split `n` into non-overlapping ranges where each range is assigned to a worker.
 * It should then call `f` within each worker on the corresponding range of jobs.
 * For example:
 *
 * ```cpp
 * template<class Function>
 * void parallelize(size_t n, Function f, size_t nworkers) {
 *     size_t jobs_per_worker = std::ceil(static_cast<double>(n) / nworkers);
 *     size_t start = 0;
 *     std::vector<std::thread> jobs;
 *
 *     for (size_t w = 0; w < nworkers; ++w) {
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
 * Worker methods (i.e., `compute()`, `add()`, `finish()`) for a single worker instance do not have to be thread-safe.
 * That is, methods for a single worker will only ever be called in a serial fashion, as each worker is local to a thread.
 *
 * However, worker methods should be thread-safe across multiple worker instances.
 * For example, `apply()` may simultaneously call the `add()` method of different workers from different threads.
 * This is most relevant when results are being written back to the shared memory in `factory`,
 * and is usually easy to achieve as long as the memory spaces for the results of two target vectors do not overlap.
 *
 * The direct calculations are trivially parallelizable - it is assumed that they can be called independently for different target indices.
 * It is also assumed that each `compute()` call will only write to the shared memory at its supplied target index.
 *
 * To support parallelization in the running calculations, we require some overloads to the `dense_running()` and `sparse_running()` methods:
 *
 * - `dense_running()` now accepts two arguments: namely, the indices of the first and one-past-the-last target vectors to be processed.
 *   This should return a worker with an `add()` method that accepts a pointer to the subinterval of the non-target vector, corresponding to the target vectors to be processed.
 *   It may also have a `finish()` method, to be called to finalize any calculations after all non-target vectors are supplied.
 * - `sparse_running()` now accepts two arguments: namely, the indices of the first and one-past-the-last target vectors to be processed.
 *   The returned worker should have an `add()` method that accepts a `SparseRange` object, specifying the non-zero elements in the subinterval of the non-target vector.
 *   It may also have a `finish()` method, to be called to finalize any calculations after all non-target vectors are supplied.
 *
 * These overloads are optional and the function will fall back to serial processing if they are not supplied (and the function decides perform a running calculation).
 *
 * Calls to `dense_direct()` and counterparts should be thread-safe for any single instance of a factory class.
 */
template<int MARGIN, typename T, typename IDX, class Factory>
void apply(const Matrix<T, IDX>* p, Factory& factory, int threads = 1) {
    size_t NR = p->nrow(), NC = p->ncol();

    /* One might question why we use MARGIN in the template if we just convert
     * it to ROW here. This is because 'MARGIN' is used when we're doing
     * something to the matrix; 'ROW' is used for the representation of the
     * matrix itself. I'm keeping the distinction clear in the interface.
     */
    constexpr bool ROW = MARGIN == 0;

    const size_t dim = (ROW ? NR : NC);
    const size_t otherdim = (ROW ? NC : NR);

    /* If we support running calculations AND the preference 
     * is not consistent with the margin, we give it a shot.
     */
    if constexpr(stats::has_sparse_running<Factory>::value || stats::has_dense_running<Factory>::value) {
        if (p->prefer_rows() != ROW){
            if constexpr(stats::has_sparse_running<Factory>::value) {
                if (p->sparse()) {
                    stats::apply_sparse_running<MARGIN == 0>(dim, otherdim, p, factory, threads);
                    return;
                }
            }

            stats::apply_dense_running<MARGIN == 0>(dim, otherdim, p, factory, threads);
            return;
        }
    }

    if constexpr(stats::has_sparse_direct<Factory>::value) {
        if (p->sparse()) {
            stats::apply_sparse_direct<MARGIN == 0>(dim, otherdim, p, factory, threads);
            return;
        }
    }

    stats::apply_dense_direct<MARGIN == 0>(dim, otherdim, p, factory, threads);
    return; 
}

}

#endif
