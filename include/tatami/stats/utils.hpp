#ifndef TATAMI_STATS_UTILS_HPP
#define TATAMI_STATS_UTILS_HPP

#include "../base/Matrix.hpp"
#include "../base/utils.hpp"
#include "../utils/Oracles.hpp"

#include <cmath>

/**
 * @file utils.hpp
 *
 * @brief Utilities for computing matrix statistics.
 */

namespace tatami {

/**
 * Apply a function to a set of tasks, distributing them to threads via OpenMP if enabled.
 * Callers can specify a custom parallelization scheme by defining a `TATAMI_CUSTOM_PARALLEL` function-like macro, 
 * which should accept the `fun`, `tasks` and `threads` arguments as below.
 *
 * @tparam parallel_ Whether the tasks should be run in parallel.
 * If `false`, no parallelization is performed and all tasks are run on the current thread.
 * @tparam Function_ Function to be applied for a contiguous range of tasks.
 * This should accept three arguments:
 * - `thread`, the thread number executing this task range.
 *   This will be passed as a `size_t`.
 * - `task_start`, the start index of the task range.
 *   This will be passed as a `Index_`.
 * - `task_length`, the number of tasks in the task range.
 *   This will be passed as a `Index_`.
 * @tparam Index_ Integer type for the number of tasks.
 *
 * @param fun Function that executes a contiguous range of tasks.
 * @param tasks Number of tasks.
 * @param threads Number of threads.
 */
template<bool parallel_ = true, class Function_, typename Index_>
void parallelize(Function_ fun, Index_ tasks, size_t threads) {
#if defined(_OPENMP) || defined(TATAMI_CUSTOM_PARALLEL)
    if constexpr(parallel_) {

        if (threads > 1) {
#ifndef TATAMI_CUSTOM_PARALLEL
            Index_ worker_size = (tasks / threads) + (tasks % threads > 0); // Ceiling of an integer division.
            threads = (tasks / worker_size) + (tasks % worker_size > 0); // Set the actual number of required threads.

            #pragma omp parallel for num_threads(threads)
            for (size_t t = 0; t < threads; ++t) {
                Index_ start = worker_size * t; // Will not overflow due to the above recomputation of 'threads'.
                Index_ remaining = tasks - start; // Must be positive, as otherwise 'tasks % worker_size = 0' and the iteration wouldn't even get here.
                fun(t, start, std::min(remaining, worker_size)); // Use 'remaining' to avoid potential overflow from computing 'end = start + worker_size'.
            }
#else
            TATAMI_CUSTOM_PARALLEL(std::move(fun), tasks, threads);
#endif
            return;
        }
    }
#endif

    fun(0, 0, tasks);
    return;
}

/**
 * @tparam row_ Whether to perform extraction on rows.
 * @tparam sparse_ Whether to perform sparse retrieval.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column index.
 * @tparam Args_ Types of further arguments to pass to `Matrix::dense_row` or `Matrix::dense_column`.
 *
 * @param mat Matrix to iterate over.
 * @param iter_start Index of the first row/column of the iteration range.
 * @param iter_length Number of rows/columns in the iteration range.
 * @param args Further arguments to pass to `Matrix::dense_row` or `Matrix::dense_column`.
 *
 * @return An `Extractor` object for iteration over consecutive rows/columns in `[iter_start, iter_start + iter_length)`.
 *
 * This function is equivalent to `new_extractor()` but additionally calls `Extractor::set_oracle()` with a `ConsecutiveOracle` instance.
 * `Matrix` implementations that are oracle-aware can then perform pre-fetching of future accesses for greater performance.
 * Of course, this assumes that the iteration over the target dimension does actually involve consecutive elements from `iter_start` to `iter_start + iter_length`.
 */
template<bool row_, bool sparse_, typename Value_, typename Index_, typename ... Args_>
auto consecutive_extractor(const Matrix<Value_, Index_>* mat, Index_ iter_start, Index_ iter_length, Args_&&... args) {
    auto ext = new_extractor<row_, sparse_>(mat, std::forward<Args_>(args)...);
    if (mat->uses_oracle(row_)) {
        ext->set_oracle(std::make_unique<ConsecutiveOracle<Index_> >(iter_start, iter_length));
    }
    return ext;
}

}

#endif
