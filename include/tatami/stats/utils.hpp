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
 * Callers can specify a custom parallelization scheme by defining a `TATAMI_CUSTOM_PARALLEL` function-like macro, which should accept four arguments:
 * `fun`, `tasks` and `threads` as below, as well as a `worker_size` argument that specifies the size of each task range
 * (the last is provided for convenience only, to avoid the need to re-compute it inside the macro).
 *
 * @tparam parallel_ Whether the tasks should be run in parallel.
 * If `false`, no parallelization is performed and all tasks are run on the current thread.
 * @tparam Function_ Function to be applied for a contiguous range of tasks.
 * This should accept three arguments:
 * - `thread`, the thread number executing this task range.
 * - `start`, the start index of the task range.
 * - `end`, the first index past the end of the task range.
 *
 * @param fun Function that executes a contiguous range of tasks.
 * @param tasks Number of tasks.
 * @param threads Number of threads.
 */
template<bool parallel_ = true, class Function_>
void parallelize(Function_ fun, size_t tasks, size_t threads) {
#if defined(_OPENMP) || defined(TATAMI_CUSTOM_PARALLEL)
    if constexpr(parallel_) {
        size_t worker_size = std::ceil(static_cast<double>(tasks) / static_cast<double>(threads));

        if (threads > 1) {
#ifndef TATAMI_CUSTOM_PARALLEL
            #pragma omp parallel for num_threads(threads)
            for (size_t t = 0; t < threads; ++t) {
                size_t start = worker_size * t, end = std::min(tasks, start + worker_size);
                if (start < end) {
                    fun(t, start, end);
                }
            }
#else
            TATAMI_CUSTOM_PARALLEL(std::move(fun), tasks, threads, worker_size);
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
 * @param from Index of the first row/column of the iteration range.
 * @param to Index of the last row/column of the iteration range.
 * @param args Further arguments to pass to `Matrix::dense_row` or `Matrix::dense_column`.
 *
 * @return An `Extractor` object for iteration over consecutive rows/columns in `[from, to)`.
 *
 * This function is equivalent to `new_extractor()` but additionally calls `Extractor::set_oracle()` with a `ConsecutiveOracle` instance.
 * `Matrix` implementations that are oracle-aware can then perform pre-fetching of future accesses for greater performance.
 * Of course, this assumes that the iteration over the target dimension does actually involve consecutive elements from `from` to `to`.
 */
template<bool row_, bool sparse_, typename Value_, typename Index_, typename ... Args_>
auto consecutive_extractor(const Matrix<Value_, Index_>* mat, Index_ from, Index_ to, Args_&&... args) {
    auto ext = new_extractor<row_, sparse_>(mat, std::forward<Args_>(args)...);
    if (mat->uses_oracle(row_)) {
        ext->set_oracle(std::make_unique<ConsecutiveOracle<Index_> >(from, to - from));
    }
    return ext;
}

}

#endif
