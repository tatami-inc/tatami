#ifndef TATAMI_STATS_PARALLELIZE_H
#define TATAMI_STATS_PARALLELIZE_H

#include "../base/Matrix.hpp"
#include "../base/utils.hpp"
#include "../utils/Oracles.hpp"

#include <cmath>

/**
 * @file parallelize.hpp
 *
 * @brief Parallelize tasks across threads.
 */

namespace tatami {

/**
 * Apply a function to a set of tasks, distributing them to threads via OpenMP if enabled.
 * Callers can specify a custom parallelization scheme by defining a `TATAMI_CUSTOM_PARALLEL` function-like macro, which should accept four arguments:
 * `fun`, `num_tasks` and `num_threads` as below, as well as a `worker_size` argument that specifies the size of each task range
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
 * @param num_tasks Number of tasks.
 * @param num_threads Number of threads.
 */
template<bool parallel_ = true, class Function_>
void parallelize(Function_ fun, size_t num_tasks, size_t num_threads) {
#if defined(_OPENMP) || defined(TATAMI_CUSTOM_PARALLEL)
    if constexpr(parallel_) {
        size_t worker_size = std::ceil(static_cast<double>(num_tasks) / static_cast<double>(num_threads));

        if (num_threads > 1) {
#ifndef TATAMI_CUSTOM_PARALLEL
            #pragma omp parallel for num_threads(threads)
            for (size_t t = 0; t < threads; ++t) {
                size_t start = worker_size * t, end = std::min(dim, start + worker_size);
                if (start < end) {
                    fun(t, start, end);
                }
            }
#else
            TATAMI_CUSTOM_PARALLEL(std::move(fun), num_tasks, num_threads, worker_size);
#endif
            return;
        }
    }
#endif

    fun(0, 0, num_tasks);
    return;
}

/**
 * @cond
 */
template<bool row_, bool sparse_, typename Value_, typename Index_>
auto direct_extractor(const Matrix<Value_, Index_>* mat, bool uses_oracle, Index_ from, Index_ to, const Options& opt) {
    if constexpr(sparse_) {
        auto ext = (row_ ? mat->sparse_row(opt) : mat->sparse_column(opt));
        if (uses_oracle) {
            ext->set_oracle(std::make_unique<ConsecutiveOracle<Index_> >(from, to - from));
        }
        return ext;
    } else {
        auto ext = (row_ ? mat->dense_row(opt) : mat->dense_column(opt));
        if (uses_oracle) {
            ext->set_oracle(std::make_unique<ConsecutiveOracle<Index_> >(from, to - from));
        }
        return ext;
    }
}

template<bool row_, bool sparse_, typename Value_, typename Index_>
auto running_extractor(const Matrix<Value_, Index_>* mat, Index_ start, Index_ len, bool uses_oracle, Index_ otherdim, const Options& opt) {
    if constexpr(sparse_) {
        auto ext = (row_ ? mat->sparse_column(start, len, opt) : mat->sparse_row(start, len, opt));
        if (uses_oracle) {
            ext->set_oracle(std::make_unique<ConsecutiveOracle<Index_> >(0, otherdim));
        }
        return ext;
    } else {
        auto ext = (row_ ? mat->dense_column(start, len, opt) : mat->dense_row(start, len, opt));
        if (uses_oracle) {
            ext->set_oracle(std::make_unique<ConsecutiveOracle<Index_> >(0, otherdim));
        }
        return ext;
    }
}
/**
 * @endcond
 */

}

#endif
