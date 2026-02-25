#ifndef TATAMI_PARALLELIZE_HPP
#define TATAMI_PARALLELIZE_HPP

#include <vector>
#include <cmath>

#ifndef TATAMI_CUSTOM_PARALLEL
#include "subpar/subpar.hpp"
#endif

/**
 * @file parallelize.hpp
 *
 * @brief Parallelized iteration over a `tatami::Matrix`.
 */

namespace tatami {

/**
 * Apply a function to a set of tasks in parallel, usually for iterating over a dimension of a `Matrix`.
 * By default, this uses `subpar::parallelize_range()` internally, which splits the integer sequence `[0, tasks)` into `K` non-overlapping ranges where `K <= workers`.
 * Each range of tasks is passed to `fun()` for parallel execution using OpenMP if available and `<thread>` otherwise.
 * Task splitting is deterministic, i.e., given the same `tasks` and `workers`, `fun()` will be called with the same combinations of worker IDs and task ranges.
 *
 * Advanced users can override the default parallelization mechanism by defining `TATAMI_CUSTOM_PARALLEL`.
 * This should be a function-like macro or the name of a function that accepts the `fun`, `tasks` and `workers` arguments (see below),
 * deterministically splits tasks into a set of ranges, calls `fun()` on each range with any parallelization scheme, and then returns the number of used workers.
 * Once `TATAMI_CUSTOM_PARALLEL` is defined (and if `parallel_ = true`), any call to `parallelize()` will invoke the user-defined scheme instead.
 *
 * @tparam parallel_ Whether the tasks should be run in parallel.
 * If `false`, no parallelization is performed and all tasks are run on the current worker.
 * @tparam Function_ Function to be applied for a contiguous range of tasks.
 * This should accept three arguments:
 * - `worker`, the worker ID executing this task range.
 *   This will be passed as an `int` in `[0, workers)`.
 * - `task_start`, the start index of the task range.
 *   This will be passed as an `Index_` in `[0, tasks)`.
 * - `task_length`, the number of tasks in the task range.
 *   This will be passed as an `Index_` in `(0, tasks)`, i.e., it is always positive.
 * @tparam Index_ Integer type for the number of tasks.
 *
 * @param fun Function that executes a contiguous range of tasks.
 * This will be called no more than once in each worker with a different non-overlapping range, where the union of all ranges will cover `[0, tasks)`. 
 * @param tasks Number of tasks.
 * This should be non-negative.
 * @param workers Number of workers.
 * This should be positive.
 *
 * @return The number of workers (`K`) that were actually used.
 * `K` is guaranteed to be no greater than `workers`.
 * `fun()` will have been called once for each of the worker IDs `[0, ..., K - 1]`.
 */
template<bool parallel_ = true, class Function_, typename Index_>
int parallelize(Function_ fun, const Index_ tasks, const int workers) {
    if constexpr(parallel_) {
#ifdef TATAMI_CUSTOM_PARALLEL
        return TATAMI_CUSTOM_PARALLEL(fun, tasks, workers);
#else
        return subpar::parallelize_range(workers, tasks, std::move(fun));
#endif
    } else {
        fun(0, 0, tasks);
        return 1;
    }
}

}

#endif
