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
 * By default, this uses `subpar::parallelize_range()` internally, which uses OpenMP if available and `<thread>` otherwise.
 * Advanced users can override the default parallelization mechanism by defining a `TATAMI_CUSTOM_PARALLEL` function-like macro. 
 * The macro should accept the `fun`, `tasks` and `threads` arguments as described below.
 *
 * @tparam parallel_ Whether the tasks should be run in parallel.
 * If `false`, no parallelization is performed and all tasks are run on the current thread.
 * @tparam Function_ Function to be applied for a contiguous range of tasks.
 * This should accept three arguments:
 * - `thread`, the thread number executing this task range.
 *   This will be passed as an `int`.
 * - `task_start`, the start index of the task range.
 *   This will be passed as an `Index_`.
 * - `task_length`, the number of tasks in the task range.
 *   This will be passed as an `Index_`.
 * @tparam Index_ Integer type for the number of tasks.
 *
 * @param fun Function that executes a contiguous range of tasks.
 * @param tasks Number of tasks.
 * @param threads Number of threads.
 */
template<bool parallel_ = true, class Function_, typename Index_>
void parallelize(Function_ fun, const Index_ tasks, const int threads) {
    if constexpr(parallel_) {
#ifdef TATAMI_CUSTOM_PARALLEL
        TATAMI_CUSTOM_PARALLEL(fun, tasks, threads);
#else
        subpar::parallelize_range(threads, tasks, std::move(fun));
#endif
    } else {
        fun(0, 0, tasks);
    }
}

}

#endif
