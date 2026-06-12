#ifndef CUSTOM_PARALLEL_H
#define CUSTOM_PARALLEL_H

// Put this before any tatami imports.
#ifdef CUSTOM_PARALLEL_TEST
#include <thread>
#include <vector>

#include "sanisizer/sanisizer.hpp"

// To test the override capabilities, we create a new function where the first worker gets the last range.
// This is because some of our conversion functions have unusual work-sharing arrangements when filling the output arrays.
// We need to check that we are not making assumptions about the first range being assigned to thread 0,
template<class Function_, typename Task_>
int weird_parallelize(Function_ run_task_range, const Task_ num_tasks, const int num_workers) {
    if (num_tasks <= 0) {
        return 0;
    }

    if (num_workers <= 1 || num_tasks == 1) {
        run_task_range(0, 0, num_tasks);
        return 1;
    }

    // All workers with indices below 'remainder' get an extra task to fill up the remainder.
    Task_ tasks_per_worker = 1;
    int remainder = 0;
    if (sanisizer::is_greater_than_or_equal(num_workers, num_tasks)) {
        num_workers = num_tasks;
    } else {
        tasks_per_worker = num_tasks / num_workers;
        remainder = num_tasks % num_workers;
    }

    std::vector<std::thread> workers;
    sanisizer::reserve(workers, num_workers); // preallocate to ensure we don't get alloc errors during emplace_back().

    // Worker 0 has the last range, etc., and the last worker gets the first range.
    for (int w = 0; w < num_workers; ++w) {
        const auto range_id = num_workers - w - 1;
        const Task_ start = range_id * tasks_per_worker + (range_id < remainder ? range_id : remainder);
        const Task_ length = tasks_per_worker + (range_id < remainder); 
        workers.emplace_back(run_task_range, w, start, length);
    }

    for (auto& wrk : workers) {
        wrk.join();
    }

    return num_workers;
}

// Overriding the default tatami parallelization.
#define TATAMI_CUSTOM_PARALLELIZE(fun, tasks, workers) weird_parallelize(fun, tasks, workers)
#endif

#endif
