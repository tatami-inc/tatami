#ifndef CUSTOM_PARALLEL_H
#define CUSTOM_PARALLEL_H

#include <vector>
#include <thread>
#include <mutex>
#include <cmath>

template<class Function_>
static void custom_parallelize(Function_ f, size_t ntasks, size_t nworkers) {
    size_t start = 0;
    std::vector<std::thread> jobs;
    jobs.reserve(nworkers);
    size_t jobs_per_worker = std::ceil(static_cast<double>(ntasks) / static_cast<double>(nworkers));

    for (size_t w = 0; w < nworkers; ++w) {
        size_t end = std::min(ntasks, start + jobs_per_worker);
        if (start >= end) {
            break;
        }
        jobs.emplace_back(f, w, start, end - start);
        start += jobs_per_worker;
    }

    for (auto& job : jobs) {
        job.join();
    }
}

#define TATAMI_CUSTOM_PARALLEL custom_parallelize

#endif
