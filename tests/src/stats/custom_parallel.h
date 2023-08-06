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
    std::vector<std::string> errors(nworkers);

    for (size_t w = 0; w < nworkers; ++w) {
        size_t end = std::min(ntasks, start + jobs_per_worker);
        if (start >= end) {
            break;
        }
        jobs.emplace_back([&f,&errors](size_t t, size_t start, size_t len) -> void {
            try {
                f(t, start, len);
            } catch (std::exception& e) {
                errors[t] = e.what();
            } catch (...) {
                errors[t] = "unknown error";
            }
        }, w, start, end - start);
        start += jobs_per_worker;
    }

    for (auto& job : jobs) {
        job.join();
    }

    for (auto& e : errors) {
        if (!e.empty()) {
            throw std::runtime_error(e);
        }
    }
}

#define TATAMI_CUSTOM_PARALLEL custom_parallelize

#endif
