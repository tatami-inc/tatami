#ifndef SIMULATE_DENSE_H
#define SIMULATE_DENSE_H

#include <random>
#include <vector>

template<typename T>
std::vector<T> simulate_dense_vector(size_t length, double lower = 0, double upper = 100, size_t seed = 1234567890) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<> unif(lower, upper);

    std::vector<T> values(length);
    for (auto& v : values) {
        v = unif(rng);
    }
    return values;
}

#endif
