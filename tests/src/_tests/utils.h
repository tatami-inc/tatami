#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <cmath>
#include "tatami/base/SparseRange.hpp"

inline std::pair<size_t, size_t> wrap_intervals(size_t first, size_t last, size_t max) {
    size_t diff = last - first;
    first %= max;
    last = std::min(max, first + diff);
    return std::make_pair(first, last);
}

template<typename T, typename IDX>
std::vector<T> expand(const tatami::SparseRangeCopy<T, IDX>& sparse, size_t dim) {
    std::vector<T> output(dim);
    for (size_t i = 0; i < sparse.index.size(); ++i) {
        output[sparse.index[i]] = sparse.value[i];
    }
    return output;
}

template<typename T, typename IDX>
std::vector<T> expand(const tatami::SparseRangeCopy<T, IDX>& sparse, size_t start, size_t end) {
    std::vector<T> output(end - start);
    for (size_t i = 0; i < sparse.index.size(); ++i) {
        output[sparse.index[i] - start] = sparse.value[i];
    }
    return output;
}

#endif
