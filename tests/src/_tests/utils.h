#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <cmath>
#include "tatami/base/SparseRange.hpp"
#include "tatami/base/Options.hpp"

typedef tatami::IterationOptions<int> IterationOptions;
typedef tatami::ExtractionOptions<int> ExtractionOptions;

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

template<typename IDX>
bool is_increasing(const std::vector<IDX>& indices) {
    for (size_t i = 1; i < indices.size(); ++i) {
        if (indices[i] <= indices[i-1]) {
            return false;
        }
    }
    return true;
}

inline void set_access_order(IterationOptions& opt, std::vector<int>& indices, int length, bool forward, int jump) {
    if (jump == 1) {
        if (forward) {
            opt.access_order = tatami::AccessOrder::CONSECUTIVE;
        } else {
            opt.access_order = tatami::AccessOrder::SEQUENCE;
            indices.resize(length);
            std::iota(indices.rbegin(), indices.rend(), 0);
            opt.sequence_start = indices.data();
            opt.sequence_length = indices.size();
        }
    } else {
        if (forward) {
            for (size_t i = 0; i < length; i += jump) {
                indices.push_back(i);
            }
            opt.access_order = tatami::AccessOrder::SEQUENCE;
            opt.sequence_start = indices.data();
            opt.sequence_length = indices.size();
        } else {
            opt.access_order = tatami::AccessOrder::RANDOM;
        }
    }
}

#endif
