#ifndef TATAMI_TEST_UTILS_HPP
#define TATAMI_TEST_UTILS_HPP

#include <cmath>
#include <vector>
#include <gtest/gtest.h>

#include "../tatami/base/SparseRange.hpp"
#include "../tatami/base/Options.hpp"

namespace tatami_test {

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

template<bool sanitize, typename T>
void sanitize_nan(std::vector<T>& values, T replacement = 1234567890) {
    if constexpr(sanitize) {
        for (auto& x : values) {
            if (std::isnan(x)) {
                x = replacement;
            }
        }
    }
}

template<class Function_>
void throws_error(Function_ fun, const std::string& msg) {
    try {
        fun();
        FAIL() << "expected error message '" << msg << "', got no error";
    } catch (std::exception& e) {
        std::string observed(e.what());
        if (observed.find(msg) == std::string::npos) {
            FAIL() << "expected error message '" << msg << "', got '" << observed << "'";
        }
    }
}

}

#endif
