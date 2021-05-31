#ifndef TEST_CORE_H 
#define TEST_CORE_H 

#include <gtest/gtest.h>

#include <vector>
#include <algorithm>

class TestCore : public ::testing::Test {
protected:
    std::vector<double> expected;
    std::vector<double> output;

    std::vector<int> outidx;
    std::vector<double> outval;

    size_t first = 0, last = 0;

protected:
    size_t set_sizes(size_t f, size_t l) {
        first = f;
        last = l;
        size_t n = l - f;
        output.resize(n);
        outidx.resize(n);
        outval.resize(n);
        expected.resize(n);
        return n;
    }

    // Buffer handling utilities.
    void wipe_output() {
        std::fill(output.begin(), output.end(), 123);
        return;
    }

    void wipe_expected() {
        std::fill(expected.begin(), expected.end(), 123);
        return;
    }

    void fill_expected(const double* ptr) {
        if (ptr!=expected.data()){ 
            std::copy(ptr, ptr + last - first, expected.begin());
        }
        return;
    }

    void fill_output(const double* ptr) {
        if (ptr!=output.data()) {
            std::copy(ptr, ptr + last - first , output.begin());
        }
        return;
    }

    // Sparse utilities.
    void wipe_sparse_buffers() {
        std::fill(outval.begin(), outval.end(), 123);
        std::fill(outidx.begin(), outidx.end(), 456);
        return;
    }

    void fill_sparse_output(const tatami::sparse_range<double, int>& info) {
        std::fill(output.begin(), output.begin() + last - first, 0);
        for (size_t i = 0; i < info.number; ++i) {
            output[info.index[i] - first] = info.value[i];
        }
        return;
    }
};

#endif
