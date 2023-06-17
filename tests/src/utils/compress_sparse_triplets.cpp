#include <gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <random>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/utils/compress_sparse_triplets.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

template<class V, class U>
void permuter(const U& values, const V& rows, const V& cols, U& values2, V& rows2, V& cols2) {
    // Creating a shuffling permutator.
    std::vector<int> permutation(rows.size());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::mt19937_64 engine;
    std::shuffle(permutation.begin(), permutation.end(), engine);

    rows2.reserve(rows.size());
    cols2.reserve(rows.size());
    values2.reserve(rows.size());

    for (auto p : permutation) {
        values2.push_back(values[p]);
        rows2.push_back(rows[p]);
        cols2.push_back(cols[p]);
    }

    return;
}

TEST(compress_sparse_triplets, CompressionByColumn) {
    size_t NR = 100, NC = 50;
    auto simulated = tatami_test::simulate_sparse_compressed<double>(NC, NR, 0.1);
    const auto& rows = simulated.index;
    const auto& values = simulated.value;

    std::vector<int> cols;
    for (size_t c = 0; c < NC; ++c) {
        size_t n = simulated.ptr[c + 1] - simulated.ptr[c];
        cols.insert(cols.end(), n, c);
    }

    // Permuting.
    std::vector<int> rows2, cols2;
    std::vector<double> values2;
    permuter(values, rows, cols, values2, rows2, cols2);

    EXPECT_NE(rows2, rows);
    EXPECT_NE(cols2, cols);
    EXPECT_NE(values2, values);

    // Unpermuting them.
    auto output = tatami::compress_sparse_triplets<false>(NR, NC, values2, rows2, cols2);
    EXPECT_EQ(rows2, rows);
    EXPECT_EQ(cols2, cols);
    EXPECT_EQ(values2, values);
    EXPECT_EQ(output, simulated.ptr);
}

TEST(compress_sparse_triplets, CompressionByRow) {
    size_t NR = 80, NC = 60;
    auto simulated = tatami_test::simulate_sparse_compressed<double>(NR, NC, 0.1);
    const auto& cols = simulated.index;
    const auto& values = simulated.value;

    std::vector<int> rows;
    for (size_t r = 0; r < NR; ++r) {
        size_t n = simulated.ptr[r + 1] - simulated.ptr[r];
        rows.insert(rows.end(), n, r);
    }

    // Permuting.
    std::vector<int> rows2, cols2;
    std::vector<double> values2;
    permuter(values, rows, cols, values2, rows2, cols2);

    EXPECT_NE(rows2, rows);
    EXPECT_NE(cols2, cols);
    EXPECT_NE(values2, values);

    // Unpermuting them.
    auto output = tatami::compress_sparse_triplets<true>(NR, NC, values2, rows2, cols2);
    EXPECT_EQ(rows2, rows);
    EXPECT_EQ(cols2, cols);
    EXPECT_EQ(values2, values);
    EXPECT_EQ(output, simulated.ptr);
}
