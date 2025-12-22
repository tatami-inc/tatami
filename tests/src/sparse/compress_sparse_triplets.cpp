#include <gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <random>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/compress_sparse_triplets.hpp"

#include "tatami_test/tatami_test.hpp"

template<class V, class U>
void permuter(const U& values, const V& rows, const V& cols, U& values2, V& rows2, V& cols2) {
    // Creating a shuffling permutator.
    std::vector<int> permutation(rows.size());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::mt19937_64 engine;
    std::shuffle(permutation.begin(), permutation.end(), engine);

    rows2.clear();
    cols2.clear();
    values2.clear();

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

template<class V, class U, class S>
void partial_permuter(const S& ptrs, const U& values, const V& secondary, U& values2, V& primary2, V& secondary2) {
    // Creating a shuffling permutator.
    std::vector<int> permutation;
    std::mt19937_64 engine;

    primary2.clear();
    secondary2.clear();
    values2.clear();

    primary2.reserve(values.size());
    secondary2.reserve(values.size());
    values2.reserve(values.size());

    auto n = ptrs.size() - 1;
    for (decltype(n) i = 0; i < n; ++i) {
        permutation.resize(ptrs[i + 1] - ptrs[i]);
        std::iota(permutation.begin(), permutation.end(), ptrs[i]);
        std::shuffle(permutation.begin(), permutation.end(), engine);
        for (auto p : permutation) {
            values2.push_back(values[p]);
            secondary2.push_back(secondary[p]);
        }
        primary2.insert(primary2.end(), permutation.size(), i);
    }
}

TEST(compress_sparse_triplets, CompressionByColumn) {
    const int NR = 100, NC = 50;
    auto simulated = tatami_test::simulate_compressed_sparse<double, int>(NC, NR, []{
        tatami_test::SimulateCompressedSparseOptions opt;
        opt.density = 0.1;
        opt.seed = 65537;
        return opt;
    }());

    const auto& rows = simulated.index;
    const auto& values = simulated.data;

    std::vector<int> cols;
    for (int c = 0; c < NC; ++c) {
        auto n = simulated.indptr[c + 1] - simulated.indptr[c];
        cols.insert(cols.end(), n, c);
    }

    std::vector<int> rows2, cols2;
    std::vector<double> values2;

    // Permuting and unpermuting.
    permuter(values, rows, cols, values2, rows2, cols2);

    EXPECT_NE(rows2, rows);
    EXPECT_NE(cols2, cols);
    EXPECT_NE(values2, values);

    auto output = tatami::compress_sparse_triplets<false>(NR, NC, values2, rows2, cols2);
    EXPECT_EQ(rows2, rows);
    EXPECT_EQ(cols2, cols);
    EXPECT_EQ(values2, values);
    EXPECT_EQ(output, simulated.indptr);

    // Same result if it was already sorted.
    values2 = values;
    rows2 = rows;
    cols2 = cols;

    auto output_sorted = tatami::compress_sparse_triplets<false>(NR, NC, values2, rows2, cols2);
    EXPECT_EQ(rows2, rows);
    EXPECT_EQ(cols2, cols);
    EXPECT_EQ(values2, values);
    EXPECT_EQ(output_sorted, simulated.indptr);

    // What if it was only partially sorted?
    partial_permuter(simulated.indptr, values, rows, values2, cols2, rows2);

    EXPECT_NE(rows2, rows);
    EXPECT_EQ(cols2, cols);
    EXPECT_NE(values2, values);

    auto output_partial = tatami::compress_sparse_triplets<false>(NR, NC, values2, rows2, cols2);
    EXPECT_EQ(rows2, rows);
    EXPECT_EQ(cols2, cols);
    EXPECT_EQ(values2, values);
    EXPECT_EQ(output_sorted, simulated.indptr);
}

TEST(compress_sparse_triplets, CompressionByRow) {
    const int NR = 80, NC = 60;
    auto simulated = tatami_test::simulate_compressed_sparse<double, int>(NR, NC, []{
        tatami_test::SimulateCompressedSparseOptions opt;
        opt.density = 0.1;
        opt.seed = 5040;
        return opt;
    }());

    const auto& cols = simulated.index;
    const auto& values = simulated.data;

    std::vector<int> rows;
    for (int r = 0; r < NR; ++r) {
        const auto n = simulated.indptr[r + 1] - simulated.indptr[r];
        rows.insert(rows.end(), n, r);
    }

    std::vector<int> rows2, cols2;
    std::vector<double> values2;

    // Permuting and unpermuting them.
    permuter(values, rows, cols, values2, rows2, cols2);

    EXPECT_NE(rows2, rows);
    EXPECT_NE(cols2, cols);
    EXPECT_NE(values2, values);

    auto output = tatami::compress_sparse_triplets<true>(NR, NC, values2, rows2, cols2);
    EXPECT_EQ(rows2, rows);
    EXPECT_EQ(cols2, cols);
    EXPECT_EQ(values2, values);
    EXPECT_EQ(output, simulated.indptr);

    // Same result if it was already sorted.
    values2 = values;
    rows2 = rows;
    cols2 = cols;

    auto output_sorted = tatami::compress_sparse_triplets<true>(NR, NC, values2, rows2, cols2);
    EXPECT_EQ(rows2, rows);
    EXPECT_EQ(cols2, cols);
    EXPECT_EQ(values2, values);
    EXPECT_EQ(output_sorted, simulated.indptr);

    // What if it was only partially sorted?
    partial_permuter(simulated.indptr, values, cols, values2, rows2, cols2);

    EXPECT_EQ(rows2, rows);
    EXPECT_NE(cols2, cols);
    EXPECT_NE(values2, values);

    auto output_partial = tatami::compress_sparse_triplets<true>(NR, NC, values2, rows2, cols2);
    EXPECT_EQ(rows2, rows);
    EXPECT_EQ(cols2, cols);
    EXPECT_EQ(values2, values);
    EXPECT_EQ(output_sorted, simulated.indptr);
}
