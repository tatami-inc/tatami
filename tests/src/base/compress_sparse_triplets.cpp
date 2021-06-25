#include <gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <random>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/utils/compress_sparse_triplets.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../data/data.h"

template<class V, class U>
void permuter(U& values, V& rows, V& cols, U& values2, V& rows2, V& cols2) {
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
    tatami::DenseRowMatrix<double> dense(sparse_nrow, sparse_ncol, sparse_matrix);
    auto sparse = tatami::convert_to_sparse(&dense, false);

    std::vector<int> rows, cols;
    std::vector<double> values;
    std::vector<int> ibuffer(sparse->nrow());
    std::vector<double> vbuffer(sparse->nrow());

    for (size_t c = 0; c < sparse->ncol(); ++c) {
        auto range = sparse->sparse_column(c, vbuffer.data(), ibuffer.data()); 
        rows.insert(rows.end(), range.index, range.index + range.number);
        values.insert(values.end(), range.value, range.value + range.number);
        cols.insert(cols.end(), range.number, c);
    }

    // Permuting.
    std::vector<int> rows2, cols2;
    std::vector<double> values2;
    permuter(values, rows, cols, values2, rows2, cols2);

    EXPECT_NE(rows2, rows);
    EXPECT_NE(cols2, cols);
    EXPECT_NE(values2, values);

    // Unpermuting them.
    auto output = tatami::compress_sparse_triplets<false>(sparse_nrow, sparse_ncol, values2, rows2, cols2);
    EXPECT_EQ(rows2, rows);
    EXPECT_EQ(cols2, cols);
    EXPECT_EQ(values2, values);
}

TEST(compress_sparse_triplets, CompressionByRow) {
    tatami::DenseRowMatrix<double> dense(sparse_nrow, sparse_ncol, sparse_matrix);
    auto sparse = tatami::convert_to_sparse(&dense, true);

    std::vector<int> rows, cols;
    std::vector<double> values;
    std::vector<int> ibuffer(sparse->ncol());
    std::vector<double> vbuffer(sparse->ncol());

    for (size_t r = 0; r < sparse->nrow(); ++r) {
        auto range = sparse->sparse_row(r, vbuffer.data(), ibuffer.data()); 
        cols.insert(cols.end(), range.index, range.index + range.number);
        values.insert(values.end(), range.value, range.value + range.number);
        rows.insert(rows.end(), range.number, r);
    }

    // Permuting.
    std::vector<int> rows2, cols2;
    std::vector<double> values2;
    permuter(values, rows, cols, values2, rows2, cols2);

    EXPECT_NE(rows2, rows);
    EXPECT_NE(cols2, cols);
    EXPECT_NE(values2, values);

    // Unpermuting them.
    auto output = tatami::compress_sparse_triplets<true>(sparse_nrow, sparse_ncol, values2, rows2, cols2);
    EXPECT_EQ(rows2, rows);
    EXPECT_EQ(cols2, cols);
    EXPECT_EQ(values2, values);
}
