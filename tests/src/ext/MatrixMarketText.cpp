#include <gtest/gtest.h>

#include "tatami/ext/MatrixMarket.hpp"
#include "tatami/ext/MatrixMarket_layered.hpp"
#include "write_matrix_market.h"

#include <limits>
#include <fstream>
#include <string>
#include <vector>
#include <random>

class MatrixMarketTextTest : public ::testing::TestWithParam<int> {
protected:
    size_t NR = 2000, NC = 1000;
    std::vector<int> vals, rows, cols;
    std::string path;
    int bufsize;

protected:
    auto extra_assemble() {
        std::mt19937_64 rng(1234567890);
        std::uniform_real_distribution<> unif(0.0, 1.0);
        for (size_t i = 0; i < NC; ++i) {
            for (size_t j = 0; j < NR; ++j) {
                if (unif(rng) < 0.05) {
                    rows.push_back(j);
                    cols.push_back(i);
                    vals.push_back(unif(rng) * 100);
                }
            }
        }

        path = temp_file_path("tatami-tests-ext-MatrixMarket.mtx");
        write_matrix_market(path, NR, NC, vals, rows, cols);
        return path;
    }
};

TEST_P(MatrixMarketTextTest, Simple) {
    auto filepath = extra_assemble();
    auto bufsize = GetParam();
    auto out = tatami::MatrixMarket::load_sparse_matrix(filepath.c_str(), bufsize);

    // Checking against a reference.
    auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
    typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(vals), decltype(rows), decltype(indptrs)> SparseMat; 
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new SparseMat(NR, NC, std::move(vals), std::move(rows), std::move(indptrs))); 

    for (size_t i = 0; i < NC; ++i) {
        auto stuff = out->column(i);
        EXPECT_EQ(stuff, ref->column(i));
    }
}

TEST_P(MatrixMarketTextTest, Layered) {
    auto filepath = extra_assemble();
    auto bufsize = GetParam();

    // Checking against a reference.
    auto loaded = tatami::MatrixMarket::load_layered_sparse_matrix(filepath.c_str(), bufsize);
    const auto& out = loaded.matrix;

    auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
    typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(vals), decltype(rows), decltype(indptrs)> SparseMat; 
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new SparseMat(NR, NC, std::move(vals), std::move(rows), std::move(indptrs))); 

    EXPECT_EQ(out->nrow(), ref->nrow());
    EXPECT_EQ(out->ncol(), ref->ncol());
    EXPECT_TRUE(out->sparse());
    EXPECT_FALSE(out->prefer_rows());

    for (size_t i = 0; i < NR; ++i) {
        int adjusted = loaded.permutation[i];
        auto stuff = out->row(adjusted);
        EXPECT_EQ(stuff, ref->row(i));
    }
}

INSTANTIATE_TEST_CASE_P(
    MatrixMarket,
    MatrixMarketTextTest,
    ::testing::Values(50, 100, 1000, 10000)
);
