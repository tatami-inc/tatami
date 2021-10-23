#include <gtest/gtest.h>

#include "tatami/ext/MatrixMarket.hpp"

#include <cstdio>
#include "zlib.h"

#include <limits>
#include <random>
#include "write_matrix_market.h"

class MatrixMarketGzipTest : public ::testing::Test {
protected:
    size_t NR = 2000, NC = 1000;
    std::vector<int> vals, rows, cols;
    std::string path, gzname;

protected:
    void SetUp() {
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

        path = temp_file_path("tatami-tests-ext-MatrixMarketGzip.mtx");
        write_matrix_market(path, NR, NC, vals, rows, cols);

        // Now applying compression.
        FILE * ihandle = std::fopen(path.c_str(), "rb");
        gzname = path + ".gz";
        gzFile ohandle = gzopen(gzname.c_str(), "w");
        constexpr int CHUNK = 16384;
        std::vector<unsigned char> in(CHUNK);
        size_t nread = 0;
        while ((nread = std::fread(in.data(), 1, CHUNK, ihandle))) {
            gzwrite(ohandle, in.data(), nread);
        }
        std::fclose(ihandle);
        gzclose(ohandle);

        return;
    }
};

TEST_F(MatrixMarketGzipTest, SimpleTest) {
    auto loaded = tatami::MatrixMarket::load_sparse_matrix_gzip(gzname.c_str());

    auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
    typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(vals), decltype(rows), decltype(indptrs)> SparseMat; 
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new SparseMat(NR, NC, std::move(vals), std::move(rows), std::move(indptrs)));

    for (size_t i = 0; i < NC; ++i) {
        auto stuff = loaded->column(i);
        EXPECT_EQ(stuff, ref->column(i));
    }
}

TEST_F(MatrixMarketGzipTest, LayeredTest) {
    // Loading everyone from file.
    auto loaded = tatami::MatrixMarket::load_layered_sparse_matrix(path.c_str());
    auto loaded_gz = tatami::MatrixMarket::load_layered_sparse_matrix_gzip(gzname.c_str());

    auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
    typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(vals), decltype(rows), decltype(indptrs)> SparseMat; 
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new SparseMat(NR, NC, std::move(vals), std::move(rows), std::move(indptrs)));

    EXPECT_EQ(ref->nrow(), loaded.matrix->nrow());
    EXPECT_EQ(ref->nrow(), loaded_gz.matrix->nrow());
    EXPECT_EQ(ref->ncol(), loaded.matrix->ncol());
    EXPECT_EQ(ref->ncol(), loaded_gz.matrix->ncol());
    EXPECT_TRUE(loaded.matrix->sparse());
    EXPECT_TRUE(loaded_gz.matrix->sparse());
    EXPECT_FALSE(loaded.matrix->prefer_rows());
    EXPECT_FALSE(loaded_gz.matrix->prefer_rows());

    for (size_t i = 0; i < NC; ++i) {
        int adjusted = loaded.permutation[i];
        EXPECT_EQ(adjusted, loaded_gz.permutation[i]);

        auto stuff = loaded.matrix->column(adjusted);
        EXPECT_EQ(stuff, loaded_gz.matrix->column(adjusted));
        EXPECT_EQ(stuff, ref->column(i));
    }
}
