#include <gtest/gtest.h>

#include "tatami/ext/MatrixMarket.hpp"
#include "tatami/ext/MatrixMarket_layered.hpp"

#include <cstdio>
#include "zlib.h"

#include <limits>
#include <random>

#include "temp_file_path.h"
#include "write_matrix_market.h"

class MatrixMarketGzipTest : public ::testing::TestWithParam<int> {
protected:
    size_t NR = 2000, NC = 1000;
    std::vector<int> vals, rows, cols;
    std::vector<unsigned char> contents;
    std::string gzname;

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

        std::stringstream ss;
        write_matrix_market(ss, NR, NC, vals, rows, cols);
        ss >> std::noskipws;
        contents.insert(contents.end(), std::istream_iterator<unsigned char>{ss}, std::istream_iterator<unsigned char>());

        gzname = temp_file_path("tatami-tests-ext-MatrixMarketGzip.mtx.gz");
        gzFile ohandle = gzopen(gzname.c_str(), "w");
        gzwrite(ohandle, contents.data(), contents.size());
        gzclose(ohandle);

        return;
    }
};

TEST_P(MatrixMarketGzipTest, SimpleTest) {
    auto bufsize = GetParam();

    auto loaded = tatami::MatrixMarket::load_sparse_matrix_from_file(gzname.c_str(), 1, bufsize);
    EXPECT_EQ(loaded->nrow(), NR);
    EXPECT_EQ(loaded->ncol(), NC);

    auto auto_loaded = tatami::MatrixMarket::load_sparse_matrix_from_file(gzname.c_str(), -1, bufsize); // auto-detect.
    EXPECT_EQ(auto_loaded->nrow(), NR);
    EXPECT_EQ(auto_loaded->ncol(), NC);

    auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
    typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(vals), decltype(rows), decltype(indptrs)> SparseMat; 
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new SparseMat(NR, NC, std::move(vals), std::move(rows), std::move(indptrs)));

    for (size_t i = 0; i < NC; ++i) {
        auto stuff = loaded->column(i);
        EXPECT_EQ(stuff, ref->column(i));
        auto auto_stuff = auto_loaded->column(i);
        EXPECT_EQ(auto_stuff, ref->column(i));
    }
}

TEST_P(MatrixMarketGzipTest, SimpleBufferTest) {
    auto bufsize = GetParam();

    auto ref = tatami::MatrixMarket::load_sparse_matrix_from_buffer(contents.data(), contents.size());
    EXPECT_EQ(ref->nrow(), NR);
    EXPECT_EQ(ref->ncol(), NC);

    std::ifstream handle(gzname, std::ios_base::binary);
    handle >> std::noskipws;
    std::vector<unsigned char> gzcontents(std::istream_iterator<unsigned char>{handle}, std::istream_iterator<unsigned char>());
    EXPECT_TRUE(gzcontents.size() < contents.size()); // compression should have an effect.

    auto obs = tatami::MatrixMarket::load_sparse_matrix_from_buffer(gzcontents.data(), gzcontents.size(), 1, bufsize);
    EXPECT_EQ(obs->nrow(), NR);
    EXPECT_EQ(obs->ncol(), NC);

    auto obs_auto = tatami::MatrixMarket::load_sparse_matrix_from_buffer(gzcontents.data(), gzcontents.size(), -1, bufsize); // now with autodetection.
    EXPECT_EQ(obs_auto->nrow(), NR);
    EXPECT_EQ(obs_auto->ncol(), NC);

    for (size_t i = 0; i < NC; ++i) {
        auto stuff = obs->column(i);
        EXPECT_EQ(stuff, ref->column(i));
        auto auto_stuff = obs_auto->column(i);
        EXPECT_EQ(auto_stuff, ref->column(i));
    }
}

TEST_P(MatrixMarketGzipTest, LayeredTest) {
    auto bufsize = GetParam();

    auto loaded = tatami::MatrixMarket::load_layered_sparse_matrix_from_buffer(contents.data(), contents.size());
    auto loaded_gz = tatami::MatrixMarket::load_layered_sparse_matrix_from_file(gzname.c_str(), 1, bufsize);

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

    // Auto-detection also works.
    auto auto_loaded = tatami::MatrixMarket::load_layered_sparse_matrix_from_file(gzname.c_str(), -1, bufsize); 
    EXPECT_EQ(auto_loaded.matrix->nrow(), NR);
    EXPECT_EQ(auto_loaded.matrix->ncol(), NC);
}

TEST_P(MatrixMarketGzipTest, LayeredBufferTest) {
    auto bufsize = GetParam();

    auto ref = tatami::MatrixMarket::load_layered_sparse_matrix_from_buffer(contents.data(), contents.size());
    EXPECT_EQ(ref.matrix->nrow(), NR);
    EXPECT_EQ(ref.matrix->ncol(), NC);

    std::ifstream handle(gzname, std::ios_base::binary);
    handle >> std::noskipws;
    std::vector<unsigned char> gzcontents(std::istream_iterator<unsigned char>{handle}, std::istream_iterator<unsigned char>());
    EXPECT_TRUE(gzcontents.size() < contents.size()); // compression should have an effect.

    auto obs = tatami::MatrixMarket::load_layered_sparse_matrix_from_buffer(gzcontents.data(), gzcontents.size(), 1, bufsize);
    EXPECT_EQ(obs.matrix->nrow(), NR);
    EXPECT_EQ(obs.matrix->ncol(), NC);

    for (size_t i = 0; i < NC; ++i) {
        auto stuff = obs.matrix->column(i);
        EXPECT_EQ(stuff, ref.matrix->column(i));
    }

    EXPECT_EQ(obs.permutation, ref.permutation);

    // Auto-detection also works.
    auto auto_obs = tatami::MatrixMarket::load_layered_sparse_matrix_from_buffer(gzcontents.data(), gzcontents.size(), -1, bufsize);
    EXPECT_EQ(auto_obs.matrix->nrow(), NR);
    EXPECT_EQ(auto_obs.matrix->ncol(), NC);
}

INSTANTIATE_TEST_CASE_P(
    MatrixMarket,
    MatrixMarketGzipTest,
    ::testing::Values(50, 100, 1000, 10000)
);
