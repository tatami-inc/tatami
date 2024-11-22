#include <gtest/gtest.h>

#include <vector>

#include "tatami_test/tatami_test.hpp"
#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/dense/transpose.hpp"

class TransposeDenseTest : public ::testing::TestWithParam<std::tuple<int, int> > {};

TEST_P(TransposeDenseTest, Basic) {
    auto params = GetParam();
    auto NR = std::get<0>(params);
    auto NC = std::get<1>(params);

    auto values = tatami_test::simulate_vector<double>(NR * NC, [&]{
        tatami_test::SimulateVectorOptions opt;
        opt.lower = 0;
        opt.upper = 100;
        opt.seed = NR * 10 + NC;
        return opt;
    }());

    std::vector<double> buffer(NR * NC);
    tatami::transpose(values.data(), NR, NC, buffer.data());

    tatami::DenseRowMatrix<double, int> original(NR, NC, std::move(values));
    tatami::DenseColumnMatrix<double, int> flipped(NR, NC, std::move(buffer));

    auto oext = original.dense_row();
    auto fext = flipped.dense_row();
    for (int r = 0; r < NR; ++r) {
        auto oout = tatami_test::fetch(*oext, r, NC);
        auto fout = tatami_test::fetch(*fext, r, NC);
        ASSERT_EQ(oout, fout);
    }
}

TEST_P(TransposeDenseTest, Strided) {
    auto params = GetParam();
    auto NR = std::get<0>(params);
    auto NC = std::get<1>(params);
    auto stride_nr = NR + 17;
    auto stride_nc = NC + 13;

    auto values = tatami_test::simulate_vector<double>(NR * stride_nc, [&]{
        tatami_test::SimulateVectorOptions opt;
        opt.lower = 0;
        opt.upper = 100;
        opt.seed = NR * 20 + NC;
        return opt;
    }());

    std::vector<double> buffer(stride_nr * NC);
    tatami::transpose(values.data(), NR, NC, stride_nc, buffer.data(), stride_nr);

    tatami::DenseRowMatrix<double, int> original(NR, stride_nc, std::move(values));
    tatami::DenseColumnMatrix<double, int> flipped(stride_nr, NC, std::move(buffer));

    auto oext = original.dense_row();
    auto fext = flipped.dense_row();
    for (int r = 0; r < NR; ++r) {
        auto oout = tatami_test::fetch(*oext, r, NC);
        auto fout = tatami_test::fetch(*fext, r, NC);
        ASSERT_EQ(oout, fout);
    }
}

INSTANTIATE_TEST_SUITE_P(
    DenseMatrix,
    TransposeDenseTest,
    ::testing::Combine(
        ::testing::Values(1, 10, 20, 40, 80), // number of rows
        ::testing::Values(1, 10, 20, 40, 80)  // number of columns
    )
);
