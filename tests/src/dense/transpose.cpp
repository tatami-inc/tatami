#include <gtest/gtest.h>

#include <vector>

#include "tatami_test/tatami_test.hpp"
#include "tatami/dense/transpose.hpp"

class TransposeDenseTest : public ::testing::TestWithParam<std::tuple<int, int> > {};

TEST_P(TransposeDenseTest, Basic) {
    auto params = GetParam();
    auto NR = std::get<0>(params);
    auto NC = std::get<1>(params);

    auto values = tatami_test::simulate_dense_vector<double>(NR * NC, 0, 100, /* seed = */ NR * 10 + NC);
    std::vector<double> buffer(NR * NC);
    tatami::transpose(values.data(), buffer.data(), NR, NC);

    tatami::DenseRowMatrix<double, int> original(NR, NC, std::move(values));
    tatami::DenseColumnMatrix<double, int> flipped(NR, NC, std::move(buffer));

    auto oext = original.dense_row();
    auto fext = flipped.dense_row();
    for (int r = 0; r < NR; ++r) {
        auto oout = tatami_test::fetch(oext.get(), r, NC);
        auto fout = tatami_test::fetch(fext.get(), r, NC);
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
