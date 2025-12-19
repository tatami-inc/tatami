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

    const auto values = tatami_test::simulate_vector<double>(NR, NC, [&]{
        tatami_test::SimulateVectorOptions opt;
        opt.lower = 0;
        opt.upper = 100;
        opt.seed = NR * 10 + NC;
        return opt;
    }());

    std::vector<double> buffer(values.size());
    tatami::transpose(values.data(), static_cast<std::size_t>(NR), static_cast<std::size_t>(NC), buffer.data());

    tatami::DenseRowMatrix<double, int> original(NR, NC, std::move(values));
    tatami::DenseColumnMatrix<double, int> flipped(NR, NC, std::move(buffer));

    auto oext = original.dense_row();
    auto fext = flipped.dense_row();
    sanisizer::as_size_type<std::vector<double> >(NC);
    std::vector<double> observed(NC), expected(NC);

    for (int r = 0; r < NR; ++r) {
        auto eptr = oext->fetch(r, expected.data());
        tatami::copy_n(eptr, NC, expected.data());
        auto optr = fext->fetch(r, observed.data());
        tatami::copy_n(optr, NC, observed.data());
        ASSERT_EQ(expected, observed);
    }
}

TEST_P(TransposeDenseTest, Strided) {
    auto params = GetParam();
    auto NR = std::get<0>(params);
    auto NC = std::get<1>(params);
    auto input_stride = NC + 13;
    auto output_stride = NR + 17;

    const auto values = tatami_test::simulate_vector<double>(NR, input_stride, [&]{
        tatami_test::SimulateVectorOptions opt;
        opt.lower = 0;
        opt.upper = 100;
        opt.seed = NR * 20 + NC;
        return opt;
    }());

    constexpr double placeholder = -123;
    const auto full_output_size = sanisizer::product<typename std::vector<double>::size_type>(output_stride, NC);
    std::vector<double> buffer(full_output_size, placeholder);
    tatami::transpose(
        values.data(),
        static_cast<std::size_t>(NR),
        static_cast<std::size_t>(NC),
        static_cast<std::size_t>(input_stride),
        buffer.data(),
        static_cast<std::size_t>(output_stride)
    );

    tatami::DenseRowMatrix<double, int> original(NR, input_stride, std::move(values));
    tatami::DenseColumnMatrix<double, int> flipped(output_stride, NC, std::move(buffer));

    auto oext = original.dense_row(0, NC);
    auto fext = flipped.dense_row();
    sanisizer::as_size_type<std::vector<double> >(input_stride);
    std::vector<double> observed(input_stride), expected(input_stride);

    for (int r = 0; r < NR; ++r) {
        auto optr = oext->fetch(r, expected.data());
        tatami::copy_n(optr, NC, expected.data());
        auto fptr = fext->fetch(r, observed.data());
        tatami::copy_n(fptr, NC, observed.data());
        ASSERT_EQ(expected, observed);
    }

    for (int r = NR; r < output_stride; ++r) {
        auto fptr = fext->fetch(r, observed.data());
        for (int c = 0; c < NC; ++c) {
            EXPECT_EQ(fptr[c], placeholder);
        }
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
