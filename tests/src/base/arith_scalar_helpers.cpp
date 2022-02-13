#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedIsometricOp.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../data/data.h"
#include "TestCore.h"

template<class PARAM> 
class ArithScalarTest : public TestCore<::testing::TestWithParam<PARAM> > {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense, sparse;
protected:
    void SetUp() {
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
        sparse = tatami::convert_to_sparse<false>(dense.get()); // column major.
        return;
    }
};

/****************************
 ********* ADDITION *********
 ****************************/

class ArithScalarAdditionTest : public ArithScalarTest<double> {};

TEST_P(ArithScalarAdditionTest, ColumnAccess) {
    double val = GetParam();
    tatami::DelayedAddScalarHelper<double> op(val);
    auto dense_mod = tatami::make_DelayedIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    // Most tests covered by the vector version, so we'll just go light here.
    for (size_t i = 0; i < dense->ncol(); ++i) {
        auto expected = extract_dense<true>(sparse.get(), i);
        for (auto& e : expected) { e += val; }

        auto output = extract_dense<true>(sparse_mod.get(), i);
        EXPECT_EQ(output, expected);

        auto output2 = extract_dense<true>(dense_mod.get(), i);
        EXPECT_EQ(output2, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithScalar,
    ArithScalarAdditionTest,
    ::testing::Values(5, 0.1, -0.7)
);

/*******************************
 ********* SUBTRACTION *********
 *******************************/

class ArithScalarSubtractionTest : public ArithScalarTest<std::tuple<double, bool> > {};

TEST_P(ArithScalarSubtractionTest, ColumnAccess) {
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        tatami::DelayedSubtractScalarHelper<true> op(val);
        dense_mod = tatami::make_DelayedIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);
    } else {
        tatami::DelayedSubtractScalarHelper<false> op(val);
        dense_mod = tatami::make_DelayedIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);
    }

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        auto expected = extract_dense<true>(sparse.get(), i);
        for (auto& e : expected) {
            if (on_right) {
                e -= val;
            } else {
                e = val - e;
            }
        }

        auto output = extract_dense<true>(sparse_mod.get(), i);
        EXPECT_EQ(output, expected);

        auto output2 = extract_dense<true>(dense_mod.get(), i);
        EXPECT_EQ(output2, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithScalar,
    ArithScalarSubtractionTest,
    ::testing::Combine(
        ::testing::Values(5, 0.1, -0.7),
        ::testing::Values(true, false)
    )
);

/**********************************
 ********* MULTIPLICATION *********
 **********************************/

class ArithScalarMultiplicationTest : public ArithScalarTest<double> {};

TEST_P(ArithScalarMultiplicationTest, ColumnAccess) {
    double val = GetParam();
    tatami::DelayedMultiplyScalarHelper<double> op(val);
    auto dense_mod = tatami::make_DelayedIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);

    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        auto expected = extract_dense<true>(sparse.get(), i);
        for (auto& e : expected) { e *= val; }

        auto output = extract_dense<true>(sparse_mod.get(), i);
        EXPECT_EQ(output, expected);

        auto output2 = extract_dense<true>(dense_mod.get(), i);
        EXPECT_EQ(output2, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithScalar,
    ArithScalarMultiplicationTest,
    ::testing::Values(5, 0.1, -0.7)
);

/****************************
 ********* DIVISION *********
 ****************************/

class ArithScalarDivisionTest : public ArithScalarTest<std::tuple<double, bool> > {};

TEST_P(ArithScalarDivisionTest, ColumnAccess) {
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        tatami::DelayedDivideScalarHelper<true> op(val);
        dense_mod = tatami::make_DelayedIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);
    } else {
        tatami::DelayedDivideScalarHelper<false> op(val);
        dense_mod = tatami::make_DelayedIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);
    }

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        auto expected = extract_dense<true>(sparse.get(), i);
        for (auto& e : expected) {
            if (on_right) {
                e /= val;
            } else {
                e = val/e;
            }
        }

        auto output = extract_dense<true>(sparse_mod.get(), i);
        EXPECT_EQ(output, expected);

        auto output2 = extract_dense<true>(dense_mod.get(), i);
        EXPECT_EQ(output2, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithScalar,
    ArithScalarDivisionTest,
    ::testing::Combine(
        ::testing::Values(5, 0.1, -0.7),
        ::testing::Values(true, false)
    )
);
