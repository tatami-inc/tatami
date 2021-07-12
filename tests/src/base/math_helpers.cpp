#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedIsometricOp.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../data/data.h"
#include "TestCore.h"

class MathTest : public TestCore<::testing::Test> {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense, sparse;
protected:
    void SetUp() {
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
        sparse = tatami::convert_to_sparse(dense.get(), false); // column major.
        return;
    }
};

TEST_F(MathTest, Abs) {
    tatami::DelayedAbsHelper op;
    auto dense_mod = tatami::make_DelayedIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    // Most DelayedIsometricOp checks are checked in arith_vector_helpers,
    // so we don't bother to do anything more complex here.
    for (size_t i = 0; i < dense->ncol(); ++i) {
        auto expected = extract_dense<false>(dense.get(), i);
        for (auto& e : expected) { e = std::abs(e); }

        auto outputD = extract_dense<false>(dense_mod.get(), i);
        EXPECT_EQ(outputD, expected);

        auto outputS = extract_sparse<false>(sparse_mod.get(), i);
        EXPECT_EQ(outputS, expected);
    }
}

TEST_F(MathTest, SqrtByColumn) {
    tatami::DelayedAbsHelper op0;
    auto dense_mod0 = tatami::make_DelayedIsometricOp(dense, op0);
    auto sparse_mod0 = tatami::make_DelayedIsometricOp(sparse, op0);

    tatami::DelayedSqrtHelper op;
    auto dense_mod = tatami::make_DelayedIsometricOp(dense_mod0, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse_mod0, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        auto expected = extract_dense<false>(dense.get(), i);
        for (auto& e : expected) { e = std::sqrt(std::abs(e)); }

        auto outputD = extract_dense<false>(dense_mod.get(), i);
        EXPECT_EQ(outputD, expected);

        auto outputS = extract_sparse<false>(sparse_mod.get(), i);
        EXPECT_EQ(outputS, expected);
    }
}

TEST_F(MathTest, LogByColumn) {
    tatami::DelayedAbsHelper op0;
    auto dense_mod0 = tatami::make_DelayedIsometricOp(dense, op0);
    auto sparse_mod0 = tatami::make_DelayedIsometricOp(sparse, op0);
    
    double CONSTANT = 5;
    tatami::DelayedAddScalarHelper<double> op1(CONSTANT);
    auto dense_mod1 = tatami::make_DelayedIsometricOp(dense_mod0, op1);
    auto sparse_mod1 = tatami::make_DelayedIsometricOp(sparse_mod0, op1);

    tatami::DelayedLogHelper op;
    auto dense_mod = tatami::make_DelayedIsometricOp(dense_mod1, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse_mod1, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        auto expected = extract_dense<false>(dense.get(), i);
        for (auto& e : expected) { e = std::log(std::abs(e) + CONSTANT); }

        auto outputD = extract_dense<false>(dense_mod.get(), i);
        EXPECT_EQ(outputD, expected);

        auto outputS = extract_sparse<false>(sparse_mod.get(), i);
        EXPECT_EQ(outputS, expected);
    }

    // Trying with another base.
    tatami::DelayedLogHelper op2(2);
    auto dense_mod2 = tatami::make_DelayedIsometricOp(dense_mod1, op2);
    auto sparse_mod2 = tatami::make_DelayedIsometricOp(sparse_mod1, op2);

    for (size_t i = 0; i < dense->ncol(); ++i) {
        auto expected = extract_dense<false>(dense.get(), i);
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::log(std::abs(expected[j]) + CONSTANT)/std::log(2);
        }

        auto outputD = extract_dense<false>(dense_mod2.get(), i);
        EXPECT_EQ(outputD, expected);

        auto outputS = extract_sparse<false>(sparse_mod2.get(), i);
        EXPECT_EQ(outputS, expected);
    }
}

TEST_F(MathTest, Log1pByColumn) {
    tatami::DelayedAbsHelper op0;
    auto dense_mod0 = tatami::make_DelayedIsometricOp(dense, op0);
    auto sparse_mod0 = tatami::make_DelayedIsometricOp(sparse, op0);

    tatami::DelayedLog1pHelper op;
    auto dense_mod = tatami::make_DelayedIsometricOp(dense_mod0, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse_mod0, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        auto expected = extract_dense<false>(dense.get(), i);
        for (auto& e : expected) { e = std::log1p(std::abs(e)); }

        auto outputD = extract_dense<false>(dense_mod.get(), i);
        EXPECT_EQ(outputD, expected);

        auto outputS = extract_sparse<false>(sparse_mod.get(), i);
        EXPECT_EQ(outputS, expected);
    }

    // Trying with another base.
    tatami::DelayedLog1pHelper op2(2);
    auto dense_mod2 = tatami::make_DelayedIsometricOp(dense_mod0, op2);
    auto sparse_mod2 = tatami::make_DelayedIsometricOp(sparse_mod0, op2);

    for (size_t i = 0; i < dense->ncol(); ++i) {
        auto expected = extract_dense<false>(dense.get(), i);
        for (auto& e : expected) { e = std::log1p(std::abs(e))/std::log(2); }

        auto outputD = extract_dense<false>(dense_mod2.get(), i);
        EXPECT_EQ(outputD, expected);

        auto outputS = extract_sparse<false>(sparse_mod2.get(), i);
        EXPECT_EQ(outputS, expected);
    }
}

TEST_F(MathTest, ExpByColumn) {
    tatami::DelayedExpHelper op;
    auto dense_mod = tatami::make_DelayedIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        auto expected = extract_dense<false>(dense.get(), i);
        for (auto& e : expected) { e = std::exp(e); }

        auto outputD = extract_dense<false>(dense_mod.get(), i);
        EXPECT_EQ(outputD, expected);

        auto outputS = extract_sparse<false>(sparse_mod.get(), i);
        EXPECT_EQ(outputS, expected);
    }
}

TEST_F(MathTest, RoundByColumn) {
    tatami::DelayedRoundHelper op;
    auto dense_mod = tatami::make_DelayedIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        auto expected = extract_dense<false>(dense.get(), i);
        for (auto& e : expected) { e = std::round(e); }

        auto outputD = extract_dense<false>(dense_mod.get(), i);
        EXPECT_EQ(outputD, expected);

        auto outputS = extract_sparse<false>(sparse_mod.get(), i);
        EXPECT_EQ(outputS, expected);
    }
}
