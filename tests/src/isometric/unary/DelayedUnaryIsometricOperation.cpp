#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOperation.hpp"
#include "tatami/isometric/unary/arithmetic_helpers.hpp"
#include "tatami/isometric/unary/helper_interface.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

TEST(DelayedUnaryIsometricOperation, ConstOverload) {
    int nrow = 23, ncol = 42;
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::vector<double>(nrow * ncol)));

    auto vec = std::vector<double>(nrow);
    auto op = tatami::make_DelayedUnaryIsometricAddVector(vec, true);
    auto mat = tatami::make_DelayedUnaryIsometricOperation(dense, std::move(op));

    // cursory checks.
    EXPECT_EQ(mat->nrow(), dense->nrow());
    EXPECT_EQ(mat->ncol(), dense->ncol());
}

class UnaryMockMissing {
public:
    static constexpr bool is_basic = false;

    bool is_sparse() const { return false; }
};

TEST(DelayedUnaryIsometricOperation, DependsChecks) {
    auto optr = std::make_shared<const tatami::ConsecutiveOracle<int> >(10, 20);

    {
        // Check that the oracle is correctly ignored.
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, tatami::DelayedUnaryIsometricMockAdvanced, int> oracle1(optr, {}, true);
        EXPECT_EQ(oracle1.get(0), 0);
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, tatami::DelayedUnaryIsometricMockAdvanced, int> oracle2(optr, {}, false);
        EXPECT_EQ(oracle2.get(0), 0);

        EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(tatami::DelayedUnaryIsometricMockAdvanced(), true));
        EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(tatami::DelayedUnaryIsometricMockAdvanced(), false));
        EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::needs_sparse_indices(tatami::DelayedUnaryIsometricMockAdvanced(), true));
        EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::needs_sparse_indices(tatami::DelayedUnaryIsometricMockAdvanced(), false));
    }

    {
        // Check that the oracle is correctly ignored.
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockMissing, int> oracle1(optr, {}, true);
        EXPECT_EQ(oracle1.get(0), 0);
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockMissing, int> oracle2(optr, {}, false);
        EXPECT_EQ(oracle2.get(0), 0);

        EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(UnaryMockMissing(), true));
        EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(UnaryMockMissing(), false));
        EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::needs_sparse_indices(UnaryMockMissing(), true));
        EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::needs_sparse_indices(UnaryMockMissing(), false));
    }

    { 
        // Check that the oracle is correctly used for rows.
        {
            auto mock = tatami::make_DelayedUnaryIsometricAddVector(std::vector<double>{1.0, 2.0}, true);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, decltype(mock), int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 10);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, decltype(mock), int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            auto mock = tatami::make_DelayedUnaryIsometricAddVector(std::vector<double>{0.0, 0.0}, true);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, decltype(mock), int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 10);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, decltype(mock), int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
            EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::needs_sparse_indices(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::needs_sparse_indices(mock, false));
        }
    }

    { 
        // Check that the oracle is correctly used for columns.
        {
            auto mock = tatami::make_DelayedUnaryIsometricAddVector(std::vector<double>{1.0, 2.0}, false);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, decltype(mock), int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, decltype(mock), int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 10);

            EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            auto mock = tatami::make_DelayedUnaryIsometricAddVector(std::vector<double>{0.0, 0.0}, false);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, decltype(mock), int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, decltype(mock), int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 10);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::needs_sparse_indices(mock, true));
            EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::needs_sparse_indices(mock, false));
        }
    }

    {
        // Check that the oracle is respected in the basic case.
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, tatami::DelayedUnaryIsometricMockBasic, int> oracle1(optr, {}, true);
        EXPECT_EQ(oracle1.get(0), 10);
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, tatami::DelayedUnaryIsometricMockBasic, int> oracle2(optr, {}, false);
        EXPECT_EQ(oracle2.get(0), 10);
    }
}

class DelayedUnaryIsometricOperationTest : public ::testing::TestWithParam<std::tuple<bool, bool> > {
protected:
    inline static int nrow = 57, ncol = 37;
    inline static std::vector<double> simulated;
    inline static std::shared_ptr<const tatami::NumericMatrix> dense, sparse, udense, usparse, ref;
    inline static std::shared_ptr<const tatami::Matrix<uint8_t, int> > i_udense, i_usparse, i_ref;

    static void SetUpTestSuite() {
        simulated = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.1;
            opt.seed = 777222777;
            return opt;
        }());

        dense.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(simulated)));
        sparse = tatami::convert_to_compressed_sparse<false, double, int>(dense.get()); 

        udense = tatami::make_DelayedUnaryIsometricOperation(dense, tatami::DelayedUnaryIsometricMockBasic());
        usparse = tatami::make_DelayedUnaryIsometricOperation(sparse, tatami::DelayedUnaryIsometricMockAdvanced());
        ref.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::vector<double>(nrow * ncol)));

        i_udense = tatami::make_DelayedUnaryIsometricOperation<uint8_t>(dense, tatami::DelayedUnaryIsometricMockBasic());
        i_usparse = tatami::make_DelayedUnaryIsometricOperation<uint8_t>(sparse, tatami::DelayedUnaryIsometricMockAdvanced());
        i_ref.reset(new tatami::DenseRowMatrix<uint8_t, int>(nrow, ncol, std::vector<uint8_t>(nrow * ncol)));
    }
};

TEST_P(DelayedUnaryIsometricOperationTest, Mock) {
    EXPECT_FALSE(udense->is_sparse());
    EXPECT_EQ(udense->is_sparse_proportion(), 0);

    EXPECT_TRUE(usparse->is_sparse());
    EXPECT_EQ(usparse->is_sparse_proportion(), 1);

    // Spamming a whole stack of tests.
    tatami_test::TestAccessOptions opts;
    auto tparam = GetParam();
    opts.use_row = std::get<0>(tparam);
    opts.use_oracle = std::get<1>(tparam);

    tatami_test::test_full_access(*udense, *ref, opts);
    tatami_test::test_block_access(*udense, *ref, 0.25, 0.34, opts);
    tatami_test::test_indexed_access(*udense, *ref, 0.3, 0.5, opts);

    tatami_test::test_full_access(*usparse, *ref, opts);
    tatami_test::test_block_access(*usparse, *ref, 0.5, 0.35, opts);
    tatami_test::test_indexed_access(*usparse, *ref, 0.2, 0.4, opts);
}

TEST_P(DelayedUnaryIsometricOperationTest, NewType) {
    auto tparam = GetParam();
    tatami_test::TestAccessOptions opts;
    opts.use_row = std::get<0>(tparam);
    opts.use_oracle = std::get<1>(tparam);

    tatami_test::test_full_access(*i_udense, *i_ref, opts);
    tatami_test::test_block_access(*i_udense, *i_ref, 0.3, 0.55, opts);
    tatami_test::test_indexed_access(*i_udense, *i_ref, 0.3, 0.55, opts);

    tatami_test::test_full_access(*i_usparse, *i_ref, opts);
    tatami_test::test_block_access(*i_usparse, *i_ref, 0.25, 0.3, opts);
    tatami_test::test_indexed_access(*i_usparse, *i_ref, 0.2, 0.3, opts);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricOperation,
    DelayedUnaryIsometricOperationTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row access
        ::testing::Values(true, false)  // oracle usage
    )
);
