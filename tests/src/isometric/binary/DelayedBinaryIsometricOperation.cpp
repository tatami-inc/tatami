#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/binary/DelayedBinaryIsometricOperation.hpp"
#include "tatami/isometric/binary/mock_helpers.hpp"
#include "tatami/isometric/binary/arithmetic_helpers.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

TEST(DelayedBinaryIsometricOperation, ConstOverload) {
    int nrow = 23, ncol = 42;
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::vector<double>(nrow * ncol)));

    auto op = tatami::make_DelayedBinaryIsometricAdd();
    auto mat = tatami::make_DelayedBinaryIsometricOperation(dense, dense, std::move(op));

    // cursory checks.
    EXPECT_EQ(mat->nrow(), nrow);
    EXPECT_EQ(mat->ncol(), ncol);
}

TEST(DelayedBinaryIsometricOperation, Misshappen) {
    std::vector<double> src(200);
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(10, 20, src));
    auto dense2 = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(20, 10, src));
    tatami_test::throws_error([&]() {
        tatami::make_DelayedBinaryIsometricOperation(dense, dense2, tatami::DelayedBinaryIsometricMockBasic());
    }, "should be the same");
}

class BinaryMockMissing {
public:
    static constexpr bool is_basic = false;

    bool is_sparse() const { return false; }
};

class BinaryMockDerived : public tatami::DelayedBinaryIsometricMockAdvanced {
public:
    BinaryMockDerived(bool sparse = false, bool zero_row = true, bool zero_col = true, bool non_zero_row = true, bool non_zero_col = true) : 
        my_sparse(sparse), my_zero_row(zero_row), my_zero_col(zero_col), my_non_zero_row(non_zero_row), my_non_zero_col(non_zero_col) {}

private:
    bool my_sparse, my_zero_row, my_zero_col, my_non_zero_row, my_non_zero_col;

public:
    static constexpr bool is_basic = false;
    bool is_sparse() const { return my_sparse; }
    bool zero_depends_on_row() const { return my_zero_row; }
    bool zero_depends_on_column() const { return my_zero_col; }
    bool non_zero_depends_on_row() const { return my_non_zero_row; }
    bool non_zero_depends_on_column() const { return my_non_zero_col; }
};

TEST(DelayedBinaryIsometricOperation, DependsChecks) {
    auto optr = std::make_shared<const tatami::ConsecutiveOracle<int> >(10, 20);

    {
        // Check that the oracle is correctly ignored.
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, tatami::DelayedBinaryIsometricMockAdvanced, int> oracle1(optr, {}, true);
        EXPECT_EQ(oracle1.get(0), 0);
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, tatami::DelayedBinaryIsometricMockAdvanced, int> oracle2(optr, {}, false);
        EXPECT_EQ(oracle2.get(0), 0);

        EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(tatami::DelayedBinaryIsometricMockAdvanced(), true));
        EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(tatami::DelayedBinaryIsometricMockAdvanced(), false));
    }

    {
        // Check that the oracle is correctly ignored.
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockMissing, int> oracle1(optr, {}, true);
        EXPECT_EQ(oracle1.get(0), 0);
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockMissing, int> oracle2(optr, {}, false);
        EXPECT_EQ(oracle2.get(0), 0);

        EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(BinaryMockMissing(), true));
        EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(BinaryMockMissing(), false));
    }

    { 
        // Check that the oracle is correctly used.
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived, int> oracle1(optr, {}, true);
        EXPECT_EQ(oracle1.get(0), 10);
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived, int> oracle2(optr, {}, false);
        EXPECT_EQ(oracle2.get(0), 10);

        EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(BinaryMockDerived(), true));
        EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(BinaryMockDerived(), false));
    }

    { 
        // Check that the oracle is correctly used for rows.
        {
            BinaryMockDerived mock(false, true, false, false, false); // zero-dependency for non-sparse operation.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 10);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            BinaryMockDerived mock(false, false, false, true, false); // non-zero dependency.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 10);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            BinaryMockDerived mock(true, true, false, false, false); // zero dependency ignored for sparse operations.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }
    }

    { 
        // Check that the oracle is correctly used for columns.
        {
            BinaryMockDerived mock(false, false, true, false, false); // zero-dependency for non-sparse operation.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 10);

            EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            BinaryMockDerived mock(false, false, false, false, true); // non-zero dependency.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 10);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            BinaryMockDerived mock(true, false, true, false, false); // zero dependency ignored for sparse operations.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }
    }

    {
        // Check that the oracle is respected in the basic case.
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, tatami::DelayedBinaryIsometricMockBasic, int> oracle1(optr, {}, true);
        EXPECT_EQ(oracle1.get(0), 10);
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, tatami::DelayedBinaryIsometricMockBasic, int> oracle2(optr, {}, false);
        EXPECT_EQ(oracle2.get(0), 10);
    }
}

class DelayedBinaryIsometricOperationTest : public ::testing::TestWithParam<std::tuple<bool, bool> > {
protected:
    inline static int nrow = 23, ncol = 42;
    inline static std::vector<double> simulated;
    inline static std::shared_ptr<const tatami::NumericMatrix> dense, sparse, bdense, bsparse, ref;
    inline static std::shared_ptr<tatami::Matrix<float, int> > fbdense, fbsparse, fref;

    static void SetUpTestSuite() {
        simulated = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.2;
            opt.seed = 918273;
            return opt;
        }());

        dense.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(simulated)));
        sparse = tatami::convert_to_compressed_sparse<false, double, int>(dense.get()); 

        bdense = tatami::make_DelayedBinaryIsometricOperation(dense, dense, tatami::DelayedBinaryIsometricMockBasic());
        bsparse = tatami::make_DelayedBinaryIsometricOperation(sparse, sparse, tatami::DelayedBinaryIsometricMockAdvanced());
        ref.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::vector<double>(nrow * ncol)));

        fbdense = tatami::make_DelayedBinaryIsometricOperation<float>(dense, dense, tatami::DelayedBinaryIsometricMockBasic());
        fbsparse = tatami::make_DelayedBinaryIsometricOperation<float>(sparse, sparse, tatami::DelayedBinaryIsometricMockAdvanced());
        fref.reset(new tatami::DenseRowMatrix<float, int>(nrow, ncol, std::vector<float>(nrow * ncol)));
    }
};

TEST_P(DelayedBinaryIsometricOperationTest, Mock) {
    EXPECT_FALSE(bdense->is_sparse());
    EXPECT_EQ(bdense->is_sparse_proportion(), 0);
    EXPECT_TRUE(bsparse->is_sparse());
    EXPECT_EQ(bsparse->is_sparse_proportion(), 1);

    // Spamming a whole stack of tests.
    tatami_test::TestAccessOptions opts;
    auto tparam = GetParam();
    opts.use_row = std::get<0>(tparam);
    opts.use_oracle = std::get<1>(tparam);

    tatami_test::test_full_access(*bdense, *ref, opts);
    tatami_test::test_block_access(*bdense, *ref, 0.1, 0.7, opts);
    tatami_test::test_indexed_access(*bdense, *ref, 0.23, 0.5, opts);

    tatami_test::test_full_access(*bsparse, *ref, opts);
    tatami_test::test_block_access(*bsparse, *ref, 0.13, 0.5, opts);
    tatami_test::test_indexed_access(*bsparse, *ref, 0.23, 0.4, opts);
}

TEST_P(DelayedBinaryIsometricOperationTest, NewType) {
    tatami_test::TestAccessOptions opts;
    auto tparam = GetParam();
    opts.use_row = std::get<0>(tparam);
    opts.use_oracle = std::get<1>(tparam);

    tatami_test::test_full_access(*fbdense, *fref, opts);
    tatami_test::test_block_access(*fbdense, *fref, 0.5, 0.5, opts);
    tatami_test::test_indexed_access(*fbdense, *fref, 0.3, 0.5, opts);

    tatami_test::test_full_access(*fbsparse, *fref, opts);
    tatami_test::test_block_access(*fbsparse, *fref, 0.2, 0.6, opts);
    tatami_test::test_indexed_access(*fbsparse, *fref, 0.2, 0.4, opts);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBinaryIsometricOperation,
    DelayedBinaryIsometricOperationTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row access
        ::testing::Values(true, false)  // oracle usage
    )
);
