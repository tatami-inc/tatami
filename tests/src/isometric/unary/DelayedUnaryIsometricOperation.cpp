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

TEST(DelayedUnaryIsometricOperation, BackCompatibility) {
    int nrow = 23, ncol = 42;
    auto dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::vector<double>(nrow * ncol)));

    auto vec = std::vector<double>(nrow);
    auto mat = tatami::make_DelayedUnaryIsometricOperation(dense, std::make_shared<tatami::DelayedUnaryIsometricAddScalarHelper<double, double, int, double> >(2.0));
    EXPECT_EQ(mat->nrow(), nrow);
    EXPECT_EQ(mat->ncol(), ncol);

    std::shared_ptr<const tatami::NumericMatrix> cdense(dense);
    auto cmat = tatami::make_DelayedUnaryIsometricOperation(cdense, std::make_shared<const tatami::DelayedUnaryIsometricAddScalarHelper<double, double, int, double> >(2.0));
    EXPECT_EQ(cmat->nrow(), nrow);
    EXPECT_EQ(cmat->ncol(), ncol);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
class UnaryMockDerived : public tatami::DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    UnaryMockDerived(bool sparse = false, bool zero_row = true, bool zero_col = true, bool non_zero_row = true, bool non_zero_col = true) : 
        my_sparse(sparse), my_zero_row(zero_row), my_zero_col(zero_col), my_non_zero_row(non_zero_row), my_non_zero_col(non_zero_col) {}

private:
    bool my_sparse, my_zero_row, my_zero_col, my_non_zero_row, my_non_zero_col;

public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool is_sparse() const { return my_sparse; }
    bool zero_depends_on_row() const { return my_zero_row; }
    bool zero_depends_on_column() const { return my_zero_col; }
    bool non_zero_depends_on_row() const { return my_non_zero_row; }
    bool non_zero_depends_on_column() const { return my_non_zero_col; }

public:
    OutputValue_ fill(bool, int) const {
        return 0;
    }

    void dense(bool, Index_, Index_, Index_ length, const InputValue_*, OutputValue_* output_buffer) const {
        std::fill_n(output_buffer, length, 0);
    }

    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_*, OutputValue_* output_buffer) const {
        std::fill_n(output_buffer, indices.size(), 0);
    }

    void sparse(bool, Index_, Index_ number, const InputValue_*, const Index_*, OutputValue_* output_buffer) const {
        std::fill_n(output_buffer, number, 0);
    }
};

TEST(DelayedUnaryIsometricOperation, DependsChecks) {
    auto optr = std::make_shared<const tatami::ConsecutiveOracle<int> >(10, 20);

    { 
        // Check that the oracle is correctly used for both rows and columns.
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockDerived<>, int> oracle1(optr, {}, true);
        EXPECT_EQ(oracle1.get(0), 10);
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockDerived<>, int> oracle2(optr, {}, false);
        EXPECT_EQ(oracle2.get(0), 10);

        EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(UnaryMockDerived<>(), true));
        EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(UnaryMockDerived<>(), false));
    }

    { 
        // Check that the oracle is correctly used for rows only.
        {
            UnaryMockDerived<> mock(false, true, false, false, false); // zero-dependency for non-sparse operation.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockDerived<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 10);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockDerived<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            UnaryMockDerived<> mock(false, false, false, true, false); // non-zero dependency.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockDerived<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 10);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockDerived<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            UnaryMockDerived<> mock(true, true, false, false, false); // zero dependency ignored for sparse operations.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockDerived<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockDerived<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }
    }

    { 
        // Check that the oracle is correctly used for columns only.
        {
            UnaryMockDerived<> mock(false, false, true, false, false); // zero-dependency for non-sparse operation.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockDerived<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockDerived<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 10);

            EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            UnaryMockDerived<> mock(false, false, false, false, true); // non-zero dependency.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockDerived<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockDerived<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 10);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            UnaryMockDerived<> mock(true, false, true, false, false); // zero dependency ignored for sparse operations.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockDerived<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMockDerived<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }
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

        dense.reset(new tatami::DenseMatrix<double, int, decltype(simulated)>(nrow, ncol, std::move(simulated), true)); // row major
        sparse = tatami::convert_to_compressed_sparse<double, int>(*dense, false, {});  // column major

        udense.reset(new tatami::DelayedUnaryIsometricOperation<double, double, int>(dense, std::make_shared<UnaryMockDerived<> >()));
        usparse.reset(new tatami::DelayedUnaryIsometricOperation<double, double, int>(sparse, std::make_shared<UnaryMockDerived<> >(true)));
        ref.reset(new tatami::DenseMatrix<double, int, std::vector<double> >(nrow, ncol, std::vector<double>(nrow * ncol), true));

        i_udense.reset(new tatami::DelayedUnaryIsometricOperation<uint8_t, double, int>(dense, std::make_shared<UnaryMockDerived<uint8_t> >()));
        i_usparse.reset(new tatami::DelayedUnaryIsometricOperation<uint8_t, double, int>(sparse, std::make_shared<UnaryMockDerived<uint8_t> >(true)));
        i_ref.reset(new tatami::DenseMatrix<uint8_t, int, std::vector<uint8_t> >(nrow, ncol, std::vector<uint8_t>(nrow * ncol), true));
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

TEST(DelayedUnaryIsometricOperation, DimMismatch) {
    auto dense = std::make_shared<tatami::DenseMatrix<double, int, std::vector<double> > >(10, 20, std::vector<double>(200), true);
    tatami_test::throws_error([&]{
        std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(
            dense,
            std::make_shared<tatami::DelayedUnaryIsometricAddVectorHelper<double, double, int, std::vector<double> > >(std::vector<double>(10), false)
        );
    }, "number of columns");
    tatami_test::throws_error([&]{
        std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(
            dense,
            std::make_shared<tatami::DelayedUnaryIsometricAddVectorHelper<double, double, int, std::vector<double> > >(std::vector<double>(20), true)
        );
    }, "number of rows");
}
