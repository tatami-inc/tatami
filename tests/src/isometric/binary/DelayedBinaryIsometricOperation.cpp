#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/binary/DelayedBinaryIsometricOperation.hpp"
#include "tatami/isometric/binary/helper_interface.hpp"
#include "tatami/isometric/binary/arithmetic_helpers.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

TEST(DelayedBinaryIsometricOperation, BackCompatibility) {
    int nrow = 23, ncol = 42;

    auto dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::vector<double>(nrow * ncol)));
    auto mat = tatami::make_DelayedBinaryIsometricOperation(dense, dense, std::make_shared<tatami::DelayedBinaryIsometricAddHelper<double, double, int> >());
    EXPECT_EQ(mat->nrow(), nrow);
    EXPECT_EQ(mat->ncol(), ncol);

    std::shared_ptr<const tatami::NumericMatrix> cdense(dense);
    auto cmat = tatami::make_DelayedBinaryIsometricOperation(cdense, cdense, std::make_shared<const tatami::DelayedBinaryIsometricAddHelper<double, double, int> >());
    EXPECT_EQ(cmat->nrow(), nrow);
    EXPECT_EQ(cmat->ncol(), ncol);
}

TEST(DelayedBinaryIsometricOperation, Misshappen) {
    std::vector<double> src(200);
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(10, 20, src));
    auto dense2 = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(20, 10, src));
    tatami_test::throws_error([&]() {
        tatami::DelayedBinaryIsometricOperation<double, double, int> mat(dense, dense2, std::make_shared<tatami::DelayedBinaryIsometricAddHelper<double, double, int> >());
    }, "should be the same");
}

TEST(DelayedBinaryIsometricOperation, MixedSparse) {
    int nrow = 123;
    int ncol = 45;
    auto simulated = tatami_test::simulate_vector<double>(nrow * ncol, []{
        tatami_test::SimulateVectorOptions opt;
        opt.density = 0.2;
        opt.seed = 918273;
        return opt;
    }());

    auto dense = std::make_shared<tatami::DenseMatrix<double, int, decltype(simulated)> >(nrow, ncol, std::move(simulated), true);
    auto sparse = tatami::convert_to_compressed_sparse<double, int>(*dense, false, {});

    // Sparsity depends on the underlying matrices if the operation is sparsity-preserving.
    {
        auto add_helper = std::make_shared<tatami::DelayedBinaryIsometricAddHelper<double, double, int> >();
        tatami::DelayedBinaryIsometricOperation<double, double, int> mat1(dense, dense, add_helper);
        EXPECT_FALSE(mat1.is_sparse());
        tatami::DelayedBinaryIsometricOperation<double, double, int> mat2(dense, sparse, add_helper);
        EXPECT_FALSE(mat2.is_sparse());
        tatami::DelayedBinaryIsometricOperation<double, double, int> mat3(sparse, sparse, add_helper);
        EXPECT_TRUE(mat3.is_sparse());
    }

    // Always dense if the operation is not sparsity-preserving.
    {
        auto pow_helper = std::make_shared<tatami::DelayedBinaryIsometricPowerHelper<double, double, int> >();
        tatami::DelayedBinaryIsometricOperation<double, double, int> mat1(dense, dense, pow_helper);
        EXPECT_FALSE(mat1.is_sparse());
        tatami::DelayedBinaryIsometricOperation<double, double, int> mat2(dense, sparse, pow_helper);
        EXPECT_FALSE(mat2.is_sparse());
        tatami::DelayedBinaryIsometricOperation<double, double, int> mat3(sparse, sparse, pow_helper);
        EXPECT_FALSE(mat3.is_sparse());
    }
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
class BinaryMockDerived : public tatami::DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    BinaryMockDerived(bool sparse = false, bool zero_row = true, bool zero_col = true, bool non_zero_row = true, bool non_zero_col = true) : 
        my_sparse(sparse), my_zero_row(zero_row), my_zero_col(zero_col), my_non_zero_row(non_zero_row), my_non_zero_col(non_zero_col) {}

private:
    bool my_sparse, my_zero_row, my_zero_col, my_non_zero_row, my_non_zero_col;

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

    void dense(bool, Index_, Index_, Index_ length, const InputValue_*, const InputValue_*, OutputValue_* output_buffer) const {
        std::fill_n(output_buffer, length, 0);
    }

    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_*, const InputValue_*, OutputValue_* output_buffer) const {
        std::fill_n(output_buffer, indices.size(), 0);
    }

    Index_ sparse(bool, Index_, const tatami::SparseRange<InputValue_, Index_>&, const tatami::SparseRange<InputValue_, Index_>&, OutputValue_*, Index_*, bool, bool) const {
        return 0;
    }
};

TEST(DelayedBinaryIsometricOperation, DependsChecks) {
    auto optr = std::make_shared<const tatami::ConsecutiveOracle<int> >(10, 20);

    { 
        // Check that the oracle is correctly used for both rows and columns.
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived<>, int> oracle1(optr, {}, true);
        EXPECT_EQ(oracle1.get(0), 10);
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived<>, int> oracle2(optr, {}, false);
        EXPECT_EQ(oracle2.get(0), 10);

        EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(BinaryMockDerived<>(), true));
        EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(BinaryMockDerived<>(), false));
    }

    { 
        // Check that the oracle is correctly used for rows only.
        {
            BinaryMockDerived<> mock(false, true, false, false, false); // zero-dependency for non-sparse operation.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 10);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            BinaryMockDerived<> mock(false, false, false, true, false); // non-zero dependency.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 10);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            BinaryMockDerived<> mock(true, true, false, false, false); // zero dependency ignored for sparse operations.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }
    }

    { 
        // Check that the oracle is correctly used for columns only.
        {
            BinaryMockDerived<> mock(false, false, true, false, false); // zero-dependency for non-sparse operation.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 10);

            EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            BinaryMockDerived<> mock(false, false, false, false, true); // non-zero dependency.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 10);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            BinaryMockDerived<> mock(true, false, true, false, false); // zero dependency ignored for sparse operations.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMockDerived<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }
    }
}

class DelayedBinaryIsometricOperationTest : public ::testing::TestWithParam<std::tuple<bool, bool> > {
protected:
    inline static int nrow = 23, ncol = 42;
    inline static std::vector<double> simulated;
    inline static std::shared_ptr<tatami::Matrix<double, int> > dense, sparse, bdense, bsparse, ref;
    inline static std::shared_ptr<tatami::Matrix<float, int> > fbdense, fbsparse, fref;

    static void SetUpTestSuite() {
        simulated = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.2;
            opt.seed = 918273;
            return opt;
        }());

        dense.reset(new tatami::DenseMatrix<double, int, decltype(simulated)>(nrow, ncol, simulated, true));
        sparse = tatami::convert_to_compressed_sparse<double, int>(*dense, false, {}); 

        bdense.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(dense, dense, std::make_shared<BinaryMockDerived<> >(false)));
        bsparse.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(sparse, sparse, std::make_shared<BinaryMockDerived<> >(true)));
        ref.reset(new tatami::DenseMatrix<double, int, std::vector<double> >(nrow, ncol, std::vector<double>(nrow * ncol), true));

        fbdense.reset(new tatami::DelayedBinaryIsometricOperation<float, double, int>(dense, dense, std::make_shared<BinaryMockDerived<float> >(false)));
        fbsparse.reset(new tatami::DelayedBinaryIsometricOperation<float, double, int>(sparse, sparse, std::make_shared<BinaryMockDerived<float> >(true)));
        fref.reset(new tatami::DenseMatrix<float, int, std::vector<float> >(nrow, ncol, std::vector<float>(nrow * ncol), true));
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
