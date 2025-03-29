#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOperation.hpp"
#include "tatami/isometric/unary/helper_interface.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

struct UnaryMockParams {
    UnaryMockParams(bool sparse = false, bool zero_row = true, bool zero_col = true, bool non_zero_row = true, bool non_zero_col = true) :
        sparse(sparse), zero_row(zero_row), zero_col(zero_col), non_zero_row(non_zero_row), non_zero_col(non_zero_col) {}
    bool sparse, zero_row, zero_col, non_zero_row, non_zero_col;
};

template<typename Index_, typename Value_>
static Value_ mock_operation(Index_ row, Index_ col, Value_ val, const UnaryMockParams& optype) {
    if (val == 0) {
        if (!optype.sparse) {
            if (optype.zero_row) {
                val += row;
            }
            if (optype.zero_col) {
                val += col * 13;
            }
        }
        return val;
    } else {
        if (optype.non_zero_row) {
            val += row;
        }
        if (optype.non_zero_col) {
            val += col * 13;
        }
        return val;
    }
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
class UnaryMock : public tatami::DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    UnaryMock() = default;
    UnaryMock(UnaryMockParams params) : UnaryMock(std::move(params), std::nullopt, std::nullopt) {}
    UnaryMock(UnaryMockParams params, std::optional<Index_> nrow, std::optional<Index_> ncol) : my_params(std::move(params)), my_nrow(nrow), my_ncol(ncol) {}

private:
    UnaryMockParams my_params;
    std::optional<Index_> my_nrow, my_ncol;

public:
    std::optional<Index_> nrow() const {
        return my_nrow;
    }

    std::optional<Index_> ncol() const {
        return my_ncol;
    }

public:
    bool is_sparse() const { return my_params.sparse; }
    bool zero_depends_on_row() const { return my_params.zero_row; }
    bool zero_depends_on_column() const { return my_params.zero_col; }
    bool non_zero_depends_on_row() const { return my_params.non_zero_row; }
    bool non_zero_depends_on_column() const { return my_params.non_zero_col; }

public:
    OutputValue_ fill(bool row, int i) const {
        return mock_operation<Index_, InputValue_>(
            row ? i : 0,
            row ? 0 : i,
            0,
            my_params
        );
    }

    void dense(bool row, Index_ i, Index_ start, Index_ length, const InputValue_* input, OutputValue_* output) const {
        for (Index_ s = 0; s < length; ++s) {
            output[s] = mock_operation(
                row ? i : s + start,
                row ? s + start : i,
                input[s],
                my_params
            );
        }
    }

    void dense(bool row, Index_ i, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        for (Index_ s = 0, end = indices.size(); s < end; ++s) {
            output[s] = mock_operation(
                row ? i : indices[s],
                row ? indices[s] : i,
                input[s],
                my_params
            );
        }
    }

    void sparse(bool row, Index_ i, Index_ number, const InputValue_* input_value, const Index_* indices, OutputValue_* output_value) const {
        for (Index_ s = 0; s < number; ++s) {
            Index_ curindex = ((!row && non_zero_depends_on_row()) || (row && non_zero_depends_on_column())) ? indices[s] : 0; // as non-target indices might not be extracted.
            output_value[s] = mock_operation(
                row ? i : curindex,
                row ? curindex : i,
                input_value[s],
                my_params
            );
        }
    }
};

TEST(DelayedUnaryIsometricOperation, BackCompatibility) {
    int nrow = 23, ncol = 42;
    auto dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::vector<double>(nrow * ncol)));

    auto vec = std::vector<double>(nrow);
    auto mat = tatami::make_DelayedUnaryIsometricOperation(dense, std::make_shared<UnaryMock<> >());
    EXPECT_EQ(mat->nrow(), nrow);
    EXPECT_EQ(mat->ncol(), ncol);

    std::shared_ptr<const tatami::NumericMatrix> cdense(dense);
    auto cmat = tatami::make_DelayedUnaryIsometricOperation(cdense, std::make_shared<const UnaryMock<> >());
    EXPECT_EQ(cmat->nrow(), nrow);
    EXPECT_EQ(cmat->ncol(), ncol);
}

TEST(DelayedUnaryIsometricOperation, HelperDimMismatch) {
    auto dense = std::make_shared<tatami::DenseMatrix<double, int, std::vector<double> > >(10, 20, std::vector<double>(200), true);

    // Matches don't cause errors.
    std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(
        dense,
        std::make_shared<UnaryMock<> >(UnaryMockParams(), std::optional<int>{ 10 }, std::optional<int>{ 20 })
    );

    tatami_test::throws_error([&]{
        std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(
            dense,
            std::make_shared<UnaryMock<> >(UnaryMockParams(), std::nullopt, std::optional<int>{ 10 })
        );
    }, "number of columns");

    tatami_test::throws_error([&]{
        std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(
            dense,
            std::make_shared<UnaryMock<> >(UnaryMockParams(), std::optional<int>{ 5 }, std::nullopt)
        );
    }, "number of rows");
}

TEST(DelayedUnaryIsometricOperation, DependsChecks) {
    auto optr = std::make_shared<const tatami::ConsecutiveOracle<int> >(10, 20);

    { 
        // Check that the oracle is correctly used for both rows and columns.
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMock<>, int> oracle1(optr, {}, true);
        EXPECT_EQ(oracle1.get(0), 10);
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMock<>, int> oracle2(optr, {}, false);
        EXPECT_EQ(oracle2.get(0), 10);

        EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(UnaryMock<>(), true));
        EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(UnaryMock<>(), false));
    }

    { 
        // Check that the oracle is correctly used for rows only.
        {
            UnaryMock<> mock({false, true, false, false, false}); // zero-dependency for non-sparse operation.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMock<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 10);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMock<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            UnaryMock<> mock({false, false, false, true, false}); // non-zero dependency.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMock<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 10);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMock<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            UnaryMock<> mock({true, true, false, false, false}); // zero dependency ignored for sparse operations.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMock<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMock<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }
    }

    { 
        // Check that the oracle is correctly used for columns only.
        {
            UnaryMock<> mock({false, false, true, false, false}); // zero-dependency for non-sparse operation.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMock<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMock<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 10);

            EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            UnaryMock<> mock({false, false, false, false, true}); // non-zero dependency.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMock<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMock<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 10);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            UnaryMock<> mock({true, false, true, false, false}); // zero dependency ignored for sparse operations.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMock<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, UnaryMock<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }
    }
}

class DelayedUnaryIsometricOperationTest : public ::testing::TestWithParam<std::tuple<bool, bool, UnaryMockParams> > {
protected:
    inline static int nrow = 57, ncol = 37;
    inline static std::vector<double> simulated;
    inline static std::shared_ptr<const tatami::NumericMatrix> dense, sparse;

    static void SetUpTestSuite() {
        simulated = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.1;
            opt.seed = 777222777;
            return opt;
        }());

        dense.reset(new tatami::DenseMatrix<double, int, decltype(simulated)>(nrow, ncol, simulated, true)); // row major
        sparse = tatami::convert_to_compressed_sparse<double, int>(*dense, false, {});  // column major
    }
};

TEST_P(DelayedUnaryIsometricOperationTest, Mock) {
    // Spamming a whole stack of tests.
    tatami_test::TestAccessOptions opts;
    auto tparam = GetParam();
    opts.use_row = std::get<0>(tparam);
    opts.use_oracle = std::get<1>(tparam);

    auto mockparams = std::get<2>(tparam);
    auto mockop = std::make_shared<UnaryMock<> >(mockparams);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, mockop);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, mockop);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.is_sparse_proportion(), 0);
    if (mockparams.sparse) {
        EXPECT_TRUE(sparse_mod.is_sparse());
        EXPECT_EQ(sparse_mod.is_sparse_proportion(), 1);
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
        EXPECT_EQ(sparse_mod.is_sparse_proportion(), 0);
    }

    auto refvec = simulated;
    size_t counter = 0;
    for (int r = 0; r < nrow; ++r) {
        for (int c = 0; c < ncol; ++c) {
            refvec[counter] = mock_operation(r, c, refvec[counter], mockparams);
            ++counter;
        }
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    tatami_test::test_full_access(dense_mod, ref, opts);
    tatami_test::test_block_access(dense_mod, ref, 0.25, 0.34, opts);
    tatami_test::test_indexed_access(dense_mod, ref, 0.3, 0.5, opts);

    tatami_test::test_full_access(sparse_mod, ref, opts);
    tatami_test::test_block_access(sparse_mod, ref, 0.5, 0.35, opts);
    tatami_test::test_indexed_access(sparse_mod, ref, 0.2, 0.4, opts);

    // Using a different type.
    {
        auto i_mockop = std::make_shared<UnaryMock<int> >(mockparams);
        tatami::DelayedUnaryIsometricOperation<int, double, int> i_dense_mod(dense, i_mockop);
        tatami::DelayedUnaryIsometricOperation<int, double, int> i_sparse_mod(sparse, i_mockop);
        tatami::DenseMatrix<int, int, std::vector<int> > i_ref(nrow, ncol, std::vector<int>(refvec.begin(), refvec.end()), true);

        tatami_test::test_full_access(i_dense_mod, i_ref, opts);
        tatami_test::test_block_access(i_dense_mod, i_ref, 0.3, 0.55, opts);
        tatami_test::test_indexed_access(i_dense_mod, i_ref, 0.3, 0.55, opts);

        tatami_test::test_full_access(i_sparse_mod, i_ref, opts);
        tatami_test::test_block_access(i_sparse_mod, i_ref, 0.25, 0.3, opts);
        tatami_test::test_indexed_access(i_sparse_mod, i_ref, 0.2, 0.3, opts);
    }
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricOperation,
    DelayedUnaryIsometricOperationTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row access
        ::testing::Values(true, false), // oracle usage
        ::testing::Values(
            UnaryMockParams({ false, true, true, true, true }),
            UnaryMockParams({ false, true, true, false, false }),
            UnaryMockParams({ false, false, false, true, true }),
            UnaryMockParams({ false, false, true, false, true }),
            UnaryMockParams({ true, false, false, false, true }),
            UnaryMockParams({ true, false, true, false, false }),
            UnaryMockParams({ true, false, true, false, true })
        )
    )
);
