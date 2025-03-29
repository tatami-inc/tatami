#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/binary/DelayedBinaryIsometricOperation.hpp"
#include "tatami/isometric/binary/helper_interface.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

struct BinaryMockParams {
    BinaryMockParams(bool sparse = false, bool zero_row = true, bool zero_col = true, bool non_zero_row = true, bool non_zero_col = true) :
        sparse(sparse), zero_row(zero_row), zero_col(zero_col), non_zero_row(non_zero_row), non_zero_col(non_zero_col) {}
    bool sparse, zero_row, zero_col, non_zero_row, non_zero_col;
};

template<typename Index_, typename Value_>
static Value_ mock_operation(Index_ row, Index_ col, Value_ left_val, Value_ right_val, const BinaryMockParams& optype) {
    if (left_val == 0 && right_val == 0) {
        if (optype.sparse) {
            return 0;
        }
        Value_ output = 0;
        if (optype.zero_row) {
            output += row;
        }
        if (optype.zero_col) {
            output += col * 13;
        }
        return output;
    } else {
        Value_ output = left_val + right_val * 2;
        if (optype.non_zero_row) {
            output += row;
        }
        if (optype.non_zero_col) {
            output += col * 13;
        }
        return output;
    }
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
class BinaryMock : public tatami::DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    BinaryMock() = default;
    BinaryMock(BinaryMockParams params) : BinaryMock(std::move(params), std::nullopt, std::nullopt) {}
    BinaryMock(BinaryMockParams params, std::optional<Index_> nrow, std::optional<Index_> ncol) : my_params(std::move(params)), my_nrow(nrow), my_ncol(ncol) {}

private:
    BinaryMockParams my_params;
    std::optional<Index_> my_nrow, my_ncol;

public:
    bool is_sparse() const { return my_params.sparse; }
    bool zero_depends_on_row() const { return my_params.zero_row; }
    bool zero_depends_on_column() const { return my_params.zero_col; }
    bool non_zero_depends_on_row() const { return my_params.non_zero_row; }
    bool non_zero_depends_on_column() const { return my_params.non_zero_col; }
    std::optional<Index_> nrow() const { return my_nrow; }
    std::optional<Index_> ncol() const { return my_ncol; }

public:
    OutputValue_ fill(bool row, int i) const {
        return mock_operation(
            row ? i : 0,
            row ? 0 : i,
            0,
            0,
            my_params
        );
    }

    void dense(bool row, Index_ i, Index_ start, Index_ length, const InputValue_* left, const InputValue_* right, OutputValue_* output) const {
        for (Index_ s = 0; s < length; ++s) {
            output[s] = mock_operation(
                row ? i : s + start,
                row ? s + start : i,
                left[s],
                right[s],
                my_params
            );
        }
    }

    void dense(bool row, Index_ i, const std::vector<Index_>& indices, const InputValue_* left, const InputValue_* right, OutputValue_* output) const {
        for (Index_ s = 0, end = indices.size(); s < end; ++s) {
            output[s] = mock_operation(
                row ? i : indices[s],
                row ? indices[s] : i,
                left[s],
                right[s],
                my_params
            );
        }
    }

    Index_ sparse(
        bool row,
        Index_ i,
        const tatami::SparseRange<InputValue_, Index_>& left,
        const tatami::SparseRange<InputValue_, Index_>& right,
        OutputValue_* output_value,
        Index_* output_index,
        bool needs_value,
        bool needs_index)
    const {
        Index_ lcount = 0, rcount = 0, output = 0;

        auto fun = [&](InputValue_ l, InputValue_ r, Index_ curindex) -> OutputValue_ {
            return mock_operation(
                row ? i : curindex,
                row ? curindex : i,
                l,
                r,
                my_params
            );
        };

        auto advance_left = [&]() -> void {
            if (needs_value) {
                output_value[output] = fun(left.value[lcount], 0, left.index[lcount]);
            }
            if (needs_index) {
                output_index[output] = left.index[lcount];
            }
            ++output;
            ++lcount;
        };

        auto advance_right = [&]() -> void {
            if (needs_value) {
                output_value[output] = fun(0, right.value[rcount], right.index[rcount]);
            }
            if (needs_index) {
                output_index[output] = right.index[rcount];
            }
            ++rcount;
            ++output;
        };

        while (lcount < left.number && rcount < right.number) {
            if (left.index[lcount] < right.index[rcount]) {
                advance_left();
            } else if (left.index[lcount] > right.index[rcount]) {
                advance_right();
            } else {
                if (needs_value) {
                    output_value[output] = fun(left.value[lcount], right.value[rcount], left.index[lcount]);
                }
                if (needs_index) {
                    output_index[output] = right.index[rcount];
                }
                ++lcount;
                ++rcount;
                ++output;
            }
        }

        while (lcount < left.number) {
            advance_left();
        }

        while (rcount < right.number) {
            advance_right();
        }

        return output;
    }
};

TEST(DelayedBinaryIsometricOperation, BackCompatibility) {
    int nrow = 23, ncol = 42;

    auto dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::vector<double>(nrow * ncol)));
    auto mat = tatami::make_DelayedBinaryIsometricOperation(dense, dense, std::make_shared<BinaryMock<> >());
    EXPECT_EQ(mat->nrow(), nrow);
    EXPECT_EQ(mat->ncol(), ncol);

    std::shared_ptr<const tatami::NumericMatrix> cdense(dense);
    auto cmat = tatami::make_DelayedBinaryIsometricOperation(cdense, cdense, std::make_shared<const BinaryMock<> >());
    EXPECT_EQ(cmat->nrow(), nrow);
    EXPECT_EQ(cmat->ncol(), ncol);
}

TEST(DelayedBinaryIsometricOperation, MatrixDimMismatch) {
    std::vector<double> src(200);
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(10, 20, src));
    auto dense2 = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(20, 10, src));
    tatami_test::throws_error([&]() {
        tatami::DelayedBinaryIsometricOperation<double, double, int> mat(dense, dense2, std::make_shared<BinaryMock<> >());
    }, "should be the same");
}

TEST(DelayedBinaryIsometricOperation, HelperDimMismatch) {
    auto dense = std::make_shared<tatami::DenseMatrix<double, int, std::vector<double> > >(10, 20, std::vector<double>(200), true);

    // No error when both match up.
    std::make_shared<tatami::DelayedBinaryIsometricOperation<double, double, int> >(
        dense,
        dense,
        std::make_shared<BinaryMock<> >(BinaryMockParams(), std::optional<int>{ 10 }, std::optional<int>{ 20 })
    );

    tatami_test::throws_error([&]{
        std::make_shared<tatami::DelayedBinaryIsometricOperation<double, double, int> >(
            dense,
            dense,
            std::make_shared<BinaryMock<> >(BinaryMockParams(), std::nullopt, std::optional<int>{ 10 })
        );
    }, "number of matrix columns");

    tatami_test::throws_error([&]{
        std::make_shared<tatami::DelayedBinaryIsometricOperation<double, double, int> >(
            dense,
            dense,
            std::make_shared<BinaryMock<> >(BinaryMockParams(), std::optional<int>{ 5 }, std::nullopt) 
        );
    }, "number of matrix rows");
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
        auto helper = std::make_shared<BinaryMock<> >(BinaryMockParams(/* sparse = */ true), std::nullopt, std::nullopt);
        tatami::DelayedBinaryIsometricOperation<double, double, int> mat1(dense, dense, helper);
        EXPECT_FALSE(mat1.is_sparse());
        tatami::DelayedBinaryIsometricOperation<double, double, int> mat2(dense, sparse, helper);
        EXPECT_FALSE(mat2.is_sparse());
        tatami::DelayedBinaryIsometricOperation<double, double, int> mat3(sparse, sparse, helper);
        EXPECT_TRUE(mat3.is_sparse());
    }

    // Always dense if the operation is not sparsity-preserving.
    {
        auto helper = std::make_shared<BinaryMock<> >(BinaryMockParams(/* sparse = */ false), std::nullopt, std::nullopt);
        tatami::DelayedBinaryIsometricOperation<double, double, int> mat1(dense, dense, helper);
        EXPECT_FALSE(mat1.is_sparse());
        tatami::DelayedBinaryIsometricOperation<double, double, int> mat2(dense, sparse, helper);
        EXPECT_FALSE(mat2.is_sparse());
        tatami::DelayedBinaryIsometricOperation<double, double, int> mat3(sparse, sparse, helper);
        EXPECT_FALSE(mat3.is_sparse());
    }
}


TEST(DelayedBinaryIsometricOperation, DependsChecks) {
    auto optr = std::make_shared<const tatami::ConsecutiveOracle<int> >(10, 20);

    { 
        // Check that the oracle is correctly used for both rows and columns.
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMock<>, int> oracle1(optr, {}, true);
        EXPECT_EQ(oracle1.get(0), 10);
        tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMock<>, int> oracle2(optr, {}, false);
        EXPECT_EQ(oracle2.get(0), 10);

        EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(BinaryMock<>(), true));
        EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(BinaryMock<>(), false));
    }

    { 
        // Check that the oracle is correctly used for rows only.
        {
            BinaryMock<> mock({ false, true, false, false, false }); // zero-dependency for non-sparse operation.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMock<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 10);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMock<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            BinaryMock<> mock({ false, false, false, true, false }); // non-zero dependency.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMock<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 10);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMock<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            BinaryMock<> mock({ true, true, false, false, false }); // zero dependency ignored for sparse operations.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMock<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMock<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }
    }

    { 
        // Check that the oracle is correctly used for columns only.
        {
            BinaryMock<> mock({ false, false, true, false, false }); // zero-dependency for non-sparse operation.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMock<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMock<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 10);

            EXPECT_FALSE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            BinaryMock<> mock({ false, false, false, false, true }); // non-zero dependency.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMock<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMock<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 10);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }

        {
            BinaryMock<> mock({ true, false, true, false, false }); // zero dependency ignored for sparse operations.
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMock<>, int> oracle1(optr, mock, true);
            EXPECT_EQ(oracle1.get(0), 0);
            tatami::DelayedIsometricOperation_internal::MaybeOracleDepends<true, BinaryMock<>, int> oracle2(optr, mock, false);
            EXPECT_EQ(oracle2.get(0), 0);

            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, true));
            EXPECT_TRUE(tatami::DelayedIsometricOperation_internal::can_dense_expand(mock, false));
        }
    }
}

class DelayedBinaryIsometricOperationTest : public ::testing::TestWithParam<std::tuple<bool, bool, BinaryMockParams> > {
protected:
    inline static int nrow = 23, ncol = 42;
    inline static std::vector<double> lsimulated, rsimulated;
    inline static std::shared_ptr<tatami::Matrix<double, int> > ldense, lsparse, rdense, rsparse;

    static void SetUpTestSuite() {
        lsimulated = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.2;
            opt.seed = 918273;
            return opt;
        }());
        ldense.reset(new tatami::DenseMatrix<double, int, decltype(lsimulated)>(nrow, ncol, lsimulated, true));
        lsparse = tatami::convert_to_compressed_sparse<double, int>(*ldense, false, {}); 

        rsimulated = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.15;
            opt.seed = 1238712;
            return opt;
        }());
        rdense.reset(new tatami::DenseMatrix<double, int, decltype(rsimulated)>(nrow, ncol, rsimulated, true));
        rsparse = tatami::convert_to_compressed_sparse<double, int>(*rdense, false, {}); 

    }
};

TEST_P(DelayedBinaryIsometricOperationTest, Mock) {
    // Spamming a whole stack of tests.
    tatami_test::TestAccessOptions opts;
    auto tparam = GetParam();
    opts.use_row = std::get<0>(tparam);
    opts.use_oracle = std::get<1>(tparam);

    auto mockparams = std::get<2>(tparam);
    auto mockop = std::make_shared<BinaryMock<> >(mockparams);
    tatami::DelayedBinaryIsometricOperation<double, double, int> dense_mod(ldense, rdense, mockop);
    tatami::DelayedBinaryIsometricOperation<double, double, int> sparse_mod(lsparse, rsparse, mockop);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.is_sparse_proportion(), 0);
    if (mockparams.sparse) {
        EXPECT_TRUE(sparse_mod.is_sparse());
        EXPECT_EQ(sparse_mod.is_sparse_proportion(), 1);
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
        EXPECT_EQ(sparse_mod.is_sparse_proportion(), 0);
    }

    std::vector<double> refvec(lsimulated.size());
    size_t counter = 0;
    for (int r = 0; r < nrow; ++r) {
        for (int c = 0; c < ncol; ++c) {
            refvec[counter] = mock_operation(r, c, lsimulated[counter], rsimulated[counter], mockparams);
            ++counter;
        }
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    tatami_test::test_full_access(dense_mod, ref, opts);
    tatami_test::test_block_access(dense_mod, ref, 0.1, 0.7, opts);
    tatami_test::test_indexed_access(dense_mod, ref, 0.23, 0.5, opts);

    tatami_test::test_full_access(sparse_mod, ref, opts);
    tatami_test::test_block_access(sparse_mod, ref, 0.13, 0.5, opts);
    tatami_test::test_indexed_access(sparse_mod, ref, 0.23, 0.4, opts);

    // Using a different type.
    {
        auto f_mockop = std::make_shared<BinaryMock<float> >(mockparams);
        tatami::DelayedBinaryIsometricOperation<float, double, int> f_dense_mod(ldense, rdense, f_mockop);
        tatami::DelayedBinaryIsometricOperation<float, double, int> f_sparse_mod(lsparse, rsparse, f_mockop);
        tatami::DenseMatrix<float, int, std::vector<float> > f_ref(nrow, ncol, std::vector<float>(refvec.begin(), refvec.end()), true);

        tatami_test::test_full_access(f_dense_mod, f_ref, opts);
        tatami_test::test_block_access(f_dense_mod, f_ref, 0.5, 0.5, opts);
        tatami_test::test_indexed_access(f_dense_mod, f_ref, 0.3, 0.5, opts);

        tatami_test::test_full_access(f_sparse_mod, f_ref, opts);
        tatami_test::test_block_access(f_sparse_mod, f_ref, 0.2, 0.6, opts);
        tatami_test::test_indexed_access(f_sparse_mod, f_ref, 0.2, 0.4, opts);
    }
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBinaryIsometricOperation,
    DelayedBinaryIsometricOperationTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row access
        ::testing::Values(true, false), // oracle usage
        ::testing::Values(
            BinaryMockParams({ false, true, true, true, true }),
            BinaryMockParams({ false, true, true, false, false }),
            BinaryMockParams({ false, false, false, true, true }),
            BinaryMockParams({ false, false, true, false, true }),
            BinaryMockParams({ true, false, false, false, true }),
            BinaryMockParams({ true, false, true, false, false }),
            BinaryMockParams({ true, false, true, false, true })
        )
    )
);
