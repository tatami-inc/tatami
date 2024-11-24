#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/other/DelayedCast.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

class CastUtils {
protected:
    inline static size_t nrow = 99, ncol = 179;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, sparse;

    inline static std::shared_ptr<tatami::Matrix<float, size_t> > fdense, fsparse;
    inline static std::shared_ptr<tatami::NumericMatrix> fdense_ref, fsparse_ref;

    inline static std::shared_ptr<tatami::Matrix<float, int> > fsparse_value;
    inline static std::shared_ptr<tatami::Matrix<double, size_t> > sparse_index;

    inline static std::shared_ptr<tatami::NumericMatrix> cast_dense, cast_fdense;
    inline static std::shared_ptr<tatami::NumericMatrix> cast_sparse, cast_fsparse, cast_fsparse_value, cast_sparse_index;

    inline static std::shared_ptr<tatami::NumericMatrix> uns_cast_fsparse;

    static void assemble() {
        if (dense) {
            return;
        }

        auto sparse_matrix = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.08;
            opt.seed = 8764823;
            return opt;
        }());

        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, sparse_matrix));
        sparse = tatami::convert_to_compressed_sparse<false, double, int>(dense.get()); // column-major.

        // Both the value and indices are changed in type.
        std::vector<float> fsparse_matrix(sparse_matrix.begin(), sparse_matrix.end());
        fdense = std::shared_ptr<tatami::Matrix<float, size_t> >(new tatami::DenseRowMatrix<float, size_t>(nrow, ncol, fsparse_matrix));
        fsparse = tatami::convert_to_compressed_sparse<false, float, size_t>(fdense.get()); // column-major.

        // Reference with reduced precision, for comparison with double->float->double casts.
        {
            std::vector<double> dsparse_matrix(fsparse_matrix.begin(), fsparse_matrix.end());
            fdense_ref = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(dsparse_matrix)));
            fsparse_ref = tatami::convert_to_compressed_sparse<false, double, int>(fdense_ref.get()); // column-major.
        }

        // Only the value is changed in type.
        {
            auto refdense = std::shared_ptr<tatami::Matrix<float, int> >(new tatami::DenseRowMatrix<float, int>(nrow, ncol, fsparse_matrix));
            fsparse_value = tatami::convert_to_compressed_sparse<false, float, int>(refdense.get()); 
        }

        // Only the index is changed in type.
        {
            auto redense = std::shared_ptr<tatami::Matrix<double, size_t> >(new tatami::DenseRowMatrix<double, size_t>(nrow, ncol, sparse_matrix));
            sparse_index = tatami::convert_to_compressed_sparse<false, double, size_t>(redense.get()); 
        }

        cast_dense = tatami::make_DelayedCast<double, int>(dense);
        cast_fdense = tatami::make_DelayedCast<double, int>(fdense);
        cast_sparse = tatami::make_DelayedCast<double, int>(sparse);
        cast_fsparse = tatami::make_DelayedCast<double, int>(fsparse);
        cast_fsparse_value = tatami::make_DelayedCast<double, int>(fsparse_value);
        cast_sparse_index = tatami::make_DelayedCast<double, int>(sparse_index);

        // Testing unsorted access.
        uns_cast_fsparse = tatami::make_DelayedCast<double, int>(
            std::shared_ptr<tatami::Matrix<float, size_t> >(
                new tatami_test::ReversedIndicesWrapper<float, size_t>(fsparse)
            )
        );
    }
};

class DelayedCastTest : public ::testing::Test, public CastUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_F(DelayedCastTest, Dense) {
    EXPECT_EQ(cast_dense->nrow(), nrow);
    EXPECT_EQ(cast_dense->ncol(), ncol);
    EXPECT_EQ(cast_dense->is_sparse(), dense->is_sparse());
    EXPECT_EQ(cast_dense->is_sparse_proportion(), dense->is_sparse_proportion());
    EXPECT_EQ(cast_dense->prefer_rows(), dense->prefer_rows());
    EXPECT_EQ(cast_dense->prefer_rows_proportion(), dense->prefer_rows_proportion());

    EXPECT_FALSE(cast_dense->uses_oracle(true));
}

TEST_F(DelayedCastTest, Sparse) {
    EXPECT_EQ(cast_sparse->nrow(), nrow);
    EXPECT_EQ(cast_sparse->ncol(), ncol);
    EXPECT_EQ(cast_sparse->is_sparse(), sparse->is_sparse());
    EXPECT_EQ(cast_sparse->is_sparse_proportion(), sparse->is_sparse_proportion());
    EXPECT_EQ(cast_sparse->prefer_rows(), sparse->prefer_rows());
    EXPECT_EQ(cast_sparse->prefer_rows_proportion(), sparse->prefer_rows_proportion());

    EXPECT_FALSE(cast_sparse->uses_oracle(true));
}

TEST_F(DelayedCastTest, ConstOverload) {
    std::shared_ptr<const tatami::NumericMatrix> const_dense = dense;
    auto tdense = tatami::make_DelayedCast<float, size_t>(const_dense);

    // Cursory checks.
    EXPECT_EQ(dense->nrow(), tdense->nrow());
    EXPECT_EQ(dense->ncol(), tdense->ncol());
}

TEST(DelayedCastMisc, CastOracle) {
    auto out = std::make_shared<tatami::ConsecutiveOracle<int> >(10, 100);
    auto casted = tatami::DelayedCast_internal::CastOracle<size_t, int>(out);
    EXPECT_EQ(casted.total(), 100);
}

/****************************************************
 ****************************************************/

class DelayedCastFullAccessTest : 
    public ::testing::TestWithParam<tatami_test::StandardTestAccessOptions>,
    public CastUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedCastFullAccessTest, Dense) {
    auto tparam = GetParam(); 
    auto options = tatami_test::convert_test_access_options(tparam);
    tatami_test::test_full_access(*cast_dense, *dense, options);
    tatami_test::test_full_access(*cast_fdense, *fdense_ref, options);
}

TEST_P(DelayedCastFullAccessTest, Sparse) {
    auto tparam = GetParam(); 
    auto options = tatami_test::convert_test_access_options(tparam);
    tatami_test::test_full_access(*cast_sparse, *sparse, options);
    tatami_test::test_full_access(*cast_fsparse, *fsparse_ref, options);
    tatami_test::test_full_access(*cast_fsparse_value, *fsparse_ref, options);
    tatami_test::test_full_access(*cast_sparse_index, *sparse, options);
    tatami_test::test_unsorted_full_access(*uns_cast_fsparse, options);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedCast,
    DelayedCastFullAccessTest,
    tatami_test::standard_test_access_options_combinations()
);

/****************************************************
 ****************************************************/

class DelayedCastBlockAccessTest : 
    public ::testing::TestWithParam<std::tuple<tatami_test::StandardTestAccessOptions, std::pair<double, double> > >, 
    public CastUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedCastBlockAccessTest, Dense) {
    auto tparam = GetParam(); 
    auto options = tatami_test::convert_test_access_options(std::get<0>(tparam));
    auto interval_info = std::get<1>(tparam);
    tatami_test::test_block_access(*cast_dense, *dense, interval_info.first, interval_info.second, options);
    tatami_test::test_block_access(*cast_fdense, *fdense_ref, interval_info.first, interval_info.second, options);
}

TEST_P(DelayedCastBlockAccessTest, Sparse) {
    auto tparam = GetParam();
    auto options = tatami_test::convert_test_access_options(std::get<0>(tparam));
    auto interval_info = std::get<1>(tparam);
    tatami_test::test_block_access(*cast_sparse, *sparse, interval_info.first, interval_info.second, options);
    tatami_test::test_block_access(*cast_fsparse, *fsparse_ref, interval_info.first, interval_info.second, options);
    tatami_test::test_block_access(*cast_fsparse_value, *fsparse_ref, interval_info.first, interval_info.second, options);
    tatami_test::test_block_access(*cast_sparse_index, *sparse, interval_info.first, interval_info.second, options);
    tatami_test::test_unsorted_block_access(*uns_cast_fsparse, interval_info.first, interval_info.second, options);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedCast,
    DelayedCastBlockAccessTest,
    ::testing::Combine(
        tatami_test::standard_test_access_options_combinations(),
        ::testing::Values(
            std::make_pair(0.0, 0.6), 
            std::make_pair(0.5, 0.25),
            std::make_pair(0.24, 0.76)
        )
    )
);

/****************************************************
 ****************************************************/

class DelayedCastIndexAccessTest : 
    public ::testing::TestWithParam<std::tuple<tatami_test::StandardTestAccessOptions, std::pair<double, double> > >, 
    public CastUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedCastIndexAccessTest, Dense) {
    auto tparam = GetParam(); 
    auto options = tatami_test::convert_test_access_options(std::get<0>(tparam));
    auto interval_info = std::get<1>(tparam);
    tatami_test::test_indexed_access(*cast_dense, *dense, interval_info.first, interval_info.second, options);
    tatami_test::test_indexed_access(*cast_fdense, *fdense_ref, interval_info.first, interval_info.second, options);
}

TEST_P(DelayedCastIndexAccessTest, Sparse) {
    auto tparam = GetParam();
    auto options = tatami_test::convert_test_access_options(std::get<0>(tparam));
    auto interval_info = std::get<1>(tparam);
    tatami_test::test_indexed_access(*cast_sparse, *sparse, interval_info.first, interval_info.second, options);
    tatami_test::test_indexed_access(*cast_fsparse, *fsparse_ref, interval_info.first, interval_info.second, options);
    tatami_test::test_indexed_access(*cast_fsparse_value, *fsparse_ref, interval_info.first, interval_info.second, options);
    tatami_test::test_indexed_access(*cast_sparse_index, *sparse, interval_info.first, interval_info.second, options);
    tatami_test::test_unsorted_indexed_access(*uns_cast_fsparse, interval_info.first, interval_info.second, options);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedCast,
    DelayedCastIndexAccessTest,
    ::testing::Combine(
        tatami_test::standard_test_access_options_combinations(),
        ::testing::Values(
            std::make_pair(0.0, 0.3), 
            std::make_pair(0.5, 0.15),
            std::make_pair(0.666, 0.2)
        )
    )
);
