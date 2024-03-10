#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/other/DelayedCast.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

template<class PARAM> 
class CastTest : public ::testing::TestWithParam<PARAM> {
protected:
    size_t nrow = 99, ncol = 179;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    std::shared_ptr<tatami::Matrix<float, size_t> > fdense, fsparse;
    std::shared_ptr<tatami::NumericMatrix> fdense_ref, fsparse_ref;

    std::shared_ptr<tatami::Matrix<float, int> > fsparse_value;
    std::shared_ptr<tatami::Matrix<double, size_t> > sparse_index;
protected:
    void SetUp() {
        auto sparse_matrix = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.08);
        {
            dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, sparse_matrix));
            sparse = tatami::convert_to_compressed_sparse<false>(dense.get()); // column-major.
        }

        // Both the value and indices are changed in type.
        std::vector<float> fsparse_matrix(sparse_matrix.begin(), sparse_matrix.end());
        {
            fdense = std::shared_ptr<tatami::Matrix<float, size_t> >(new tatami::DenseRowMatrix<float, size_t>(nrow, ncol, fsparse_matrix));
            fsparse = tatami::convert_to_compressed_sparse<false, float, size_t>(fdense.get()); // column-major.
        }

        // Reference with reduced precision, for comparison with double->float->double casts.
        {
            std::vector<double> dsparse_matrix(fsparse_matrix.begin(), fsparse_matrix.end());
            fdense_ref = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(dsparse_matrix)));
            fsparse_ref = tatami::convert_to_compressed_sparse<false>(fdense_ref.get()); // column-major.
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

        return;
    }
};

/****************************************************
 ****************************************************/

class DelayedCastFullAccess : public CastTest<std::tuple<bool, bool, tatami_test::TestAccessOrder, size_t> > {};

TEST_P(DelayedCastFullAccess, Dense) {
    auto tparam = GetParam(); 

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<0>(tparam);
    params.use_oracle = std::get<1>(tparam);
    params.order = std::get<2>(tparam);
    params.jump = std::get<3>(tparam);

    auto cast_dense = tatami::make_DelayedCast<double, int>(dense);
    EXPECT_EQ(cast_dense->nrow(), nrow);
    EXPECT_EQ(cast_dense->ncol(), ncol);
    EXPECT_EQ(cast_dense->sparse(), dense->sparse());
    EXPECT_EQ(cast_dense->sparse_proportion(), dense->sparse_proportion());
    EXPECT_EQ(cast_dense->prefer_rows(), dense->prefer_rows());
    EXPECT_EQ(cast_dense->prefer_rows_proportion(), dense->prefer_rows_proportion());
    tatami_test::test_full_access(params, cast_dense.get(), dense.get());

    auto cast_fdense = tatami::make_DelayedCast<double, int>(fdense);
    tatami_test::test_full_access(params, cast_fdense.get(), fdense_ref.get());
}

TEST_P(DelayedCastFullAccess, Sparse) {
    auto tparam = GetParam(); 

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<0>(tparam);
    params.use_oracle = std::get<1>(tparam);
    params.order = std::get<2>(tparam);
    params.jump = std::get<3>(tparam);

    auto cast_sparse = tatami::make_DelayedCast<double, int>(sparse);
    EXPECT_EQ(cast_sparse->nrow(), nrow);
    EXPECT_EQ(cast_sparse->ncol(), ncol);
    EXPECT_EQ(cast_sparse->sparse(), sparse->sparse());
    EXPECT_EQ(cast_sparse->sparse_proportion(), sparse->sparse_proportion());
    EXPECT_EQ(cast_sparse->prefer_rows(), sparse->prefer_rows());
    EXPECT_EQ(cast_sparse->prefer_rows_proportion(), sparse->prefer_rows_proportion());
    tatami_test::test_full_access(params, cast_sparse.get(), sparse.get());

    auto cast_fsparse = tatami::make_DelayedCast<double, int>(fsparse);
    tatami_test::test_full_access(params, cast_fsparse.get(), fsparse_ref.get());

    auto cast_fsparse_value = tatami::make_DelayedCast<double, int>(fsparse_value);
    tatami_test::test_full_access(params, cast_fsparse_value.get(), fsparse_ref.get());

    auto cast_sparse_index = tatami::make_DelayedCast<double, int>(sparse_index);
    tatami_test::test_full_access(params, cast_sparse_index.get(), sparse.get());
}

INSTANTIATE_TEST_SUITE_P(
    DelayedCast,
    DelayedCastFullAccess,
    ::testing::Combine(
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3) // jump, to test the workspace memory.
    )
);

/****************************************************
 ****************************************************/

class DelayedCastBlockAccess : public CastTest<std::tuple<bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > > {};

TEST_P(DelayedCastBlockAccess, Dense) {
    auto tparam = GetParam(); 

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<0>(tparam);
    params.use_oracle = std::get<1>(tparam);
    params.order = std::get<2>(tparam);
    params.jump = std::get<3>(tparam);

    auto interval_info = std::get<4>(tparam);
    auto len = (params.use_row ? ncol : nrow);
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;

    auto cast_dense = tatami::make_DelayedCast<double, int>(dense);
    tatami_test::test_block_access(params, cast_dense.get(), dense.get(), FIRST, LAST);

    auto cast_fdense = tatami::make_DelayedCast<double, int>(fdense);
    tatami_test::test_block_access(params, cast_fdense.get(), fdense_ref.get(), FIRST, LAST);
}

TEST_P(DelayedCastBlockAccess, Sparse) {
    auto tparam = GetParam();

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<0>(tparam);
    params.use_oracle = std::get<1>(tparam);
    params.order = std::get<2>(tparam);
    params.jump = std::get<3>(tparam);

    auto interval_info = std::get<4>(tparam);
    auto len = (params.use_row ? ncol : nrow);
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;

    auto cast_sparse = tatami::make_DelayedCast<double, int>(sparse);
    tatami_test::test_block_access(params, cast_sparse.get(), sparse.get(), FIRST, LAST);

    auto cast_fsparse = tatami::make_DelayedCast<double, int>(fsparse);
    tatami_test::test_block_access(params, cast_fsparse.get(), fsparse_ref.get(), FIRST, LAST);

    auto cast_fsparse_value = tatami::make_DelayedCast<double, int>(fsparse_value);
    tatami_test::test_block_access(params, cast_fsparse_value.get(), fsparse_ref.get(), FIRST, LAST);

    auto cast_sparse_index = tatami::make_DelayedCast<double, int>(sparse_index);
    tatami_test::test_block_access(params, cast_sparse_index.get(), sparse.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedCast,
    DelayedCastBlockAccess,
    ::testing::Combine(
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to check the workspace memory
        ::testing::Values(
            std::make_pair(0, 0.6), 
            std::make_pair(0.5, 0.75),
            std::make_pair(0.24, 1)
        )
    )
);

/****************************************************
 ****************************************************/

class DelayedCastIndexAccess : public CastTest<std::tuple<bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > > {};

TEST_P(DelayedCastIndexAccess, Dense) {
    auto tparam = GetParam(); 

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<0>(tparam);
    params.use_oracle = std::get<1>(tparam);
    params.order = std::get<2>(tparam);
    params.jump = std::get<3>(tparam);

    auto interval_info = std::get<4>(tparam);
    auto len = (params.use_row ? ncol : nrow);
    size_t FIRST = interval_info.first * len, STEP = interval_info.second * len;

    auto cast_dense = tatami::make_DelayedCast<double, int>(dense);
    tatami_test::test_indexed_access(params, cast_dense.get(), dense.get(), FIRST, STEP);

    auto cast_fdense = tatami::make_DelayedCast<double, int>(fdense);
    tatami_test::test_indexed_access(params, cast_fdense.get(), fdense_ref.get(), FIRST, STEP);
}

TEST_P(DelayedCastIndexAccess, Sparse) {
    auto tparam = GetParam();

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<0>(tparam);
    params.use_oracle = std::get<1>(tparam);
    params.order = std::get<2>(tparam);
    params.jump = std::get<3>(tparam);

    auto interval_info = std::get<4>(tparam);
    auto len = (params.use_row ? ncol : nrow);
    size_t FIRST = interval_info.first * len, STEP = interval_info.second * len;

    auto cast_sparse = tatami::make_DelayedCast<double, int>(sparse);
    tatami_test::test_indexed_access(params, cast_sparse.get(), sparse.get(), FIRST, STEP);

    auto cast_fsparse = tatami::make_DelayedCast<double, int>(fsparse);
    tatami_test::test_indexed_access(params, cast_fsparse.get(), fsparse_ref.get(), FIRST, STEP);

    auto cast_fsparse_value = tatami::make_DelayedCast<double, int>(fsparse_value);
    tatami_test::test_indexed_access(params, cast_fsparse_value.get(), fsparse_ref.get(), FIRST, STEP);

    auto cast_sparse_index = tatami::make_DelayedCast<double, int>(sparse_index);
    tatami_test::test_indexed_access(params, cast_sparse_index.get(), sparse.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedCast,
    DelayedCastIndexAccess,
    ::testing::Combine(
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to check the workspace memory
        ::testing::Values(
            std::make_pair(0, 0.015), 
            std::make_pair(0.5, 0.025),
            std::make_pair(0.666, 0.055)
        )
    )
);

/****************************************************
 ****************************************************/

TEST(DelayedCastTest, ConstOverload) {
    int nrow = 10, ncol = 15;
    auto simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.05);
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated));
    auto tdense = tatami::make_DelayedCast<float, size_t>(dense);

    // Cursory checks.
    EXPECT_EQ(dense->nrow(), tdense->nrow());
    EXPECT_EQ(dense->ncol(), tdense->ncol());
}

