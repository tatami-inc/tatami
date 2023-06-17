#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/other/DelayedCast.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

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
            sparse = tatami::convert_to_sparse<false>(dense.get()); // column-major.
        }

        // Both the value and indices are changed in type.
        std::vector<float> fsparse_matrix(sparse_matrix.begin(), sparse_matrix.end());
        {
            fdense = std::shared_ptr<tatami::Matrix<float, size_t> >(new tatami::DenseRowMatrix<float, size_t>(nrow, ncol, fsparse_matrix));
            fsparse = tatami::convert_to_sparse<false, float, size_t>(fdense.get()); // column-major.
        }

        // Reference with reduced precision, for comparison with double->float->double casts.
        {
            std::vector<double> dsparse_matrix(fsparse_matrix.begin(), fsparse_matrix.end());
            fdense_ref = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(dsparse_matrix)));
            fsparse_ref = tatami::convert_to_sparse<false>(fdense_ref.get()); // column-major.
        }

        // Only the value is changed in type.
        {
            auto refdense = std::shared_ptr<tatami::Matrix<float, int> >(new tatami::DenseRowMatrix<float, int>(nrow, ncol, fsparse_matrix));
            fsparse_value = tatami::convert_to_sparse<false, float, int>(refdense.get()); 
        }

        // Only the index is changed in type.
        {
            auto redense = std::shared_ptr<tatami::Matrix<double, size_t> >(new tatami::DenseRowMatrix<double, size_t>(nrow, ncol, sparse_matrix));
            sparse_index = tatami::convert_to_sparse<false, double, size_t>(redense.get()); 
        }

        return;
    }
};

/****************************************************
 ****************************************************/

class DelayedCastFullAccess : public CastTest<std::tuple<bool, size_t> > {};

TEST_P(DelayedCastFullAccess, Dense) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

    auto cast_dense = tatami::make_DelayedCast<double, int>(dense);
    EXPECT_EQ(cast_dense->nrow(), nrow);
    EXPECT_EQ(cast_dense->ncol(), ncol);
    EXPECT_EQ(cast_dense->sparse(), dense->sparse());
    EXPECT_EQ(cast_dense->sparse_proportion(), dense->sparse_proportion());
    EXPECT_EQ(cast_dense->prefer_rows(), dense->prefer_rows());
    EXPECT_EQ(cast_dense->prefer_rows_proportion(), dense->prefer_rows_proportion());

    tatami_test::test_simple_row_access(cast_dense.get(), dense.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(cast_dense.get(), dense.get(), FORWARD, JUMP);

    auto cast_fdense = tatami::make_DelayedCast<double, int>(fdense);
    tatami_test::test_simple_row_access(cast_fdense.get(), fdense_ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(cast_fdense.get(), fdense_ref.get(), FORWARD, JUMP);
}

TEST_P(DelayedCastFullAccess, Sparse) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

    auto cast_sparse = tatami::make_DelayedCast<double, int>(sparse);
    EXPECT_EQ(cast_sparse->nrow(), nrow);
    EXPECT_EQ(cast_sparse->ncol(), ncol);
    EXPECT_EQ(cast_sparse->sparse(), sparse->sparse());
    EXPECT_EQ(cast_sparse->sparse_proportion(), sparse->sparse_proportion());
    EXPECT_EQ(cast_sparse->prefer_rows(), sparse->prefer_rows());
    EXPECT_EQ(cast_sparse->prefer_rows_proportion(), sparse->prefer_rows_proportion());

    tatami_test::test_simple_row_access(cast_sparse.get(), sparse.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(cast_sparse.get(), sparse.get(), FORWARD, JUMP);

    auto cast_fsparse = tatami::make_DelayedCast<double, int>(fsparse);
    tatami_test::test_simple_row_access(cast_fsparse.get(), fsparse_ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(cast_fsparse.get(), fsparse_ref.get(), FORWARD, JUMP);

    auto cast_fsparse_value = tatami::make_DelayedCast<double, int>(fsparse_value);
    tatami_test::test_simple_row_access(cast_fsparse_value.get(), fsparse_ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(cast_fsparse_value.get(), fsparse_ref.get(), FORWARD, JUMP);

    auto cast_sparse_index = tatami::make_DelayedCast<double, int>(sparse_index);
    tatami_test::test_simple_row_access(cast_sparse_index.get(), sparse.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(cast_sparse_index.get(), sparse.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    DelayedCast,
    DelayedCastFullAccess,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 3) // jump, to test the workspace memory.
    )
);

/****************************************************
 ****************************************************/

class DelayedCastBlockAccess : public CastTest<std::tuple<bool, size_t, std::vector<double> > > {};

TEST_P(DelayedCastBlockAccess, Dense) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t RFIRST = interval_info[0] * nrow, RLAST = interval_info[1] * nrow;
    size_t CFIRST = interval_info[0] * ncol, CLAST = interval_info[1] * ncol;

    auto cast_dense = tatami::make_DelayedCast<double, int>(dense);
    tatami_test::test_sliced_row_access(cast_dense.get(), dense.get(), true, JUMP, CFIRST, CLAST);
    tatami_test::test_sliced_column_access(cast_dense.get(), dense.get(), true, JUMP, RFIRST, RLAST);

    auto cast_fdense = tatami::make_DelayedCast<double, int>(fdense);
    tatami_test::test_sliced_row_access(cast_fdense.get(), fdense_ref.get(), true, JUMP, CFIRST, CLAST);
    tatami_test::test_sliced_column_access(cast_fdense.get(), fdense_ref.get(), true, JUMP, RFIRST, RLAST);
}

TEST_P(DelayedCastBlockAccess, Sparse) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t RFIRST = interval_info[0] * nrow, RLAST = interval_info[1] * nrow;
    size_t CFIRST = interval_info[0] * ncol, CLAST = interval_info[1] * ncol;

    auto cast_sparse = tatami::make_DelayedCast<double, int>(sparse);
    tatami_test::test_sliced_row_access(cast_sparse.get(), sparse.get(), true, JUMP, CFIRST, CLAST);
    tatami_test::test_sliced_column_access(cast_sparse.get(), sparse.get(), true, JUMP, RFIRST, RLAST);

    auto cast_fsparse = tatami::make_DelayedCast<double, int>(fsparse);
    tatami_test::test_sliced_row_access(cast_fsparse.get(), fsparse_ref.get(), true, JUMP, CFIRST, CLAST);
    tatami_test::test_sliced_column_access(cast_fsparse.get(), fsparse_ref.get(), true, JUMP, RFIRST, RLAST);

    auto cast_fsparse_value = tatami::make_DelayedCast<double, int>(fsparse_value);
    tatami_test::test_sliced_row_access(cast_fsparse_value.get(), fsparse_ref.get(), true, JUMP, CFIRST, CLAST);
    tatami_test::test_sliced_column_access(cast_fsparse_value.get(), fsparse_ref.get(), true, JUMP, RFIRST, RLAST);

    auto cast_sparse_index = tatami::make_DelayedCast<double, int>(sparse_index);
    tatami_test::test_sliced_row_access(cast_sparse_index.get(), sparse.get(), true, JUMP, CFIRST, CLAST);
    tatami_test::test_sliced_column_access(cast_sparse_index.get(), sparse.get(), true, JUMP, RFIRST, RLAST);
}

INSTANTIATE_TEST_CASE_P(
    DelayedCast,
    DelayedCastBlockAccess,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 3), // jump, to check the workspace memory
        ::testing::Values(
            std::vector<double>({ 0, 0.6 }), 
            std::vector<double>({ 0.5, 0.75 }),
            std::vector<double>({ 0.24, 1 })
        )
    )
);

/****************************************************
 ****************************************************/

class DelayedCastIndexAccess : public CastTest<std::tuple<bool, size_t, std::vector<double> > > {};

TEST_P(DelayedCastIndexAccess, Dense) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t RFIRST = interval_info[0] * nrow, RSTEP = interval_info[1] * nrow;
    size_t CFIRST = interval_info[0] * ncol, CSTEP = interval_info[1] * ncol;

    auto cast_dense = tatami::make_DelayedCast<double, int>(dense);
    tatami_test::test_indexed_row_access(cast_dense.get(), dense.get(), true, JUMP, CFIRST, CSTEP);
    tatami_test::test_indexed_column_access(cast_dense.get(), dense.get(), true, JUMP, RFIRST, RSTEP);

    auto cast_fdense = tatami::make_DelayedCast<double, int>(fdense);
    tatami_test::test_indexed_row_access(cast_fdense.get(), fdense_ref.get(), true, JUMP, CFIRST, CSTEP);
    tatami_test::test_indexed_column_access(cast_fdense.get(), fdense_ref.get(), true, JUMP, RFIRST, RSTEP);
}

TEST_P(DelayedCastIndexAccess, Sparse) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t RFIRST = interval_info[0] * nrow, RSTEP = interval_info[1] * nrow;
    size_t CFIRST = interval_info[0] * ncol, CSTEP = interval_info[1] * ncol;

    auto cast_sparse = tatami::make_DelayedCast<double, int>(sparse);
    tatami_test::test_indexed_row_access(cast_sparse.get(), sparse.get(), true, JUMP, CFIRST, CSTEP);
    tatami_test::test_indexed_column_access(cast_sparse.get(), sparse.get(), true, JUMP, RFIRST, RSTEP);

    auto cast_fsparse = tatami::make_DelayedCast<double, int>(fsparse);
    tatami_test::test_indexed_row_access(cast_fsparse.get(), fsparse_ref.get(), true, JUMP, CFIRST, CSTEP);
    tatami_test::test_indexed_column_access(cast_fsparse.get(), fsparse_ref.get(), true, JUMP, RFIRST, RSTEP);

    auto cast_fsparse_value = tatami::make_DelayedCast<double, int>(fsparse_value);
    tatami_test::test_indexed_row_access(cast_fsparse_value.get(), fsparse_ref.get(), true, JUMP, CFIRST, CSTEP);
    tatami_test::test_indexed_column_access(cast_fsparse_value.get(), fsparse_ref.get(), true, JUMP, RFIRST, RSTEP);

    auto cast_sparse_index = tatami::make_DelayedCast<double, int>(sparse_index);
    tatami_test::test_indexed_row_access(cast_sparse_index.get(), sparse.get(), true, JUMP, CFIRST, CSTEP);
    tatami_test::test_indexed_column_access(cast_sparse_index.get(), sparse.get(), true, JUMP, RFIRST, RSTEP);
}

INSTANTIATE_TEST_CASE_P(
    DelayedCast,
    DelayedCastIndexAccess,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 3), // jump, to check the workspace memory
        ::testing::Values(
            std::vector<double>({ 0, 0.015 }), 
            std::vector<double>({ 0.5, 0.025 }),
            std::vector<double>({ 0.666, 0.055 })
        )
    )
);

/****************************************************
 ****************************************************/

class DelayedCastOracleTest : public ::testing::TestWithParam<bool> {
protected:
    size_t NR = 109, NC = 88;
    std::shared_ptr<tatami::NumericMatrix> cast_dense, cast_sparse, wrapped_cast_dense, wrapped_cast_sparse;

protected:
    void assemble() {
        auto simulated = tatami_test::simulate_sparse_vector<float>(NR * NC, 0.2);
        auto raw_dense = std::shared_ptr<tatami::Matrix<float, size_t> >(new tatami::DenseRowMatrix<float, size_t>(NR, NC, simulated));
        auto raw_sparse = tatami::convert_to_sparse<false, float, size_t>(raw_dense.get()); 

        cast_dense = tatami::make_DelayedCast<double, int>(raw_dense);
        cast_sparse = tatami::make_DelayedCast<double, int>(raw_sparse);

        wrapped_cast_dense = tatami::make_DelayedCast<double, int>(tatami_test::make_CrankyMatrix<float, size_t>(raw_dense));
        wrapped_cast_sparse = tatami::make_DelayedCast<double, int>(tatami_test::make_CrankyMatrix<float, size_t>(raw_sparse));
    }
};

TEST_P(DelayedCastOracleTest, Validate) {
    assemble();
    auto random = GetParam();

    EXPECT_FALSE(cast_dense->uses_oracle(true));
    EXPECT_TRUE(wrapped_cast_dense->uses_oracle(true));

    tatami_test::test_oracle_column_access(wrapped_cast_dense.get(), cast_dense.get(), random);
    tatami_test::test_oracle_column_access(wrapped_cast_sparse.get(), cast_sparse.get(), random);

    tatami_test::test_oracle_row_access(wrapped_cast_dense.get(), cast_dense.get(), random);
    tatami_test::test_oracle_row_access(wrapped_cast_sparse.get(), cast_sparse.get(), random);
}

INSTANTIATE_TEST_CASE_P(
    DelayedCast,
    DelayedCastOracleTest,
    ::testing::Values(true, false)  // use random or consecutive oracle.
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

