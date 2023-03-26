#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedCast.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../_tests/test_row_access.h"
#include "../_tests/test_column_access.h"
#include "../_tests/simulate_vector.h"

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
        auto sparse_matrix = simulate_sparse_vector<double>(nrow * ncol, 0.08);
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
    EXPECT_EQ(cast_dense->prefer_rows(), dense->prefer_rows());
    EXPECT_EQ(cast_dense->dimension_preference(), dense->dimension_preference());

    test_simple_row_access(cast_dense.get(), dense.get(), FORWARD, JUMP);
    test_simple_column_access(cast_dense.get(), dense.get(), FORWARD, JUMP);

    auto cast_fdense = tatami::make_DelayedCast<double, int>(fdense);
    test_simple_row_access(cast_fdense.get(), fdense_ref.get(), FORWARD, JUMP);
    test_simple_column_access(cast_fdense.get(), fdense_ref.get(), FORWARD, JUMP);
}

TEST_P(DelayedCastFullAccess, Sparse) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

    auto cast_sparse = tatami::make_DelayedCast<double, int>(sparse);
    test_simple_row_access(cast_sparse.get(), sparse.get(), FORWARD, JUMP);
    test_simple_column_access(cast_sparse.get(), sparse.get(), FORWARD, JUMP);

    auto cast_fsparse = tatami::make_DelayedCast<double, int>(fsparse);
    test_simple_row_access(cast_fsparse.get(), fsparse_ref.get(), FORWARD, JUMP);
    test_simple_column_access(cast_fsparse.get(), fsparse_ref.get(), FORWARD, JUMP);

    auto cast_fsparse_value = tatami::make_DelayedCast<double, int>(fsparse_value);
    test_simple_row_access(cast_fsparse_value.get(), fsparse_ref.get(), FORWARD, JUMP);
    test_simple_column_access(cast_fsparse_value.get(), fsparse_ref.get(), FORWARD, JUMP);

    auto cast_sparse_index = tatami::make_DelayedCast<double, int>(sparse_index);
    test_simple_row_access(cast_sparse_index.get(), sparse.get(), FORWARD, JUMP);
    test_simple_column_access(cast_sparse_index.get(), sparse.get(), FORWARD, JUMP);
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
    test_sliced_row_access(cast_dense.get(), dense.get(), true, JUMP, CFIRST, CLAST);
    test_sliced_column_access(cast_dense.get(), dense.get(), true, JUMP, RFIRST, RLAST);

    auto cast_fdense = tatami::make_DelayedCast<double, int>(fdense);
    test_sliced_row_access(cast_fdense.get(), fdense_ref.get(), true, JUMP, CFIRST, CLAST);
    test_sliced_column_access(cast_fdense.get(), fdense_ref.get(), true, JUMP, RFIRST, RLAST);
}

TEST_P(DelayedCastBlockAccess, Sparse) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t RFIRST = interval_info[0] * nrow, RLAST = interval_info[1] * nrow;
    size_t CFIRST = interval_info[0] * ncol, CLAST = interval_info[1] * ncol;

    auto cast_sparse = tatami::make_DelayedCast<double, int>(sparse);
    test_sliced_row_access(cast_sparse.get(), sparse.get(), true, JUMP, CFIRST, CLAST);
    test_sliced_column_access(cast_sparse.get(), sparse.get(), true, JUMP, RFIRST, RLAST);

    auto cast_fsparse = tatami::make_DelayedCast<double, int>(fsparse);
    test_sliced_row_access(cast_fsparse.get(), fsparse_ref.get(), true, JUMP, CFIRST, CLAST);
    test_sliced_column_access(cast_fsparse.get(), fsparse_ref.get(), true, JUMP, RFIRST, RLAST);

    auto cast_fsparse_value = tatami::make_DelayedCast<double, int>(fsparse_value);
    test_sliced_row_access(cast_fsparse_value.get(), fsparse_ref.get(), true, JUMP, CFIRST, CLAST);
    test_sliced_column_access(cast_fsparse_value.get(), fsparse_ref.get(), true, JUMP, RFIRST, RLAST);

    auto cast_sparse_index = tatami::make_DelayedCast<double, int>(sparse_index);
    test_sliced_row_access(cast_sparse_index.get(), sparse.get(), true, JUMP, CFIRST, CLAST);
    test_sliced_column_access(cast_sparse_index.get(), sparse.get(), true, JUMP, RFIRST, RLAST);
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
    test_indexed_row_access(cast_dense.get(), dense.get(), true, JUMP, CFIRST, CSTEP);
    test_indexed_column_access(cast_dense.get(), dense.get(), true, JUMP, RFIRST, RSTEP);

    auto cast_fdense = tatami::make_DelayedCast<double, int>(fdense);
    test_indexed_row_access(cast_fdense.get(), fdense_ref.get(), true, JUMP, CFIRST, CSTEP);
    test_indexed_column_access(cast_fdense.get(), fdense_ref.get(), true, JUMP, RFIRST, RSTEP);
}

TEST_P(DelayedCastIndexAccess, Sparse) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t RFIRST = interval_info[0] * nrow, RSTEP = interval_info[1] * nrow;
    size_t CFIRST = interval_info[0] * ncol, CSTEP = interval_info[1] * ncol;

    auto cast_sparse = tatami::make_DelayedCast<double, int>(sparse);
    test_indexed_row_access(cast_sparse.get(), sparse.get(), true, JUMP, CFIRST, CSTEP);
    test_indexed_column_access(cast_sparse.get(), sparse.get(), true, JUMP, RFIRST, RSTEP);

    auto cast_fsparse = tatami::make_DelayedCast<double, int>(fsparse);
    test_indexed_row_access(cast_fsparse.get(), fsparse_ref.get(), true, JUMP, CFIRST, CSTEP);
    test_indexed_column_access(cast_fsparse.get(), fsparse_ref.get(), true, JUMP, RFIRST, RSTEP);

    auto cast_fsparse_value = tatami::make_DelayedCast<double, int>(fsparse_value);
    test_indexed_row_access(cast_fsparse_value.get(), fsparse_ref.get(), true, JUMP, CFIRST, CSTEP);
    test_indexed_column_access(cast_fsparse_value.get(), fsparse_ref.get(), true, JUMP, RFIRST, RSTEP);

    auto cast_sparse_index = tatami::make_DelayedCast<double, int>(sparse_index);
    test_indexed_row_access(cast_sparse_index.get(), sparse.get(), true, JUMP, CFIRST, CSTEP);
    test_indexed_column_access(cast_sparse_index.get(), sparse.get(), true, JUMP, RFIRST, RSTEP);
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


