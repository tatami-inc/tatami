#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedCast.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../_tests/test_row_access.h"
#include "../_tests/test_column_access.h"
#include "../data/data.h"
#include "TestCore.h"

template<class PARAM> 
class CastTest : public TestCore<::testing::TestWithParam<PARAM> > {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    std::shared_ptr<tatami::Matrix<float, size_t> > fdense, fsparse;
    std::shared_ptr<tatami::NumericMatrix> fdense_ref, fsparse_ref;
protected:
    void SetUp() {
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
        sparse = tatami::convert_to_sparse<false>(dense.get()); // column-major.

        std::vector<float> fsparse_matrix(sparse_matrix.begin(), sparse_matrix.end());
        fdense = std::shared_ptr<tatami::Matrix<float, size_t> >(new tatami::DenseRowMatrix<float, size_t>(sparse_nrow, sparse_ncol, fsparse_matrix));
        fsparse = tatami::convert_to_sparse<false, float, size_t>(fdense.get()); // column-major.

        std::vector<double> dsparse_matrix(fsparse_matrix.begin(), fsparse_matrix.end());
        fdense_ref = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, dsparse_matrix));
        fsparse_ref = tatami::convert_to_sparse<false>(fdense_ref.get()); // column-major.
        return;
    }
};

/****************************************************
 ****************************************************/

class DelayedCastFullAccess : public CastTest<size_t> {};

TEST_P(DelayedCastFullAccess, Dense) {
    size_t JUMP = GetParam();

    auto cast_dense = tatami::make_DelayedCast<double, int>(dense);
    test_simple_row_access(cast_dense.get(), dense.get(), true, JUMP);
    test_simple_column_access(cast_dense.get(), dense.get(), true, JUMP);

    auto cast_fdense = tatami::make_DelayedCast<double, int>(fdense);
    test_simple_row_access(cast_fdense.get(), fdense_ref.get(), true, JUMP);
    test_simple_column_access(cast_fdense.get(), fdense_ref.get(), true, JUMP);
}

TEST_P(DelayedCastFullAccess, Sparse) {
    size_t JUMP = GetParam();

    auto cast_sparse = tatami::make_DelayedCast<double, int>(sparse);
    test_simple_row_access(cast_sparse.get(), sparse.get(), true, JUMP);
    test_simple_column_access(cast_sparse.get(), sparse.get(), true, JUMP);

    auto cast_fsparse = tatami::make_DelayedCast<double, int>(fsparse);
    test_simple_row_access(cast_fsparse.get(), fsparse_ref.get(), true, JUMP);
    test_simple_column_access(cast_fsparse.get(), fsparse_ref.get(), true, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    DelayedCast,
    DelayedCastFullAccess,
    ::testing::Values(1, 3) // jump, to test the workspace memory.
);

/****************************************************
 ****************************************************/

class DelayedCastSlicedAccess : public CastTest<std::tuple<size_t, std::vector<size_t> > > {};

TEST_P(DelayedCastSlicedAccess, Dense) {
    auto param = GetParam();
    size_t JUMP = std::get<0>(param);
    auto slice = std::get<1>(param);
    size_t FIRST = slice[0], LEN = slice[1], SHIFT = slice[2];

    auto cast_dense = tatami::make_DelayedCast<double, int>(dense);
    test_sliced_row_access(cast_dense.get(), dense.get(), true, JUMP, FIRST, LEN, SHIFT);
    test_sliced_column_access(cast_dense.get(), dense.get(), true, JUMP, FIRST, LEN, SHIFT);

    auto cast_fdense = tatami::make_DelayedCast<double, int>(fdense);
    test_sliced_row_access(cast_fdense.get(), fdense_ref.get(), true, JUMP, FIRST, LEN, SHIFT);
    test_sliced_column_access(cast_fdense.get(), fdense_ref.get(), true, JUMP, FIRST, LEN, SHIFT);
}

TEST_P(DelayedCastSlicedAccess, Sparse) {
    auto param = GetParam();
    size_t JUMP = std::get<0>(param);
    auto slice = std::get<1>(param);
    size_t FIRST = slice[0], LEN = slice[1], SHIFT = slice[2];

    auto cast_sparse = tatami::make_DelayedCast<double, int>(sparse);
    test_sliced_row_access(cast_sparse.get(), sparse.get(), true, JUMP, FIRST, LEN, SHIFT);
    test_sliced_column_access(cast_sparse.get(), sparse.get(), true, JUMP, FIRST, LEN, SHIFT);

    auto cast_fsparse = tatami::make_DelayedCast<double, int>(fsparse);
    test_sliced_row_access(cast_fsparse.get(), fsparse_ref.get(), true, JUMP, FIRST, LEN, SHIFT);
    test_sliced_column_access(cast_fsparse.get(), fsparse_ref.get(), true, JUMP, FIRST, LEN, SHIFT);
}

INSTANTIATE_TEST_CASE_P(
    DelayedCast,
    DelayedCastSlicedAccess,
    ::testing::Combine(
        ::testing::Values(1, 3), // jump, to check the workspace memory
        ::testing::Values(
            std::vector<size_t>({ 0, 6, 1 }), // start, length, shift
            std::vector<size_t>({ 5, 5, 2 }),
            std::vector<size_t>({ 3, 7, 0 })
        )
    )
);

