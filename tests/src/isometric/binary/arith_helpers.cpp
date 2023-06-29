#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/binary/DelayedBinaryIsometricOp.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOp.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class BinaryArithUtils {
protected:
    size_t nrow = 91, ncol = 121;
    std::shared_ptr<tatami::NumericMatrix> dense_left, sparse_left, dense_right, sparse_right;
    std::vector<double> simulated_left, simulated_right;
protected:
    void assemble() {
        simulated_left = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1, /* lower = */ -5, /* upper = */ 5, /* seed */ 12345);
        dense_left = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated_left));
        sparse_left = tatami::convert_to_sparse<false>(dense_left.get()); // column major.

        simulated_right = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1, /* lower = */ -5, /* upper = */ 5, /* seed */ 67890);
        dense_right = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated_right));
        sparse_right = tatami::convert_to_sparse<false>(dense_right.get()); // column major.
        return;
    }
};

/****************************
 ********* ADDITION *********
 ****************************/

template<class PARAM>
class BinaryArithAdditionTest : public ::testing::TestWithParam<PARAM>, public BinaryArithUtils {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, mixed_mod, ref;

    void SetUp() {
        assemble();

        auto op = tatami::make_DelayedBinaryAddHelper();
        dense_mod = tatami::make_DelayedBinaryIsometricOp(this->dense_left, this->dense_right, op);
        sparse_mod = tatami::make_DelayedBinaryIsometricOp(this->sparse_left, this->sparse_right, op);
        mixed_mod = tatami::make_DelayedBinaryIsometricOp(this->sparse_left, this->dense_right, op);

        auto refvec = this->simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] += this->simulated_right[i];
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));
    }
};

using BinaryArithAdditionFullTest = BinaryArithAdditionTest<std::tuple<bool, size_t> >;

TEST_P(BinaryArithAdditionFullTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    int JUMP = std::get<1>(param);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense_mod->sparse_proportion(), 0);
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(sparse_mod->sparse_proportion(), 1);
    EXPECT_EQ(nrow, dense_mod->nrow());
    EXPECT_EQ(ncol, dense_mod->ncol());
    EXPECT_EQ(mixed_mod->sparse_proportion(), 0.5);

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_EQ(dense_mod->prefer_rows_proportion(), 1);
    EXPECT_FALSE(sparse_mod->prefer_rows());
    EXPECT_EQ(sparse_mod->prefer_rows_proportion(), 0);
    EXPECT_EQ(mixed_mod->prefer_rows_proportion(), 0.5);

    EXPECT_FALSE(dense_mod->uses_oracle(true));
    EXPECT_FALSE(dense_mod->uses_oracle(false));

    tatami_test::test_simple_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);

    tatami_test::test_simple_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);

    tatami_test::test_simple_column_access(mixed_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(mixed_mod.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithAdditionFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using BinaryArithAdditionBlockTest = BinaryArithAdditionTest<std::tuple<bool, size_t, std::vector<double> > >;

TEST_P(BinaryArithAdditionBlockTest, Basic) {
    auto param = GetParam();

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    {
        size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;
        tatami_test::test_sliced_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }

    {
        size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;
        tatami_test::test_sliced_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithAdditionBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.35 }),
            std::vector<double>({ 0.38, 0.61 }), 
            std::vector<double>({ 0.777, 1 })

        )
    )
);

using BinaryArithAdditionIndexTest = BinaryArithAdditionTest<std::tuple<bool, size_t, std::vector<double> > >;

TEST_P(BinaryArithAdditionIndexTest, Basic) {
    auto param = GetParam();

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    {
        size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1];
        tatami_test::test_indexed_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }

    {
        size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1];
        tatami_test::test_indexed_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithAdditionIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 6 }),
            std::vector<double>({ 0.21, 4 }), 
            std::vector<double>({ 0.56, 2 })
        )
    )
);

/*******************************
 ********* SUBTRACTION *********
 *******************************/

template<class PARAM>
class BinaryArithSubtractionTest : public ::testing::TestWithParam<PARAM>, public BinaryArithUtils {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, mixed_mod, ref;

    void SetUp() {
        assemble();

        auto op = tatami::make_DelayedBinarySubtractHelper();
        dense_mod = tatami::make_DelayedBinaryIsometricOp(this->dense_left, this->dense_right, op);
        sparse_mod = tatami::make_DelayedBinaryIsometricOp(this->sparse_left, this->sparse_right, op);
        mixed_mod = tatami::make_DelayedBinaryIsometricOp(this->sparse_left, this->dense_right, op);

        auto refvec = this->simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] -= this->simulated_right[i];
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));
    }
};

using BinaryArithSubtractionFullTest = BinaryArithSubtractionTest<std::tuple<bool, size_t> >;

TEST_P(BinaryArithSubtractionFullTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    int JUMP = std::get<1>(param);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_FALSE(mixed_mod->sparse());

    tatami_test::test_simple_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(mixed_mod.get(), ref.get(), FORWARD, JUMP);

    tatami_test::test_simple_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(mixed_mod.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithSubtractionFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using BinaryArithSubtractionBlockTest = BinaryArithSubtractionTest<std::tuple<bool, size_t, std::vector<double> > >;

TEST_P(BinaryArithSubtractionBlockTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    {
        size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;
        tatami_test::test_sliced_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }

    {
        size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;
        tatami_test::test_sliced_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithSubtractionBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.541 }), 
            std::vector<double>({ 0.111, 0.999 }), 
            std::vector<double>({ 0.42, 1 })
        )
    )
);

using BinaryArithSubtractionIndexTest = BinaryArithSubtractionTest<std::tuple<bool, size_t, std::vector<double> > >;

TEST_P(BinaryArithSubtractionIndexTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    {
        size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1];
        tatami_test::test_indexed_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }

    {
        size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1];
        tatami_test::test_indexed_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithSubtractionIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 10 }), 
            std::vector<double>({ 0.21, 5 }), 
            std::vector<double>({ 0.56, 2 })
        )
    )
);

/**********************************
 ********* MULTIPLICATION *********
 **********************************/

template<class PARAM>
class BinaryArithMultiplicationTest : public ::testing::TestWithParam<PARAM>, public BinaryArithUtils {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, mixed_mod, ref;

    void SetUp() {
        assemble();

        auto op = tatami::make_DelayedBinaryMultiplyHelper();
        dense_mod = tatami::make_DelayedBinaryIsometricOp(this->dense_left, this->dense_right, op);
        sparse_mod = tatami::make_DelayedBinaryIsometricOp(this->sparse_left, this->sparse_right, op);
        mixed_mod = tatami::make_DelayedBinaryIsometricOp(this->sparse_left, this->dense_right, op);

        auto refvec = this->simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] *= this->simulated_right[i];
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));
    }
};

using BinaryArithMultiplicationFullTest = BinaryArithMultiplicationTest<std::tuple<bool, size_t> >;

TEST_P(BinaryArithMultiplicationFullTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    int JUMP = std::get<1>(param);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());

    tatami_test::test_simple_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(mixed_mod.get(), ref.get(), FORWARD, JUMP);

    tatami_test::test_simple_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(mixed_mod.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithMultiplicationFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using BinaryArithMultiplicationBlockTest = BinaryArithMultiplicationTest<std::tuple<bool, size_t, std::vector<double> > >;

TEST_P(BinaryArithMultiplicationBlockTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    {
        size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;
        tatami_test::test_sliced_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }

    {
        size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;
        tatami_test::test_sliced_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithMultiplicationBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.35 }),
            std::vector<double>({ 0.38, 0.61 }), 
            std::vector<double>({ 0.777, 1 })

        )
    )
);

using BinaryArithMultiplicationIndexTest = BinaryArithMultiplicationTest<std::tuple<bool, size_t, std::vector<double> > >;

TEST_P(BinaryArithMultiplicationIndexTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    {
        size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1];
        tatami_test::test_indexed_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }

    {
        size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1];
        tatami_test::test_indexed_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithMultiplicationIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 11 }),
            std::vector<double>({ 0.21, 7 }), 
            std::vector<double>({ 0.56, 5 })
        )
    )
);

/****************************
 ********* DIVISION *********
 ****************************/

template<class PARAM>
class BinaryArithDivisionTest : public ::testing::TestWithParam<PARAM>, public BinaryArithUtils {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, mixed_mod, ref;

    void SetUp() {
        assemble();

        auto op = tatami::make_DelayedBinaryDivideHelper();
        dense_mod = tatami::make_DelayedBinaryIsometricOp(this->dense_left, this->dense_right, op);
        sparse_mod = tatami::make_DelayedBinaryIsometricOp(this->sparse_left, this->sparse_right, op);
        mixed_mod = tatami::make_DelayedBinaryIsometricOp(this->sparse_left, this->dense_right, op);

        auto refvec = this->simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] = careful_division(refvec[i], this->simulated_right[i]);
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));
    }
};

using BinaryArithDivisionFullTest = BinaryArithDivisionTest<std::tuple<bool, size_t> >;

TEST_P(BinaryArithDivisionFullTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    int JUMP = std::get<1>(param);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());

    tatami_test::test_simple_column_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP);

    tatami_test::test_simple_row_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithDivisionFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using BinaryArithDivisionBlockTest = BinaryArithDivisionTest<std::tuple<bool, size_t, std::vector<double> > >;

TEST_P(BinaryArithDivisionBlockTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    {
        size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;
        tatami_test::test_sliced_column_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }

    {
        size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;
        tatami_test::test_sliced_row_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithDivisionBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.541 }), 
            std::vector<double>({ 0.111, 0.999 }), 
            std::vector<double>({ 0.42, 1 })
        )
    )
);

using BinaryArithDivisionIndexTest = BinaryArithDivisionTest<std::tuple<bool, size_t, std::vector<double> > >;

TEST_P(BinaryArithDivisionIndexTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    {
        size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1];
        tatami_test::test_indexed_column_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }

    {
        size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1];
        tatami_test::test_indexed_row_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithDivisionIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 7 }),
            std::vector<double>({ 0.21, 9 }), 
            std::vector<double>({ 0.56, 10 })
        )
    )
);

/*******************************
 ************ POWER ************
 *******************************/

template<class PARAM>
class BinaryArithPowerTest : public ::testing::TestWithParam<PARAM>, public BinaryArithUtils {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, mixed_mod, ref;

    void SetUp() {
        assemble();

        tatami::DelayedAbsHelper op0;
        auto dense_left = tatami::make_DelayedUnaryIsometricOp(this->dense_left, op0);
        auto sparse_left = tatami::make_DelayedUnaryIsometricOp(this->sparse_left, op0);
        auto dense_right = tatami::make_DelayedUnaryIsometricOp(this->dense_right, op0);
        auto sparse_right = tatami::make_DelayedUnaryIsometricOp(this->sparse_right, op0);

        auto op = tatami::make_DelayedBinaryPowerHelper();
        dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
        sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);
        mixed_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, dense_right, op);

        auto refvec = this->simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] = std::pow(std::abs(refvec[i]), std::abs(this->simulated_right[i]));
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));
    }
};

using BinaryArithPowerFullTest = BinaryArithPowerTest<std::tuple<bool, size_t> >;

TEST_P(BinaryArithPowerFullTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    int JUMP = std::get<1>(param);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_FALSE(mixed_mod->sparse());

    tatami_test::test_simple_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(mixed_mod.get(), ref.get(), FORWARD, JUMP);

    tatami_test::test_simple_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(mixed_mod.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithPowerFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using BinaryArithPowerBlockTest = BinaryArithPowerTest<std::tuple<bool, size_t, std::vector<double> > >;

TEST_P(BinaryArithPowerBlockTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    {
        size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;
        tatami_test::test_sliced_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }

    {
        size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;
        tatami_test::test_sliced_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithPowerBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.35 }),
            std::vector<double>({ 0.38, 0.61 }),
            std::vector<double>({ 0.777, 1 })
        )
    )
);

using BinaryArithPowerIndexTest = BinaryArithPowerTest<std::tuple<bool, size_t, std::vector<double> > >;

TEST_P(BinaryArithPowerIndexTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    {
        size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1];
        tatami_test::test_indexed_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }

    {
        size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1];
        tatami_test::test_indexed_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithPowerIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 11 }),
            std::vector<double>({ 0.21, 7 }),
            std::vector<double>({ 0.56, 5 })
        )
    )
);

/****************************
 ********** MODULO **********
 ****************************/

template<class PARAM>
class BinaryArithModuloTest : public ::testing::TestWithParam<PARAM>, public BinaryArithUtils {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, mixed_mod, ref;

    void SetUp() {
        assemble();

        auto op = tatami::make_DelayedBinaryModuloHelper();
        dense_mod = tatami::make_DelayedBinaryIsometricOp(this->dense_left, this->dense_right, op);
        sparse_mod = tatami::make_DelayedBinaryIsometricOp(this->sparse_left, this->sparse_right, op);
        mixed_mod = tatami::make_DelayedBinaryIsometricOp(this->sparse_left, this->dense_right, op);

        auto refvec = this->simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] = std::fmod(refvec[i], this->simulated_right[i]);
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));
    }
};

using BinaryArithModuloFullTest = BinaryArithModuloTest<std::tuple<bool, size_t> >;

TEST_P(BinaryArithModuloFullTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    int JUMP = std::get<1>(param);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());

    tatami_test::test_simple_column_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP);

    tatami_test::test_simple_row_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithModuloFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using BinaryArithModuloBlockTest = BinaryArithModuloTest<std::tuple<bool, size_t, std::vector<double> > >;

TEST_P(BinaryArithModuloBlockTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    {
        size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;
        tatami_test::test_sliced_column_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }

    {
        size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;
        tatami_test::test_sliced_row_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithModuloBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.541 }),
            std::vector<double>({ 0.111, 0.999 }),
            std::vector<double>({ 0.42, 1 })
        )
    )
);

using BinaryArithModuloIndexTest = BinaryArithModuloTest<std::tuple<bool, size_t, std::vector<double> > >;

TEST_P(BinaryArithModuloIndexTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    {
        size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1];
        tatami_test::test_indexed_column_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }

    {
        size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1];
        tatami_test::test_indexed_row_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithModuloIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 7 }),
            std::vector<double>({ 0.21, 9 }),
            std::vector<double>({ 0.56, 10 })
        )
    )
);

/****************************
 ***** INTEGER DIVISION *****
 ****************************/

template<class PARAM>
class BinaryArithIntegerDivisionTest : public ::testing::TestWithParam<PARAM>, public BinaryArithUtils {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, mixed_mod, ref;

    void SetUp() {
        assemble();

        auto op = tatami::make_DelayedBinaryIntegerDivideHelper();
        dense_mod = tatami::make_DelayedBinaryIsometricOp(this->dense_left, this->dense_right, op);
        sparse_mod = tatami::make_DelayedBinaryIsometricOp(this->sparse_left, this->sparse_right, op);
        mixed_mod = tatami::make_DelayedBinaryIsometricOp(this->sparse_left, this->dense_right, op);

        auto refvec = this->simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            // x == (x %% y) + y * (x %/% y)
            refvec[i] = std::floor(refvec[i] / this->simulated_right[i]);
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));
    }
};

using BinaryArithIntegerDivisionFullTest = BinaryArithIntegerDivisionTest<std::tuple<bool, size_t> >;

TEST_P(BinaryArithIntegerDivisionFullTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    int JUMP = std::get<1>(param);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());

    tatami_test::test_simple_column_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP);

    tatami_test::test_simple_row_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithIntegerDivisionFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using BinaryArithIntegerDivisionBlockTest = BinaryArithIntegerDivisionTest<std::tuple<bool, size_t, std::vector<double> > >;

TEST_P(BinaryArithIntegerDivisionBlockTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    {
        size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;
        tatami_test::test_sliced_column_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }

    {
        size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;
        tatami_test::test_sliced_row_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithIntegerDivisionBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.541 }),
            std::vector<double>({ 0.111, 0.999 }),
            std::vector<double>({ 0.42, 1 })
        )
    )
);

using BinaryArithIntegerDivisionIndexTest = BinaryArithIntegerDivisionTest<std::tuple<bool, size_t, std::vector<double> > >;

TEST_P(BinaryArithIntegerDivisionIndexTest, Basic) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    {
        size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1];
        tatami_test::test_indexed_column_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }

    {
        size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1];
        tatami_test::test_indexed_row_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access<true>(mixed_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithIntegerDivisionIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 7 }),
            std::vector<double>({ 0.21, 9 }),
            std::vector<double>({ 0.56, 10 })
        )
    )
);

/**************************
 ********* ORACLE *********
 **************************/

class BinaryArithOracleTest : public ::testing::TestWithParam<bool>, public BinaryArithUtils {
protected:
    std::shared_ptr<tatami::NumericMatrix> wrapped_both, wrapped_left, wrapped_right, ref;

    void SetUp() {
        assemble();

        auto op = tatami::make_DelayedBinaryAddHelper();
        auto cranky_left = tatami_test::make_CrankyMatrix(this->dense_left);
        auto cranky_right = tatami_test::make_CrankyMatrix(this->dense_right);

        wrapped_both = tatami::make_DelayedBinaryIsometricOp(cranky_left, cranky_right, op);
        wrapped_left = tatami::make_DelayedBinaryIsometricOp(cranky_left, this->dense_right, op);
        wrapped_right = tatami::make_DelayedBinaryIsometricOp(this->dense_left, cranky_right, op);

        auto refvec = this->simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] += this->simulated_right[i];
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));

    }
};

TEST_P(BinaryArithOracleTest, Validate) {
    auto randomized = GetParam();
    EXPECT_TRUE(wrapped_both->uses_oracle(true));
    EXPECT_TRUE(wrapped_left->uses_oracle(true));
    EXPECT_TRUE(wrapped_right->uses_oracle(true));

    tatami_test::test_oracle_column_access(wrapped_both.get(), ref.get(), randomized);
    tatami_test::test_oracle_column_access(wrapped_left.get(), ref.get(), randomized);
    tatami_test::test_oracle_column_access(wrapped_right.get(), ref.get(), randomized);

    tatami_test::test_oracle_row_access(wrapped_both.get(), ref.get(), randomized);
    tatami_test::test_oracle_row_access(wrapped_left.get(), ref.get(), randomized);
    tatami_test::test_oracle_row_access(wrapped_right.get(), ref.get(), randomized);
}

INSTANTIATE_TEST_CASE_P(
    BinaryArith,
    BinaryArithOracleTest,
    ::testing::Values(true, false)  // use random or consecutive oracle.
);

class BinaryArithOracleTest2 : public ::testing::Test, public BinaryArithUtils {};

TEST_F(BinaryArithOracleTest2, Elongated) {
    assemble();

    // Adding some variety in the prediction numbers.
    auto cranky_left = tatami_test::make_CrankyMatrix(this->sparse_left, 23);
    auto cranky_right = tatami_test::make_CrankyMatrix(this->sparse_right, 37);

    auto op = tatami::make_DelayedBinarySubtractHelper();
    auto wrapped_both = tatami::make_DelayedBinaryIsometricOp(cranky_left, cranky_right, op);
    auto ref = tatami::make_DelayedBinaryIsometricOp(this->sparse_left, this->sparse_right, op);

    // Use a very long simulated sequence.
    // This checks that the collection of expired predictions works correctly
    // in the DelayedBinaryIsometricOp::ParentOracle class.

    std::mt19937_64 rng(45987324); 
    std::vector<int> fixed(50000);
    for (auto& x : fixed) {
        x = rng() % this->ncol;
    }

    auto swork = ref->sparse_column();
    auto swork_o = wrapped_both->sparse_column();
    swork_o->set_oracle(std::make_unique<tatami::FixedOracle<int> >(fixed.data(), fixed.size()));

    for (auto i : fixed) {
        auto sexpected = swork->fetch(i);
        auto sobserved = swork_o->fetch(i);
        EXPECT_EQ(sexpected.index, sobserved.index);
        EXPECT_EQ(sexpected.value, sobserved.value);
    }
}

/***********************************
 ********* CONST OVERLOADS *********
 ***********************************/

TEST(BinaryArith, ConstOverload) {
    int nrow = 23, ncol = 42;
    auto simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1);
    auto dense_left = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated));
    auto dense_right = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated));

    auto op = tatami::make_DelayedBinaryAddHelper();
    auto mat = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, std::move(op));

    // cursory checks.
    EXPECT_EQ(mat->nrow(), nrow);
    EXPECT_EQ(mat->ncol(), ncol);
}
