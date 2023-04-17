#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/base/dense/DenseMatrix.hpp"
#include "tatami/base/isometric/DelayedIsometricOp.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../_tests/test_column_access.h"
#include "../_tests/test_row_access.h"
#include "../_tests/simulate_vector.h"

template<class PARAM> 
class ArithScalarTest : public ::testing::TestWithParam<PARAM> {
protected:
    size_t nrow = 123, ncol = 89;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    std::vector<double> simulated;
protected:
    void SetUp() {
        simulated = simulate_sparse_vector<double>(nrow * ncol, 0.1);
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated));
        sparse = tatami::convert_to_sparse<false>(dense.get()); // column major.
        return;
    }
};

/****************************
 ********* ADDITION *********
 ****************************/

class ArithScalarAdditionTest : public ArithScalarTest<double> {
protected:
    tatami::DenseRowMatrix<double> reference(double val) {
        auto refvec = simulated;
        for (auto& r : refvec) {
            r += val;
        }
        return tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec));
    }
};

TEST_P(ArithScalarAdditionTest, Basic) {
    double val = GetParam();
    tatami::DelayedAddScalarHelper<double> op(val);

    auto dense_mod = tatami::make_DelayedIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    bool FORWARD = true;
    int JUMP = 1;
    auto ref = reference(val);

    test_simple_column_access(dense_mod.get(), &ref, FORWARD, JUMP);
    test_simple_column_access(sparse_mod.get(), &ref, FORWARD, JUMP);

    test_simple_row_access(dense_mod.get(), &ref, FORWARD, JUMP);
    test_simple_row_access(sparse_mod.get(), &ref, FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    ArithScalar,
    ArithScalarAdditionTest,
    ::testing::Values(5, 0.1, -0.7)
);

/*******************************
 ********* SUBTRACTION *********
 *******************************/

class ArithScalarSubtractionTest : public ArithScalarTest<std::tuple<double, bool> > {
protected:
    tatami::DenseRowMatrix<double> reference(double val, bool on_right) {
        auto refvec = simulated;
        for (auto& r : refvec) {
            if (on_right) {
                r -= val;
            } else {
                r = val - r;
            }
        }
        return tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec));
    }
};

TEST_P(ArithScalarSubtractionTest, ColumnAccess) {
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        tatami::DelayedSubtractScalarHelper<true> op(val);
        dense_mod = tatami::make_DelayedIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);
    } else {
        tatami::DelayedSubtractScalarHelper<false> op(val);
        dense_mod = tatami::make_DelayedIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);
    }

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    // Again, doing some light tests.
    bool FORWARD = true;
    int JUMP = 1;
    auto ref = reference(val, on_right);

    test_simple_column_access(dense_mod.get(), &ref, FORWARD, JUMP);
    test_simple_column_access(sparse_mod.get(), &ref, FORWARD, JUMP);

    test_simple_row_access(dense_mod.get(), &ref, FORWARD, JUMP);
    test_simple_row_access(sparse_mod.get(), &ref, FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    ArithScalar,
    ArithScalarSubtractionTest,
    ::testing::Combine(
        ::testing::Values(5, 0.1, -0.7),
        ::testing::Values(true, false)
    )
);

/**********************************
 ********* MULTIPLICATION *********
 **********************************/

class ArithScalarMultiplicationTest : public ArithScalarTest<double> {
protected:
    tatami::DenseRowMatrix<double> reference(double val) {
        auto refvec = simulated;
        for (auto& r : refvec) {
            r *= val;
        }
        return tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec));
    }
};

TEST_P(ArithScalarMultiplicationTest, ColumnAccess) {
    double val = GetParam();
    tatami::DelayedMultiplyScalarHelper<double> op(val);
    auto dense_mod = tatami::make_DelayedIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);

    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());

    // Again, doing some light tests.
    bool FORWARD = true;
    int JUMP = 1;
    auto ref = reference(val);

    test_simple_column_access(dense_mod.get(), &ref, FORWARD, JUMP);
    test_simple_column_access(sparse_mod.get(), &ref, FORWARD, JUMP);

    test_simple_row_access(dense_mod.get(), &ref, FORWARD, JUMP);
    test_simple_row_access(sparse_mod.get(), &ref, FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    ArithScalar,
    ArithScalarMultiplicationTest,
    ::testing::Values(5, 0.1, -0.7)
);

/****************************
 ********* DIVISION *********
 ****************************/

class ArithScalarDivisionTest : public ArithScalarTest<std::tuple<double, bool> > {
protected:
    tatami::DenseRowMatrix<double> reference(double val, bool on_right) {
        auto refvec = simulated;
        for (auto& r : refvec) {
            if (on_right) {
                r /= val;
            } else {
                if (r) {
                    r = val / r;
                } else {
                    r = std::numeric_limits<double>::infinity();
                }
            }
        }
        return tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec));
    }
};

TEST_P(ArithScalarDivisionTest, ColumnAccess) {
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        tatami::DelayedDivideScalarHelper<true> op(val);
        dense_mod = tatami::make_DelayedIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);
    } else {
        tatami::DelayedDivideScalarHelper<false> op(val);
        dense_mod = tatami::make_DelayedIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);
    }

    EXPECT_FALSE(dense_mod->sparse());
    if (on_right) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    // Again, doing some light tests.
    bool FORWARD = true;
    int JUMP = 1;
    auto ref = reference(val, on_right);

    test_simple_column_access(dense_mod.get(), &ref, FORWARD, JUMP);
    test_simple_column_access(sparse_mod.get(), &ref, FORWARD, JUMP);

    test_simple_row_access(dense_mod.get(), &ref, FORWARD, JUMP);
    test_simple_row_access(sparse_mod.get(), &ref, FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    ArithScalar,
    ArithScalarDivisionTest,
    ::testing::Combine(
        ::testing::Values(5, 0.1, -0.7),
        ::testing::Values(true, false)
    )
);
