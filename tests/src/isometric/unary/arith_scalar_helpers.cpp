#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/base/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOp.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "utils.h"

class ArithScalarUtils { 
protected:
    size_t nrow = 123, ncol = 89;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    std::vector<double> simulated;
protected:
    void assemble() {
        simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1);
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated));
        sparse = tatami::convert_to_sparse<false>(dense.get()); // column major.
        return;
    }
};

template<class PARAM> 
class ArithScalarTest : public ::testing::TestWithParam<PARAM>, public ArithScalarUtils {
protected:    
    void SetUp() {
        assemble();
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
    auto op = tatami::make_DelayedAddScalarHelper(val);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (val) {
        EXPECT_FALSE(sparse_mod->sparse());
    } else {
        EXPECT_TRUE(sparse_mod->sparse());
    }

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto ref = reference(val);

    // Full access.
    tatami_test::test_simple_column_access(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_column_access(sparse_mod.get(), &ref, true, 1); 

    tatami_test::test_simple_row_access(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_row_access(sparse_mod.get(), &ref, true, 1);

    // Block access.
    tatami_test::test_sliced_column_access(dense_mod.get(), &ref, true, 1, nrow * 0.25, nrow * 0.9);
    tatami_test::test_sliced_column_access(sparse_mod.get(), &ref, true, 1, nrow * 0.2, nrow * 0.8); 

    tatami_test::test_sliced_row_access(dense_mod.get(), &ref, true, 1, ncol * 0.1, ncol * 0.5);
    tatami_test::test_sliced_row_access(sparse_mod.get(), &ref, true, 1, ncol * 0.4, ncol * 0.9);

    // Indexed access.
    tatami_test::test_indexed_column_access(dense_mod.get(), &ref, true, 1, nrow * 0.25, 10);
    tatami_test::test_indexed_column_access(sparse_mod.get(), &ref, true, 1, nrow * 0.2, 5);

    tatami_test::test_indexed_row_access(dense_mod.get(), &ref, true, 1, ncol * 0.1, 7);
    tatami_test::test_indexed_row_access(sparse_mod.get(), &ref, true, 1, ncol * 0.4, 11);
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

TEST_P(ArithScalarSubtractionTest, Basic) {
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        auto op = tatami::make_DelayedSubtractScalarHelper<true>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    } else {
        auto op = tatami::make_DelayedSubtractScalarHelper<false>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    }

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (val) {
        EXPECT_FALSE(sparse_mod->sparse());
    } else {
        EXPECT_TRUE(sparse_mod->sparse());
    }

    // Again, doing some light tests.
    auto ref = reference(val, on_right);

    // Full access.
    tatami_test::test_simple_column_access(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_column_access(sparse_mod.get(), &ref, true, 1); 

    tatami_test::test_simple_row_access(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_row_access(sparse_mod.get(), &ref, true, 1);

    // Block access.
    tatami_test::test_sliced_column_access(dense_mod.get(), &ref, true, 1, nrow * 0.25, nrow * 0.9);
    tatami_test::test_sliced_column_access(sparse_mod.get(), &ref, true, 1, nrow * 0.2, nrow * 0.8); 

    tatami_test::test_sliced_row_access(dense_mod.get(), &ref, true, 1, ncol * 0.1, ncol * 0.5);
    tatami_test::test_sliced_row_access(sparse_mod.get(), &ref, true, 1, ncol * 0.4, ncol * 0.9);

    // Indexed access.
    tatami_test::test_indexed_column_access(dense_mod.get(), &ref, true, 1, nrow * 0.25, 10);
    tatami_test::test_indexed_column_access(sparse_mod.get(), &ref, true, 1, nrow * 0.2, 5);

    tatami_test::test_indexed_row_access(dense_mod.get(), &ref, true, 1, ncol * 0.1, 7);
    tatami_test::test_indexed_row_access(sparse_mod.get(), &ref, true, 1, ncol * 0.4, 11);
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

TEST_P(ArithScalarMultiplicationTest, Basic) {
    double val = GetParam();
    auto op = tatami::make_DelayedMultiplyScalarHelper(val);
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());

    // Again, doing some light tests.
    auto ref = reference(val);

    // Full access.
    tatami_test::test_simple_column_access(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_column_access(sparse_mod.get(), &ref, true, 1); 

    tatami_test::test_simple_row_access(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_row_access(sparse_mod.get(), &ref, true, 1);

    // Block access.
    tatami_test::test_sliced_column_access(dense_mod.get(), &ref, true, 1, nrow * 0.25, nrow * 0.9);
    tatami_test::test_sliced_column_access(sparse_mod.get(), &ref, true, 1, nrow * 0.2, nrow * 0.8); 

    tatami_test::test_sliced_row_access(dense_mod.get(), &ref, true, 1, ncol * 0.1, ncol * 0.5);
    tatami_test::test_sliced_row_access(sparse_mod.get(), &ref, true, 1, ncol * 0.4, ncol * 0.9);

    // Indexed access.
    tatami_test::test_indexed_column_access(dense_mod.get(), &ref, true, 1, nrow * 0.25, 10);
    tatami_test::test_indexed_column_access(sparse_mod.get(), &ref, true, 1, nrow * 0.2, 5);

    tatami_test::test_indexed_row_access(dense_mod.get(), &ref, true, 1, ncol * 0.1, 7);
    tatami_test::test_indexed_row_access(sparse_mod.get(), &ref, true, 1, ncol * 0.4, 11);
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
                r = careful_division(r, val);
            } else {
                r = careful_division(val, r);
            }
        }
        return tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec));
    }
};

TEST_P(ArithScalarDivisionTest, Basic) {
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        auto op = tatami::make_DelayedDivideScalarHelper<true>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    } else {
        auto op = tatami::make_DelayedDivideScalarHelper<false>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
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
    auto ref = reference(val, on_right);

    // Full access.
    tatami_test::test_simple_column_access(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_column_access(sparse_mod.get(), &ref, true, 1); 

    tatami_test::test_simple_row_access(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_row_access(sparse_mod.get(), &ref, true, 1);

    // Block access.
    tatami_test::test_sliced_column_access(dense_mod.get(), &ref, true, 1, nrow * 0.25, nrow * 0.9);
    tatami_test::test_sliced_column_access(sparse_mod.get(), &ref, true, 1, nrow * 0.2, nrow * 0.8); 

    tatami_test::test_sliced_row_access(dense_mod.get(), &ref, true, 1, ncol * 0.1, ncol * 0.5);
    tatami_test::test_sliced_row_access(sparse_mod.get(), &ref, true, 1, ncol * 0.4, ncol * 0.9);

    // Indexed access.
    tatami_test::test_indexed_column_access(dense_mod.get(), &ref, true, 1, nrow * 0.25, 10);
    tatami_test::test_indexed_column_access(sparse_mod.get(), &ref, true, 1, nrow * 0.2, 5);

    tatami_test::test_indexed_row_access(dense_mod.get(), &ref, true, 1, ncol * 0.1, 7);
    tatami_test::test_indexed_row_access(sparse_mod.get(), &ref, true, 1, ncol * 0.4, 11);
}

INSTANTIATE_TEST_CASE_P(
    ArithScalar,
    ArithScalarDivisionTest,
    ::testing::Combine(
        ::testing::Values(5, 0.1, -0.7),
        ::testing::Values(true, false)
    )
);

/**************************
 ********* ZEROED *********
 **************************/

class ArithScalarZeroedTest : public ::testing::Test, public ArithScalarUtils {
protected:
    void SetUp() {
        assemble();
    }
};

TEST_F(ArithScalarZeroedTest, Addition) {
    auto op = tatami::make_DelayedAddScalarHelper(0);
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());

    tatami_test::test_simple_column_access(dense_mod.get(), dense.get(), true, 1);
    tatami_test::test_simple_column_access(sparse_mod.get(), sparse.get(), true, 1); 

    tatami_test::test_simple_row_access(dense_mod.get(), dense.get(), true, 1);
    tatami_test::test_simple_row_access(sparse_mod.get(), sparse.get(), true, 1);
}

TEST_F(ArithScalarZeroedTest, Subtraction) {
    {
        auto op = tatami::make_DelayedSubtractScalarHelper<true>(0);

        auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        EXPECT_FALSE(dense_mod->sparse());
        tatami_test::test_simple_column_access(dense_mod.get(), dense.get(), true, 1);
        tatami_test::test_simple_row_access(dense_mod.get(), dense.get(), true, 1);

        auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
        EXPECT_TRUE(sparse_mod->sparse());
        tatami_test::test_simple_row_access(sparse_mod.get(), sparse.get(), true, 1);
        tatami_test::test_simple_column_access(sparse_mod.get(), sparse.get(), true, 1); 
    }

    {
        auto op = tatami::make_DelayedSubtractScalarHelper<false>(0);

        auto copy = simulated;
        for (auto& x : copy) {
            x *= -1;
        }
        tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

        auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        EXPECT_FALSE(dense_mod->sparse());
        tatami_test::test_simple_column_access(dense_mod.get(), &ref, true, 1);
        tatami_test::test_simple_row_access(dense_mod.get(), &ref, true, 1);

        auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
        EXPECT_TRUE(sparse_mod->sparse());
        tatami_test::test_simple_row_access(sparse_mod.get(), &ref, true, 1);
        tatami_test::test_simple_column_access(sparse_mod.get(), &ref, true, 1); 
    }
}

TEST_F(ArithScalarZeroedTest, Multiplication) {
    auto op = tatami::make_DelayedMultiplyScalarHelper(0);
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::vector<double>(nrow * ncol));

    tatami_test::test_simple_column_access(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_column_access(sparse_mod.get(), &ref, true, 1); 

    tatami_test::test_simple_row_access(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_row_access(sparse_mod.get(), &ref, true, 1);
}

TEST_F(ArithScalarZeroedTest, Division) {
    auto op = tatami::make_DelayedDivideScalarHelper<true>(0);

    auto copy = simulated;
    for (auto& x : copy) {
        x = careful_division(x, 0.0);
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    EXPECT_FALSE(dense_mod->sparse());
    test_nan_access(&ref, dense_mod.get());

    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    EXPECT_FALSE(sparse_mod->sparse());
    test_nan_access(&ref, sparse_mod.get());
}
