#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOperation.hpp"
#include "tatami/isometric/unary/boolean_helpers.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class DelayedUnaryIsometricBooleanScalarUtils {
protected:
    inline static size_t nrow = 83, ncol = 111;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    inline static std::vector<double> simulated;

    static void assemble() {
        if (dense) {
            return;
        }

        simulated = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.1;
            opt.lower = -2;
            opt.upper = 2;
            opt.seed = 387126837;
            return opt;
        }());

        dense.reset(new tatami::DenseMatrix<double, int, decltype(simulated)>(nrow, ncol, simulated, true)); // row major.
        sparse = tatami::convert_to_compressed_sparse<double, int>(*dense, false, {}); // column major.
    }
};

class DelayedUnaryIsometricBooleanScalarTest : public ::testing::TestWithParam<bool>, public DelayedUnaryIsometricBooleanScalarUtils { 
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedUnaryIsometricBooleanScalarTest, AND) {
    bool other = GetParam();
    auto op = std::make_shared<tatami::DelayedUnaryIsometricBooleanAndScalarHelper<double, double, int> >(other);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);
    EXPECT_TRUE(sparse_mod.is_sparse());

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r && other;
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricBooleanScalarTest, OR) {
    bool other = GetParam();
    auto op = std::make_shared<tatami::DelayedUnaryIsometricBooleanOrScalarHelper<double, double, int> >(other);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);

    if (other) {
        EXPECT_FALSE(sparse_mod.is_sparse());
    } else {
        EXPECT_TRUE(sparse_mod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r || other;
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricBooleanScalarTest, XOR) {
    bool other = GetParam();
    auto op = std::make_shared<tatami::DelayedUnaryIsometricBooleanXorScalarHelper<double, double, int> >(other);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);

    if (other) {
        EXPECT_FALSE(sparse_mod.is_sparse());
    } else {
        EXPECT_TRUE(sparse_mod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = static_cast<bool>(r) != other;
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricBooleanScalarTest, EQUAL) {
    bool other = GetParam();
    auto op = std::make_shared<tatami::DelayedUnaryIsometricBooleanEqualScalarHelper<double, double, int> >(other);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);

    if (other) {
        EXPECT_TRUE(sparse_mod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = static_cast<bool>(r) == other;
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricBooleanScalarTest, NewType) {
    bool other = GetParam();
    auto op = std::make_shared<tatami::DelayedUnaryIsometricBooleanAndScalarHelper<uint8_t, double, int> >(other);
    tatami::DelayedUnaryIsometricOperation<uint8_t, double, int> dense_umod(dense, op);
    tatami::DelayedUnaryIsometricOperation<uint8_t, double, int> sparse_umod(sparse, op);

    EXPECT_FALSE(dense_umod.is_sparse());
    EXPECT_TRUE(sparse_umod.is_sparse());

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    std::vector<uint8_t> urefvec(simulated.size());
    for (size_t i = 0; i < simulated.size(); ++i) {
        urefvec[i] = (simulated[i] && other);
    }
    tatami::DenseMatrix<uint8_t, int, decltype(urefvec)> uref(nrow, ncol, std::move(urefvec), true);

    quick_test_all<uint8_t, int>(dense_umod, uref);
    quick_test_all<uint8_t, int>(sparse_umod, uref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricBooleanScalar,
    DelayedUnaryIsometricBooleanScalarTest,
    ::testing::Values(true, false)
);

/*******************************************
 *******************************************/

class DelayedUnaryIsometricBooleanNotTest : public ::testing::Test, public DelayedUnaryIsometricBooleanScalarUtils { 
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_F(DelayedUnaryIsometricBooleanNotTest, Basic) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricBooleanNotHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);
    EXPECT_FALSE(sparse_mod.is_sparse());

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = !static_cast<bool>(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // For coverage purposes.
    EXPECT_FALSE(op->zero_depends_on_row());
    EXPECT_FALSE(op->non_zero_depends_on_row());
    EXPECT_FALSE(op->zero_depends_on_column());
    EXPECT_FALSE(op->non_zero_depends_on_column());
}

TEST_F(DelayedUnaryIsometricBooleanNotTest, NewType) {
    auto uop = std::make_shared<tatami::DelayedUnaryIsometricBooleanNotHelper<uint8_t, double, int> >();
    tatami::DelayedUnaryIsometricOperation<uint8_t, double, int> dense_umod(dense, uop);
    tatami::DelayedUnaryIsometricOperation<uint8_t, double, int> sparse_umod(sparse, uop);

    EXPECT_FALSE(dense_umod.is_sparse());
    EXPECT_FALSE(sparse_umod.is_sparse());

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    std::vector<uint8_t> urefvec(simulated.size());
    for (size_t i = 0; i < simulated.size(); ++i) {
        urefvec[i] = !(simulated[i]);
    }
    tatami::DenseMatrix<uint8_t, int, decltype(urefvec)> uref(nrow, ncol, std::move(urefvec), true);

    quick_test_all<uint8_t, int>(dense_umod, uref);
    quick_test_all<uint8_t, int>(sparse_umod, uref);
}

/*******************************************
 *******************************************/

class DelayedUnaryIsometricBooleanCastTest : public ::testing::Test, public DelayedUnaryIsometricBooleanScalarUtils { 
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_F(DelayedUnaryIsometricBooleanCastTest, Basic) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricBooleanCastHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);
    EXPECT_TRUE(sparse_mod.is_sparse());

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = static_cast<bool>(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // For coverage purposes.
    EXPECT_FALSE(op->zero_depends_on_row());
    EXPECT_FALSE(op->non_zero_depends_on_row());
    EXPECT_FALSE(op->zero_depends_on_column());
    EXPECT_FALSE(op->non_zero_depends_on_column());
}

TEST_F(DelayedUnaryIsometricBooleanCastTest, NewType) {
    auto uop = std::make_shared<tatami::DelayedUnaryIsometricBooleanCastHelper<uint8_t, double, int> >();
    tatami::DelayedUnaryIsometricOperation<uint8_t, double, int> dense_umod(dense, uop);
    tatami::DelayedUnaryIsometricOperation<uint8_t, double, int> sparse_umod(sparse, uop);

    EXPECT_FALSE(dense_umod.is_sparse());
    EXPECT_TRUE(sparse_umod.is_sparse());

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    std::vector<uint8_t> urefvec(simulated.size());
    for (size_t i = 0; i < simulated.size(); ++i) {
        urefvec[i] = static_cast<bool>(simulated[i]);
    }
    tatami::DenseMatrix<uint8_t, int, decltype(urefvec)> uref(nrow, ncol, std::move(urefvec), true);

    quick_test_all<uint8_t, int>(dense_umod, uref);
    quick_test_all<uint8_t, int>(sparse_umod, uref);
}

TEST(DelayedUnaryIsometricBooleanScalar, BackCompatibility) {
    auto eqop = tatami::make_DelayedUnaryIsometricBooleanEqualScalar(false);
    EXPECT_FALSE(eqop->is_sparse());
    auto andop = tatami::make_DelayedUnaryIsometricBooleanAndScalar(false);
    EXPECT_TRUE(andop->is_sparse());
    auto orop = tatami::make_DelayedUnaryIsometricBooleanOrScalar(false);
    EXPECT_TRUE(orop->is_sparse());
    auto xorop = tatami::make_DelayedUnaryIsometricBooleanXorScalar(false);
    EXPECT_TRUE(xorop->is_sparse());
    tatami::DelayedUnaryIsometricBooleanNot notop;
    EXPECT_FALSE(notop.is_sparse());
    tatami::DelayedUnaryIsometricBooleanCast castop;
    EXPECT_TRUE(castop.is_sparse());
}
