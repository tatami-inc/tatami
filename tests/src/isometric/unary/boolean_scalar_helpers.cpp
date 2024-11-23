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

        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulated));
        sparse = tatami::convert_to_compressed_sparse<false, double, int>(dense.get()); // column major.
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
    auto op = tatami::make_DelayedUnaryIsometricBooleanAndScalar(other);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense_mod->nrow(), nrow);
    EXPECT_EQ(dense_mod->ncol(), ncol);
    EXPECT_TRUE(sparse_mod->is_sparse());

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r && other;
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricBooleanScalarTest, OR) {
    bool other = GetParam();
    auto op = tatami::make_DelayedUnaryIsometricBooleanOrScalar(other);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (other) {
        EXPECT_FALSE(sparse_mod->is_sparse());
    } else {
        EXPECT_TRUE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r || other;
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricBooleanScalarTest, XOR) {
    bool other = GetParam();
    auto op = tatami::make_DelayedUnaryIsometricBooleanXorScalar(other);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (other) {
        EXPECT_FALSE(sparse_mod->is_sparse());
    } else {
        EXPECT_TRUE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = static_cast<bool>(r) != other;
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricBooleanScalarTest, EQUAL) {
    bool other = GetParam();
    auto op = tatami::make_DelayedUnaryIsometricBooleanEqualScalar(other);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (other) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = static_cast<bool>(r) == other;
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricBooleanScalarTest, NewType) {
    bool other = GetParam();

    auto op = tatami::make_DelayedUnaryIsometricBooleanAndScalar(other);
    auto dense_umod = tatami::make_DelayedUnaryIsometricOperation<uint8_t>(dense, op);
    auto sparse_umod = tatami::make_DelayedUnaryIsometricOperation<uint8_t>(sparse, op);

    EXPECT_FALSE(dense_umod->is_sparse());
    EXPECT_TRUE(sparse_umod->is_sparse());

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    std::vector<uint8_t> urefvec(simulated.size());
    for (size_t i = 0; i < simulated.size(); ++i) {
        urefvec[i] = (simulated[i] && other);
    }
    tatami::DenseRowMatrix<uint8_t, int> uref(nrow, ncol, std::move(urefvec));

    quick_test_all<uint8_t, int>(*dense_umod, uref);
    quick_test_all<uint8_t, int>(*sparse_umod, uref);
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
    auto op = tatami::make_DelayedUnaryIsometricBooleanNot();
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense_mod->nrow(), nrow);
    EXPECT_EQ(dense_mod->ncol(), ncol);
    EXPECT_FALSE(sparse_mod->is_sparse());

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = !static_cast<bool>(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_F(DelayedUnaryIsometricBooleanNotTest, NewType) {
    auto op = tatami::make_DelayedUnaryIsometricBooleanNot();
    auto dense_umod = tatami::make_DelayedUnaryIsometricOperation<uint8_t>(dense, op);
    auto sparse_umod = tatami::make_DelayedUnaryIsometricOperation<uint8_t>(sparse, op);

    EXPECT_FALSE(dense_umod->is_sparse());
    EXPECT_FALSE(sparse_umod->is_sparse());

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    std::vector<uint8_t> urefvec(simulated.size());
    for (size_t i = 0; i < simulated.size(); ++i) {
        urefvec[i] = !(simulated[i]);
    }
    tatami::DenseRowMatrix<uint8_t, int> uref(nrow, ncol, std::move(urefvec));

    quick_test_all<uint8_t, int>(*dense_umod, uref);
    quick_test_all<uint8_t, int>(*sparse_umod, uref);
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
    auto op = tatami::DelayedUnaryIsometricBooleanCast();
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense_mod->nrow(), nrow);
    EXPECT_EQ(dense_mod->ncol(), ncol);
    EXPECT_TRUE(sparse_mod->is_sparse());

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = static_cast<bool>(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_F(DelayedUnaryIsometricBooleanCastTest, NewType) {
    auto op = tatami::DelayedUnaryIsometricBooleanCast();
    auto dense_umod = tatami::make_DelayedUnaryIsometricOperation<uint8_t>(dense, op);
    auto sparse_umod = tatami::make_DelayedUnaryIsometricOperation<uint8_t>(sparse, op);

    EXPECT_FALSE(dense_umod->is_sparse());
    EXPECT_TRUE(sparse_umod->is_sparse());

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    std::vector<uint8_t> urefvec(simulated.size());
    for (size_t i = 0; i < simulated.size(); ++i) {
        urefvec[i] = static_cast<bool>(simulated[i]);
    }
    tatami::DenseRowMatrix<uint8_t, int> uref(nrow, ncol, std::move(urefvec));

    quick_test_all<uint8_t, int>(*dense_umod, uref);
    quick_test_all<uint8_t, int>(*sparse_umod, uref);
}
