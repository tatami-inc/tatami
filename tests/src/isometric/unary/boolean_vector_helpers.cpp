#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOperation.hpp"
#include "tatami/isometric/unary/boolean_helpers.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class DelayedUnaryIsometricBooleanVectorTest : public ::testing::TestWithParam<std::tuple<bool, bool> > {
protected:
    inline static size_t nrow = 191, ncol = 88;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    inline static std::vector<double> simulated;

    static void SetUpTestSuite() {
        simulated = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.1;
            opt.lower = -3;
            opt.upper = 3;
            opt.seed = 817236584;
            return opt;
        }());

        dense.reset(new tatami::DenseMatrix<double, int, decltype(simulated)>(nrow, ncol, simulated, true)); // row major.
        sparse = tatami::convert_to_compressed_sparse<false, double, int>(dense.get()); // column major.
    }

    static void fill_default_vector(std::vector<char>& vec) {
        bool val = true;
        for (auto& x : vec) {
            x = val;
            val = !val;
        }
    }
};

TEST_P(DelayedUnaryIsometricBooleanVectorTest, AND) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<char> vec(row ? nrow : ncol);
    if (is_sparse) {
        std::fill(vec.begin(), vec.end(), 1);
    } else {
        fill_default_vector(vec);
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricBooleanAndVectorHelper<double, double, int, decltype(vec)> >(vec, row);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());
    EXPECT_TRUE(sparse_mod.is_sparse());

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = simulated;
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            auto& x = refvec[r * ncol + c];
            x = (x && vec[row ? r : c]);
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricBooleanVectorTest, OR) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<char> vec(row ? nrow : ncol);
    if (!is_sparse) {
        fill_default_vector(vec);
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricBooleanOrVectorHelper<double, double, int, decltype(vec)> >(vec, row);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());
    if (is_sparse) {
        EXPECT_TRUE(sparse_mod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = simulated;
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            auto& x = refvec[r * ncol + c];
            x = (x || vec[row ? r : c]);
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricBooleanVectorTest, XOR) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<char> vec(row ? nrow : ncol);
    if (!is_sparse) {
        fill_default_vector(vec);
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricBooleanXorVectorHelper<double, double, int, decltype(vec)> >(vec, row);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());
    if (is_sparse) {
        EXPECT_TRUE(sparse_mod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = simulated;
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            auto& x = refvec[r * ncol + c];
            x = (static_cast<bool>(x) != static_cast<bool>(vec[row ? r : c]));
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricBooleanVectorTest, EQUAL) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<char> vec(row ? nrow : ncol);
    if (is_sparse) {
        std::fill(vec.begin(), vec.end(), 1);
    } else {
        fill_default_vector(vec);
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricBooleanEqualVectorHelper<double, double, int, decltype(vec)> >(vec, row);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());
    if (is_sparse) {
        EXPECT_TRUE(sparse_mod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = simulated;
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            auto& x = refvec[r * ncol + c];
            x = (static_cast<bool>(x) == static_cast<bool>(vec[row ? r : c]));
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricBooleanVectorTest, NewType) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<char> vec(row ? nrow : ncol);
    if (is_sparse) {
        std::fill(vec.begin(), vec.end(), 1);
    } else {
        fill_default_vector(vec);
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricBooleanAndVectorHelper<uint8_t, double, int, decltype(vec)> >(vec, row);
    tatami::DelayedUnaryIsometricOperation<uint8_t, double, int> dense_umod(dense, op);
    tatami::DelayedUnaryIsometricOperation<uint8_t, double, int> sparse_umod(sparse, op);

    EXPECT_FALSE(dense_umod.is_sparse());
    EXPECT_TRUE(sparse_umod.is_sparse());

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<uint8_t> urefvec(simulated.size());
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            size_t offset = r * ncol + c;
            urefvec[offset] = (simulated[offset] && vec[row ? r : c]);
        }
    }
    
    tatami::DenseMatrix<uint8_t, int, decltype(urefvec)> uref(nrow, ncol, std::move(urefvec), true);
    quick_test_all<uint8_t, int>(dense_umod, uref);
    quick_test_all<uint8_t, int>(sparse_umod, uref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricBooleanVector,
    DelayedUnaryIsometricBooleanVectorTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or by column
        ::testing::Values(true, false) // check sparse case
    )
);
