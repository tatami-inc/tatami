#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOperation.hpp"
#include "tatami/isometric/unary/compare_helpers.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class DelayedUnaryIsometricCompareVectorTest : public ::testing::TestWithParam<std::tuple<bool, bool> > {
protected:
    inline static size_t nrow = 291, ncol = 188;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    inline static std::vector<double> simulated;

    static void SetUpTestSuite() {
        simulated = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.1;
            opt.lower = -3;
            opt.upper = 3;
            opt.seed = 128376111;
            return opt;
        }());

        for (auto& x : simulated) {
            if (x) {
                // Rounding for easier tests of exact equality.
                x = std::round(x);
                if (x == 0) {
                    x = 1;
                }
            }
        }

        dense.reset(new tatami::DenseMatrix<double, int, decltype(simulated)>(nrow, ncol, simulated, true)); // row major.
        sparse = tatami::convert_to_compressed_sparse<double, int>(*dense, false, {}); // column major.
    }

    static void fill_default_vector(std::vector<double>& vec) {
        int val = 1;
        for (auto& x : vec) {
            x = (val % 3) - 1;
            ++val;
        }
    }
};

TEST_P(DelayedUnaryIsometricCompareVectorTest, Equal) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<double> vec(row ? nrow : ncol);
    if (is_sparse) {
        int val = 1;
        for (auto& x : vec) {
            // i.e., no zero equality here.
            x = (val % 2 ? 1 : -1);
            ++val;
        }
    } else {
        fill_default_vector(vec);
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricEqualVectorHelper<double, double, int, decltype(vec)> >(vec, row);
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
            x = (x == vec[row ? r : c]);
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareVectorTest, GreaterThan) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<double> vec(row ? nrow : ncol);
    if (is_sparse) {
        int val = 1;
        for (auto& x : vec) {
            x = (val % 3);
            ++val;
        }
    } else {
        fill_default_vector(vec);
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricGreaterThanVectorHelper<double, double, int, decltype(vec)> >(vec, row);
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
            x = (x > vec[row ? r : c]);
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareVectorTest, LessThan) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<double> vec(row ? nrow : ncol);
    if (is_sparse) {
        int val = 1;
        for (auto& x : vec) {
            x = -(val % 3);
            ++val;
        }
    } else {
        fill_default_vector(vec);
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricLessThanVectorHelper<double, double, int, decltype(vec)> >(vec, row);
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
            x = (x < vec[row ? r : c]);
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareVectorTest, GreaterThanOrEqual) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<double> vec(row ? nrow : ncol);
    if (is_sparse) {
        int val = 1;
        for (auto& x : vec) {
            x = (val % 2 ? 1.0 : 1.5);
            ++val;
        }
    } else {
        fill_default_vector(vec);
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricGreaterThanOrEqualVectorHelper<double, double, int, decltype(vec)> >(vec, row);
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
            x = (x >= vec[row ? r : c]);
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareVectorTest, LessThanOrEqual) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<double> vec(row ? nrow : ncol);
    if (is_sparse) {
        int val = 1;
        for (auto& x : vec) {
            x = (val % 2 ? -1.0 : -2.0);
            ++val;
        }
    } else {
        fill_default_vector(vec);
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricLessThanOrEqualVectorHelper<double, double, int, decltype(vec)> >(vec, row);
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
            x = (x <= vec[row ? r : c]);
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareVectorTest, NotEqual) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<double> vec(row ? nrow : ncol);
    if (!is_sparse) {
        fill_default_vector(vec);
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricNotEqualVectorHelper<double, double, int, decltype(vec)> >(vec, row);
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
            x = (x != vec[row ? r : c]);
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareVectorTest, NewType) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<double> vec(row ? nrow : ncol);
    if (is_sparse) {
        int val = 1;
        for (auto& x : vec) {
            // i.e., no zero equality here.
            x = (val % 2 ? 1 : -1);
            ++val;
        }
    } else {
        fill_default_vector(vec);
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricEqualVectorHelper<uint8_t, double, int, decltype(vec)> >(vec, row);
    tatami::DelayedUnaryIsometricOperation<uint8_t, double, int> dense_umod(dense, op);
    tatami::DelayedUnaryIsometricOperation<uint8_t, double, int> sparse_umod(sparse, op);

    EXPECT_FALSE(dense_umod.is_sparse());
    if (is_sparse) {
        EXPECT_TRUE(sparse_umod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_umod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<uint8_t> urefvec(simulated.size());
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            size_t offset = r * ncol + c;
            urefvec[offset] = (simulated[offset] == vec[row ? r : c]);
        }
    }

    tatami::DenseMatrix<uint8_t, int, decltype(urefvec)> uref(nrow, ncol, std::move(urefvec), true);
    quick_test_all<uint8_t, int>(dense_umod, uref);
    quick_test_all<uint8_t, int>(sparse_umod, uref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricCompareVector,
    DelayedUnaryIsometricCompareVectorTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or by column
        ::testing::Values(true, false) // check sparse case
    )
);
