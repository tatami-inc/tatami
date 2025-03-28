#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOperation.hpp"
#include "tatami/isometric/unary/substitute_helpers.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class DelayedUnaryIsometricSubstituteVectorTest : public ::testing::TestWithParam<std::tuple<bool, bool> > {
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
            opt.seed = 1237961238;
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

    static std::vector<double> create_comparison_vector(size_t len) {
        std::vector<double> vec(len);
        int val = 1;
        for (auto& x : vec) {
            x = (val % 3) - 1;
            ++val;
        }
        return vec;
    }

    static std::vector<double> create_replacement_vector(size_t len, uint64_t seed) {
        return tatami_test::simulate_vector<double>(len, [&]{
            tatami_test::SimulateVectorOptions opt;
            opt.lower = -10;
            opt.upper = 10;
            opt.seed = seed;
            return opt;
        }());
    }
};

TEST_P(DelayedUnaryIsometricSubstituteVectorTest, Equal) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    auto comp = create_comparison_vector(row ? nrow : ncol);
    auto sub = create_replacement_vector(comp.size(), static_cast<int>(row) * 1000 + static_cast<int>(is_sparse) * 13);
    if (is_sparse) {
        for (size_t i = 0, end = comp.size(); i < end; ++i){
            if (comp[i] == 0) {
                sub[i] = 0;
            }
        }
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricSubstituteEqualVectorHelper<double, double, int, decltype(comp), decltype(sub) > >(comp, sub, row);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());
    EXPECT_EQ(sparse_mod.is_sparse(), is_sparse);

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = simulated;
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            auto& x = refvec[r * ncol + c];
            auto i = (row ? r : c);
            x = (x == comp[i] ? sub[i] : x);
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSubstituteVectorTest, GreaterThan) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    auto comp = create_comparison_vector(row ? nrow : ncol);
    auto sub = create_replacement_vector(comp.size(), static_cast<int>(row) * 1000 + static_cast<int>(is_sparse) * 17);
    if (is_sparse) {
        for (size_t i = 0, end = comp.size(); i < end; ++i){
            if (comp[i] < 0) {
                sub[i] = 0;
            }
        }
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricSubstituteGreaterThanVectorHelper<double, double, int, decltype(comp), decltype(sub)> >(comp, sub, row);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());
    EXPECT_EQ(sparse_mod.is_sparse(), is_sparse);

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = simulated;
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            auto& x = refvec[r * ncol + c];
            auto i = (row ? r : c);
            x = (x > comp[i] ? sub[i] : x);
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSubstituteVectorTest, LessThan) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    auto comp = create_comparison_vector(row ? nrow : ncol);
    auto sub = create_replacement_vector(comp.size(), static_cast<int>(row) * 1000 + static_cast<int>(is_sparse) * 31);
    if (is_sparse) {
        for (size_t i = 0, end = comp.size(); i < end; ++i){
            if (comp[i] > 0) {
                sub[i] = 0;
            }
        }
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricSubstituteLessThanVectorHelper<double, double, int, decltype(comp), decltype(sub)> >(comp, sub, row);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());
    EXPECT_EQ(sparse_mod.is_sparse(), is_sparse);

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = simulated;
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            auto& x = refvec[r * ncol + c];
            auto i = (row ? r : c);
            x = (x < comp[i] ? sub[i] : x);
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSubstituteVectorTest, GreaterThanOrEqual) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    auto comp = create_comparison_vector(row ? nrow : ncol);
    auto sub = create_replacement_vector(comp.size(), static_cast<int>(row) * 1000 + static_cast<int>(is_sparse) * 19);
    if (is_sparse) {
        for (size_t i = 0, end = comp.size(); i < end; ++i){
            if (comp[i] <= 0) {
                sub[i] = 0;
            }
        }
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricSubstituteGreaterThanOrEqualVectorHelper<double, double, int, decltype(comp), decltype(sub)> >(comp, sub, row);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());
    EXPECT_EQ(sparse_mod.is_sparse(), is_sparse);

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = simulated;
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            auto& x = refvec[r * ncol + c];
            auto i = (row ? r : c);
            x = (x >= comp[i] ? sub[i] : x);
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSubstituteVectorTest, LessThanOrEqual) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    auto comp = create_comparison_vector(row ? nrow : ncol);
    auto sub = create_replacement_vector(comp.size(), static_cast<int>(row) * 1000 + static_cast<int>(is_sparse) * 23);
    if (is_sparse) {
        for (size_t i = 0, end = comp.size(); i < end; ++i){
            if (comp[i] >= 0) {
                sub[i] = 0;
            }
        }
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricSubstituteLessThanOrEqualVectorHelper<double, double, int, decltype(comp), decltype(sub)> >(comp, sub, row);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());
    EXPECT_EQ(sparse_mod.is_sparse(), is_sparse);

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = simulated;
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            auto& x = refvec[r * ncol + c];
            auto i = (row ? r : c);
            x = (x <= comp[i] ? sub[i] : x);
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSubstituteVectorTest, NotEqual) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    auto comp = create_comparison_vector(row ? nrow : ncol);
    auto sub = create_replacement_vector(comp.size(), static_cast<int>(row) * 1000 + static_cast<int>(is_sparse) * 29);
    if (is_sparse) {
        for (size_t i = 0, end = comp.size(); i < end; ++i){
            if (comp[i] != 0) {
                sub[i] = 0;
            }
        }
    }

    auto op = std::make_shared<tatami::DelayedUnaryIsometricSubstituteNotEqualVectorHelper<double, double, int, decltype(comp), decltype(sub)> >(comp, sub, row);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());
    EXPECT_EQ(sparse_mod.is_sparse(), is_sparse);

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = simulated;
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            auto& x = refvec[r * ncol + c];
            auto i = (row ? r : c);
            x = (x != comp[i] ? sub[i] : x);
        }
    }

    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricSubstituteVector,
    DelayedUnaryIsometricSubstituteVectorTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or by column
        ::testing::Values(true, false) // check sparse case
    )
);
