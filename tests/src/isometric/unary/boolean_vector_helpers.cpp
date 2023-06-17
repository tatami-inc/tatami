#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOp.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class BooleanVectorTest : public ::testing::TestWithParam<std::tuple<bool, bool> > {
protected:
    size_t nrow = 191, ncol = 88;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    std::vector<double> simulated;
protected:
    void SetUp() {
        simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1, -3, 3);
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated));
        sparse = tatami::convert_to_sparse<false>(dense.get()); // column major.
        return;
    }

    static void fill_default_vector(std::vector<char>& vec) {
        bool val = true;
        for (auto& x : vec) {
            x = val;
            val = !val;
        }
    }
};

TEST_P(BooleanVectorTest, AND) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<char> vec(row ? nrow : ncol);
    if (is_sparse) {
        std::fill(vec.begin(), vec.end(), 1);
    } else {
        fill_default_vector(vec);
    }

    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;
    if (row) {
        auto op = tatami::make_DelayedBooleanAndVectorHelper<0>(vec);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
    } else {
        auto op = tatami::make_DelayedBooleanAndVectorHelper<1>(vec);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
    }

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    EXPECT_TRUE(sparse_mod->sparse());

    auto refvec = this->simulated;
    for (size_t r = 0; r < this->nrow; ++r) {
        for (size_t c = 0; c < this->ncol; ++c) {
            auto& x = refvec[r * this->ncol + c];
            x = (x && vec[row ? r : c]);
        }
    }
    
    tatami::DenseRowMatrix<double> ref(this->nrow, this->ncol, std::move(refvec));
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(BooleanVectorTest, OR) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<char> vec(row ? nrow : ncol);
    if (!is_sparse) {
        fill_default_vector(vec);
    }

    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;
    if (row) {
        auto op = tatami::make_DelayedBooleanOrVectorHelper<0>(vec);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
    } else {
        auto op = tatami::make_DelayedBooleanOrVectorHelper<1>(vec);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
    }

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    if (is_sparse) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }

    auto refvec = this->simulated;
    for (size_t r = 0; r < this->nrow; ++r) {
        for (size_t c = 0; c < this->ncol; ++c) {
            auto& x = refvec[r * this->ncol + c];
            x = (x || vec[row ? r : c]);
        }
    }
    
    tatami::DenseRowMatrix<double> ref(this->nrow, this->ncol, std::move(refvec));
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(BooleanVectorTest, XOR) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<char> vec(row ? nrow : ncol);
    if (!is_sparse) {
        fill_default_vector(vec);
    }

    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;
    if (row) {
        auto op = tatami::make_DelayedBooleanXorVectorHelper<0>(vec);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
    } else {
        auto op = tatami::make_DelayedBooleanXorVectorHelper<1>(vec);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
    }

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    if (is_sparse) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }

    auto refvec = this->simulated;
    for (size_t r = 0; r < this->nrow; ++r) {
        for (size_t c = 0; c < this->ncol; ++c) {
            auto& x = refvec[r * this->ncol + c];
            x = (static_cast<bool>(x) != static_cast<bool>(vec[row ? r : c]));
        }
    }
    
    tatami::DenseRowMatrix<double> ref(this->nrow, this->ncol, std::move(refvec));
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(BooleanVectorTest, EQUAL) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<char> vec(row ? nrow : ncol);
    if (is_sparse) {
        std::fill(vec.begin(), vec.end(), 1);
    } else {
        fill_default_vector(vec);
    }

    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;
    if (row) {
        auto op = tatami::make_DelayedBooleanEqualVectorHelper<0>(vec);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
    } else {
        auto op = tatami::make_DelayedBooleanEqualVectorHelper<1>(vec);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
    }

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    if (is_sparse) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }

    auto refvec = this->simulated;
    for (size_t r = 0; r < this->nrow; ++r) {
        for (size_t c = 0; c < this->ncol; ++c) {
            auto& x = refvec[r * this->ncol + c];
            x = (static_cast<bool>(x) == static_cast<bool>(vec[row ? r : c]));
        }
    }
    
    tatami::DenseRowMatrix<double> ref(this->nrow, this->ncol, std::move(refvec));
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

INSTANTIATE_TEST_CASE_P(
    BooleanVector,
    BooleanVectorTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or by column
        ::testing::Values(true, false) // check sparse case
    )
);
