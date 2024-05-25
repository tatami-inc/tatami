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
        simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1, -3, 3);
        for (auto& x : simulated) {
            if (x) {
                // Rounding for easier tests of exact equality.
                x = std::round(x);
                if (x == 0) {
                    x = 1;
                }
            }
        }

        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulated));
        sparse = tatami::convert_to_compressed_sparse<false, double, int>(dense.get()); // column major.
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

    auto op = tatami::make_DelayedUnaryIsometricEqualVector(vec, row);
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(this->dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(this->sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    if (is_sparse) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = this->simulated;
    for (size_t r = 0; r < this->nrow; ++r) {
        for (size_t c = 0; c < this->ncol; ++c) {
            auto& x = refvec[r * this->ncol + c];
            x = (x == vec[row ? r : c]);
        }
    }
    
    tatami::DenseRowMatrix<double, int> ref(this->nrow, this->ncol, std::move(refvec));
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
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

    auto op = tatami::make_DelayedUnaryIsometricGreaterThanVector(vec, row);
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(this->dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(this->sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    if (is_sparse) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = this->simulated;
    for (size_t r = 0; r < this->nrow; ++r) {
        for (size_t c = 0; c < this->ncol; ++c) {
            auto& x = refvec[r * this->ncol + c];
            x = (x > vec[row ? r : c]);
        }
    }
    
    tatami::DenseRowMatrix<double, int> ref(this->nrow, this->ncol, std::move(refvec));
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
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

    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;
    auto op = tatami::make_DelayedUnaryIsometricLessThanVector(vec, row);
    dense_mod = tatami::make_DelayedUnaryIsometricOperation(this->dense, op);
    sparse_mod = tatami::make_DelayedUnaryIsometricOperation(this->sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    if (is_sparse) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = this->simulated;
    for (size_t r = 0; r < this->nrow; ++r) {
        for (size_t c = 0; c < this->ncol; ++c) {
            auto& x = refvec[r * this->ncol + c];
            x = (x < vec[row ? r : c]);
        }
    }
    
    tatami::DenseRowMatrix<double, int> ref(this->nrow, this->ncol, std::move(refvec));
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
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

    auto op = tatami::make_DelayedUnaryIsometricGreaterThanOrEqualVector(vec, row);
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(this->dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(this->sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    if (is_sparse) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = this->simulated;
    for (size_t r = 0; r < this->nrow; ++r) {
        for (size_t c = 0; c < this->ncol; ++c) {
            auto& x = refvec[r * this->ncol + c];
            x = (x >= vec[row ? r : c]);
        }
    }
    
    tatami::DenseRowMatrix<double, int> ref(this->nrow, this->ncol, std::move(refvec));
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
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

    auto op = tatami::make_DelayedUnaryIsometricLessThanOrEqualVector(vec, row);
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(this->dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(this->sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    if (is_sparse) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = this->simulated;
    for (size_t r = 0; r < this->nrow; ++r) {
        for (size_t c = 0; c < this->ncol; ++c) {
            auto& x = refvec[r * this->ncol + c];
            x = (x <= vec[row ? r : c]);
        }
    }
    
    tatami::DenseRowMatrix<double, int> ref(this->nrow, this->ncol, std::move(refvec));
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(DelayedUnaryIsometricCompareVectorTest, NotEqual) {
    auto param = GetParam();
    bool row = std::get<0>(param);
    bool is_sparse = std::get<1>(param);

    std::vector<double> vec(row ? nrow : ncol);
    if (!is_sparse) {
        fill_default_vector(vec);
    }

    auto op = tatami::make_DelayedUnaryIsometricNotEqualVector(vec, row);
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(this->dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(this->sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    if (is_sparse) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    auto refvec = this->simulated;
    for (size_t r = 0; r < this->nrow; ++r) {
        for (size_t c = 0; c < this->ncol; ++c) {
            auto& x = refvec[r * this->ncol + c];
            x = (x != vec[row ? r : c]);
        }
    }
    
    tatami::DenseRowMatrix<double, int> ref(this->nrow, this->ncol, std::move(refvec));
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricCompareVector,
    DelayedUnaryIsometricCompareVectorTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or by column
        ::testing::Values(true, false) // check sparse case
    )
);
