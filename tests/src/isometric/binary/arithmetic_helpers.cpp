#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/binary/DelayedBinaryIsometricOperation.hpp"
#include "tatami/isometric/binary/arithmetic_helpers.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOperation.hpp"
#include "tatami/isometric/unary/math_helpers.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class DelayedBinaryIsometricArithmeticUtils {
protected:
    inline static size_t nrow = 91, ncol = 121;
    inline static std::shared_ptr<tatami::NumericMatrix> dense_left, sparse_left, dense_right, sparse_right;
    inline static std::vector<double> simulated_left, simulated_right;
    
    static void assemble() {
        if (dense_left) {
            return;
        }

        simulated_left = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.12;
            opt.lower = -5;
            opt.upper = 5;
            opt.seed = 12345;
            return opt;
        }());
        dense_left.reset(new tatami::DenseMatrix<double, int, decltype(simulated_left)>(nrow, ncol, simulated_left, true)); // row major.
        sparse_left = tatami::convert_to_compressed_sparse<double, int>(*dense_left, false, {}); // column major.

        simulated_right = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.15;
            opt.lower = -5;
            opt.upper = 5;
            opt.seed = 67890;
            return opt;
        }());
        dense_right.reset(new tatami::DenseMatrix<double, int, decltype(simulated_right)>(nrow, ncol, simulated_right, true)); // row major.
        sparse_right = tatami::convert_to_compressed_sparse<double, int>(*dense_right, false, {}); // column major.
        return;
    }
};

#define BINARY_ARITH_BASIC_SETUP(name, base) \
    class name : \
        public ::testing::Test, \
        public base { \
    protected: \
        static void SetUpTestSuite() { \
            assemble(); \
        } \
    };

#define BINARY_ARITH_FULL_TEST(name, base) \
    class name : \
        public ::testing::TestWithParam<tatami_test::StandardTestAccessOptions>, \
        public base { \
    protected: \
        static void SetUpTestSuite() { \
            assemble(); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto opts = tatami_test::convert_test_access_options(GetParam()); \
        tatami_test::test_full_access(*dense_mod, *ref, opts); \
        tatami_test::test_full_access(*sparse_mod, *ref, opts); \
        tatami_test::test_unsorted_full_access(*sparse_uns, opts); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        DelayedBinaryIsometricArithmetic, \
        name, \
        tatami_test::standard_test_access_options_combinations() \
    );

#define BINARY_ARITH_BLOCK_TEST(name, base) \
    class name : \
        public ::testing::TestWithParam<std::tuple<tatami_test::StandardTestAccessOptions, std::pair<double, double> > >, \
        public base { \
    protected: \
        static void SetUpTestSuite() { \
            assemble(); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto tparam = GetParam(); \
        auto opts = tatami_test::convert_test_access_options(std::get<0>(tparam)); \
        auto interval_info = std::get<1>(tparam); \
        tatami_test::test_block_access(*dense_mod, *ref, interval_info.first, interval_info.second, opts); \
        tatami_test::test_block_access(*sparse_mod, *ref, interval_info.first, interval_info.second, opts); \
        tatami_test::test_unsorted_block_access(*sparse_uns, interval_info.first, interval_info.second, opts); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        DelayedBinaryIsometricArithmetic, \
        name, \
        ::testing::Combine( \
            tatami_test::standard_test_access_options_combinations(), \
            ::testing::Values( \
                std::make_pair(0.0, 0.35), \
                std::make_pair(0.27, 0.6), \
                std::make_pair(0.67, 0.33) \
            ) \
        ) \
    );

#define BINARY_ARITH_INDEX_TEST(name, base) \
    class name : \
        public ::testing::TestWithParam<std::tuple<tatami_test::StandardTestAccessOptions, std::pair<double, double> > >, \
        public base { \
    protected: \
        static void SetUpTestSuite() { \
            assemble(); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto tparam = GetParam(); \
        auto opts = tatami_test::convert_test_access_options(std::get<0>(tparam)); \
        auto interval_info = std::get<1>(tparam); \
        tatami_test::test_indexed_access(*dense_mod, *ref, interval_info.first, interval_info.second, opts); \
        tatami_test::test_indexed_access(*sparse_mod, *ref, interval_info.first, interval_info.second, opts); \
        tatami_test::test_unsorted_indexed_access(*sparse_uns, interval_info.first, interval_info.second, opts); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        DelayedBinaryIsometricArithmetic, \
        name, \
        ::testing::Combine( \
            tatami_test::standard_test_access_options_combinations(), \
            ::testing::Values( \
                std::make_pair(0.0, 0.4), \
                std::make_pair(0.21, 0.15), \
                std::make_pair(0.56, 0.2) \
            ) \
        ) \
    );

/****************************
 ********* ADDITION *********
 ****************************/

class DelayedBinaryIsometricAddUtils : public DelayedBinaryIsometricArithmeticUtils {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;

    static void assemble() {
        if (ref) {
            return;
        }

        DelayedBinaryIsometricArithmeticUtils::assemble();

        auto op = std::make_shared<tatami::DelayedBinaryIsometricAddHelper<double, double, int> >();
        dense_mod.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(dense_left, dense_right, op));
        sparse_mod.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(sparse_left, sparse_right, op));
        sparse_uns.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::ReversedIndicesWrapper<double, int>(sparse_left)),
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::ReversedIndicesWrapper<double, int>(sparse_right)),
            op
        ));

        auto refvec = simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] += simulated_right[i];
        }
        ref.reset(new tatami::DenseMatrix<double, int, decltype(refvec)>(nrow, ncol, std::move(refvec), true));
    }
};

BINARY_ARITH_BASIC_SETUP(DelayedBinaryIsometricAddTest, DelayedBinaryIsometricAddUtils)
TEST_F(DelayedBinaryIsometricAddTest, Basic) {
    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense_mod->is_sparse_proportion(), 0);
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(sparse_mod->is_sparse_proportion(), 1);
    EXPECT_EQ(nrow, dense_mod->nrow());
    EXPECT_EQ(ncol, dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_EQ(dense_mod->prefer_rows_proportion(), 1);
    EXPECT_FALSE(sparse_mod->prefer_rows());
    EXPECT_EQ(sparse_mod->prefer_rows_proportion(), 0);

    tatami::DelayedBinaryIsometricOperation<double, double, int> mixed_mod(
        sparse_left,
        dense_right,
        std::make_shared<tatami::DelayedBinaryIsometricAddHelper<double, double, int> >()
    );
    EXPECT_FALSE(mixed_mod.is_sparse());
    EXPECT_EQ(mixed_mod.prefer_rows_proportion(), 0.5);
    EXPECT_EQ(mixed_mod.is_sparse_proportion(), 0.5);
}

BINARY_ARITH_FULL_TEST(DelayedBinaryIsometricAddFullTest, DelayedBinaryIsometricAddUtils)
BINARY_ARITH_BLOCK_TEST(DelayedBinaryIsometricAddBlockTest, DelayedBinaryIsometricAddUtils)
BINARY_ARITH_INDEX_TEST(DelayedBinaryIsometricAddIndexTest, DelayedBinaryIsometricAddUtils)

TEST_F(DelayedBinaryIsometricAddTest, NewType) {
    auto op = std::make_shared<tatami::DelayedBinaryIsometricAddHelper<float, double, int> >();
    tatami::DelayedBinaryIsometricOperation<float, double, int> dense_fmod(dense_left, dense_right, op);
    tatami::DelayedBinaryIsometricOperation<float, double, int> sparse_fmod(sparse_left, sparse_right, op);

    std::vector<float> frefvec(simulated_left.size());
    for (size_t i = 0; i < frefvec.size(); ++i) {
        frefvec[i] = simulated_left[i] + simulated_right[i];
    }
    tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

    quick_test_all<float, int>(dense_fmod, fref);
    quick_test_all<float, int>(sparse_fmod, fref);
}

/*******************************
 ********* SUBTRACTION *********
 *******************************/

class DelayedBinaryIsometricSubtractUtils : public DelayedBinaryIsometricArithmeticUtils {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;

    static void assemble() {
        if (ref) {
            return;
        }

        DelayedBinaryIsometricArithmeticUtils::assemble();

        auto op = std::make_shared<tatami::DelayedBinaryIsometricSubtractHelper<double, double, int> >();
        dense_mod.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(dense_left, dense_right, op));
        sparse_mod.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(sparse_left, sparse_right, op));
        sparse_uns.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::ReversedIndicesWrapper<double, int>(sparse_left)),
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::ReversedIndicesWrapper<double, int>(sparse_right)), 
            op
        ));

        auto refvec = simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] -= simulated_right[i];
        }
        ref.reset(new tatami::DenseMatrix<double, int, decltype(refvec)>(nrow, ncol, std::move(refvec), true));
    }
};

BINARY_ARITH_BASIC_SETUP(DelayedBinaryIsometricSubtractTest, DelayedBinaryIsometricSubtractUtils)
TEST_F(DelayedBinaryIsometricSubtractTest, Basic) {
    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());

    auto op = std::make_shared<tatami::DelayedBinaryIsometricSubtractHelper<double, double, int> >();
    tatami::DelayedBinaryIsometricOperation<double, double, int> mixed_mod(sparse_left, dense_right, std::move(op));
    EXPECT_FALSE(mixed_mod.is_sparse());
}

BINARY_ARITH_FULL_TEST(DelayedBinaryIsometricSubtractFullTest, DelayedBinaryIsometricSubtractUtils)
BINARY_ARITH_BLOCK_TEST(DelayedBinaryIsometricSubtractBlockTest, DelayedBinaryIsometricSubtractUtils)
BINARY_ARITH_INDEX_TEST(DelayedBinaryIsometricSubtractIndexTest, DelayedBinaryIsometricSubtractUtils)

/**********************************
 ********* MULTIPLICATION *********
 **********************************/

class DelayedBinaryIsometricMultiplyUtils : public DelayedBinaryIsometricArithmeticUtils {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;

    static void assemble() {
        if (ref) {
            return;
        }

        DelayedBinaryIsometricArithmeticUtils::assemble();

        auto op = std::make_shared<tatami::DelayedBinaryIsometricMultiplyHelper<double, double, int> >();
        dense_mod.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(dense_left, dense_right, op));
        sparse_mod.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(sparse_left, sparse_right, op));
        sparse_uns.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::ReversedIndicesWrapper<double, int>(sparse_left)),
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::ReversedIndicesWrapper<double, int>(sparse_right)), 
            op
        ));

        auto refvec = simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] *= simulated_right[i];
        }
        ref.reset(new tatami::DenseMatrix<double, int, decltype(refvec)>(nrow, ncol, std::move(refvec), true));
    }
};

BINARY_ARITH_BASIC_SETUP(DelayedBinaryIsometricMultiplyTest, DelayedBinaryIsometricMultiplyUtils)
TEST_F(DelayedBinaryIsometricMultiplyTest, Basic) {
    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
}

BINARY_ARITH_FULL_TEST(DelayedBinaryIsometricMultiplyFullTest, DelayedBinaryIsometricMultiplyUtils)
BINARY_ARITH_BLOCK_TEST(DelayedBinaryIsometricMultiplyBlockTest, DelayedBinaryIsometricMultiplyUtils)
BINARY_ARITH_INDEX_TEST(DelayedBinaryIsometricMultiplyIndexTest, DelayedBinaryIsometricMultiplyUtils)

/****************************
 ********* DIVISION *********
 ****************************/

class DelayedBinaryIsometricDivideUtils : public DelayedBinaryIsometricArithmeticUtils {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;

    static void assemble() {
        if (ref) {
            return;
        }

        DelayedBinaryIsometricArithmeticUtils::assemble();

        auto op = std::make_shared<tatami::DelayedBinaryIsometricDivideHelper<double, double, int> >();
        dense_mod.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(dense_left, dense_right, op));
        sparse_mod.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(sparse_left, sparse_right, op));
        sparse_uns.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::ReversedIndicesWrapper<double, int>(sparse_left)),
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::ReversedIndicesWrapper<double, int>(sparse_right)), 
            op
        ));

        auto refvec = simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] = careful_division(refvec[i], simulated_right[i]);
        }
        ref.reset(new tatami::DenseMatrix<double, int, decltype(refvec)>(nrow, ncol, std::move(refvec), true));
    }
};

BINARY_ARITH_BASIC_SETUP(DelayedBinaryIsometricDivideTest, DelayedBinaryIsometricDivideUtils)
TEST_F(DelayedBinaryIsometricDivideTest, Basic) {
    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());
    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

BINARY_ARITH_FULL_TEST(DelayedBinaryIsometricDivideFullTest, DelayedBinaryIsometricDivideUtils)
BINARY_ARITH_BLOCK_TEST(DelayedBinaryIsometricDivideBlockTest, DelayedBinaryIsometricDivideUtils)
BINARY_ARITH_INDEX_TEST(DelayedBinaryIsometricDivideIndexTest, DelayedBinaryIsometricDivideUtils)

TEST_F(DelayedBinaryIsometricDivideTest, NewType) {
    auto op = std::make_shared<tatami::DelayedBinaryIsometricDivideHelper<float, double, int> >();
    tatami::DelayedBinaryIsometricOperation<float, double, int> dense_fmod(dense_left, dense_right, op);
    tatami::DelayedBinaryIsometricOperation<float, double, int> sparse_fmod(sparse_left, sparse_right, op);

    std::vector<float> frefvec(simulated_left.size());
    for (size_t i = 0; i < frefvec.size(); ++i) {
        frefvec[i] = simulated_left[i] / simulated_right[i];
    }
    tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

    quick_test_all<float, int>(dense_fmod, fref);
    quick_test_all<float, int>(sparse_fmod, fref);
}

/*******************************
 ************ POWER ************
 *******************************/

class DelayedBinaryIsometricPowerUtils : public DelayedBinaryIsometricArithmeticUtils {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;

    static void assemble() {
        if (ref) {
            return;
        }

        DelayedBinaryIsometricArithmeticUtils::assemble();

        auto abs_left = simulated_left;
        for (auto& x : abs_left) {
            x = std::abs(x);
        }
        auto dense_left0 = std::make_shared<tatami::DenseMatrix<double, int, decltype(abs_left)> >(nrow, ncol, std::move(abs_left), true); // row major.
        auto sparse_left0 = tatami::convert_to_compressed_sparse<double, int>(*dense_left0, false, {}); // column major.

        auto abs_right = simulated_right;
        for (auto& x : abs_right) {
            x = std::abs(x);
        }
        auto dense_right0 = std::make_shared<tatami::DenseMatrix<double, int, decltype(abs_right)> >(nrow, ncol, std::move(abs_right), true); // row major.
        auto sparse_right0 = tatami::convert_to_compressed_sparse<double, int>(*dense_right0, false, {}); // column major.

        auto op = std::make_shared<tatami::DelayedBinaryIsometricPowerHelper<double, double, int> >();
        dense_mod.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(dense_left0, dense_right0, op));
        sparse_mod.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(sparse_left0, sparse_right0, op));
        sparse_uns.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::ReversedIndicesWrapper<double, int>(sparse_left0)),
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::ReversedIndicesWrapper<double, int>(sparse_right0)), 
            op
        ));

        auto refvec = simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] = std::pow(std::abs(refvec[i]), std::abs(simulated_right[i]));
        }
        ref.reset(new tatami::DenseMatrix<double, int, decltype(refvec)>(nrow, ncol, std::move(refvec), true));
    }
};

BINARY_ARITH_BASIC_SETUP(DelayedBinaryIsometricPowerTest, DelayedBinaryIsometricPowerUtils)
TEST_F(DelayedBinaryIsometricPowerTest, Basic) {
    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());
}

BINARY_ARITH_FULL_TEST(DelayedBinaryIsometricPowerFullTest, DelayedBinaryIsometricPowerUtils)
BINARY_ARITH_BLOCK_TEST(DelayedBinaryIsometricPowerBlockTest, DelayedBinaryIsometricPowerUtils)
BINARY_ARITH_INDEX_TEST(DelayedBinaryIsometricPowerIndexTest, DelayedBinaryIsometricPowerUtils)

/****************************
 ********** MODULO **********
 ****************************/

class DelayedBinaryIsometricModuloUtils : public DelayedBinaryIsometricArithmeticUtils {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;

    static void assemble() {
        if (ref) {
            return;
        }

        DelayedBinaryIsometricArithmeticUtils::assemble();

        auto op = std::make_shared<tatami::DelayedBinaryIsometricModuloHelper<double, double, int> >();
        dense_mod.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(dense_left, dense_right, op));
        sparse_mod.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(sparse_left, sparse_right, op));
        sparse_uns.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::ReversedIndicesWrapper<double, int>(sparse_left)),
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::ReversedIndicesWrapper<double, int>(sparse_right)), 
            op
        ));

        auto refvec = simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] = careful_modulo(refvec[i], simulated_right[i]);
        }
        ref.reset(new tatami::DenseMatrix<double, int, decltype(refvec)>(nrow, ncol, std::move(refvec), true));
    }
};

BINARY_ARITH_BASIC_SETUP(DelayedBinaryIsometricModuloTest, DelayedBinaryIsometricModuloUtils)
TEST_F(DelayedBinaryIsometricModuloTest, Basic) {
    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());
    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

BINARY_ARITH_FULL_TEST(DelayedBinaryIsometricModuloFullTest, DelayedBinaryIsometricModuloUtils)
BINARY_ARITH_BLOCK_TEST(DelayedBinaryIsometricModuloBlockTest, DelayedBinaryIsometricModuloUtils)
BINARY_ARITH_INDEX_TEST(DelayedBinaryIsometricModuloIndexTest, DelayedBinaryIsometricModuloUtils)

/****************************
 ***** INTEGER DIVISION *****
 ****************************/

class DelayedBinaryIsometricIntegerDivideUtils : public DelayedBinaryIsometricArithmeticUtils {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;

    static void assemble() {
        if (ref) {
            return;
        }

        DelayedBinaryIsometricArithmeticUtils::assemble();

        auto op = std::make_shared<tatami::DelayedBinaryIsometricIntegerDivideHelper<double, double, int> >();
        dense_mod.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(dense_left, dense_right, op));
        sparse_mod.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(sparse_left, sparse_right, op));
        sparse_uns.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::ReversedIndicesWrapper<double, int>(sparse_left)),
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::ReversedIndicesWrapper<double, int>(sparse_right)), 
            op
        ));

        auto refvec = simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            // x == (x %% y) + y * (x %/% y)
            refvec[i] = std::floor(refvec[i] / simulated_right[i]);
        }
        ref.reset(new tatami::DenseMatrix<double, int, decltype(refvec)>(nrow, ncol, std::move(refvec), true));
    }
};

BINARY_ARITH_BASIC_SETUP(DelayedBinaryIsometricIntegerDivideTest, DelayedBinaryIsometricIntegerDivideUtils)
TEST_F(DelayedBinaryIsometricIntegerDivideTest, Basic) {
    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());
    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

BINARY_ARITH_FULL_TEST(DelayedBinaryIsometricIntegerDivideFullTest, DelayedBinaryIsometricIntegerDivideUtils)
BINARY_ARITH_BLOCK_TEST(DelayedBinaryIsometricIntegerDivideBlockTest, DelayedBinaryIsometricIntegerDivideUtils)
BINARY_ARITH_INDEX_TEST(DelayedBinaryIsometricIntegerDivideIndexTest, DelayedBinaryIsometricIntegerDivideUtils)

TEST(DelayedBinaryIsometricArithmetic, BackCompatibility) {
    auto add = tatami::make_DelayedBinaryIsometricAdd();
    EXPECT_TRUE(add->is_sparse());
    auto sub = tatami::make_DelayedBinaryIsometricSubtract();
    EXPECT_TRUE(sub->is_sparse());
    auto mult = tatami::make_DelayedBinaryIsometricMultiply();
    EXPECT_TRUE(mult->is_sparse());
    auto div = tatami::make_DelayedBinaryIsometricDivide();
    EXPECT_FALSE(div->is_sparse());
    auto pow = tatami::make_DelayedBinaryIsometricPower();
    EXPECT_FALSE(pow->is_sparse());
    auto mod = tatami::make_DelayedBinaryIsometricModulo();
    EXPECT_FALSE(mod->is_sparse());
    auto intdiv = tatami::make_DelayedBinaryIsometricIntegerDivide();
    EXPECT_FALSE(intdiv->is_sparse());
}
