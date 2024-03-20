#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/binary/DelayedBinaryIsometricOp.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOp.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class BinaryArithUtils {
protected:
    inline static size_t nrow = 91, ncol = 121;
    inline static std::shared_ptr<tatami::NumericMatrix> dense_left, sparse_left, dense_right, sparse_right;
    inline static std::vector<double> simulated_left, simulated_right;
    
    static void assemble() {
        if (dense_left) {
            return;
        }

        simulated_left = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1, /* lower = */ -5, /* upper = */ 5, /* seed */ 12345);
        dense_left = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated_left));
        sparse_left = tatami::convert_to_compressed_sparse<false>(dense_left.get()); // column major.

        simulated_right = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1, /* lower = */ -5, /* upper = */ 5, /* seed */ 67890);
        dense_right = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated_right));
        sparse_right = tatami::convert_to_compressed_sparse<false>(dense_right.get()); // column major.
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
        public ::testing::TestWithParam<tatami_test::StandardTestAccessParameters>, \
        public base { \
    protected: \
        static void SetUpTestSuite() { \
            assemble(); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto params = tatami_test::convert_access_parameters(GetParam()); \
        tatami_test::test_full_access(params, dense_mod.get(), ref.get()); \
        tatami_test::test_full_access(params, sparse_mod.get(), ref.get()); \
        tatami_test::test_full_access(params, sparse_uns.get(), ref.get()); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        BinaryArith, \
        name, \
        tatami_test::standard_test_access_parameter_combinations() \
    );

#define BINARY_ARITH_BLOCK_TEST(name, base) \
    class name : \
        public ::testing::TestWithParam<std::tuple<tatami_test::StandardTestAccessParameters, std::pair<double, double> > >, \
        public base { \
    protected: \
        static void SetUpTestSuite() { \
            assemble(); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto tparam = GetParam(); \
        auto params = tatami_test::convert_access_parameters(std::get<0>(tparam)); \
        auto interval_info = std::get<1>(tparam); \
        auto len = (params.use_row ? ref->ncol() : ref->nrow()); \
        size_t FIRST = interval_info.first * len, LAST = interval_info.second * len; \
        tatami_test::test_block_access(params, dense_mod.get(), ref.get(), FIRST, LAST); \
        tatami_test::test_block_access(params, sparse_mod.get(), ref.get(), FIRST, LAST); \
        tatami_test::test_block_access(params, sparse_uns.get(), ref.get(), FIRST, LAST); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        BinaryArith, \
        name, \
        ::testing::Combine( \
            tatami_test::standard_test_access_parameter_combinations(), \
            ::testing::Values( \
                std::make_pair(0, 0.35), \
                std::make_pair(0.27, 0.87), \
                std::make_pair(0.67, 1.0) \
            ) \
        ) \
    );

#define BINARY_ARITH_INDEX_TEST(name, base) \
    class name : \
        public ::testing::TestWithParam<std::tuple<tatami_test::StandardTestAccessParameters, std::pair<double, int> > >, \
        public base { \
    protected: \
        static void SetUpTestSuite() { \
            assemble(); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto tparam = GetParam(); \
        auto params = tatami_test::convert_access_parameters(std::get<0>(tparam)); \
        auto interval_info = std::get<1>(tparam); \
        auto len = (params.use_row ? ref->ncol() : ref->nrow()); \
        size_t FIRST = interval_info.first * len, STEP = interval_info.second; \
        tatami_test::test_indexed_access(params, dense_mod.get(), ref.get(), FIRST, STEP); \
        tatami_test::test_indexed_access(params, sparse_mod.get(), ref.get(), FIRST, STEP); \
        tatami_test::test_indexed_access(params, sparse_uns.get(), ref.get(), FIRST, STEP); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        BinaryArith, \
        name, \
        ::testing::Combine( \
            tatami_test::standard_test_access_parameter_combinations(), \
            ::testing::Values( \
                std::make_pair(0.0, 4), \
                std::make_pair(0.21, 9), \
                std::make_pair(0.56, 10) \
            ) \
        ) \
    );

/****************************
 ********* ADDITION *********
 ****************************/

class BinaryArithAdditionUtils : public BinaryArithUtils {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;

    static void assemble() {
        if (ref) {
            return;
        }

        BinaryArithUtils::assemble();
        auto op = tatami::make_DelayedBinaryAddHelper();
        dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
        sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);
        sparse_uns = tatami::make_DelayedBinaryIsometricOp(
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::UnsortedWrapper<double, int>(sparse_left)),
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::UnsortedWrapper<double, int>(sparse_right)),
            op
        );

        auto refvec = simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] += simulated_right[i];
        }
        ref.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec)));
    }
};

BINARY_ARITH_BASIC_SETUP(BinaryArithAdditionTest, BinaryArithAdditionUtils)
TEST_F(BinaryArithAdditionTest, Basic) {
    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense_mod->sparse_proportion(), 0);
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(sparse_mod->sparse_proportion(), 1);
    EXPECT_EQ(nrow, dense_mod->nrow());
    EXPECT_EQ(ncol, dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_EQ(dense_mod->prefer_rows_proportion(), 1);
    EXPECT_FALSE(sparse_mod->prefer_rows());
    EXPECT_EQ(sparse_mod->prefer_rows_proportion(), 0);

    auto mixed_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, dense_right, tatami::make_DelayedBinaryAddHelper());
    EXPECT_FALSE(mixed_mod->sparse());
    EXPECT_EQ(mixed_mod->prefer_rows_proportion(), 0.5);
    EXPECT_EQ(mixed_mod->sparse_proportion(), 0.5);
}

BINARY_ARITH_FULL_TEST(BinaryArithAdditionFullTest, BinaryArithAdditionUtils)
BINARY_ARITH_BLOCK_TEST(BinaryArithAdditionBlockTest, BinaryArithAdditionUtils)
BINARY_ARITH_INDEX_TEST(BinaryArithAdditionIndexTest, BinaryArithAdditionUtils)

/*******************************
 ********* SUBTRACTION *********
 *******************************/

class BinaryArithSubtractionUtils : public BinaryArithUtils {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;

    static void assemble() {
        if (ref) {
            return;
        }

        BinaryArithUtils::assemble();
        auto op = tatami::make_DelayedBinarySubtractHelper();
        dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
        sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);
        sparse_uns = tatami::make_DelayedBinaryIsometricOp(
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::UnsortedWrapper<double, int>(sparse_left)),
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::UnsortedWrapper<double, int>(sparse_right)), 
            op
        );

        auto refvec = simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] -= simulated_right[i];
        }
        ref.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec)));
    }
};

BINARY_ARITH_BASIC_SETUP(BinaryArithSubtractionTest, BinaryArithSubtractionUtils)
TEST_F(BinaryArithSubtractionTest, Basic) {
    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());

    auto mixed_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, dense_right, tatami::make_DelayedBinarySubtractHelper());
    EXPECT_FALSE(mixed_mod->sparse());
}

BINARY_ARITH_FULL_TEST(BinaryArithSubtractionFullTest, BinaryArithSubtractionUtils)
BINARY_ARITH_BLOCK_TEST(BinaryArithSubtractionBlockTest, BinaryArithSubtractionUtils)
BINARY_ARITH_INDEX_TEST(BinaryArithSubtractionIndexTest, BinaryArithSubtractionUtils)

/**********************************
 ********* MULTIPLICATION *********
 **********************************/

class BinaryArithMultiplicationUtils : public BinaryArithUtils {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;

    static void assemble() {
        if (ref) {
            return;
        }

        BinaryArithUtils::assemble();

        auto op = tatami::make_DelayedBinaryMultiplyHelper();
        dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
        sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);
        sparse_uns = tatami::make_DelayedBinaryIsometricOp(
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::UnsortedWrapper<double, int>(sparse_left)),
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::UnsortedWrapper<double, int>(sparse_right)), 
            op
        );

        auto refvec = simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] *= simulated_right[i];
        }
        ref.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec)));
    }
};

BINARY_ARITH_BASIC_SETUP(BinaryArithMultiplicationTest, BinaryArithMultiplicationUtils)
TEST_F(BinaryArithMultiplicationTest, Basic) {
    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
}

BINARY_ARITH_FULL_TEST(BinaryArithMultiplicationFullTest, BinaryArithMultiplicationUtils)
BINARY_ARITH_BLOCK_TEST(BinaryArithMultiplicationBlockTest, BinaryArithMultiplicationUtils)
BINARY_ARITH_INDEX_TEST(BinaryArithMultiplicationIndexTest, BinaryArithMultiplicationUtils)

/****************************
 ********* DIVISION *********
 ****************************/

#define BINARY_ARITH_FULL_TEST_WITH_NAN(name, base) \
    class name : \
        public ::testing::TestWithParam<tatami_test::StandardTestAccessParameters>, \
        public base { \
    protected: \
        static void SetUpTestSuite() { \
            assemble(); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto params = tatami_test::convert_access_parameters(GetParam()); \
        params.has_nan = true; \
        tatami_test::test_full_access(params, dense_mod.get(), ref.get()); \
        tatami_test::test_full_access(params, sparse_mod.get(), ref.get()); \
        tatami_test::test_full_access(params, sparse_uns.get(), ref.get()); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        BinaryArith, \
        name, \
        tatami_test::standard_test_access_parameter_combinations() \
    );

#define BINARY_ARITH_BLOCK_TEST_WITH_NAN(name, base) \
    class name : \
        public ::testing::TestWithParam<std::tuple<tatami_test::StandardTestAccessParameters, std::pair<double, double> > >, \
        public base { \
    protected: \
        static void SetUpTestSuite() { \
            assemble(); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto tparam = GetParam(); \
        auto params = tatami_test::convert_access_parameters(std::get<0>(tparam)); \
        params.has_nan = true; \
        auto interval_info = std::get<1>(tparam); \
        auto len = (params.use_row ? ref->ncol() : ref->nrow()); \
        size_t FIRST = interval_info.first * len, LAST = interval_info.second * len; \
        tatami_test::test_block_access(params, dense_mod.get(), ref.get(), FIRST, LAST); \
        tatami_test::test_block_access(params, sparse_mod.get(), ref.get(), FIRST, LAST); \
        tatami_test::test_block_access(params, sparse_uns.get(), ref.get(), FIRST, LAST); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        BinaryArith, \
        name, \
        ::testing::Combine( \
            tatami_test::standard_test_access_parameter_combinations(), \
            ::testing::Values( \
                std::make_pair(0, 0.35), \
                std::make_pair(0.27, 0.87), \
                std::make_pair(0.67, 1.0) \
            ) \
        ) \
    );

#define BINARY_ARITH_INDEX_TEST_WITH_NAN(name, base) \
    class name : \
        public ::testing::TestWithParam<std::tuple<tatami_test::StandardTestAccessParameters, std::pair<double, int> > >, \
        public base { \
    protected: \
        static void SetUpTestSuite() { \
            assemble(); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto tparam = GetParam(); \
        auto params = tatami_test::convert_access_parameters(std::get<0>(tparam)); \
        params.has_nan = true; \
        auto interval_info = std::get<1>(tparam); \
        auto len = (params.use_row ? ref->ncol() : ref->nrow()); \
        size_t FIRST = interval_info.first * len, STEP = interval_info.second; \
        tatami_test::test_indexed_access(params, dense_mod.get(), ref.get(), FIRST, STEP); \
        tatami_test::test_indexed_access(params, sparse_mod.get(), ref.get(), FIRST, STEP); \
        tatami_test::test_indexed_access(params, sparse_uns.get(), ref.get(), FIRST, STEP); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        BinaryArith, \
        name, \
        ::testing::Combine( \
            tatami_test::standard_test_access_parameter_combinations(), \
            ::testing::Values( \
                std::make_pair(0.0, 4), \
                std::make_pair(0.21, 9), \
                std::make_pair(0.56, 10) \
            ) \
        ) \
    );


class BinaryArithDivisionUtils : public BinaryArithUtils {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;

    static void assemble() {
        if (ref) {
            return;
        }

        BinaryArithUtils::assemble();

        auto op = tatami::make_DelayedBinaryDivideHelper();
        dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
        sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);
        sparse_uns = tatami::make_DelayedBinaryIsometricOp(
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::UnsortedWrapper<double, int>(sparse_left)),
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::UnsortedWrapper<double, int>(sparse_right)), 
            op
        );

        auto refvec = simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] = careful_division(refvec[i], simulated_right[i]);
        }
        ref.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec)));
    }
};

BINARY_ARITH_BASIC_SETUP(BinaryArithDivisionTest, BinaryArithDivisionUtils)
TEST_F(BinaryArithDivisionTest, Basic) {
    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

BINARY_ARITH_FULL_TEST_WITH_NAN(BinaryArithDivisionFullTest, BinaryArithDivisionUtils)
BINARY_ARITH_BLOCK_TEST_WITH_NAN(BinaryArithDivisionBlockTest, BinaryArithDivisionUtils)
BINARY_ARITH_INDEX_TEST_WITH_NAN(BinaryArithDivisionIndexTest, BinaryArithDivisionUtils)

/*******************************
 ************ POWER ************
 *******************************/

class BinaryArithPowerUtils : public BinaryArithUtils {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;

    static void assemble() {
        if (ref) {
            return;
        }

        BinaryArithUtils::assemble();

        tatami::DelayedAbsHelper op0;
        auto dense_left0 = tatami::make_DelayedUnaryIsometricOp(dense_left, op0);
        auto sparse_left0 = tatami::make_DelayedUnaryIsometricOp(sparse_left, op0);
        auto dense_right0 = tatami::make_DelayedUnaryIsometricOp(dense_right, op0);
        auto sparse_right0 = tatami::make_DelayedUnaryIsometricOp(sparse_right, op0);

        auto op = tatami::make_DelayedBinaryPowerHelper();
        dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left0, dense_right0, op);
        sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left0, sparse_right0, op);
        sparse_uns = tatami::make_DelayedBinaryIsometricOp(
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::UnsortedWrapper<double, int>(sparse_left0)),
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::UnsortedWrapper<double, int>(sparse_right0)), 
            op
        );

        auto refvec = simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] = std::pow(std::abs(refvec[i]), std::abs(simulated_right[i]));
        }
        ref.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec)));
    }
};

BINARY_ARITH_BASIC_SETUP(BinaryArithPowerTest, BinaryArithPowerUtils)
TEST_F(BinaryArithPowerTest, Basic) {
    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
}

BINARY_ARITH_FULL_TEST(BinaryArithPowerFullTest, BinaryArithPowerUtils)
BINARY_ARITH_BLOCK_TEST(BinaryArithPowerBlockTest, BinaryArithPowerUtils)
BINARY_ARITH_INDEX_TEST(BinaryArithPowerIndexTest, BinaryArithPowerUtils)

/****************************
 ********** MODULO **********
 ****************************/

class BinaryArithModuloUtils : public BinaryArithUtils {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;

    static void assemble() {
        if (ref) {
            return;
        }

        BinaryArithUtils::assemble();

        auto op = tatami::make_DelayedBinaryModuloHelper();
        dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
        sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);
        sparse_uns = tatami::make_DelayedBinaryIsometricOp(
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::UnsortedWrapper<double, int>(sparse_left)),
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::UnsortedWrapper<double, int>(sparse_right)), 
            op
        );

        auto refvec = simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            refvec[i] = std::fmod(refvec[i], simulated_right[i]);
        }
        ref.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec)));
    }
};

BINARY_ARITH_BASIC_SETUP(BinaryArithModuloTest, BinaryArithModuloUtils)
TEST_F(BinaryArithModuloTest, Basic) {
    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

BINARY_ARITH_FULL_TEST_WITH_NAN(BinaryArithModuloFullTest, BinaryArithModuloUtils)
BINARY_ARITH_BLOCK_TEST_WITH_NAN(BinaryArithModuloBlockTest, BinaryArithModuloUtils)
BINARY_ARITH_INDEX_TEST_WITH_NAN(BinaryArithModuloIndexTest, BinaryArithModuloUtils)

/****************************
 ***** INTEGER DIVISION *****
 ****************************/

class BinaryArithIntegerDivisionUtils : public BinaryArithUtils {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;

    static void assemble() {
        if (ref) {
            return;
        }

        BinaryArithUtils::assemble();

        auto op = tatami::make_DelayedBinaryIntegerDivideHelper();
        dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
        sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);
        sparse_uns = tatami::make_DelayedBinaryIsometricOp(
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::UnsortedWrapper<double, int>(sparse_left)),
            std::shared_ptr<tatami::NumericMatrix>(new tatami_test::UnsortedWrapper<double, int>(sparse_right)), 
            op
        );

        auto refvec = simulated_left;
        for (size_t i = 0; i < refvec.size(); ++i) {
            // x == (x %% y) + y * (x %/% y)
            refvec[i] = std::floor(refvec[i] / simulated_right[i]);
        }
        ref.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec)));
    }
};

BINARY_ARITH_BASIC_SETUP(BinaryArithIntegerDivisionTest, BinaryArithIntegerDivisionUtils)
TEST_F(BinaryArithIntegerDivisionTest, Basic) {
    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

BINARY_ARITH_FULL_TEST_WITH_NAN(BinaryArithIntegerDivisionFullTest, BinaryArithIntegerDivisionUtils)
BINARY_ARITH_BLOCK_TEST_WITH_NAN(BinaryArithIntegerDivisionBlockTest, BinaryArithIntegerDivisionUtils)
BINARY_ARITH_INDEX_TEST_WITH_NAN(BinaryArithIntegerDivisionIndexTest, BinaryArithIntegerDivisionUtils)
