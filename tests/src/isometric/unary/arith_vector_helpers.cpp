#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOp.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class ArithVectorUtils {
protected:
    inline static size_t nrow = 291, ncol = 188;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    inline static std::vector<double> simulated;

    static void assemble() {
        if (dense) {
            return;
        }
        simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1);
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated));
        sparse = tatami::convert_to_compressed_sparse<false>(dense.get()); // column major.
        return;
    }

    static std::vector<double> create_vector(size_t n, double starter = 0, double jump = 1) { 
        std::vector<double> output(n, starter);
        for (size_t i = 1; i < n; ++i) {
            output[i] = output[i-1] + jump;
        }
        return output;
    }

    static void test_simple_row_access_wt_nan(const tatami::NumericMatrix* test, const tatami::NumericMatrix* ref) {
        tatami_test::TestAccessParameters params;
        params.has_nan = true;
        params.use_row = true;
        tatami_test::test_full_access(params, test, ref);
    }

    static void test_simple_column_access_wt_nan(const tatami::NumericMatrix* test, const tatami::NumericMatrix* ref) {
        tatami_test::TestAccessParameters params;
        params.has_nan = true;
        params.use_row = false;
        tatami_test::test_full_access(params, test, ref);
    }
};

#define ARITH_VECTOR_CONFIGURE(y) \
    if (ref && y == last_params) { \
        return; \
    } \
    last_params = y;

#define ARITH_VECTOR_BASIC_SETUP(name, base) \
    class name : \
        public ::testing::TestWithParam<typename base::SimulationParameters>, \
        public base { \
    protected: \
        void SetUp() { \
            assemble(GetParam()); \
        } \
    }; \
    \
    INSTANTIATE_TEST_SUITE_P( \
        ArithVector, \
        name, \
        base::simulation_parameter_combinations() \
    );

#define ARITH_VECTOR_FULL_TEST(name, base) \
    class name : \
        public ::testing::TestWithParam<std::tuple<typename base::SimulationParameters, tatami_test::StandardTestAccessParameters> >, \
        public base { \
    protected: \
        void SetUp() { \
            assemble(std::get<0>(GetParam())); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto params = tatami_test::convert_access_parameters(std::get<1>(GetParam())); \
        tatami_test::test_full_access(params, dense_mod.get(), ref.get()); \
        tatami_test::test_full_access(params, sparse_mod.get(), ref.get()); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        ArithVector, \
        name, \
        ::testing::Combine( \
            name::simulation_parameter_combinations(), \
            tatami_test::standard_test_access_parameter_combinations() \
        ) \
    );

#define ARITH_VECTOR_BLOCK_TEST(name, base) \
    class name : \
        public ::testing::TestWithParam<std::tuple<typename base::SimulationParameters, tatami_test::StandardTestAccessParameters, std::pair<double, double> > >, \
        public base { \
    protected: \
        void SetUp() { \
            assemble(std::get<0>(GetParam())); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto tparam = GetParam(); \
        auto params = tatami_test::convert_access_parameters(std::get<1>(tparam)); \
        auto interval_info = std::get<2>(tparam); \
        auto len = (params.use_row ? ref->ncol() : ref->nrow()); \
        size_t FIRST = interval_info.first * len, LAST = interval_info.second * len; \
        tatami_test::test_block_access(params, dense_mod.get(), ref.get(), FIRST, LAST); \
        tatami_test::test_block_access(params, sparse_mod.get(), ref.get(), FIRST, LAST); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        ArithVector, \
        name, \
        ::testing::Combine( \
            name::simulation_parameter_combinations(), \
            tatami_test::standard_test_access_parameter_combinations(), \
            ::testing::Values( \
                std::make_pair(0.0, 0.35), \
                std::make_pair(0.38, 0.61), \
                std::make_pair(0.777, 1.0) \
            ) \
        ) \
    );

#define ARITH_VECTOR_INDEX_TEST(name, base) \
    class name : \
        public ::testing::TestWithParam<std::tuple<typename base::SimulationParameters, tatami_test::StandardTestAccessParameters, std::pair<double, int> > >, \
        public base { \
    protected: \
        void SetUp() { \
            assemble(std::get<0>(GetParam())); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto tparam = GetParam(); \
        auto params = tatami_test::convert_access_parameters(std::get<1>(tparam)); \
        auto interval_info = std::get<2>(tparam); \
        auto len = (params.use_row ? ref->ncol() : ref->nrow()); \
        size_t FIRST = interval_info.first * len, STEP = interval_info.second; \
        tatami_test::test_indexed_access(params, dense_mod.get(), ref.get(), FIRST, STEP); \
        tatami_test::test_indexed_access(params, sparse_mod.get(), ref.get(), FIRST, STEP); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        ArithVector, \
        name, \
        ::testing::Combine( \
            name::simulation_parameter_combinations(), \
            tatami_test::standard_test_access_parameter_combinations(), \
            ::testing::Values( \
                std::make_pair(0.0, 7), \
                std::make_pair(0.21, 5), \
                std::make_pair(0.56, 3) \
            ) \
        ) \
    );

/****************************
 ********* ADDITION *********
 ****************************/

class ArithVectorAdditionUtils : public ArithVectorUtils {
public:
    typedef std::tuple<bool> SimulationParameters;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false) // add by row (or column).
        );
    }

    static void apply_operation(
        bool row,
        const std::vector<double>& vec,
        std::shared_ptr<tatami::NumericMatrix>& dense_ptr, 
        std::shared_ptr<tatami::NumericMatrix>& sparse_ptr) 
    {
        if (row) {
            auto op = tatami::make_DelayedAddVectorHelper<0>(vec);
            dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
            sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
        } else {
            auto op = tatami::make_DelayedAddVectorHelper<1>(vec);
            dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
            sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
        }
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    inline static SimulationParameters last_params;

    static void assemble(SimulationParameters sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        ArithVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, 5, 0.5);
        apply_operation(row, vec, dense_mod, sparse_mod);

        auto refvec = simulated;
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                refvec[r * ncol + c] += vec[row ? r : c];
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec)));
    }
};

ARITH_VECTOR_BASIC_SETUP(ArithVectorAdditionTest, ArithVectorAdditionUtils)
TEST_P(ArithVectorAdditionTest, Basic) {
    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense_mod->sparse_proportion(), 0);
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(sparse_mod->sparse_proportion(), 0);
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_EQ(dense_mod->prefer_rows_proportion(), 1);
    EXPECT_FALSE(sparse_mod->prefer_rows());
    EXPECT_EQ(sparse_mod->prefer_rows_proportion(), 0);
}

ARITH_VECTOR_FULL_TEST(ArithVectorAdditionFullTest, ArithVectorAdditionUtils)
ARITH_VECTOR_BLOCK_TEST(ArithVectorAdditionBlockTest, ArithVectorAdditionUtils)
ARITH_VECTOR_INDEX_TEST(ArithVectorAdditionIndexTest, ArithVectorAdditionUtils)

class ArithVectorAdditionZeroedTest : public ::testing::TestWithParam<typename ArithVectorAdditionUtils::SimulationParameters>, public ArithVectorUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(ArithVectorAdditionZeroedTest, Basic) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    ArithVectorAdditionUtils::apply_operation(row, zeroed, dense_z, sparse_z);

    EXPECT_FALSE(dense_z->sparse());
    EXPECT_TRUE(sparse_z->sparse());

    tatami_test::test_simple_column_access(dense_z.get(), dense.get());
    tatami_test::test_simple_column_access(sparse_z.get(), sparse.get()); 

    tatami_test::test_simple_row_access(dense_z.get(), dense.get());
    tatami_test::test_simple_row_access(sparse_z.get(), sparse.get());
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector, 
    ArithVectorAdditionZeroedTest, 
    ArithVectorAdditionUtils::simulation_parameter_combinations()
);

/*******************************
 ********* SUBTRACTION *********
 *******************************/

class ArithVectorSubtractionUtils : public ArithVectorUtils {
public:
    typedef std::tuple<bool, bool> SimulationParameters;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false), // by row or by column
            ::testing::Values(true, false)  // on the right or left
        );
    }

    static void apply_operation(
        bool row, 
        bool right, 
        const std::vector<double>& vec, 
        std::shared_ptr<tatami::NumericMatrix>& dense_ptr, 
        std::shared_ptr<tatami::NumericMatrix>& sparse_ptr)
    {
        if (row) {
            if (right) {
                auto op = tatami::make_DelayedSubtractVectorHelper<true, 0>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            } else {
                auto op = tatami::make_DelayedSubtractVectorHelper<false, 0>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            }
        } else {
            if (right) {
                auto op = tatami::make_DelayedSubtractVectorHelper<true, 1>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            } else {
                auto op = tatami::make_DelayedSubtractVectorHelper<false, 1>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            }
        }
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    inline static SimulationParameters last_params;

    static void assemble(SimulationParameters sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        ArithVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto right = std::get<1>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, -10, 2.222);
        apply_operation(row, right, vec, dense_mod, sparse_mod);

        auto refvec = simulated;
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                auto& x = refvec[r * ncol + c];
                if (right) {
                    x -= vec[row ? r : c];
                } else {
                    x = vec[row ? r : c] - x;
                }
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec)));
    }
};

ARITH_VECTOR_BASIC_SETUP(ArithVectorSubtractionTest, ArithVectorSubtractionUtils)
TEST_P(ArithVectorSubtractionTest, Basic) {
    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

ARITH_VECTOR_FULL_TEST(ArithVectorSubtractionFullTest, ArithVectorSubtractionUtils)
ARITH_VECTOR_BLOCK_TEST(ArithVectorSubtractionBlockTest, ArithVectorSubtractionUtils)
ARITH_VECTOR_INDEX_TEST(ArithVectorSubtractionIndexTest, ArithVectorSubtractionUtils)

class ArithVectorSubtractionZeroedTest : public ::testing::TestWithParam<typename ArithVectorSubtractionUtils::SimulationParameters>, public ArithVectorUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(ArithVectorSubtractionZeroedTest, Basic) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    ArithVectorSubtractionUtils::apply_operation(row, right, zeroed, dense_z, sparse_z);

    EXPECT_FALSE(dense_z->sparse());
    EXPECT_TRUE(sparse_z->sparse());

    if (right) {
        tatami_test::test_simple_column_access(dense_z.get(), dense.get());
        tatami_test::test_simple_column_access(sparse_z.get(), sparse.get()); 

        tatami_test::test_simple_row_access(dense_z.get(), dense.get());
        tatami_test::test_simple_row_access(sparse_z.get(), sparse.get());
    } else {
        auto copy = simulated;
        for (auto& x : copy) {
            x *= -1;
        }
        tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

        tatami_test::test_simple_column_access(dense_z.get(), &ref);
        tatami_test::test_simple_column_access(sparse_z.get(), &ref);

        tatami_test::test_simple_row_access(dense_z.get(), &ref);
        tatami_test::test_simple_row_access(sparse_z.get(), &ref);
    }
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector, 
    ArithVectorSubtractionZeroedTest, 
    ArithVectorSubtractionUtils::simulation_parameter_combinations()
);

/**********************************
 ********* MULTIPLICATION *********
 **********************************/

class ArithVectorMultiplicationUtils : public ArithVectorUtils {
public:
    typedef std::tuple<bool> SimulationParameters;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false) // by row or by column
        );
    }

    static void apply_operation(
        bool row, 
        const std::vector<double>& vec, 
        std::shared_ptr<tatami::NumericMatrix>& dense_ptr, 
        std::shared_ptr<tatami::NumericMatrix>& sparse_ptr)
    {
        if (row) {
            auto op = tatami::make_DelayedMultiplyVectorHelper<0>(vec);
            dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
            sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
        } else {
            auto op = tatami::make_DelayedMultiplyVectorHelper<1>(vec);
            dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
            sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
        }
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    inline static SimulationParameters last_params;

    static void assemble(SimulationParameters sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        ArithVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, 99, -2.5);
        apply_operation(row, vec, dense_mod, sparse_mod);

        auto refvec = simulated;
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                refvec[r * ncol + c] *= vec[row ? r : c];
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec)));
    }
};

ARITH_VECTOR_BASIC_SETUP(ArithVectorMultiplicationTest, ArithVectorMultiplicationUtils)
TEST_P(ArithVectorMultiplicationTest, Basic) {
    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense_mod->sparse_proportion(), 0);
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(sparse_mod->sparse_proportion(), 1);
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

ARITH_VECTOR_FULL_TEST(ArithVectorMultiplicationFullTest, ArithVectorMultiplicationUtils)
ARITH_VECTOR_BLOCK_TEST(ArithVectorMultiplicationBlockTest, ArithVectorMultiplicationUtils)
ARITH_VECTOR_INDEX_TEST(ArithVectorMultiplicationIndexTest, ArithVectorMultiplicationUtils)

class ArithVectorMultiplicationZeroedTest : public ::testing::TestWithParam<typename ArithVectorMultiplicationUtils::SimulationParameters>, public ArithVectorUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(ArithVectorMultiplicationZeroedTest, Basic) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    ArithVectorMultiplicationUtils::apply_operation(row, zeroed, dense_z, sparse_z);

    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::vector<double>(nrow * ncol));

    tatami_test::test_simple_column_access(dense_z.get(), &ref);
    tatami_test::test_simple_column_access(sparse_z.get(), &ref);

    tatami_test::test_simple_row_access(dense_z.get(), &ref);
    tatami_test::test_simple_row_access(sparse_z.get(), &ref);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector, 
    ArithVectorMultiplicationZeroedTest, 
    ArithVectorMultiplicationUtils::simulation_parameter_combinations()
);

/****************************
 ********* DIVISION *********
 ****************************/

class ArithVectorDivisionUtils : public ArithVectorUtils {
public:
    typedef std::tuple<bool, bool> SimulationParameters;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false), // by row or by column
            ::testing::Values(true, false)  // on the right or left
        );
    }

    static void apply_operation(
        bool row,
        bool right,
        const std::vector<double>& vec,
        std::shared_ptr<tatami::NumericMatrix>& dense_ptr,
        std::shared_ptr<tatami::NumericMatrix>& sparse_ptr)
    {
        if (row) {
            if (right) {
                auto op = tatami::make_DelayedDivideVectorHelper<true, 0>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            } else {
                auto op = tatami::make_DelayedDivideVectorHelper<false, 0>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            }
        } else {
            if (right) {
                auto op = tatami::make_DelayedDivideVectorHelper<true, 1>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            } else {
                auto op = tatami::make_DelayedDivideVectorHelper<false, 1>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            }
        }
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    inline static SimulationParameters last_params;

    static void assemble(SimulationParameters sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        ArithVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto right = std::get<1>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, -19, 2.11);
        apply_operation(row, right, vec, dense_mod, sparse_mod);

        auto refvec = simulated;
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                auto& x = refvec[r * ncol + c];
                auto val = vec[row ? r : c];
                if (right) {
                    x /= val;
                } else {
                    if (x) {
                        x = val / x;
                    } else if (val > 0) {
                        x = std::numeric_limits<double>::infinity();
                    } else {
                        x = -std::numeric_limits<double>::infinity();
                    }
                }
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec)));
    }
};

ARITH_VECTOR_BASIC_SETUP(ArithVectorDivisionTest, ArithVectorDivisionUtils)
TEST_P(ArithVectorDivisionTest, Basic) {
    EXPECT_FALSE(dense_mod->sparse());

    auto right = std::get<1>(last_params);
    if (right) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }

    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

ARITH_VECTOR_FULL_TEST(ArithVectorDivisionFullTest, ArithVectorDivisionUtils);
ARITH_VECTOR_BLOCK_TEST(ArithVectorDivisionBlockTest, ArithVectorDivisionUtils);
ARITH_VECTOR_INDEX_TEST(ArithVectorDivisionIndexTest, ArithVectorDivisionUtils);

class ArithVectorDivisionZeroedTest : public ::testing::TestWithParam<typename ArithVectorDivisionUtils::SimulationParameters>, public ArithVectorUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(ArithVectorDivisionZeroedTest, AllZero) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    ArithVectorDivisionUtils::apply_operation(row, right, zeroed, dense_z, sparse_z);

    auto copy = simulated;
    if (right) {
        for (auto& x : copy) {
            x = careful_division(x, 0.0);
        }
    } else {
        for (auto& x : copy) {
            x = careful_division(0.0, x);
        }
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_z->sparse());
    EXPECT_FALSE(sparse_z->sparse());

    // Turning on NaN protection.
    test_simple_column_access_wt_nan(dense_z.get(), &ref);
    test_simple_column_access_wt_nan(sparse_z.get(), &ref);

    test_simple_row_access_wt_nan(dense_z.get(), &ref);
    test_simple_row_access_wt_nan(sparse_z.get(), &ref);
}

TEST_P(ArithVectorDivisionZeroedTest, OneZero) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    // But actually, even 1 zero is enough to break sparsity.
    std::vector<double> solo_zero(row ? nrow : ncol, 1);
    solo_zero[0] = 0;
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    ArithVectorDivisionUtils::apply_operation(row, right, solo_zero, dense_z, sparse_z);

    auto copy = simulated;
    if (row) {
        if (right) {
            for (size_t c = 0; c < ncol; ++c) {
                copy[c] = careful_division(copy[c], 0.0);
            }
        } else {
            int it = 0;
            for (size_t r = 0; r < nrow; ++r) {
                for (size_t c = 0; c < ncol; ++c, ++it) {
                    copy[it] = careful_division((r == 0 ? 0.0 : 1.0), copy[it]);
                }
            }
        }
    } else {
        if (right) {
            for (size_t r = 0; r < nrow; ++r) {
                copy[r * ncol] = careful_division(copy[r * ncol], 0.0);
            }
        } else {
            int it = 0;
            for (size_t r = 0; r < nrow; ++r) {
                for (size_t c = 0; c < ncol; ++c, ++it) {
                    copy[it] = careful_division((c == 0 ? 0.0 : 1.0), copy[it]);
                }
            }
        }
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_z->sparse());
    EXPECT_FALSE(sparse_z->sparse());

    // Turning on NaN protection.
    test_simple_column_access_wt_nan(dense_z.get(), &ref);
    test_simple_column_access_wt_nan(sparse_z.get(), &ref);

    test_simple_row_access_wt_nan(dense_z.get(), &ref);
    test_simple_row_access_wt_nan(sparse_z.get(), &ref);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector, 
    ArithVectorDivisionZeroedTest, 
    ArithVectorDivisionUtils::simulation_parameter_combinations()
);

/****************************
 ********** POWER ***********
 ****************************/

class ArithVectorPowerUtils : public ArithVectorUtils {
public:
    typedef std::tuple<bool, bool> SimulationParameters;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false), // by row or by column
            ::testing::Values(true, false)  // on the right or left
        );
    }
    
    static void apply_operation(
        bool row, 
        bool right, 
        const std::vector<double>& vec,
        std::shared_ptr<tatami::NumericMatrix>& dense_ptr, 
        std::shared_ptr<tatami::NumericMatrix>& sparse_ptr)
    {
        tatami::DelayedAbsHelper op0;
        auto dense_tmp = tatami::make_DelayedUnaryIsometricOp(dense, op0);
        auto sparse_tmp = tatami::make_DelayedUnaryIsometricOp(sparse, op0);

        if (row) {
            if (right) {
                auto op = tatami::make_DelayedPowerVectorHelper<true, 0>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense_tmp, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse_tmp, op);
            } else {
                auto op = tatami::make_DelayedPowerVectorHelper<false, 0>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense_tmp, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse_tmp, op);
            }
        } else {
            if (right) {
                auto op = tatami::make_DelayedPowerVectorHelper<true, 1>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense_tmp, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse_tmp, op);
            } else {
                auto op = tatami::make_DelayedPowerVectorHelper<false, 1>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense_tmp, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse_tmp, op);
            }
        }
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    inline static SimulationParameters last_params;

    static void assemble(SimulationParameters sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        ArithVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto right = std::get<1>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, 0.01, 2.11);
        apply_operation(row, right, vec, dense_mod, sparse_mod);

        auto refvec = simulated;
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                auto& x = refvec[r * ncol + c];
                auto val = vec[row ? r : c];
                if (right) {
                    x = std::pow(std::abs(x), val);
                } else {
                    x = std::pow(val, std::abs(x));
                }
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec)));
    }
};

ARITH_VECTOR_BASIC_SETUP(ArithVectorPowerTest, ArithVectorPowerUtils)
TEST_P(ArithVectorPowerTest, Basic) {
    EXPECT_FALSE(dense_mod->sparse());

    auto right = std::get<1>(last_params);
    if (right) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }

    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

ARITH_VECTOR_FULL_TEST(ArithVectorPowerFullTest, ArithVectorPowerUtils);
ARITH_VECTOR_BLOCK_TEST(ArithVectorPowerBlockTest, ArithVectorPowerUtils);
ARITH_VECTOR_INDEX_TEST(ArithVectorPowerIndexTest, ArithVectorPowerUtils);

class ArithVectorPowerZeroedTest : public ::testing::TestWithParam<typename ArithVectorPowerUtils::SimulationParameters>, public ArithVectorUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(ArithVectorPowerZeroedTest, AllZero) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    ArithVectorPowerUtils::apply_operation(row, right, zeroed, dense_z, sparse_z);

    auto copy = simulated;
    if (right) {
        for (auto& x : copy) {
            x = std::pow(std::abs(x), 0.0);
        }
    } else {
        for (auto& x : copy) {
            x = std::pow(0.0, std::abs(x));
        }
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_z->sparse());
    EXPECT_FALSE(sparse_z->sparse());

    // Turning on NaN protection.
    test_simple_column_access_wt_nan(dense_z.get(), &ref);
    test_simple_column_access_wt_nan(sparse_z.get(), &ref);

    test_simple_row_access_wt_nan(dense_z.get(), &ref);
    test_simple_row_access_wt_nan(sparse_z.get(), &ref);
}

TEST_P(ArithVectorPowerZeroedTest, OneZero) {
    // But actually, even 1 zero is enough to break sparsity.
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    std::vector<double> solo_zero(row ? nrow : ncol, 1);
    solo_zero[0] = 0;
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    ArithVectorPowerUtils::apply_operation(row, right, solo_zero, dense_z, sparse_z);

    auto copy = simulated;
    for (auto& x : copy) {
        x = std::abs(x);
    }
    if (row) {
        if (right) {
            for (size_t c = 0; c < ncol; ++c) {
                copy[c] = std::pow(copy[c], 0.0);
            }
        } else {
            int it = 0;
            for (size_t r = 0; r < nrow; ++r) {
                for (size_t c = 0; c < ncol; ++c, ++it) {
                    copy[it] = std::pow((r == 0 ? 0.0 : 1.0), copy[it]);
                }
            }
        }
    } else {
        if (right) {
            for (size_t r = 0; r < nrow; ++r) {
                copy[r * ncol] = std::pow(copy[r * ncol], 0.0);
            }
        } else {
            int it = 0;
            for (size_t r = 0; r < nrow; ++r) {
                for (size_t c = 0; c < ncol; ++c, ++it) {
                    copy[it] = std::pow((c == 0 ? 0.0 : 1.0), copy[it]);
                }
            }
        }
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_z->sparse());
    EXPECT_FALSE(sparse_z->sparse());

    // Turning on NaN protection.
    test_simple_column_access_wt_nan(dense_z.get(), &ref);
    test_simple_column_access_wt_nan(sparse_z.get(), &ref);

    test_simple_row_access_wt_nan(dense_z.get(), &ref);
    test_simple_row_access_wt_nan(sparse_z.get(), &ref);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector, 
    ArithVectorPowerZeroedTest, 
    ArithVectorPowerUtils::simulation_parameter_combinations()
);

/****************************
 ********** MODULO **********
 ****************************/

#define ARITH_VECTOR_FULL_TEST_WITH_NAN(name, base) \
    class name : \
        public ::testing::TestWithParam<std::tuple<typename base::SimulationParameters, tatami_test::StandardTestAccessParameters> >, \
        public base { \
    protected: \
        void SetUp() { \
            assemble(std::get<0>(GetParam())); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto params = tatami_test::convert_access_parameters(std::get<1>(GetParam())); \
        params.has_nan = true; \
        tatami_test::test_full_access(params, dense_mod.get(), ref.get()); \
        tatami_test::test_full_access(params, sparse_mod.get(), ref.get()); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        ArithVector, \
        name, \
        ::testing::Combine( \
            name::simulation_parameter_combinations(), \
            tatami_test::standard_test_access_parameter_combinations() \
        ) \
    ); \

#define ARITH_VECTOR_BLOCK_TEST_WITH_NAN(name, base) \
    class name : \
        public ::testing::TestWithParam<std::tuple<typename base::SimulationParameters, tatami_test::StandardTestAccessParameters, std::pair<double, double> > >, \
        public base { \
    protected: \
        void SetUp() { \
            assemble(std::get<0>(GetParam())); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto tparam = GetParam(); \
        auto params = tatami_test::convert_access_parameters(std::get<1>(tparam)); \
        params.has_nan = true; \
        auto interval_info = std::get<2>(tparam); \
        auto len = (params.use_row ? ref->ncol() : ref->nrow()); \
        size_t FIRST = interval_info.first * len, LAST = interval_info.second * len; \
        tatami_test::test_block_access(params, dense_mod.get(), ref.get(), FIRST, LAST); \
        tatami_test::test_block_access(params, sparse_mod.get(), ref.get(), FIRST, LAST); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        ArithVector, \
        name, \
        ::testing::Combine( \
            name::simulation_parameter_combinations(), \
            tatami_test::standard_test_access_parameter_combinations(), \
            ::testing::Values( \
                std::make_pair(0.0, 0.35), \
                std::make_pair(0.38, 0.61), \
                std::make_pair(0.777, 1.0) \
            ) \
        ) \
    );

#define ARITH_VECTOR_INDEX_TEST_WITH_NAN(name, base) \
    class name : \
        public ::testing::TestWithParam<std::tuple<typename base::SimulationParameters, tatami_test::StandardTestAccessParameters, std::pair<double, int> > >, \
        public base { \
    protected: \
        void SetUp() { \
            assemble(std::get<0>(GetParam())); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto tparam = GetParam(); \
        auto params = tatami_test::convert_access_parameters(std::get<1>(tparam)); \
        params.has_nan = true; \
        auto interval_info = std::get<2>(tparam); \
        auto len = (params.use_row ? ref->ncol() : ref->nrow()); \
        size_t FIRST = interval_info.first * len, STEP = interval_info.second; \
        tatami_test::test_indexed_access(params, dense_mod.get(), ref.get(), FIRST, STEP); \
        tatami_test::test_indexed_access(params, sparse_mod.get(), ref.get(), FIRST, STEP); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        ArithVector, \
        name, \
        ::testing::Combine( \
            name::simulation_parameter_combinations(), \
            tatami_test::standard_test_access_parameter_combinations(), \
            ::testing::Values( \
                std::make_pair(0.0, 7), \
                std::make_pair(0.21, 5), \
                std::make_pair(0.56, 3) \
            ) \
        ) \
    );


class ArithVectorModuloUtils : public ArithVectorUtils {
public:
    typedef std::tuple<bool, bool> SimulationParameters;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false), // by row or by column
            ::testing::Values(true, false)  // on the right or left
        );
    }

    static void apply_operation(
        bool row,
        bool right,
        const std::vector<double>& vec,
        std::shared_ptr<tatami::NumericMatrix>& dense_ptr, 
        std::shared_ptr<tatami::NumericMatrix>& sparse_ptr)
    {
        if (row) {
            if (right) {
                auto op = tatami::make_DelayedModuloVectorHelper<true, 0>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            } else {
                auto op = tatami::make_DelayedModuloVectorHelper<false, 0>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            }
        } else {
            if (right) {
                auto op = tatami::make_DelayedModuloVectorHelper<true, 1>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            } else {
                auto op = tatami::make_DelayedModuloVectorHelper<false, 1>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            }
        }
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    inline static SimulationParameters last_params;

    static void assemble(SimulationParameters sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        ArithVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto right = std::get<1>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, -19, 2.11);
        apply_operation(row, right, vec, dense_mod, sparse_mod);

        auto refvec = simulated;
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                auto& x = refvec[r * ncol + c];
                auto val = vec[row ? r : c];
                if (right) {
                    x = std::fmod(x, val);
                } else {
                    x = std::fmod(val, x);
                }
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec)));
    }
};

ARITH_VECTOR_BASIC_SETUP(ArithVectorModuloTest, ArithVectorModuloUtils);
TEST_P(ArithVectorModuloTest, Basic) {
    EXPECT_FALSE(dense_mod->sparse());

    auto right = std::get<1>(last_params);
    if (right) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }

    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

ARITH_VECTOR_FULL_TEST_WITH_NAN(ArithVectorModuloFullTest, ArithVectorModuloUtils);
ARITH_VECTOR_BLOCK_TEST_WITH_NAN(ArithVectorModuloBlockTest, ArithVectorModuloUtils);
ARITH_VECTOR_INDEX_TEST_WITH_NAN(ArithVectorModuloIndexTest, ArithVectorModuloUtils);

class ArithVectorModuloZeroedTest : public ::testing::TestWithParam<typename ArithVectorModuloUtils::SimulationParameters>, public ArithVectorUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(ArithVectorModuloZeroedTest, AllZero) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    ArithVectorModuloUtils::apply_operation(row, right, zeroed, dense_z, sparse_z);

    auto copy = simulated;
    if (right) {
        for (auto& x : copy) {
            x = std::fmod(x, 0.0);
        }
    } else {
        for (auto& x : copy) {
            x = std::fmod(0.0, x);
        }
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_z->sparse());
    EXPECT_FALSE(sparse_z->sparse());

    // Turning on NaN protection.
    test_simple_column_access_wt_nan(dense_z.get(), &ref);
    test_simple_column_access_wt_nan(sparse_z.get(), &ref);

    test_simple_row_access_wt_nan(dense_z.get(), &ref);
    test_simple_row_access_wt_nan(sparse_z.get(), &ref);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector, 
    ArithVectorModuloZeroedTest, 
    ArithVectorModuloUtils::simulation_parameter_combinations()
);

/****************************
 ***** INTEGER DIVISION *****
 ****************************/

class ArithVectorIntegerDivisionUtils : public ArithVectorUtils {
public:
    typedef std::tuple<bool, bool> SimulationParameters;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false), // by row or by column
            ::testing::Values(true, false)  // on the right or left
        );
    }

    static void apply_operation(
        bool row,
        bool right,
        const std::vector<double>& vec,
        std::shared_ptr<tatami::NumericMatrix>& dense_ptr,
        std::shared_ptr<tatami::NumericMatrix>& sparse_ptr)
    {
        if (row) {
            if (right) {
                auto op = tatami::make_DelayedIntegerDivideVectorHelper<true, 0>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            } else {
                auto op = tatami::make_DelayedIntegerDivideVectorHelper<false, 0>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            }
        } else {
            if (right) {
                auto op = tatami::make_DelayedIntegerDivideVectorHelper<true, 1>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            } else {
                auto op = tatami::make_DelayedIntegerDivideVectorHelper<false, 1>(vec);
                dense_ptr = tatami::make_DelayedUnaryIsometricOp(dense, op);
                sparse_ptr = tatami::make_DelayedUnaryIsometricOp(sparse, op);
            }
        }
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    inline static SimulationParameters last_params;

    void assemble(SimulationParameters sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        ArithVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto right = std::get<1>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, -19, 2.11);
        apply_operation(row, right, vec, dense_mod, sparse_mod);

        auto refvec = simulated;
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                auto& x = refvec[r * ncol + c];
                auto val = vec[row ? r : c];
                // x == (x %% y) + y * (x %/% y)
                if (right) {
                    x = std::floor(careful_division(x, val));
                } else {
                    x = std::floor(careful_division(val, x));
                }
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(refvec)));
    }
};

ARITH_VECTOR_BASIC_SETUP(ArithVectorIntegerDivisionTest, ArithVectorIntegerDivisionUtils);
TEST_P(ArithVectorIntegerDivisionTest, Basic) {
    EXPECT_FALSE(dense_mod->sparse());

    auto right = std::get<1>(last_params);
    if (right) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }

    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

ARITH_VECTOR_FULL_TEST_WITH_NAN(ArithVectorIntegerDivisionFullTest, ArithVectorIntegerDivisionUtils);
ARITH_VECTOR_BLOCK_TEST_WITH_NAN(ArithVectorIntegerDivisionBlockTest, ArithVectorIntegerDivisionUtils);
ARITH_VECTOR_INDEX_TEST_WITH_NAN(ArithVectorIntegerDivisionIndexTest, ArithVectorIntegerDivisionUtils);

class ArithVectorIntegerDivisionZeroedTest : public ::testing::TestWithParam<typename ArithVectorIntegerDivisionUtils::SimulationParameters>, public ArithVectorUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(ArithVectorIntegerDivisionZeroedTest, AllZero) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    ArithVectorIntegerDivisionUtils::apply_operation(row, right, zeroed, dense_z, sparse_z);

    auto copy = simulated;
    if (right) {
        for (auto& x : copy) {
            // x == (x %% y) + y * (x %/% y)
            x = std::floor(careful_division(x, 0.0));
        }
    } else {
        for (auto& x : copy) {
            x = std::floor(careful_division(0.0, x));
        }
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_z->sparse());
    EXPECT_FALSE(sparse_z->sparse());

    // Turning on NaN protection.
    test_simple_column_access_wt_nan(dense_z.get(), &ref);
    test_simple_column_access_wt_nan(sparse_z.get(), &ref);

    test_simple_row_access_wt_nan(dense_z.get(), &ref);
    test_simple_row_access_wt_nan(sparse_z.get(), &ref);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorIntegerDivisionZeroedTest,
    ArithVectorIntegerDivisionUtils::simulation_parameter_combinations()
);
