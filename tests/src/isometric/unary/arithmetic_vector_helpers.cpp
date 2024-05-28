#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOperation.hpp"
#include "tatami/isometric/unary/arithmetic_helpers.hpp"
#include "tatami/isometric/unary/math_helpers.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class DelayedUnaryIsometricArithmeticVectorUtils {
protected:
    inline static size_t nrow = 291, ncol = 188;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    inline static std::vector<double> simulated;

    static void assemble() {
        if (dense) {
            return;
        }
        simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1);
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulated));
        sparse = tatami::convert_to_compressed_sparse<false, double, int>(dense.get()); // column major.
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
        DelayedUnaryIsometricArithmeticVector, \
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
        DelayedUnaryIsometricArithmeticVector, \
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
        DelayedUnaryIsometricArithmeticVector, \
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
        DelayedUnaryIsometricArithmeticVector, \
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

class DelayedUnaryIsometricAddVectorUtils : public DelayedUnaryIsometricArithmeticVectorUtils {
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
        auto op = tatami::make_DelayedUnaryIsometricAddVector(vec, row);
        dense_ptr = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_ptr = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    inline static SimulationParameters last_params;

    static void assemble(SimulationParameters sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        DelayedUnaryIsometricArithmeticVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, 5, 0.5);
        apply_operation(row, vec, dense_mod, sparse_mod);

        auto refvec = simulated;
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                refvec[r * ncol + c] += vec[row ? r : c];
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(refvec)));
    }
};

ARITH_VECTOR_BASIC_SETUP(DelayedUnaryIsometricAddVectorTest, DelayedUnaryIsometricAddVectorUtils)
TEST_P(DelayedUnaryIsometricAddVectorTest, Basic) {
    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense_mod->is_sparse_proportion(), 0);
    EXPECT_FALSE(sparse_mod->is_sparse());
    EXPECT_EQ(sparse_mod->is_sparse_proportion(), 0);
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_EQ(dense_mod->prefer_rows_proportion(), 1);
    EXPECT_FALSE(sparse_mod->prefer_rows());
    EXPECT_EQ(sparse_mod->prefer_rows_proportion(), 0);
}

ARITH_VECTOR_FULL_TEST(DelayedUnaryIsometricAddVectorFullTest, DelayedUnaryIsometricAddVectorUtils)
ARITH_VECTOR_BLOCK_TEST(DelayedUnaryIsometricAddVectorBlockTest, DelayedUnaryIsometricAddVectorUtils)
ARITH_VECTOR_INDEX_TEST(DelayedUnaryIsometricAddVectorIndexTest, DelayedUnaryIsometricAddVectorUtils)

class DelayedUnaryIsometricAddVectorZeroedTest : public ::testing::TestWithParam<typename DelayedUnaryIsometricAddVectorUtils::SimulationParameters>, public DelayedUnaryIsometricArithmeticVectorUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedUnaryIsometricAddVectorZeroedTest, Basic) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    DelayedUnaryIsometricAddVectorUtils::apply_operation(row, zeroed, dense_z, sparse_z);

    EXPECT_FALSE(dense_z->is_sparse());
    EXPECT_TRUE(sparse_z->is_sparse());

    tatami_test::test_simple_column_access(dense_z.get(), dense.get());
    tatami_test::test_simple_column_access(sparse_z.get(), sparse.get()); 

    tatami_test::test_simple_row_access(dense_z.get(), dense.get());
    tatami_test::test_simple_row_access(sparse_z.get(), sparse.get());
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricArithmeticVector, 
    DelayedUnaryIsometricAddVectorZeroedTest, 
    DelayedUnaryIsometricAddVectorUtils::simulation_parameter_combinations()
);

class DelayedUnaryIsometricAddVectorNewTypeTest : 
    public ::testing::TestWithParam<typename DelayedUnaryIsometricAddVectorUtils::SimulationParameters>, 
    public DelayedUnaryIsometricArithmeticVectorUtils 
{
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedUnaryIsometricAddVectorNewTypeTest, Basic) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);

    auto vec = create_vector(row ? nrow : ncol, 5, 0.5);
    auto op = tatami::make_DelayedUnaryIsometricAddVector(vec, row);
    auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
    auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

    std::vector<float> frefvec(nrow * ncol);
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            auto offset = r * ncol + c;
            frefvec[offset] = simulated[offset] + vec[row ? r : c];
        }
    }
    tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

    quick_test_all(dense_fmod.get(), &fref);
    quick_test_all(sparse_fmod.get(), &fref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricArithmeticVector, 
    DelayedUnaryIsometricAddVectorNewTypeTest, 
    DelayedUnaryIsometricAddVectorUtils::simulation_parameter_combinations()
);

/*******************************
 ********* SUBTRACTION *********
 *******************************/

class DelayedUnaryIsometricSubtractVectorUtils : public DelayedUnaryIsometricArithmeticVectorUtils {
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
        if (right) {
            auto op = tatami::make_DelayedUnaryIsometricSubtractVector<true>(vec, row);
            dense_ptr = tatami::make_DelayedUnaryIsometricOperation(dense, op);
            sparse_ptr = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
        } else {
            auto op = tatami::make_DelayedUnaryIsometricSubtractVector<false>(vec, row);
            dense_ptr = tatami::make_DelayedUnaryIsometricOperation(dense, op);
            sparse_ptr = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
        }
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    inline static SimulationParameters last_params;

    static void assemble(SimulationParameters sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        DelayedUnaryIsometricArithmeticVectorUtils::assemble();

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
        ref.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(refvec)));
    }
};

ARITH_VECTOR_BASIC_SETUP(DelayedUnaryIsometricSubtractVectorTest, DelayedUnaryIsometricSubtractVectorUtils)
TEST_P(DelayedUnaryIsometricSubtractVectorTest, Basic) {
    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

ARITH_VECTOR_FULL_TEST(DelayedUnaryIsometricSubtractVectorFullTest, DelayedUnaryIsometricSubtractVectorUtils)
ARITH_VECTOR_BLOCK_TEST(DelayedUnaryIsometricSubtractVectorBlockTest, DelayedUnaryIsometricSubtractVectorUtils)
ARITH_VECTOR_INDEX_TEST(DelayedUnaryIsometricSubtractVectorIndexTest, DelayedUnaryIsometricSubtractVectorUtils)

class DelayedUnaryIsometricSubtractVectorZeroedTest : public ::testing::TestWithParam<typename DelayedUnaryIsometricSubtractVectorUtils::SimulationParameters>, public DelayedUnaryIsometricArithmeticVectorUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedUnaryIsometricSubtractVectorZeroedTest, Basic) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    DelayedUnaryIsometricSubtractVectorUtils::apply_operation(row, right, zeroed, dense_z, sparse_z);

    EXPECT_FALSE(dense_z->is_sparse());
    EXPECT_TRUE(sparse_z->is_sparse());

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
        tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(copy));

        tatami_test::test_simple_column_access(dense_z.get(), &ref);
        tatami_test::test_simple_column_access(sparse_z.get(), &ref);

        tatami_test::test_simple_row_access(dense_z.get(), &ref);
        tatami_test::test_simple_row_access(sparse_z.get(), &ref);
    }
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricArithmeticVector, 
    DelayedUnaryIsometricSubtractVectorZeroedTest, 
    DelayedUnaryIsometricSubtractVectorUtils::simulation_parameter_combinations()
);

/**********************************
 ********* MULTIPLICATION *********
 **********************************/

class DelayedUnaryIsometricMultiplyVectorUtils : public DelayedUnaryIsometricArithmeticVectorUtils {
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
        auto op = tatami::make_DelayedUnaryIsometricMultiplyVector(vec, row);
        dense_ptr = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_ptr = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    inline static SimulationParameters last_params;

    static void assemble(SimulationParameters sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        DelayedUnaryIsometricArithmeticVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, 99, -2.5);
        apply_operation(row, vec, dense_mod, sparse_mod);

        auto refvec = simulated;
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                refvec[r * ncol + c] *= vec[row ? r : c];
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(refvec)));
    }
};

ARITH_VECTOR_BASIC_SETUP(DelayedUnaryIsometricMultiplyVectorTest, DelayedUnaryIsometricMultiplyVectorUtils)
TEST_P(DelayedUnaryIsometricMultiplyVectorTest, Basic) {
    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense_mod->is_sparse_proportion(), 0);
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(sparse_mod->is_sparse_proportion(), 1);
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

ARITH_VECTOR_FULL_TEST(DelayedUnaryIsometricMultiplyVectorFullTest, DelayedUnaryIsometricMultiplyVectorUtils)
ARITH_VECTOR_BLOCK_TEST(DelayedUnaryIsometricMultiplyVectorBlockTest, DelayedUnaryIsometricMultiplyVectorUtils)
ARITH_VECTOR_INDEX_TEST(DelayedUnaryIsometricMultiplyVectorIndexTest, DelayedUnaryIsometricMultiplyVectorUtils)

class DelayedUnaryIsometricMultiplyVectorZeroedTest : public ::testing::TestWithParam<typename DelayedUnaryIsometricMultiplyVectorUtils::SimulationParameters>, public DelayedUnaryIsometricArithmeticVectorUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedUnaryIsometricMultiplyVectorZeroedTest, Basic) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    DelayedUnaryIsometricMultiplyVectorUtils::apply_operation(row, zeroed, dense_z, sparse_z);

    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::vector<double>(nrow * ncol));

    tatami_test::test_simple_column_access(dense_z.get(), &ref);
    tatami_test::test_simple_column_access(sparse_z.get(), &ref);

    tatami_test::test_simple_row_access(dense_z.get(), &ref);
    tatami_test::test_simple_row_access(sparse_z.get(), &ref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricArithmeticVector, 
    DelayedUnaryIsometricMultiplyVectorZeroedTest, 
    DelayedUnaryIsometricMultiplyVectorUtils::simulation_parameter_combinations()
);

/****************************
 ********* DIVISION *********
 ****************************/

class DelayedUnaryIsometricDivideVectorUtils : public DelayedUnaryIsometricArithmeticVectorUtils {
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
        if (right) {
            auto op = tatami::make_DelayedUnaryIsometricDivideVector<true>(vec, row);
            dense_ptr = tatami::make_DelayedUnaryIsometricOperation(dense, op);
            sparse_ptr = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
        } else {
            auto op = tatami::make_DelayedUnaryIsometricDivideVector<false>(vec, row);
            dense_ptr = tatami::make_DelayedUnaryIsometricOperation(dense, op);
            sparse_ptr = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
        }
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    inline static SimulationParameters last_params;

    static void assemble(SimulationParameters sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        DelayedUnaryIsometricArithmeticVectorUtils::assemble();

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
        ref.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(refvec)));
    }
};

ARITH_VECTOR_BASIC_SETUP(DelayedUnaryIsometricDivideVectorTest, DelayedUnaryIsometricDivideVectorUtils)
TEST_P(DelayedUnaryIsometricDivideVectorTest, Basic) {
    EXPECT_FALSE(dense_mod->is_sparse());

    auto right = std::get<1>(last_params);
    if (right) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

ARITH_VECTOR_FULL_TEST(DelayedUnaryIsometricDivideVectorFullTest, DelayedUnaryIsometricDivideVectorUtils);
ARITH_VECTOR_BLOCK_TEST(DelayedUnaryIsometricDivideVectorBlockTest, DelayedUnaryIsometricDivideVectorUtils);
ARITH_VECTOR_INDEX_TEST(DelayedUnaryIsometricDivideVectorIndexTest, DelayedUnaryIsometricDivideVectorUtils);

class DelayedUnaryIsometricDivideVectorZeroedTest : public ::testing::TestWithParam<typename DelayedUnaryIsometricDivideVectorUtils::SimulationParameters>, public DelayedUnaryIsometricArithmeticVectorUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedUnaryIsometricDivideVectorZeroedTest, AllZero) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    DelayedUnaryIsometricDivideVectorUtils::apply_operation(row, right, zeroed, dense_z, sparse_z);

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
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_z->is_sparse());
    EXPECT_FALSE(sparse_z->is_sparse());

    // Turning on NaN protection.
    test_simple_column_access_wt_nan(dense_z.get(), &ref);
    test_simple_column_access_wt_nan(sparse_z.get(), &ref);

    test_simple_row_access_wt_nan(dense_z.get(), &ref);
    test_simple_row_access_wt_nan(sparse_z.get(), &ref);
}

TEST_P(DelayedUnaryIsometricDivideVectorZeroedTest, OneZero) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    // But actually, even 1 zero is enough to break sparsity.
    std::vector<double> solo_zero(row ? nrow : ncol, 1);
    solo_zero[0] = 0;
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    DelayedUnaryIsometricDivideVectorUtils::apply_operation(row, right, solo_zero, dense_z, sparse_z);

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
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_z->is_sparse());
    EXPECT_FALSE(sparse_z->is_sparse());

    // Turning on NaN protection.
    test_simple_column_access_wt_nan(dense_z.get(), &ref);
    test_simple_column_access_wt_nan(sparse_z.get(), &ref);

    test_simple_row_access_wt_nan(dense_z.get(), &ref);
    test_simple_row_access_wt_nan(sparse_z.get(), &ref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricArithmeticVector, 
    DelayedUnaryIsometricDivideVectorZeroedTest, 
    DelayedUnaryIsometricDivideVectorUtils::simulation_parameter_combinations()
);

/****************************
 ********** POWER ***********
 ****************************/

class DelayedUnaryIsometricPowerVectorUtils : public DelayedUnaryIsometricArithmeticVectorUtils {
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
        tatami::DelayedUnaryIsometricAbs op0;
        auto dense_tmp = tatami::make_DelayedUnaryIsometricOperation(dense, op0);
        auto sparse_tmp = tatami::make_DelayedUnaryIsometricOperation(sparse, op0);

        if (right) {
            auto op = tatami::make_DelayedUnaryIsometricPowerVector<true>(vec, row);
            dense_ptr = tatami::make_DelayedUnaryIsometricOperation(dense_tmp, op);
            sparse_ptr = tatami::make_DelayedUnaryIsometricOperation(sparse_tmp, op);
        } else {
            auto op = tatami::make_DelayedUnaryIsometricPowerVector<false>(vec, row);
            dense_ptr = tatami::make_DelayedUnaryIsometricOperation(dense_tmp, op);
            sparse_ptr = tatami::make_DelayedUnaryIsometricOperation(sparse_tmp, op);
        }
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    inline static SimulationParameters last_params;

    static void assemble(SimulationParameters sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        DelayedUnaryIsometricArithmeticVectorUtils::assemble();

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
        ref.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(refvec)));
    }
};

ARITH_VECTOR_BASIC_SETUP(DelayedUnaryIsometricPowerVectorTest, DelayedUnaryIsometricPowerVectorUtils)
TEST_P(DelayedUnaryIsometricPowerVectorTest, Basic) {
    EXPECT_FALSE(dense_mod->is_sparse());

    auto right = std::get<1>(last_params);
    if (right) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

ARITH_VECTOR_FULL_TEST(DelayedUnaryIsometricPowerVectorFullTest, DelayedUnaryIsometricPowerVectorUtils);
ARITH_VECTOR_BLOCK_TEST(DelayedUnaryIsometricPowerVectorBlockTest, DelayedUnaryIsometricPowerVectorUtils);
ARITH_VECTOR_INDEX_TEST(DelayedUnaryIsometricPowerVectorIndexTest, DelayedUnaryIsometricPowerVectorUtils);

class DelayedUnaryIsometricPowerVectorZeroedTest : public ::testing::TestWithParam<typename DelayedUnaryIsometricPowerVectorUtils::SimulationParameters>, public DelayedUnaryIsometricArithmeticVectorUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedUnaryIsometricPowerVectorZeroedTest, AllZero) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    DelayedUnaryIsometricPowerVectorUtils::apply_operation(row, right, zeroed, dense_z, sparse_z);

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
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_z->is_sparse());
    EXPECT_FALSE(sparse_z->is_sparse());

    // Turning on NaN protection.
    test_simple_column_access_wt_nan(dense_z.get(), &ref);
    test_simple_column_access_wt_nan(sparse_z.get(), &ref);

    test_simple_row_access_wt_nan(dense_z.get(), &ref);
    test_simple_row_access_wt_nan(sparse_z.get(), &ref);
}

TEST_P(DelayedUnaryIsometricPowerVectorZeroedTest, OneZero) {
    // But actually, even 1 zero is enough to break sparsity.
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    std::vector<double> solo_zero(row ? nrow : ncol, 1);
    solo_zero[0] = 0;
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    DelayedUnaryIsometricPowerVectorUtils::apply_operation(row, right, solo_zero, dense_z, sparse_z);

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
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_z->is_sparse());
    EXPECT_FALSE(sparse_z->is_sparse());

    // Turning on NaN protection.
    test_simple_column_access_wt_nan(dense_z.get(), &ref);
    test_simple_column_access_wt_nan(sparse_z.get(), &ref);

    test_simple_row_access_wt_nan(dense_z.get(), &ref);
    test_simple_row_access_wt_nan(sparse_z.get(), &ref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricArithmeticVector, 
    DelayedUnaryIsometricPowerVectorZeroedTest, 
    DelayedUnaryIsometricPowerVectorUtils::simulation_parameter_combinations()
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
        DelayedUnaryIsometricArithmeticVector, \
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
        DelayedUnaryIsometricArithmeticVector, \
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
        DelayedUnaryIsometricArithmeticVector, \
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


class DelayedUnaryIsometricModuloVectorUtils : public DelayedUnaryIsometricArithmeticVectorUtils {
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
        if (right) {
            auto op = tatami::make_DelayedUnaryIsometricModuloVector<true>(vec, row);
            dense_ptr = tatami::make_DelayedUnaryIsometricOperation(dense, op);
            sparse_ptr = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
        } else {
            auto op = tatami::make_DelayedUnaryIsometricModuloVector<false>(vec, row);
            dense_ptr = tatami::make_DelayedUnaryIsometricOperation(dense, op);
            sparse_ptr = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
        }
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    inline static SimulationParameters last_params;

    static void assemble(SimulationParameters sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        DelayedUnaryIsometricArithmeticVectorUtils::assemble();

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
        ref.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(refvec)));
    }
};

ARITH_VECTOR_BASIC_SETUP(DelayedUnaryIsometricModuloVectorTest, DelayedUnaryIsometricModuloVectorUtils);
TEST_P(DelayedUnaryIsometricModuloVectorTest, Basic) {
    EXPECT_FALSE(dense_mod->is_sparse());

    auto right = std::get<1>(last_params);
    if (right) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

ARITH_VECTOR_FULL_TEST_WITH_NAN(DelayedUnaryIsometricModuloVectorFullTest, DelayedUnaryIsometricModuloVectorUtils);
ARITH_VECTOR_BLOCK_TEST_WITH_NAN(DelayedUnaryIsometricModuloVectorBlockTest, DelayedUnaryIsometricModuloVectorUtils);
ARITH_VECTOR_INDEX_TEST_WITH_NAN(DelayedUnaryIsometricModuloVectorIndexTest, DelayedUnaryIsometricModuloVectorUtils);

class DelayedUnaryIsometricModuloVectorZeroedTest : public ::testing::TestWithParam<typename DelayedUnaryIsometricModuloVectorUtils::SimulationParameters>, public DelayedUnaryIsometricArithmeticVectorUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedUnaryIsometricModuloVectorZeroedTest, AllZero) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    DelayedUnaryIsometricModuloVectorUtils::apply_operation(row, right, zeroed, dense_z, sparse_z);

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
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_z->is_sparse());
    EXPECT_FALSE(sparse_z->is_sparse());

    // Turning on NaN protection.
    test_simple_column_access_wt_nan(dense_z.get(), &ref);
    test_simple_column_access_wt_nan(sparse_z.get(), &ref);

    test_simple_row_access_wt_nan(dense_z.get(), &ref);
    test_simple_row_access_wt_nan(sparse_z.get(), &ref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricArithmeticVector, 
    DelayedUnaryIsometricModuloVectorZeroedTest, 
    DelayedUnaryIsometricModuloVectorUtils::simulation_parameter_combinations()
);

class DelayedUnaryIsometricModuloVectorNewTypeTest : 
    public ::testing::TestWithParam<typename DelayedUnaryIsometricModuloVectorUtils::SimulationParameters>, 
    public DelayedUnaryIsometricArithmeticVectorUtils 
{
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedUnaryIsometricModuloVectorNewTypeTest, Basic) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    auto vec = create_vector(row ? nrow : ncol, 5, 0.5);
    std::shared_ptr<tatami::Matrix<float, int> > dense_fmod, sparse_fmod;
    if (right) {
        auto op = tatami::make_DelayedUnaryIsometricModuloVector<true>(vec, row);
        dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);
    } else {
        auto op = tatami::make_DelayedUnaryIsometricModuloVector<false>(vec, row);
        dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);
    }

    std::vector<float> frefvec(nrow * ncol);
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            auto offset = r * ncol + c;
            auto val = vec[row ? r : c];
            if (right) {
                frefvec[offset] = std::fmod(simulated[offset], val);
            } else {
                frefvec[offset] = std::fmod(val, simulated[offset]);
            }

        }
    }
    tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

    quick_test_all(dense_fmod.get(), &fref, /* has_nan = */ true);
    quick_test_all(sparse_fmod.get(), &fref, /* has_nan = */ true);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricArithmeticVector, 
    DelayedUnaryIsometricModuloVectorNewTypeTest, 
    DelayedUnaryIsometricModuloVectorUtils::simulation_parameter_combinations()
);

/****************************
 ***** INTEGER DIVISION *****
 ****************************/

class DelayedUnaryIsometricIntegerDivideVectorUtils : public DelayedUnaryIsometricArithmeticVectorUtils {
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
        if (right) {
            auto op = tatami::make_DelayedUnaryIsometricIntegerDivideVector<true>(vec, row);
            dense_ptr = tatami::make_DelayedUnaryIsometricOperation(dense, op);
            sparse_ptr = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
        } else {
            auto op = tatami::make_DelayedUnaryIsometricIntegerDivideVector<false>(vec, row);
            dense_ptr = tatami::make_DelayedUnaryIsometricOperation(dense, op);
            sparse_ptr = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
        }
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    inline static SimulationParameters last_params;

    void assemble(SimulationParameters sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        DelayedUnaryIsometricArithmeticVectorUtils::assemble();

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
        ref.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(refvec)));
    }
};

ARITH_VECTOR_BASIC_SETUP(DelayedUnaryIsometricIntegerDivideVectorTest, DelayedUnaryIsometricIntegerDivideVectorUtils);
TEST_P(DelayedUnaryIsometricIntegerDivideVectorTest, Basic) {
    EXPECT_FALSE(dense_mod->is_sparse());

    auto right = std::get<1>(last_params);
    if (right) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());
}

ARITH_VECTOR_FULL_TEST_WITH_NAN(DelayedUnaryIsometricIntegerDivideVectorFullTest, DelayedUnaryIsometricIntegerDivideVectorUtils);
ARITH_VECTOR_BLOCK_TEST_WITH_NAN(DelayedUnaryIsometricIntegerDivideVectorBlockTest, DelayedUnaryIsometricIntegerDivideVectorUtils);
ARITH_VECTOR_INDEX_TEST_WITH_NAN(DelayedUnaryIsometricIntegerDivideVectorIndexTest, DelayedUnaryIsometricIntegerDivideVectorUtils);

class DelayedUnaryIsometricIntegerDivideVectorZeroedTest : public ::testing::TestWithParam<typename DelayedUnaryIsometricIntegerDivideVectorUtils::SimulationParameters>, public DelayedUnaryIsometricArithmeticVectorUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedUnaryIsometricIntegerDivideVectorZeroedTest, AllZero) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    std::shared_ptr<tatami::NumericMatrix> dense_z, sparse_z;
    DelayedUnaryIsometricIntegerDivideVectorUtils::apply_operation(row, right, zeroed, dense_z, sparse_z);

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
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_z->is_sparse());
    EXPECT_FALSE(sparse_z->is_sparse());

    // Turning on NaN protection.
    test_simple_column_access_wt_nan(dense_z.get(), &ref);
    test_simple_column_access_wt_nan(sparse_z.get(), &ref);

    test_simple_row_access_wt_nan(dense_z.get(), &ref);
    test_simple_row_access_wt_nan(sparse_z.get(), &ref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricArithmeticVector,
    DelayedUnaryIsometricIntegerDivideVectorZeroedTest,
    DelayedUnaryIsometricIntegerDivideVectorUtils::simulation_parameter_combinations()
);
