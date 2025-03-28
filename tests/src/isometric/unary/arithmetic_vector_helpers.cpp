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

        simulated = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.1;
            opt.seed = 42428481111;
            return opt;
        }());

        dense.reset(new tatami::DenseMatrix<double, int, decltype(simulated)>(nrow, ncol, simulated, true)); // row major
        sparse = tatami::convert_to_compressed_sparse<double, int>(*dense, false, {}); // column major.
        return;
    }

    static std::vector<double> create_vector(size_t n, double starter = 0, double jump = 1) { 
        std::vector<double> output(n, starter);
        for (size_t i = 1; i < n; ++i) {
            output[i] = output[i-1] + jump;
        }
        return output;
    }
};

#define ARITH_VECTOR_CONFIGURE(y) \
    if (ref && y == last_params) { \
        return; \
    } \
    last_params = y;

#define ARITH_VECTOR_BASIC_SETUP(name, base) \
    class name : \
        public ::testing::TestWithParam<typename base::SimulationOptions>, \
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
        public ::testing::TestWithParam<std::tuple<typename base::SimulationOptions, tatami_test::StandardTestAccessOptions> >, \
        public base { \
    protected: \
        void SetUp() { \
            assemble(std::get<0>(GetParam())); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto opts = tatami_test::convert_test_access_options(std::get<1>(GetParam())); \
        tatami_test::test_full_access<double, int>(*dense_mod, *ref, opts); \
        tatami_test::test_full_access<double, int>(*sparse_mod, *ref, opts); \
        tatami_test::test_unsorted_full_access<double, int>(*sparse_uns, opts); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        DelayedUnaryIsometricArithmeticVector, \
        name, \
        ::testing::Combine( \
            name::simulation_parameter_combinations(), \
            tatami_test::standard_test_access_options_combinations() \
        ) \
    );

#define ARITH_VECTOR_BLOCK_TEST(name, base) \
    class name : \
        public ::testing::TestWithParam<std::tuple<typename base::SimulationOptions, tatami_test::StandardTestAccessOptions, std::pair<double, double> > >, \
        public base { \
    protected: \
        void SetUp() { \
            assemble(std::get<0>(GetParam())); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto tparam = GetParam(); \
        auto opts = tatami_test::convert_test_access_options(std::get<1>(tparam)); \
        auto interval_info = std::get<2>(tparam); \
        tatami_test::test_block_access<double, int>(*dense_mod, *ref, interval_info.first, interval_info.second, opts); \
        tatami_test::test_block_access<double, int>(*sparse_mod, *ref, interval_info.first, interval_info.second, opts); \
        tatami_test::test_unsorted_block_access<double, int>(*sparse_uns, interval_info.first, interval_info.second, opts); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        DelayedUnaryIsometricArithmeticVector, \
        name, \
        ::testing::Combine( \
            name::simulation_parameter_combinations(), \
            tatami_test::standard_test_access_options_combinations(), \
            ::testing::Values( \
                std::make_pair(0.0, 0.35), \
                std::make_pair(0.38, 0.23), \
                std::make_pair(0.777, 0.223) \
            ) \
        ) \
    );

#define ARITH_VECTOR_INDEX_TEST(name, base) \
    class name : \
        public ::testing::TestWithParam<std::tuple<typename base::SimulationOptions, tatami_test::StandardTestAccessOptions, std::pair<double, double> > >, \
        public base { \
    protected: \
        void SetUp() { \
            assemble(std::get<0>(GetParam())); \
        } \
    }; \
    \
    TEST_P(name, Basic) { \
        auto tparam = GetParam(); \
        auto opts = tatami_test::convert_test_access_options(std::get<1>(tparam)); \
        auto interval_info = std::get<2>(tparam); \
        tatami_test::test_indexed_access<double, int>(*dense_mod, *ref, interval_info.first, interval_info.second, opts); \
        tatami_test::test_indexed_access<double, int>(*sparse_mod, *ref, interval_info.first, interval_info.second, opts); \
        tatami_test::test_unsorted_indexed_access<double, int>(*sparse_uns, interval_info.first, interval_info.second, opts); \
    } \
    \
    INSTANTIATE_TEST_SUITE_P( \
        DelayedUnaryIsometricArithmeticVector, \
        name, \
        ::testing::Combine( \
            name::simulation_parameter_combinations(), \
            tatami_test::standard_test_access_options_combinations(), \
            ::testing::Values( \
                std::make_pair(0.0, 0.15), \
                std::make_pair(0.21, 0.2), \
                std::make_pair(0.56, 0.3) \
            ) \
        ) \
    );

/****************************
 ********* ADDITION *********
 ****************************/

class DelayedUnaryIsometricAddVectorUtils : public DelayedUnaryIsometricArithmeticVectorUtils {
public:
    typedef std::tuple<bool> SimulationOptions;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false) // add by row (or column).
        );
    }

    static std::shared_ptr<tatami::NumericMatrix> apply_operation(bool row, std::vector<double> vec, std::shared_ptr<tatami::NumericMatrix> source) {
        auto op = std::make_shared<tatami::DelayedUnaryIsometricAddVectorHelper<double, double, int, decltype(vec)> >(std::move(vec), row);
        return std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(std::move(source), std::move(op));
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;
    inline static SimulationOptions last_params;

    static void assemble(SimulationOptions sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        DelayedUnaryIsometricArithmeticVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, 5, 0.5);

        dense_mod = apply_operation(row, vec, dense);
        sparse_mod = apply_operation(row, vec, sparse);
        sparse_uns = apply_operation(row, vec, std::make_shared<tatami_test::ReversedIndicesWrapper<double, int> >(sparse));

        auto refvec = simulated;
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                refvec[r * ncol + c] += vec[row ? r : c];
            }
        }
        ref.reset(new tatami::DenseMatrix<double, int, decltype(refvec)>(nrow, ncol, std::move(refvec), true));
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

class DelayedUnaryIsometricAddVectorZeroedTest : 
    public ::testing::TestWithParam<typename DelayedUnaryIsometricAddVectorUtils::SimulationOptions>,
    public DelayedUnaryIsometricArithmeticVectorUtils
{
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedUnaryIsometricAddVectorZeroedTest, Basic) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    auto dense_z = DelayedUnaryIsometricAddVectorUtils::apply_operation(row, zeroed, dense);
    auto sparse_z = DelayedUnaryIsometricAddVectorUtils::apply_operation(row, zeroed, sparse);

    EXPECT_FALSE(dense_z->is_sparse());
    EXPECT_TRUE(sparse_z->is_sparse());

    tatami_test::test_simple_column_access<double, int>(*dense_z, *dense);
    tatami_test::test_simple_column_access<double, int>(*sparse_z, *sparse); 

    tatami_test::test_simple_row_access<double, int>(*dense_z, *dense);
    tatami_test::test_simple_row_access<double, int>(*sparse_z, *sparse);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricArithmeticVector, 
    DelayedUnaryIsometricAddVectorZeroedTest, 
    DelayedUnaryIsometricAddVectorUtils::simulation_parameter_combinations()
);

class DelayedUnaryIsometricAddVectorNewTypeTest : 
    public ::testing::TestWithParam<typename DelayedUnaryIsometricAddVectorUtils::SimulationOptions>, 
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
    auto op = std::make_shared<tatami::DelayedUnaryIsometricAddVectorHelper<float, double, int, decltype(vec)> >(vec, row);
    tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, op);
    tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, op);

    std::vector<float> frefvec(nrow * ncol);
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            auto offset = r * ncol + c;
            frefvec[offset] = simulated[offset] + vec[row ? r : c];
        }
    }
    tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

    quick_test_all<float, int>(dense_fmod, fref);
    quick_test_all<float, int>(sparse_fmod, fref);
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
    typedef std::tuple<bool, bool> SimulationOptions;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false), // by row or by column
            ::testing::Values(true, false)  // on the right or left
        );
    }

    static std::shared_ptr<tatami::NumericMatrix> apply_operation(bool row, bool right, std::vector<double> vec, std::shared_ptr<tatami::NumericMatrix> source) {
        std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > op;
        if (right) {
            op.reset(new tatami::DelayedUnaryIsometricSubtractVectorHelper<true, double, double, int, decltype(vec)>(std::move(vec), row));
        } else {
            op.reset(new tatami::DelayedUnaryIsometricSubtractVectorHelper<false, double, double, int, decltype(vec)>(std::move(vec), row));
        }
        return std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(std::move(source), std::move(op));
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;
    inline static SimulationOptions last_params;

    static void assemble(SimulationOptions sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        DelayedUnaryIsometricArithmeticVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto right = std::get<1>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, -10, 2.222);

        dense_mod = apply_operation(row, right, vec, dense);
        sparse_mod = apply_operation(row, right, vec, sparse);
        sparse_uns = apply_operation(row, right, vec, std::make_shared<tatami_test::ReversedIndicesWrapper<double, int> >(sparse));

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
        ref.reset(new tatami::DenseMatrix<double, int, decltype(refvec)>(nrow, ncol, std::move(refvec), true));
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

class DelayedUnaryIsometricSubtractVectorZeroedTest : 
    public ::testing::TestWithParam<typename DelayedUnaryIsometricSubtractVectorUtils::SimulationOptions>,
    public DelayedUnaryIsometricArithmeticVectorUtils 
{
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
    auto dense_z = DelayedUnaryIsometricSubtractVectorUtils::apply_operation(row, right, zeroed, dense);
    auto sparse_z = DelayedUnaryIsometricSubtractVectorUtils::apply_operation(row, right, zeroed, sparse);

    EXPECT_FALSE(dense_z->is_sparse());
    EXPECT_TRUE(sparse_z->is_sparse());

    if (right) {
        tatami_test::test_simple_column_access<double, int>(*dense_z, *dense);
        tatami_test::test_simple_column_access<double, int>(*sparse_z, *sparse); 

        tatami_test::test_simple_row_access<double, int>(*dense_z, *dense);
        tatami_test::test_simple_row_access<double, int>(*sparse_z, *sparse);
    } else {
        auto copy = simulated;
        for (auto& x : copy) {
            x *= -1;
        }
        tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(copy));

        tatami_test::test_simple_column_access<double, int>(*dense_z, ref);
        tatami_test::test_simple_column_access<double, int>(*sparse_z, ref);

        tatami_test::test_simple_row_access<double, int>(*dense_z, ref);
        tatami_test::test_simple_row_access<double, int>(*sparse_z, ref);
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
    typedef std::tuple<bool> SimulationOptions;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false) // by row or by column
        );
    }

    static std::shared_ptr<tatami::NumericMatrix> apply_operation(bool row, std::vector<double> vec, std::shared_ptr<tatami::NumericMatrix> source) {
        auto op = std::make_shared<tatami::DelayedUnaryIsometricMultiplyVectorHelper<double, double, int, decltype(vec)> >(std::move(vec), row);
        return std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(std::move(source), std::move(op));
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;
    inline static SimulationOptions last_params;

    static void assemble(SimulationOptions sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        DelayedUnaryIsometricArithmeticVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, 99, -2.5);

        dense_mod = apply_operation(row, vec, dense);
        sparse_mod = apply_operation(row, vec, sparse);
        sparse_uns = apply_operation(row, vec, std::make_shared<tatami_test::ReversedIndicesWrapper<double, int> >(sparse));

        auto refvec = simulated;
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                refvec[r * ncol + c] *= vec[row ? r : c];
            }
        }
        ref.reset(new tatami::DenseMatrix<double, int, decltype(refvec)>(nrow, ncol, std::move(refvec), true));
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

class DelayedUnaryIsometricMultiplyVectorZeroedTest : 
    public ::testing::TestWithParam<typename DelayedUnaryIsometricMultiplyVectorUtils::SimulationOptions>, 
    public DelayedUnaryIsometricArithmeticVectorUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedUnaryIsometricMultiplyVectorZeroedTest, Basic) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);

    std::vector<double> zeroed(row ? nrow : ncol);
    auto dense_z = DelayedUnaryIsometricMultiplyVectorUtils::apply_operation(row, zeroed, dense);
    auto sparse_z = DelayedUnaryIsometricMultiplyVectorUtils::apply_operation(row, zeroed, sparse);

    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::vector<double>(nrow * ncol));

    tatami_test::test_simple_column_access<double, int>(*dense_z, ref);
    tatami_test::test_simple_column_access<double, int>(*sparse_z, ref);

    tatami_test::test_simple_row_access<double, int>(*dense_z, ref);
    tatami_test::test_simple_row_access<double, int>(*sparse_z, ref);
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
    typedef std::tuple<bool, bool> SimulationOptions;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false), // by row or by column
            ::testing::Values(true, false)  // on the right or left
        );
    }

    static std::shared_ptr<tatami::NumericMatrix> apply_operation(bool row, bool right, std::vector<double> vec, std::shared_ptr<tatami::NumericMatrix> source) {
        std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > op;
        if (right) {
            op.reset(new tatami::DelayedUnaryIsometricDivideVectorHelper<true, double, double, int, decltype(vec)>(std::move(vec), row));
        } else {
            op.reset(new tatami::DelayedUnaryIsometricDivideVectorHelper<false, double, double, int, decltype(vec)>(std::move(vec), row));
        }
        return std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(std::move(source), std::move(op));
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;
    inline static SimulationOptions last_params;

    static void assemble(SimulationOptions sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        DelayedUnaryIsometricArithmeticVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto right = std::get<1>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, -19, 2.11);

        dense_mod = apply_operation(row, right, vec, dense);
        sparse_mod = apply_operation(row, right, vec, sparse);
        sparse_uns = apply_operation(row, right, vec, std::make_shared<tatami_test::ReversedIndicesWrapper<double, int> >(sparse));

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

class DelayedUnaryIsometricDivideVectorZeroedTest : 
    public ::testing::TestWithParam<typename DelayedUnaryIsometricDivideVectorUtils::SimulationOptions>,
    public DelayedUnaryIsometricArithmeticVectorUtils {
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
    auto dense_z = DelayedUnaryIsometricDivideVectorUtils::apply_operation(row, right, zeroed, dense);
    auto sparse_z = DelayedUnaryIsometricDivideVectorUtils::apply_operation(row, right, zeroed, sparse);

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

    tatami_test::test_simple_column_access<double, int>(*dense_z, ref);
    tatami_test::test_simple_column_access<double, int>(*sparse_z, ref);

    tatami_test::test_simple_row_access<double, int>(*dense_z, ref);
    tatami_test::test_simple_row_access<double, int>(*sparse_z, ref);
}

TEST_P(DelayedUnaryIsometricDivideVectorZeroedTest, OneZero) {
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    // But actually, even 1 zero is enough to break sparsity.
    std::vector<double> solo_zero(row ? nrow : ncol, 1);
    solo_zero[0] = 0;
    auto dense_z = DelayedUnaryIsometricDivideVectorUtils::apply_operation(row, right, solo_zero, dense);
    auto sparse_z = DelayedUnaryIsometricDivideVectorUtils::apply_operation(row, right, solo_zero, sparse);

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

    tatami_test::test_simple_column_access<double, int>(*dense_z, ref);
    tatami_test::test_simple_column_access<double, int>(*sparse_z, ref);

    tatami_test::test_simple_row_access<double, int>(*dense_z, ref);
    tatami_test::test_simple_row_access<double, int>(*sparse_z, ref);
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
    typedef std::tuple<bool, bool> SimulationOptions;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false), // by row or by column
            ::testing::Values(true, false)  // on the right or left
        );
    }

    static std::shared_ptr<tatami::NumericMatrix> apply_operation(bool row, bool right, std::vector<double> vec, std::shared_ptr<tatami::NumericMatrix> source) {
        std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > op;
        if (right) {
            op.reset(new tatami::DelayedUnaryIsometricPowerVectorHelper<true, double, double, int, decltype(vec)>(std::move(vec), row));
        } else {
            op.reset(new tatami::DelayedUnaryIsometricPowerVectorHelper<false, double, double, int, decltype(vec)>(std::move(vec), row));
        }
        return std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(std::move(source), std::move(op));
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;
    inline static SimulationOptions last_params;

    static void assemble(SimulationOptions sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        DelayedUnaryIsometricArithmeticVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto right = std::get<1>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, 0.01, 2.11);

        auto asimulated = simulated; 
        for (auto& a : asimulated) {
            a = std::abs(a);
        }
        auto adense = std::make_shared<tatami::DenseMatrix<double, int, decltype(asimulated)> >(nrow, ncol, asimulated, true);
        auto asparse = tatami::convert_to_compressed_sparse<double, int>(*adense, false, {});

        dense_mod = apply_operation(row, right, vec, adense);
        sparse_mod = apply_operation(row, right, vec, asparse);
        sparse_uns = apply_operation(row, right, vec, std::make_shared<tatami_test::ReversedIndicesWrapper<double, int> >(asparse));

        auto refvec = asimulated;
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                auto& x = refvec[r * ncol + c];
                auto val = vec[row ? r : c];
                if (right) {
                    x = std::pow(x, val);
                } else {
                    x = std::pow(val, x);
                }
            }
        }
        ref.reset(new tatami::DenseMatrix<double, int, decltype(refvec)>(nrow, ncol, std::move(refvec), true));
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

class DelayedUnaryIsometricPowerVectorZeroedTest : 
    public ::testing::TestWithParam<typename DelayedUnaryIsometricPowerVectorUtils::SimulationOptions>,
    public DelayedUnaryIsometricArithmeticVectorUtils {
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
    auto dense_z = DelayedUnaryIsometricPowerVectorUtils::apply_operation(row, right, zeroed, dense);
    auto sparse_z = DelayedUnaryIsometricPowerVectorUtils::apply_operation(row, right, zeroed, sparse);

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

    tatami_test::test_simple_column_access<double, int>(*dense_z, ref);
    tatami_test::test_simple_column_access<double, int>(*sparse_z, ref);

    tatami_test::test_simple_row_access<double, int>(*dense_z, ref);
    tatami_test::test_simple_row_access<double, int>(*sparse_z, ref);
}

TEST_P(DelayedUnaryIsometricPowerVectorZeroedTest, OneZero) {
    // But actually, even 1 zero is enough to break sparsity.
    auto sim_params = GetParam();
    auto row = std::get<0>(sim_params);
    auto right = std::get<1>(sim_params);

    std::vector<double> solo_zero(row ? nrow : ncol, 1);
    solo_zero[0] = 0;
    auto dense_z = DelayedUnaryIsometricPowerVectorUtils::apply_operation(row, right, solo_zero, dense);
    auto sparse_z = DelayedUnaryIsometricPowerVectorUtils::apply_operation(row, right, solo_zero, sparse);

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

    tatami_test::test_simple_column_access<double, int>(*dense_z, ref);
    tatami_test::test_simple_column_access<double, int>(*sparse_z, ref);

    tatami_test::test_simple_row_access<double, int>(*dense_z, ref);
    tatami_test::test_simple_row_access<double, int>(*sparse_z, ref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricArithmeticVector, 
    DelayedUnaryIsometricPowerVectorZeroedTest, 
    DelayedUnaryIsometricPowerVectorUtils::simulation_parameter_combinations()
);

/****************************
 ********** MODULO **********
 ****************************/

class DelayedUnaryIsometricModuloVectorUtils : public DelayedUnaryIsometricArithmeticVectorUtils {
public:
    typedef std::tuple<bool, bool> SimulationOptions;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false), // by row or by column
            ::testing::Values(true, false)  // on the right or left
        );
    }

    static std::shared_ptr<tatami::NumericMatrix> apply_operation(bool row, bool right, std::vector<double> vec, std::shared_ptr<tatami::NumericMatrix> source) {
        std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > op;
        if (right) {
            op.reset(new tatami::DelayedUnaryIsometricModuloVectorHelper<true, double, double, int, decltype(vec)>(std::move(vec), row));
        } else {
            op.reset(new tatami::DelayedUnaryIsometricModuloVectorHelper<false, double, double, int, decltype(vec)>(std::move(vec), row));
        }
        return std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(std::move(source), std::move(op));
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;
    inline static SimulationOptions last_params;

    static void assemble(SimulationOptions sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        DelayedUnaryIsometricArithmeticVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto right = std::get<1>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, -19, 2.11);

        dense_mod = apply_operation(row, right, vec, dense);
        sparse_mod = apply_operation(row, right, vec, sparse);
        sparse_uns = apply_operation(row, right, vec, std::make_shared<tatami_test::ReversedIndicesWrapper<double, int> >(sparse));

        auto refvec = simulated;
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                auto& x = refvec[r * ncol + c];
                auto val = vec[row ? r : c];
                if (right) {
                    x = careful_modulo(x, val);
                } else {
                    x = careful_modulo(val, x);
                }
            }
        }
        ref.reset(new tatami::DenseMatrix<double, int, decltype(refvec)>(nrow, ncol, std::move(refvec), true));
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

ARITH_VECTOR_FULL_TEST(DelayedUnaryIsometricModuloVectorFullTest, DelayedUnaryIsometricModuloVectorUtils);
ARITH_VECTOR_BLOCK_TEST(DelayedUnaryIsometricModuloVectorBlockTest, DelayedUnaryIsometricModuloVectorUtils);
ARITH_VECTOR_INDEX_TEST(DelayedUnaryIsometricModuloVectorIndexTest, DelayedUnaryIsometricModuloVectorUtils);

class DelayedUnaryIsometricModuloVectorZeroedTest : 
    public ::testing::TestWithParam<typename DelayedUnaryIsometricModuloVectorUtils::SimulationOptions>,
    public DelayedUnaryIsometricArithmeticVectorUtils {
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
    auto dense_z = DelayedUnaryIsometricModuloVectorUtils::apply_operation(row, right, zeroed, dense);
    auto sparse_z = DelayedUnaryIsometricModuloVectorUtils::apply_operation(row, right, zeroed, sparse);

    auto copy = simulated;
    if (right) {
        for (auto& x : copy) {
            x = careful_modulo(x, 0.0);
        }
    } else {
        for (auto& x : copy) {
            x = careful_modulo(0.0, x);
        }
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_z->is_sparse());
    EXPECT_FALSE(sparse_z->is_sparse());

    tatami_test::test_simple_column_access<double, int>(*dense_z, ref);
    tatami_test::test_simple_column_access<double, int>(*sparse_z, ref);

    tatami_test::test_simple_row_access<double, int>(*dense_z, ref);
    tatami_test::test_simple_row_access<double, int>(*sparse_z, ref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricArithmeticVector, 
    DelayedUnaryIsometricModuloVectorZeroedTest, 
    DelayedUnaryIsometricModuloVectorUtils::simulation_parameter_combinations()
);

class DelayedUnaryIsometricModuloVectorNewTypeTest : 
    public ::testing::TestWithParam<typename DelayedUnaryIsometricModuloVectorUtils::SimulationOptions>, 
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
    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<float, double, int> > op;
    if (right) {
        op.reset(new tatami::DelayedUnaryIsometricModuloVectorHelper<true, float, double, int, decltype(vec)>(vec, row));
    } else {
        op.reset(new tatami::DelayedUnaryIsometricModuloVectorHelper<false, float, double, int, decltype(vec)>(vec, row));
    }
    tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, op);
    tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, op);

    std::vector<float> frefvec(nrow * ncol);
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t c = 0; c < ncol; ++c) {
            auto offset = r * ncol + c;
            auto val = vec[row ? r : c];
            if (right) {
                frefvec[offset] = careful_modulo(simulated[offset], val);
            } else {
                frefvec[offset] = careful_modulo(val, simulated[offset]);
            }

        }
    }
    tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

    quick_test_all<float, int>(dense_fmod, fref);
    quick_test_all<float, int>(sparse_fmod, fref);
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
    typedef std::tuple<bool, bool> SimulationOptions;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false), // by row or by column
            ::testing::Values(true, false)  // on the right or left
        );
    }

    static std::shared_ptr<tatami::NumericMatrix> apply_operation(bool row, bool right, const std::vector<double>& vec, std::shared_ptr<tatami::NumericMatrix> source) {
        if (right) {
            auto op = tatami::make_DelayedUnaryIsometricIntegerDivideVector<true>(vec, row);
            return tatami::make_DelayedUnaryIsometricOperation(std::move(source), std::move(op));
        } else {
            auto op = tatami::make_DelayedUnaryIsometricIntegerDivideVector<false>(vec, row);
            return tatami::make_DelayedUnaryIsometricOperation(std::move(source), std::move(op));
        }
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, sparse_uns, ref;
    inline static SimulationOptions last_params;

    void assemble(SimulationOptions sim_params) {
        ARITH_VECTOR_CONFIGURE(sim_params);
        DelayedUnaryIsometricArithmeticVectorUtils::assemble();

        auto row = std::get<0>(sim_params);
        auto right = std::get<1>(sim_params);
        auto vec = create_vector(row ? nrow : ncol, -19, 2.11);

        dense_mod = apply_operation(row, right, vec, dense);
        sparse_mod = apply_operation(row, right, vec, sparse);
        sparse_uns = apply_operation(row, right, vec, std::make_shared<tatami_test::ReversedIndicesWrapper<double, int> >(sparse));

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

ARITH_VECTOR_FULL_TEST(DelayedUnaryIsometricIntegerDivideVectorFullTest, DelayedUnaryIsometricIntegerDivideVectorUtils);
ARITH_VECTOR_BLOCK_TEST(DelayedUnaryIsometricIntegerDivideVectorBlockTest, DelayedUnaryIsometricIntegerDivideVectorUtils);
ARITH_VECTOR_INDEX_TEST(DelayedUnaryIsometricIntegerDivideVectorIndexTest, DelayedUnaryIsometricIntegerDivideVectorUtils);

class DelayedUnaryIsometricIntegerDivideVectorZeroedTest : 
    public ::testing::TestWithParam<typename DelayedUnaryIsometricIntegerDivideVectorUtils::SimulationOptions>,
    public DelayedUnaryIsometricArithmeticVectorUtils {
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
    auto dense_z = DelayedUnaryIsometricIntegerDivideVectorUtils::apply_operation(row, right, zeroed, dense);
    auto sparse_z = DelayedUnaryIsometricIntegerDivideVectorUtils::apply_operation(row, right, zeroed, sparse);

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

    tatami_test::test_simple_column_access<double, int>(*dense_z, ref);
    tatami_test::test_simple_column_access<double, int>(*sparse_z, ref);

    tatami_test::test_simple_row_access<double, int>(*dense_z, ref);
    tatami_test::test_simple_row_access<double, int>(*sparse_z, ref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricArithmeticVector,
    DelayedUnaryIsometricIntegerDivideVectorZeroedTest,
    DelayedUnaryIsometricIntegerDivideVectorUtils::simulation_parameter_combinations()
);

TEST(DelayedUnaryIsometricArithmeticVector, BackCompatibility) {
    std::vector<double> vec(10, 1);
    auto add = tatami::make_DelayedUnaryIsometricAddVector(vec, true);
    EXPECT_FALSE(add->is_sparse());
    auto sub = tatami::make_DelayedUnaryIsometricSubtractVector<true>(vec, true);
    EXPECT_FALSE(sub->is_sparse());
    auto mult = tatami::make_DelayedUnaryIsometricMultiplyVector(vec, true);
    EXPECT_TRUE(mult->is_sparse());
    auto div = tatami::make_DelayedUnaryIsometricDivideVector<true>(vec, true);
    EXPECT_TRUE(div->is_sparse());
    auto mod = tatami::make_DelayedUnaryIsometricModuloVector<true>(vec, true);
    EXPECT_TRUE(mod->is_sparse());
    auto pow = tatami::make_DelayedUnaryIsometricPowerVector<true>(vec, true);
    EXPECT_TRUE(pow->is_sparse());
    auto idiv = tatami::make_DelayedUnaryIsometricIntegerDivideVector<true>(vec, true);
    EXPECT_TRUE(idiv->is_sparse());
}
