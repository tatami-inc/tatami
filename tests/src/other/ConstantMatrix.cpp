#include <gtest/gtest.h>
#include "tatami/other/ConstantMatrix.hpp"
#include "tatami/dense/DenseMatrix.hpp"
#include "tatami_test/tatami_test.hpp"

class ConstantMatrixUtils {
public:
    typedef std::tuple<double> SimulationParameters;

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> mat, ref;
    inline static SimulationParameters last_params;

    static void assemble(const SimulationParameters& params) {
        if (ref && params == last_params) {
            return;
        }
        last_params = params;

        int NR = 20, NC = 50;
        auto val = std::get<0>(params);
        mat.reset(new tatami::ConstantMatrix<double, int>(NR, NC, val));
        ref.reset(new tatami::DenseRowMatrix<double, int>(NR, NC, std::vector<double>(NR * NC, val)));
    }
};

/**********************************
 **********************************/

class ConstantMatrixUtilsTest: public ::testing::Test, public ConstantMatrixUtils {};

TEST_F(ConstantMatrixUtilsTest, Simple) {
    assemble(std::make_tuple(1.0));
    EXPECT_EQ(mat->nrow(), 20);
    EXPECT_EQ(mat->ncol(), 50);
    EXPECT_FALSE(mat->is_sparse());
    EXPECT_EQ(mat->is_sparse_proportion(), 0);
    EXPECT_TRUE(mat->prefer_rows());
    EXPECT_EQ(mat->prefer_rows_proportion(), 1);
}

TEST_F(ConstantMatrixUtilsTest, Zeroed) {
    assemble(std::make_tuple(0.0));
    EXPECT_EQ(mat->nrow(), 20);
    EXPECT_EQ(mat->ncol(), 50);
    EXPECT_TRUE(mat->is_sparse());
    EXPECT_EQ(mat->is_sparse_proportion(), 1);
}

/**********************************
 **********************************/

class ConstantMatrixFullTest : 
    public ::testing::TestWithParam<std::tuple<ConstantMatrixUtils::SimulationParameters, tatami_test::StandardTestAccessOptions> >,
    public ConstantMatrixUtils
{
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(ConstantMatrixFullTest, Basic) {
    auto opt = tatami_test::convert_test_access_options(std::get<1>(GetParam()));
    test_full_access(*mat, *ref, opt);
}

INSTANTIATE_TEST_SUITE_P(
    ConstantMatrix,
    ConstantMatrixFullTest,
    ::testing::Combine(
        ::testing::Combine(
            ::testing::Values(0, 1.0)
        ),
        tatami_test::standard_test_access_options_combinations()
    )
);

/**********************************
 **********************************/

class ConstantMatrixBlockTest :
    public ::testing::TestWithParam<std::tuple<ConstantMatrixUtils::SimulationParameters, tatami_test::StandardTestAccessOptions, std::pair<double, double> > >,
    public ConstantMatrixUtils {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(ConstantMatrixBlockTest, Basic) {
    auto tparam = GetParam();
    auto opt = tatami_test::convert_test_access_options(std::get<1>(tparam));
    auto block = std::get<2>(tparam);
    test_block_access(*mat, *ref, block.first, block.second, opt);
}

INSTANTIATE_TEST_SUITE_P(
    ConstantMatrix,
    ConstantMatrixBlockTest,
    ::testing::Combine(
        ::testing::Combine(
            ::testing::Values(0, 1.0)
        ),
        tatami_test::standard_test_access_options_combinations(),
        ::testing::Values(
            std::pair(0.0, 0.45),
            std::pair(0.2, 0.6),
            std::pair(0.4, 0.53)
        )
    )
);

/**********************************
 **********************************/

class ConstantMatrixIndexTest :
    public ::testing::TestWithParam<std::tuple<ConstantMatrixUtils::SimulationParameters, tatami_test::StandardTestAccessOptions, std::pair<double, double> > >,
    public ConstantMatrixUtils {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(ConstantMatrixIndexTest, Basic) {
    auto tparam = GetParam();
    auto opts = tatami_test::convert_test_access_options(std::get<1>(tparam));
    auto block = std::get<2>(tparam);
    test_indexed_access(*mat, *ref, block.first, block.second, opts);
}

INSTANTIATE_TEST_SUITE_P(
    ConstantMatrix,
    ConstantMatrixIndexTest,
    ::testing::Combine(
        ::testing::Combine(
            ::testing::Values(0, 1.0)
        ),
        tatami_test::standard_test_access_options_combinations(),
        ::testing::Values(
            std::pair(0.0, 0.2),
            std::pair(0.21, 0.3),
            std::pair(0.42, 0.5) 
        )
    )
);
