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
    public ::testing::TestWithParam<std::tuple<ConstantMatrixUtils::SimulationParameters, tatami_test::StandardTestAccessParameters> >,
    public ConstantMatrixUtils
{
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(ConstantMatrixFullTest, Basic) {
    auto tparams = tatami_test::convert_access_parameters(std::get<1>(GetParam()));
    test_full_access(tparams, mat.get(), ref.get());
}

INSTANTIATE_TEST_SUITE_P(
    ConstantMatrix,
    ConstantMatrixFullTest,
    ::testing::Combine(
        ::testing::Combine(
            ::testing::Values(0, 1.0)
        ),
        tatami_test::standard_test_access_parameter_combinations()
    )
);

/**********************************
 **********************************/

class ConstantMatrixBlockTest :
    public ::testing::TestWithParam<std::tuple<ConstantMatrixUtils::SimulationParameters, tatami_test::StandardTestAccessParameters, std::pair<double, double> > >,
    public ConstantMatrixUtils {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(ConstantMatrixBlockTest, Basic) {
    auto param = GetParam();
    auto tparams = tatami_test::convert_access_parameters(std::get<1>(param));

    auto block = std::get<2>(param);
    int len = (tparams.use_row ? ref->ncol() : ref->nrow());
    int block_start = len * block.first, block_end = len * block.second;

    test_block_access(tparams, mat.get(), ref.get(), block_start, block_end);
}

INSTANTIATE_TEST_SUITE_P(
    ConstantMatrix,
    ConstantMatrixBlockTest,
    ::testing::Combine(
        ::testing::Combine(
            ::testing::Values(0, 1.0)
        ),
        tatami_test::standard_test_access_parameter_combinations(),
        ::testing::Values(
            std::pair(0.0, 0.45),
            std::pair(0.2, 0.8),
            std::pair(0.4, 0.93)
        )
    )
);

/**********************************
 **********************************/

class ConstantMatrixIndexTest :
    public ::testing::TestWithParam<std::tuple<ConstantMatrixUtils::SimulationParameters, tatami_test::StandardTestAccessParameters, std::pair<double, int> > >,
    public ConstantMatrixUtils {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(ConstantMatrixIndexTest, Basic) {
    auto param = GetParam();
    auto tparams = tatami_test::convert_access_parameters(std::get<1>(param));

    auto block = std::get<2>(param);
    int len = (tparams.use_row ? ref->ncol() : ref->nrow());
    int index_start = len * block.first, index_step = block.second;

    test_indexed_access(tparams, mat.get(), ref.get(), index_start, index_step);
}

INSTANTIATE_TEST_SUITE_P(
    ConstantMatrix,
    ConstantMatrixIndexTest,
    ::testing::Combine(
        ::testing::Combine(
            ::testing::Values(0, 1.0)
        ),
        tatami_test::standard_test_access_parameter_combinations(),
        ::testing::Values(
            std::pair(0.0, 4),
            std::pair(0.21, 3),
            std::pair(0.42, 2) 
        )
    )
);
