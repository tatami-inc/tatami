#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/ForcedDense.hpp"
#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

class ForcedDenseUtils {
protected:
    inline static int nrow = 80, ncol = 100;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, forced;

    static void assemble() {
        if (dense) {
            return;
        }

        auto vec = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.2;
            opt.seed = 59;
            return opt;
        }());

        dense.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(vec)));
        forced.reset(new tatami::ForcedDense<double, int>(tatami::convert_to_compressed_sparse<double, int>(*dense, true, {})));
    }
};

class ForcedDenseTest : public ::testing::Test, public ForcedDenseUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_F(ForcedDenseTest, Basic) {
    EXPECT_EQ(forced->ncol(), ncol);
    EXPECT_EQ(forced->nrow(), nrow);

    EXPECT_FALSE(dense->is_sparse());
    EXPECT_FALSE(forced->is_sparse());
    EXPECT_EQ(forced->is_sparse_proportion(), 0);

    EXPECT_TRUE(forced->prefer_rows());
    EXPECT_EQ(forced->prefer_rows_proportion(), 1);

    EXPECT_FALSE(forced->uses_oracle(true));
}

/*************************************
 *************************************/

class ForcedDenseFullAccessTest :
    public ::testing::TestWithParam<tatami_test::StandardTestAccessOptions>,
    public ForcedDenseUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(ForcedDenseFullAccessTest, Full) {
    auto tparam = GetParam(); 
    auto opts = tatami_test::convert_test_access_options(tparam);
    tatami_test::test_full_access(*forced, *dense, opts);
}

INSTANTIATE_TEST_SUITE_P(
    CompressedForcedDenseMatrix,
    ForcedDenseFullAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::TestAccessOrder::FORWARD, tatami_test::TestAccessOrder::REVERSE, tatami_test::TestAccessOrder::RANDOM),
        ::testing::Values(1, 4, 10, 20) // jump, to test the workspace's memory; trying out longer jumps to check secondaries.
    )
);

/*************************************
 *************************************/

class ForcedDenseBlockAccessTest :
    public ::testing::TestWithParam<std::tuple<tatami_test::StandardTestAccessOptions, std::pair<double, double> > >,
    public ForcedDenseUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(ForcedDenseBlockAccessTest, Block) {
    auto tparam = GetParam(); 
    auto opts = tatami_test::convert_test_access_options(std::get<0>(tparam));
    auto interval_info = std::get<1>(tparam);
    tatami_test::test_block_access(*forced, *dense, interval_info.first, interval_info.second, opts);
}

INSTANTIATE_TEST_SUITE_P(
    CompressedForcedDenseMatrix,
    ForcedDenseBlockAccessTest,
    ::testing::Combine(
        tatami_test::standard_test_access_options_combinations(),
        ::testing::Values(
            std::make_pair(0.0, 0.51),
            std::make_pair(0.25, 0.65), 
            std::make_pair(0.63, 0.37)
        )
    )
);

/*************************************
 *************************************/

class ForcedDenseIndexedAccessTest : 
    public ::testing::TestWithParam<std::tuple<tatami_test::StandardTestAccessOptions, std::pair<double, double> > >, 
    public ForcedDenseUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(ForcedDenseIndexedAccessTest, Indexed) {
    auto tparam = GetParam(); 
    auto opts = tatami_test::convert_test_access_options(std::get<0>(tparam));
    auto interval_info = std::get<1>(tparam);
    tatami_test::test_indexed_access(*forced, *dense, interval_info.first, interval_info.second, opts);
}

INSTANTIATE_TEST_SUITE_P(
    CompressedForcedDenseMatrix,
    ForcedDenseIndexedAccessTest,
    ::testing::Combine(
        tatami_test::standard_test_access_options_combinations(),
        ::testing::Values(
            std::make_pair(0.0, 0.2),
            std::make_pair(0.2, 0.15), 
            std::make_pair(0.7, 0.3)
        )
    )
);
