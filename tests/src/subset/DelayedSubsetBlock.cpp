#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/subset/DelayedSubsetBlock.hpp"
#include "tatami/subset/make_DelayedSubset.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

class SubsetBlockUtils {
protected:
    inline static int NR = 192, NC = 132;
    inline static std::vector<double> simulated;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, sparse;

    static void assemble() {
        if (dense) {
            return;
        }
        simulated = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.2);
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(NR, NC, simulated));
        sparse = tatami::convert_to_compressed_sparse<false, double, int>(dense.get()); // column-major.
    }

public:
    typedef std::tuple<bool, std::pair<double, double> > SimulationParameters;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false), // row or column subsetting, respectively.
            ::testing::Values( // the block dimensions.
                std::make_pair(0.0, 0.5),
                std::make_pair(0.25, 0.8),
                std::make_pair(0.4, 1)
            )
        );
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_block, sparse_block, ref;
    inline static SimulationParameters last_params;
    inline static int block_length;

    static void assemble(SimulationParameters sim_params) {
        if (ref && last_params == sim_params) {
            return;
        }
        last_params = sim_params;

        assemble();

        auto bind_rows = std::get<0>(sim_params);
        auto interval_info = std::get<1>(sim_params);

        auto full = (bind_rows ? NR : NC);
        int first = full * interval_info.first, last = full * interval_info.second;
        block_length = last - first;

        if (bind_rows) {
            std::vector<double> sub(simulated.data() + first * NC, simulated.data() + last * NC);
            ref.reset(new tatami::DenseRowMatrix<double, int>(block_length, NC, std::move(sub)));
            dense_block = tatami::make_DelayedSubsetBlock<0>(dense, first, block_length);
            sparse_block = tatami::make_DelayedSubsetBlock<0>(sparse, first, block_length);
        } else {
            std::vector<double> sub;
            sub.reserve(NR * block_length);
            for (int r = 0; r < NR; ++r) {
                auto row = simulated.data() + r * NC;
                sub.insert(sub.end(), row + first, row + last);
            }
            ref.reset(new tatami::DenseRowMatrix<double, int>(NR, block_length, std::move(sub)));
            dense_block = tatami::make_DelayedSubsetBlock<1>(dense, first, block_length);
            sparse_block = tatami::make_DelayedSubsetBlock<1>(sparse, first, block_length);
        }
    }
};

/*****************************
 *****************************/

class SubsetBlockTest : public ::testing::TestWithParam<typename SubsetBlockUtils::SimulationParameters>, public SubsetBlockUtils {
protected:
    void SetUp() {
        assemble(GetParam());
    }
};

TEST_P(SubsetBlockTest, Basic) {
    EXPECT_EQ(ref->nrow(), dense_block->nrow());
    EXPECT_EQ(ref->ncol(), dense_block->ncol());

    auto bind_rows = std::get<0>(last_params);
    if (bind_rows) {
        EXPECT_EQ(block_length, dense_block->nrow());
        EXPECT_EQ(dense->ncol(), dense_block->ncol());
    } else {
        EXPECT_EQ(dense->nrow(), dense_block->nrow());
        EXPECT_EQ(block_length, dense_block->ncol());
    }

    EXPECT_FALSE(dense_block->sparse());
    EXPECT_EQ(dense_block->sparse_proportion(), 0);
    EXPECT_TRUE(sparse_block->sparse());
    EXPECT_EQ(sparse_block->sparse_proportion(), 1);

    EXPECT_TRUE(dense_block->prefer_rows());
    EXPECT_EQ(dense_block->prefer_rows_proportion(), 1);
    EXPECT_FALSE(sparse_block->prefer_rows());
    EXPECT_EQ(sparse_block->prefer_rows_proportion(), 0);

    EXPECT_FALSE(dense_block->uses_oracle(false));
}

INSTANTIATE_TEST_SUITE_P(
    DelayedSubsetBlock,
    SubsetBlockTest,
    SubsetBlockUtils::simulation_parameter_combinations()
);

TEST(SubsetBlockMisc, SubsetOracle) {
    auto out = std::make_shared<tatami::ConsecutiveOracle<int> >(10, 100);
    auto casted = tatami::DelayedSubsetBlock_internal::SubsetOracle<int>(out, 50);
    EXPECT_EQ(casted.total(), 100);
}

/*****************************
 *****************************/

class SubsetBlockFullAccessTest : public ::testing::TestWithParam<std::tuple<SubsetBlockUtils::SimulationParameters, tatami_test::StandardTestAccessParameters> >, public SubsetBlockUtils {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(SubsetBlockFullAccessTest, Basic) {
    auto tparam = GetParam();
    auto params = tatami_test::convert_access_parameters(std::get<1>(tparam));
    tatami_test::test_full_access(params, dense_block.get(), ref.get());
    tatami_test::test_full_access(params, sparse_block.get(), ref.get());
}

INSTANTIATE_TEST_SUITE_P(
    DelayedSubsetBlock,
    SubsetBlockFullAccessTest,
    ::testing::Combine(
        SubsetBlockUtils::simulation_parameter_combinations(),
        tatami_test::standard_test_access_parameter_combinations()
    )
);

/*****************************
 *****************************/

class SubsetBlockSlicedAccessTest : public ::testing::TestWithParam<std::tuple<SubsetBlockUtils::SimulationParameters, tatami_test::StandardTestAccessParameters, std::pair<double, double> > >, public SubsetBlockUtils {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(SubsetBlockSlicedAccessTest, Sliced) {
    auto tparam = GetParam();
    auto params = tatami_test::convert_access_parameters(std::get<1>(tparam));

    auto interval_info = std::get<2>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;

    tatami_test::test_block_access(params, dense_block.get(), ref.get(), FIRST, LAST);
    tatami_test::test_block_access(params, sparse_block.get(), ref.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedSubsetBlock,
    SubsetBlockSlicedAccessTest,
    ::testing::Combine(
        SubsetBlockUtils::simulation_parameter_combinations(),
        tatami_test::standard_test_access_parameter_combinations(),
        ::testing::Values(
            std::make_pair(0.0, 0.45), 
            std::make_pair(0.33, 0.66),
            std::make_pair(0.56, 1.0)
        )
    )
);

/*****************************
 *****************************/

class SubsetBlockIndexedAccessTest : public ::testing::TestWithParam<std::tuple<SubsetBlockUtils::SimulationParameters, tatami_test::StandardTestAccessParameters, std::pair<double, int> > >, public SubsetBlockUtils {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(SubsetBlockIndexedAccessTest, Indexed) {
    auto tparam = GetParam();
    auto params = tatami_test::convert_access_parameters(std::get<1>(tparam));

    auto interval_info = std::get<2>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, STEP = interval_info.second;

    tatami_test::test_indexed_access(params, dense_block.get(), ref.get(), FIRST, STEP);
    tatami_test::test_indexed_access(params, sparse_block.get(), ref.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedSubsetBlock,
    SubsetBlockIndexedAccessTest,
    ::testing::Combine(
        SubsetBlockUtils::simulation_parameter_combinations(),
        tatami_test::standard_test_access_parameter_combinations(),
        ::testing::Values(
            std::make_pair(0.0, 12), 
            std::make_pair(0.33, 6),
            std::make_pair(0.56, 9)
        )
    )
);

/****************************************************
 ****************************************************/

TEST(DelayedSubsetBlock, ConstOverload) {
    int NR = 9, NC = 7;
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(NR, NC, tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1)));
    auto sub = tatami::make_DelayedSubsetBlock<1>(dense, static_cast<int>(5), static_cast<int>(3));
    EXPECT_EQ(sub->ncol(), 3);
    EXPECT_EQ(sub->nrow(), NR);
}

TEST(DelayedSubsetBlock, CorrectMaker) {
    int NR = 90, NC = 50;
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(NR, NC, tatami_test::simulate_dense_vector<double>(NR * NC)));

    // Checking that the make function dispatches correctly to the block subset class.
    {
        std::vector<int> indices { 5, 6, 7, 8, 9, 10 };
        auto sub = tatami::make_DelayedSubset<1>(dense, indices);
        EXPECT_EQ(sub->ncol(), 6);
        EXPECT_EQ(sub->nrow(), NR);

        auto ref = tatami::make_DelayedSubsetBlock<1>(dense, static_cast<int>(5), static_cast<int>(6));
        tatami_test::test_simple_row_access(sub.get(), ref.get());
        tatami_test::test_simple_column_access(sub.get(), ref.get());
    }

    // Checking that it behaves correctly with an empty index vector.
    {
        std::vector<int> indices;
        auto sub = tatami::make_DelayedSubset<0>(dense, indices);
        EXPECT_EQ(sub->ncol(), NC);
        EXPECT_EQ(sub->nrow(), 0);
    }
}
