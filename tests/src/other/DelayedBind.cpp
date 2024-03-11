#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/other/DelayedBind.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

class DelayedBindTestMethods {
protected:
    std::shared_ptr<tatami::NumericMatrix> bound_dense, bound_sparse, manual;

    void assemble(const std::vector<int>& lengths, int dim, bool row) {
        std::vector<double> concat;
        size_t n_total = 0;
        std::vector<std::shared_ptr<tatami::NumericMatrix> > collected_dense, collected_sparse;

        for (size_t i = 0; i < lengths.size(); ++i) {
            auto to_add = tatami_test::simulate_sparse_vector<double>(lengths[i] * dim, 0.2, /* lower = */ -10, /* upper = */ 10, /* seed = */ i * 1000);
            concat.insert(concat.end(), to_add.begin(), to_add.end());
            n_total += lengths[i];

            if (row) {
                collected_dense.emplace_back(new tatami::DenseRowMatrix<double, int>(lengths[i], dim, to_add));
            } else {
                collected_dense.emplace_back(new tatami::DenseColumnMatrix<double, int>(dim, lengths[i], to_add));
            }
            collected_sparse.push_back(tatami::convert_to_compressed_sparse<false>(collected_dense.back().get())); // always CSC
        }

        if (row) {
            bound_dense = tatami::make_DelayedBind<0>(std::move(collected_dense));
            bound_sparse = tatami::make_DelayedBind<0>(std::move(collected_sparse));
            manual.reset(new tatami::DenseRowMatrix<double, int>(n_total, dim, std::move(concat)));
        } else {
            bound_dense = tatami::make_DelayedBind<1>(std::move(collected_dense));
            bound_sparse = tatami::make_DelayedBind<1>(std::move(collected_sparse));
            manual.reset(new tatami::DenseColumnMatrix<double, int>(dim, n_total, std::move(concat)));
        }

        return;
    }
};

class DelayedBindUtilsTest : public ::testing::Test, public DelayedBindTestMethods {};

TEST_F(DelayedBindUtilsTest, ByRow) {
    assemble({ 10, 20, 5 }, 20, true);

    EXPECT_EQ(bound_dense->nrow(), 35);
    EXPECT_EQ(bound_dense->ncol(), 20);
    EXPECT_FALSE(bound_dense->sparse());
    EXPECT_EQ(bound_dense->sparse_proportion(), 0);
    EXPECT_TRUE(bound_dense->prefer_rows());
    EXPECT_EQ(bound_dense->prefer_rows_proportion(), 1);

    EXPECT_EQ(bound_sparse->nrow(), 35);
    EXPECT_EQ(bound_sparse->ncol(), 20);
    EXPECT_TRUE(bound_sparse->sparse());
    EXPECT_EQ(bound_sparse->sparse_proportion(), 1);
    EXPECT_FALSE(bound_sparse->prefer_rows());
    EXPECT_EQ(bound_sparse->prefer_rows_proportion(), 0);

    EXPECT_FALSE(bound_sparse->uses_oracle(true));
    EXPECT_FALSE(bound_sparse->uses_oracle(true));
}

TEST_F(DelayedBindUtilsTest, ByColumn) {
    assemble({ 10, 20, 5 }, 20, false);

    EXPECT_EQ(bound_dense->nrow(), 20);
    EXPECT_EQ(bound_dense->ncol(), 35);
    EXPECT_FALSE(bound_dense->sparse());
    EXPECT_FALSE(bound_dense->prefer_rows());

    EXPECT_EQ(bound_sparse->nrow(), 20);
    EXPECT_EQ(bound_sparse->ncol(), 35);
    EXPECT_TRUE(bound_sparse->sparse());
    EXPECT_FALSE(bound_sparse->prefer_rows());

    EXPECT_FALSE(bound_dense->uses_oracle(false));
    EXPECT_FALSE(bound_sparse->uses_oracle(false));
}

TEST_F(DelayedBindUtilsTest, InconsistentBinds) {
    assemble({ 10, 20, 5 }, 20, true); 

    // Bound_sparse is CSC, bound_dense is row-major.
    auto combined = tatami::make_DelayedBind<1>(std::vector<std::shared_ptr<tatami::NumericMatrix> >{ bound_sparse, bound_dense });

    EXPECT_FLOAT_EQ(combined->sparse_proportion(), 0.5);
    EXPECT_FALSE(combined->sparse());

    EXPECT_FLOAT_EQ(combined->prefer_rows_proportion(), 0.5);
    EXPECT_FALSE(combined->prefer_rows());
}

TEST_F(DelayedBindUtilsTest, EmptyBind) {
    assemble({}, 20, true); 
    EXPECT_EQ(bound_dense->nrow(), 0);
    EXPECT_EQ(bound_dense->ncol(), 0);

    // Checking that empty workspaces can be constructed.
    {
        auto rthing = bound_dense->dense_row();
        EXPECT_NE(rthing.get(), nullptr);

        auto cthing = bound_dense->dense_column();
        EXPECT_NE(cthing.get(), nullptr);
    }

    {
        auto rthing = bound_dense->dense_row(0, 0);
        EXPECT_NE(rthing.get(), nullptr);

        auto cthing = bound_dense->dense_column(0, 0);
        EXPECT_NE(cthing.get(), nullptr);
    }

    {
        auto rthing = bound_dense->dense_row(std::vector<int>());
        EXPECT_NE(rthing.get(), nullptr);

        auto cthing = bound_dense->dense_column(std::vector<int>());
        EXPECT_NE(cthing.get(), nullptr);
    }
}

TEST_F(DelayedBindUtilsTest, ConstOverloads) {
    assemble({ 10, 50 }, 20, true); 
    std::vector<std::shared_ptr<const tatami::NumericMatrix> > const_collected({ bound_dense, bound_sparse });
    auto const_combined = tatami::make_DelayedBind<0>(std::move(const_collected));

    // Some cursory checks.
    EXPECT_EQ(const_combined->nrow(), 120); // i.e., (10 + 50) * 2 
    EXPECT_EQ(const_combined->ncol(), 20);
}

/****************************
 ****************************/

class DelayedBindFullAccessTest : 
    public ::testing::TestWithParam<std::tuple<std::vector<int>, bool, bool, bool, tatami_test::TestAccessOrder, int> >, 
    public DelayedBindTestMethods {};

TEST_P(DelayedBindFullAccessTest, Basic) {
    auto tparam = GetParam();
    assemble(std::get<0>(tparam), 50, std::get<1>(tparam));

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);

    tatami_test::test_full_access(params, bound_sparse.get(), manual.get());
    tatami_test::test_full_access(params, bound_dense.get(), manual.get());
}

static auto spawn_bind_scenarios () {
    return ::testing::Values(
        std::vector<int>{ 10 },
        std::vector<int>{ 10, 20 },
        std::vector<int>{ 20, 10 },
        std::vector<int>{ 5, 2, 5 },
        std::vector<int>{ 5, 10, 20 },
        std::vector<int>{ 20, 10, 5 },
        std::vector<int>{ 5, 0, 5 }
    );
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindFullAccessTest,
    ::testing::Combine(
        spawn_bind_scenarios(),
        ::testing::Values(true, false), // bind by row or by column
        ::testing::Values(true, false), // access by row or column
        ::testing::Values(true, false), // use oracle or not.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

/****************************
 ****************************/

class DelayedBindSlicedAccessTest : 
    public ::testing::TestWithParam<std::tuple<std::vector<int>, bool, bool, bool, tatami_test::TestAccessOrder, int, std::pair<double, double> > >, 
    public DelayedBindTestMethods {};

TEST_P(DelayedBindSlicedAccessTest, Basic) {
    auto tparam = GetParam();
    assemble(std::get<0>(tparam), 50, std::get<1>(tparam));

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);

    auto interval_info = std::get<6>(tparam);
    auto len = (params.use_row ?  manual->ncol() : manual->nrow());
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;

    tatami_test::test_block_access(params, bound_sparse.get(), manual.get(), FIRST, LAST);
    tatami_test::test_block_access(params, bound_dense.get(), manual.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindSlicedAccessTest,
    ::testing::Combine(
        spawn_bind_scenarios(),
        ::testing::Values(true, false), // bind by row or by column
        ::testing::Values(true, false), // access by row or column
        ::testing::Values(true, false), // use oracle or not.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair(0, 0.6), 
            std::make_pair(0.25, 0.75), 
            std::make_pair(0.55, 1)
        )
    )
);

/****************************
 ****************************/

class DelayedBindIndexedAccessTest : 
    public ::testing::TestWithParam<std::tuple<std::vector<int>, bool, bool, bool, tatami_test::TestAccessOrder, int, std::pair<double, int> > >, 
    public DelayedBindTestMethods {};

TEST_P(DelayedBindIndexedAccessTest, Basic) {
    auto tparam = GetParam();
    assemble(std::get<0>(tparam), 50, std::get<1>(tparam));

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);

    auto interval_info = std::get<6>(tparam);
    size_t len = (params.use_row ? manual->ncol() : manual->nrow());
    size_t FIRST = interval_info.first * len, STEP = interval_info.second;

    tatami_test::test_indexed_access(params, bound_sparse.get(), manual.get(), FIRST, STEP);
    tatami_test::test_indexed_access(params, bound_dense.get(), manual.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindIndexedAccessTest,
    ::testing::Combine(
        spawn_bind_scenarios(),
        ::testing::Values(true, false), // bind by row or by column
        ::testing::Values(true, false), // access by row or column
        ::testing::Values(true, false), // use oracle or not.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::pair<double, int>(0, 5), 
            std::pair<double, int>(0.33, 3),
            std::pair<double, int>(0.5, 2)
        )
    )
);
