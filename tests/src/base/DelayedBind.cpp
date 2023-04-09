#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedBind.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../_tests/test_row_access.h"
#include "../_tests/test_column_access.h"
#include "../_tests/simulate_vector.h"

class DelayedBindTestMethods {
protected:
    std::shared_ptr<tatami::NumericMatrix> bound_dense, bound_sparse, manual;

    void assemble(const std::vector<int>& lengths, int dim, bool row) {
        std::vector<double> concat;
        size_t n_total = 0;
        std::vector<std::shared_ptr<tatami::NumericMatrix> > collected_dense, collected_sparse;

        for (size_t i = 0; i < lengths.size(); ++i) {
            auto to_add = simulate_sparse_vector<double>(lengths[i] * dim, 0.2, /* lower = */ -10, /* upper = */ 10, /* seed = */ i * 1000);
            concat.insert(concat.end(), to_add.begin(), to_add.end());
            n_total += lengths[i];

            if (row) {
                collected_dense.emplace_back(new tatami::DenseRowMatrix<double, int>(lengths[i], dim, to_add));
            } else {
                collected_dense.emplace_back(new tatami::DenseColumnMatrix<double, int>(dim, lengths[i], to_add));
            }
            collected_sparse.push_back(tatami::convert_to_sparse<false>(collected_dense.back().get())); // always CSC
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
    EXPECT_TRUE(bound_dense->prefer_rows());

    EXPECT_EQ(bound_sparse->nrow(), 35);
    EXPECT_EQ(bound_sparse->ncol(), 20);
    EXPECT_TRUE(bound_sparse->sparse());
    EXPECT_FALSE(bound_sparse->prefer_rows());
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
    EXPECT_FALSE(bound_dense->prefer_rows());
}

TEST_F(DelayedBindUtilsTest, InconsistentBinds) {
    assemble({ 10, 20, 5 }, 20, true); 

    // Bound_sparse is CSC, bound_dense is row-major.
    auto combined = tatami::make_DelayedBind<1>(std::vector<std::shared_ptr<tatami::NumericMatrix> >{ bound_sparse, bound_dense });
    auto preference = combined->dimension_preference();

    EXPECT_TRUE(preference.first > 0);
    EXPECT_TRUE(preference.second > 0);
    EXPECT_EQ(preference.first, preference.second);

    EXPECT_FALSE(combined->sparse());
}

TEST_F(DelayedBindUtilsTest, EmptyBind) {
    assemble({}, 20, true); 
    EXPECT_EQ(bound_dense->nrow(), 0);
    EXPECT_EQ(bound_dense->ncol(), 0);

    // Checking that empty workspaces can be constructed.
    {
        auto rthing = bound_dense->dense_row_workspace();
        EXPECT_NE(rthing.get(), nullptr);

        auto cthing = bound_dense->dense_column_workspace();
        EXPECT_NE(cthing.get(), nullptr);
    }

    {
        auto rthing = bound_dense->dense_row_workspace(0, 0);
        EXPECT_NE(rthing.get(), nullptr);

        auto cthing = bound_dense->dense_column_workspace(0, 0);
        EXPECT_NE(cthing.get(), nullptr);
    }

    {
        auto rthing = bound_dense->dense_row_workspace(std::vector<int>());
        EXPECT_NE(rthing.get(), nullptr);

        auto cthing = bound_dense->dense_column_workspace(std::vector<int>());
        EXPECT_NE(cthing.get(), nullptr);
    }
}

/****************************
 ****************************/

class DelayedBindFullAccessTest : public ::testing::TestWithParam<std::tuple<std::vector<int>, bool, bool, int> >, public DelayedBindTestMethods {};

TEST_P(DelayedBindFullAccessTest, Basic) {
    auto param = GetParam();
    assemble(std::get<0>(param), 50, std::get<1>(param));
    int FORWARD = std::get<2>(param);
    int JUMP = std::get<3>(param);

    test_simple_column_access(bound_sparse.get(), manual.get(), FORWARD, JUMP);
    test_simple_column_access(bound_dense.get(), manual.get(), FORWARD, JUMP);

    test_simple_row_access(bound_sparse.get(), manual.get(), FORWARD, JUMP);
    test_simple_row_access(bound_dense.get(), manual.get(), FORWARD, JUMP);
}

static auto spawn_bind_scenarios () {
    return ::testing::Values(
        std::vector<int>{ 10 },
        std::vector<int>{ 10, 20 },
        std::vector<int>{ 5, 2, 5 },
        std::vector<int>{ 5, 10, 20 },
        std::vector<int>{ 5, 0, 5 }
    );
}

INSTANTIATE_TEST_CASE_P(
    DelayedBind,
    DelayedBindFullAccessTest,
    ::testing::Combine(
        spawn_bind_scenarios(),
        ::testing::Values(true, false), // bind by row or by column
        ::testing::Values(true, false), // forward or backward traversal.
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

/****************************
 ****************************/

class DelayedBindSlicedAccessTest : public ::testing::TestWithParam<std::tuple<std::vector<int>, bool, bool, int, std::vector<double> > >, public DelayedBindTestMethods {};

TEST_P(DelayedBindSlicedAccessTest, Basic) {
    auto param = GetParam();
    assemble(std::get<0>(param), 50, std::get<1>(param));
    int FORWARD = std::get<2>(param);
    int JUMP = std::get<3>(param);

    auto interval_info = std::get<4>(param);
    size_t RFIRST = interval_info[0] * manual->nrow(), RLAST = interval_info[1] * manual->nrow();
    size_t CFIRST = interval_info[0] * manual->ncol(), CLAST = interval_info[1] * manual->ncol();

    test_sliced_column_access(bound_sparse.get(), manual.get(), FORWARD, JUMP, RFIRST, RLAST);
    test_sliced_column_access(bound_dense.get(), manual.get(), FORWARD, JUMP, RFIRST, RLAST);

    test_sliced_row_access(bound_sparse.get(), manual.get(), FORWARD, JUMP, CFIRST, CLAST);
    test_sliced_row_access(bound_dense.get(), manual.get(), FORWARD, JUMP, CFIRST, CLAST);
}

INSTANTIATE_TEST_CASE_P(
    DelayedBind,
    DelayedBindSlicedAccessTest,
    ::testing::Combine(
        spawn_bind_scenarios(),
        ::testing::Values(true, false), // bind by row or by column
        ::testing::Values(true, false), // forward or backward traversal.
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.6 }), 
            std::vector<double>({ 0.25, 0.75 }), 
            std::vector<double>({ 0.55, 1 })
        )
    )
);

/****************************
 ****************************/

class DelayedBindIndexedAccessTest : public ::testing::TestWithParam<std::tuple<std::vector<int>, bool, bool, int, std::vector<double> > >, public DelayedBindTestMethods {};

TEST_P(DelayedBindIndexedAccessTest, Basic) {
    auto param = GetParam();
    assemble(std::get<0>(param), 50, std::get<1>(param));
    int FORWARD = std::get<2>(param);
    int JUMP = std::get<3>(param);

    auto interval_info = std::get<4>(param);
    size_t RFIRST = interval_info[0] * manual->nrow(),
        CFIRST = interval_info[0] * manual->ncol(), 
        STEP = interval_info[1];

    test_indexed_column_access(bound_sparse.get(), manual.get(), FORWARD, JUMP, RFIRST, STEP);
    test_indexed_column_access(bound_dense.get(), manual.get(), FORWARD, JUMP, RFIRST, STEP);

    test_indexed_row_access(bound_sparse.get(), manual.get(), FORWARD, JUMP, CFIRST, STEP);
    test_indexed_row_access(bound_dense.get(), manual.get(), FORWARD, JUMP, CFIRST, STEP);
}

INSTANTIATE_TEST_CASE_P(
    DelayedBind,
    DelayedBindIndexedAccessTest,
    ::testing::Combine(
        spawn_bind_scenarios(),
        ::testing::Values(true, false), // bind by row or by column
        ::testing::Values(true, false), // forward or backward traversal.
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 5 }), 
            std::vector<double>({ 0.33, 3 }),
            std::vector<double>({ 0.5, 2 })
        )
    )
);
