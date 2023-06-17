#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/other/DelayedBind.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

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

class DelayedBindFullAccessTest : public ::testing::TestWithParam<std::tuple<std::vector<int>, bool, bool, int> >, public DelayedBindTestMethods {};

TEST_P(DelayedBindFullAccessTest, Basic) {
    auto param = GetParam();
    assemble(std::get<0>(param), 50, std::get<1>(param));
    int FORWARD = std::get<2>(param);
    int JUMP = std::get<3>(param);

    tatami_test::test_simple_column_access(bound_sparse.get(), manual.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(bound_dense.get(), manual.get(), FORWARD, JUMP);

    tatami_test::test_simple_row_access(bound_sparse.get(), manual.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(bound_dense.get(), manual.get(), FORWARD, JUMP);
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

    tatami_test::test_sliced_column_access(bound_sparse.get(), manual.get(), FORWARD, JUMP, RFIRST, RLAST);
    tatami_test::test_sliced_column_access(bound_dense.get(), manual.get(), FORWARD, JUMP, RFIRST, RLAST);

    tatami_test::test_sliced_row_access(bound_sparse.get(), manual.get(), FORWARD, JUMP, CFIRST, CLAST);
    tatami_test::test_sliced_row_access(bound_dense.get(), manual.get(), FORWARD, JUMP, CFIRST, CLAST);
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

    tatami_test::test_indexed_column_access(bound_sparse.get(), manual.get(), FORWARD, JUMP, RFIRST, STEP);
    tatami_test::test_indexed_column_access(bound_dense.get(), manual.get(), FORWARD, JUMP, RFIRST, STEP);

    tatami_test::test_indexed_row_access(bound_sparse.get(), manual.get(), FORWARD, JUMP, CFIRST, STEP);
    tatami_test::test_indexed_row_access(bound_dense.get(), manual.get(), FORWARD, JUMP, CFIRST, STEP);
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

/****************************
 ****************************/

class DelayedBindOracleTestCore {
protected:
    std::vector<std::shared_ptr<tatami::NumericMatrix> > collected;
    std::shared_ptr<tatami::NumericMatrix> bound;

    void assemble(const std::vector<int>& lengths, int dim, bool row) {
        for (size_t i = 0; i < lengths.size(); ++i) {
            int len = lengths[i];
            auto to_add = tatami_test::simulate_sparse_compressed<double>(len, dim, 0.2, /* lower = */ -10, /* upper = */ 10, /* seed = */ i * 99);
            if (row) {
                collected.emplace_back(new tatami::CompressedSparseRowMatrix<double, int>(len, dim, std::move(to_add.value), std::move(to_add.index), std::move(to_add.ptr)));
            } else {
                collected.emplace_back(new tatami::CompressedSparseColumnMatrix<double, int>(dim, len, std::move(to_add.value), std::move(to_add.index), std::move(to_add.ptr)));
            }
        }

        bound = combine(collected, row);
    }

    static std::shared_ptr<tatami::NumericMatrix> combine(std::vector<std::shared_ptr<tatami::NumericMatrix> > inputs, bool row) {
        if (row) {
            return tatami::make_DelayedBind<0>(std::move(inputs));
        } else {
            return tatami::make_DelayedBind<1>(std::move(inputs));
        }
    }
};

class DelayedBindOracleTest : public ::testing::TestWithParam<std::tuple<std::vector<int>, bool, bool> >, public DelayedBindOracleTestCore {};

TEST_P(DelayedBindOracleTest, AllOracular) {
    auto param = GetParam();
    assemble(std::get<0>(param), 500, std::get<1>(param));
    auto random = std::get<2>(param);

    for (size_t m = 0; m < collected.size(); ++m) {
        int step_size = (m + 1) * 10; // variable prediction number across bound matrices, for some variety.
        collected[m] = tatami_test::make_CrankyMatrix(std::move(collected[m]), step_size);
    }
    auto wrapped_bound = combine(std::move(collected), std::get<1>(param));

    EXPECT_FALSE(bound->uses_oracle(true));
    EXPECT_TRUE(wrapped_bound->uses_oracle(true));

    tatami_test::test_oracle_column_access(wrapped_bound.get(), bound.get(), random);
    tatami_test::test_oracle_row_access(wrapped_bound.get(), bound.get(), random);
}

TEST_P(DelayedBindOracleTest, FirstOracular) {
    auto param = GetParam();
    assemble(std::get<0>(param), 350, std::get<1>(param));
    auto random = std::get<2>(param);

    collected.front() = tatami_test::make_CrankyMatrix(std::move(collected.front()));
    auto wrapped_bound = combine(std::move(collected), std::get<1>(param));

    EXPECT_FALSE(bound->uses_oracle(true));
    EXPECT_TRUE(wrapped_bound->uses_oracle(true));

    tatami_test::test_oracle_column_access(wrapped_bound.get(), bound.get(), random);
    tatami_test::test_oracle_row_access(wrapped_bound.get(), bound.get(), random);
}

TEST_P(DelayedBindOracleTest, LastOracular) {
    auto param = GetParam();
    assemble(std::get<0>(param), 540, std::get<1>(param));
    auto random = std::get<2>(param);

    collected.back() = tatami_test::make_CrankyMatrix(std::move(collected.back()));
    auto wrapped_bound = combine(std::move(collected), std::get<1>(param));

    EXPECT_FALSE(bound->uses_oracle(true));
    EXPECT_TRUE(wrapped_bound->uses_oracle(true));

    tatami_test::test_oracle_column_access(wrapped_bound.get(), bound.get(), random);
    tatami_test::test_oracle_row_access(wrapped_bound.get(), bound.get(), random);
}

INSTANTIATE_TEST_CASE_P(
    DelayedBind,
    DelayedBindOracleTest,
    ::testing::Combine(
        ::testing::Values(
            std::vector<int>{ 100 },
            std::vector<int>{ 150, 100 },
            std::vector<int>{ 50, 200, 150 }
        ),
        ::testing::Values(true, false), // bind by row or by column
        ::testing::Values(true, false)  // use random or consecutive oracle.
    )
);

class DelayedBindOracleTest2 : public ::testing::TestWithParam<std::tuple<std::vector<int>, bool, bool> >, public DelayedBindOracleTestCore {};

TEST_F(DelayedBindOracleTest2, Elongated) {
    size_t NC = 200;
    assemble({ 10, 20, 30 }, NC, true);
    for (size_t m = 0; m < collected.size(); ++m) {
        collected[m] = tatami_test::make_CrankyMatrix(std::move(collected[m]), 20 - m); // again, some variety in the prediction numbers.
    }
    auto wrapped_bound = combine(std::move(collected), true); // combining by row.

    // Use a very long simulated sequence.
    // This checks that the collection of expired predictions works correctly
    // in the ParallelExtractor::ParentOracle class.

    std::mt19937_64 rng(4123123); 
    std::vector<int> fixed(50000);
    for (auto& x : fixed) {
        x = rng() % NC;
    }

    auto swork = bound->sparse_column();
    auto swork_o = wrapped_bound->sparse_column();
    swork_o->set_oracle(std::make_unique<tatami::FixedOracle<int> >(fixed.data(), fixed.size()));

    for (auto i : fixed) {
        auto sexpected = swork->fetch(i);
        auto sobserved = swork_o->fetch(i);
        EXPECT_EQ(sexpected.index, sobserved.index);
        EXPECT_EQ(sexpected.value, sobserved.value);
    }
}
