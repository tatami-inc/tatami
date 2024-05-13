#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/other/DelayedBind.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

class DelayedBindUtils {
public:
    typedef std::tuple<std::vector<int>, bool> SimulationParameters;

    static auto spawn_bind_scenarios () {
        return ::testing::Combine(
            ::testing::Values(
                std::vector<int>{ 10 },
                std::vector<int>{ 10, 20 },
                std::vector<int>{ 20, 10 },
                std::vector<int>{ 5, 2, 5 },
                std::vector<int>{ 5, 10, 20 },
                std::vector<int>{ 20, 10, 5 },
                std::vector<int>{ 5, 0, 5 }
            ),
            ::testing::Values(true, false) // bind by row or by column
        );
    }

protected:
    inline static size_t otherdim = 97;
    inline static std::shared_ptr<tatami::NumericMatrix> bound_dense, bound_sparse, manual, forced_bound_dense, forced_bound_sparse;
    inline static SimulationParameters last_params;

    static void assemble(SimulationParameters sim_params) {
        if (bound_dense && sim_params == last_params) {
            return;
        }
        last_params = sim_params;

        const auto& lengths = std::get<0>(sim_params);
        bool row = std::get<1>(sim_params);

        std::vector<double> concat;
        size_t n_total = 0;
        std::vector<std::shared_ptr<tatami::NumericMatrix> > collected_dense, collected_sparse;
        std::vector<std::shared_ptr<tatami::NumericMatrix> > forced_collected_dense, forced_collected_sparse;

        for (size_t i = 0; i < lengths.size(); ++i) {
            auto to_add = tatami_test::simulate_sparse_vector<double>(lengths[i] * otherdim, 0.2, /* lower = */ -10, /* upper = */ 10, /* seed = */ i * 1000 + lengths[i]);
            concat.insert(concat.end(), to_add.begin(), to_add.end());
            n_total += lengths[i];

            if (row) {
                collected_dense.emplace_back(new tatami::DenseRowMatrix<double, int>(lengths[i], otherdim, to_add));
            } else {
                collected_dense.emplace_back(new tatami::DenseColumnMatrix<double, int>(otherdim, lengths[i], to_add));
            }
            collected_sparse.push_back(tatami::convert_to_compressed_sparse<false, double, int>(collected_dense.back().get())); // always CSC

            forced_collected_dense.emplace_back(std::make_shared<tatami_test::ForcedOracleWrapper<double, int> >(collected_dense.back()));
            forced_collected_sparse.emplace_back(std::make_shared<tatami_test::ForcedOracleWrapper<double, int> >(collected_sparse.back()));
        }

        if (row) {
            bound_dense = tatami::make_DelayedBind<0>(std::move(collected_dense));
            bound_sparse = tatami::make_DelayedBind<0>(std::move(collected_sparse));
            manual.reset(new tatami::DenseRowMatrix<double, int>(n_total, otherdim, std::move(concat)));
            forced_bound_dense = tatami::make_DelayedBind<0>(std::move(forced_collected_dense));
            forced_bound_sparse = tatami::make_DelayedBind<0>(std::move(forced_collected_sparse));
        } else {
            bound_dense = tatami::make_DelayedBind<1>(std::move(collected_dense));
            bound_sparse = tatami::make_DelayedBind<1>(std::move(collected_sparse));
            manual.reset(new tatami::DenseColumnMatrix<double, int>(otherdim, n_total, std::move(concat)));
            forced_bound_dense = tatami::make_DelayedBind<1>(std::move(forced_collected_dense));
            forced_bound_sparse = tatami::make_DelayedBind<1>(std::move(forced_collected_sparse));
        }
    }
};

class DelayedBindUtilsTest : public ::testing::Test, public DelayedBindUtils {};

TEST_F(DelayedBindUtilsTest, ByRow) {
    assemble(SimulationParameters({ 10, 20, 5 }, true));

    EXPECT_EQ(bound_dense->nrow(), 35);
    EXPECT_EQ(bound_dense->ncol(), otherdim);
    EXPECT_FALSE(bound_dense->sparse());
    EXPECT_EQ(bound_dense->sparse_proportion(), 0);
    EXPECT_TRUE(bound_dense->prefer_rows());
    EXPECT_EQ(bound_dense->prefer_rows_proportion(), 1);

    EXPECT_EQ(bound_sparse->nrow(), 35);
    EXPECT_EQ(bound_sparse->ncol(), otherdim);
    EXPECT_TRUE(bound_sparse->sparse());
    EXPECT_EQ(bound_sparse->sparse_proportion(), 1);
    EXPECT_FALSE(bound_sparse->prefer_rows());
    EXPECT_EQ(bound_sparse->prefer_rows_proportion(), 0);

    EXPECT_FALSE(bound_sparse->uses_oracle(true));
    EXPECT_FALSE(bound_sparse->uses_oracle(true));
}

TEST_F(DelayedBindUtilsTest, ByColumn) {
    assemble(SimulationParameters({ 10, 20, 5 }, false));

    EXPECT_EQ(bound_dense->nrow(), otherdim);
    EXPECT_EQ(bound_dense->ncol(), 35);
    EXPECT_FALSE(bound_dense->sparse());
    EXPECT_FALSE(bound_dense->prefer_rows());

    EXPECT_EQ(bound_sparse->nrow(), otherdim);
    EXPECT_EQ(bound_sparse->ncol(), 35);
    EXPECT_TRUE(bound_sparse->sparse());
    EXPECT_FALSE(bound_sparse->prefer_rows());

    EXPECT_FALSE(bound_dense->uses_oracle(false));
    EXPECT_FALSE(bound_sparse->uses_oracle(false));
}

TEST_F(DelayedBindUtilsTest, InconsistentBinds) {
    assemble(SimulationParameters({ 10, 20, 5 }, true));

    // Bound_sparse is CSC, bound_dense is row-major.
    auto combined = tatami::make_DelayedBind<1>(std::vector<std::shared_ptr<tatami::NumericMatrix> >{ bound_sparse, bound_dense });

    EXPECT_FLOAT_EQ(combined->sparse_proportion(), 0.5);
    EXPECT_FALSE(combined->sparse());

    EXPECT_FLOAT_EQ(combined->prefer_rows_proportion(), 0.5);
    EXPECT_FALSE(combined->prefer_rows());
}

TEST_F(DelayedBindUtilsTest, ConstOverloads) {
    assemble(SimulationParameters({ 10, 50 }, true));

    std::vector<std::shared_ptr<const tatami::NumericMatrix> > const_collected({ bound_dense, bound_sparse });
    auto const_combined = tatami::make_DelayedBind<0>(std::move(const_collected));

    // Some cursory checks.
    EXPECT_EQ(const_combined->nrow(), 120); // i.e., (10 + 50) * 2 
    EXPECT_EQ(const_combined->ncol(), otherdim);
}

TEST(DelayedBindMisc, ErrorCheck) {
    std::vector<std::shared_ptr<tatami::NumericMatrix> > collected;
    collected.emplace_back(new tatami::DenseRowMatrix<double, int>(10, 20, std::vector<double>(200)));
    collected.emplace_back(new tatami::DenseRowMatrix<double, int>(20, 10, std::vector<double>(200)));
    tatami_test::throws_error([&]() { tatami::make_DelayedBind<0>(collected); }, "same number of columns");
    tatami_test::throws_error([&]() { tatami::make_DelayedBind<1>(collected); }, "same number of rows");
}

TEST(DelayedBindMisc, PartialOracleUsage) {
    std::vector<std::shared_ptr<tatami::NumericMatrix> > collected;
    collected.emplace_back(new tatami::DenseRowMatrix<double, int>(10, 20, std::vector<double>(200)));
    collected.emplace_back(new tatami::DenseRowMatrix<double, int>(10, 20, std::vector<double>(200)));

    {
        auto combined = tatami::make_DelayedBind<0>(collected); 
        EXPECT_FALSE(combined->uses_oracle(true));
        EXPECT_FALSE(combined->uses_oracle(false));
    }

    {
        auto p = std::make_shared<tatami_test::ForcedOracleWrapper<double, int> >(std::move(collected.back()));
        collected.back() = std::move(p);
        auto combined = tatami::make_DelayedBind<0>(collected); 
        EXPECT_TRUE(combined->uses_oracle(true));
        EXPECT_TRUE(combined->uses_oracle(false));
    }

    {
        auto p = std::make_shared<tatami_test::ForcedOracleWrapper<double, int> >(std::move(collected.front()));
        collected.front() = std::move(p);
        auto combined = tatami::make_DelayedBind<0>(collected); 
        EXPECT_TRUE(combined->uses_oracle(true));
        EXPECT_TRUE(combined->uses_oracle(false));
    }
}

TEST(DelayedBindMisc, AllEmpty) {
    auto empty = tatami::make_DelayedBind<0>(std::vector<std::shared_ptr<tatami::Matrix<double, int> > >{});
    EXPECT_EQ(empty->nrow(), 0);
    EXPECT_EQ(empty->ncol(), 0);
}

/****************************
 ****************************/

class DelayedBindEmptyAccessTest : public ::testing::TestWithParam<std::tuple<bool, bool, bool> >, public DelayedBindUtils {};

TEST_P(DelayedBindEmptyAccessTest, Empty) {
    auto tparam = GetParam();
    auto bind_rows = std::get<0>(tparam);
    assemble(SimulationParameters({0, 0}, bind_rows));

    if (bind_rows) {
        EXPECT_EQ(bound_dense->nrow(), 0);
        EXPECT_EQ(bound_dense->ncol(), otherdim);
    } else {
        EXPECT_EQ(bound_dense->nrow(), otherdim);
        EXPECT_EQ(bound_dense->ncol(), 0);
    }

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<1>(tparam);
    params.use_oracle = std::get<2>(tparam);

    // Check that we can perform all of the accesses.
    tatami_test::test_full_access(params, bound_dense.get(), manual.get());
    tatami_test::test_block_access(params, bound_dense.get(), manual.get(), 0, 0);
    tatami_test::test_indexed_access(params, bound_dense.get(), manual.get(), 0, 1);

    tatami_test::test_full_access(params, bound_sparse.get(), manual.get());
    tatami_test::test_block_access(params, bound_sparse.get(), manual.get(), 0, 0);
    tatami_test::test_indexed_access(params, bound_sparse.get(), manual.get(), 0, 1);

    if (params.use_oracle) {
        tatami_test::test_full_access(params, forced_bound_dense.get(), manual.get());
        tatami_test::test_block_access(params, forced_bound_dense.get(), manual.get(), 0, 0);
        tatami_test::test_indexed_access(params, forced_bound_dense.get(), manual.get(), 0, 1);

        tatami_test::test_full_access(params, forced_bound_sparse.get(), manual.get());
        tatami_test::test_block_access(params, forced_bound_sparse.get(), manual.get(), 0, 0);
        tatami_test::test_indexed_access(params, forced_bound_sparse.get(), manual.get(), 0, 1);
    }
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindEmptyAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // bind by row or by column
        ::testing::Values(true, false), // access by row or column
        ::testing::Values(true, false)  // whether to use an oracle
    )
);

/****************************
 ****************************/

class DelayedBindFullAccessTest : public ::testing::TestWithParam<std::tuple<typename DelayedBindUtils::SimulationParameters, tatami_test::StandardTestAccessParameters> >, public DelayedBindUtils {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(DelayedBindFullAccessTest, Basic) {
    auto tparam = GetParam();
    auto params = tatami_test::convert_access_parameters(std::get<1>(tparam));

    tatami_test::test_full_access(params, bound_sparse.get(), manual.get());
    tatami_test::test_full_access(params, bound_dense.get(), manual.get());

    if (params.use_oracle) {
        tatami_test::test_full_access(params, forced_bound_dense.get(), manual.get());
        tatami_test::test_full_access(params, forced_bound_sparse.get(), manual.get());
    }
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindFullAccessTest,
    ::testing::Combine(
        DelayedBindUtils::spawn_bind_scenarios(),
        tatami_test::standard_test_access_parameter_combinations()
    )
);

/****************************
 ****************************/

class DelayedBindSlicedAccessTest : public ::testing::TestWithParam<std::tuple<typename DelayedBindUtils::SimulationParameters, tatami_test::StandardTestAccessParameters, std::pair<double, double> > >, public DelayedBindUtils {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(DelayedBindSlicedAccessTest, Basic) {
    auto tparam = GetParam();
    auto params = tatami_test::convert_access_parameters(std::get<1>(tparam));

    auto interval_info = std::get<2>(tparam);
    auto len = (params.use_row ?  manual->ncol() : manual->nrow());
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;

    tatami_test::test_block_access(params, bound_sparse.get(), manual.get(), FIRST, LAST);
    tatami_test::test_block_access(params, bound_dense.get(), manual.get(), FIRST, LAST);

    if (params.use_oracle) {
        tatami_test::test_block_access(params, forced_bound_sparse.get(), manual.get(), FIRST, LAST);
        tatami_test::test_block_access(params, forced_bound_dense.get(), manual.get(), FIRST, LAST);
    }
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindSlicedAccessTest,
    ::testing::Combine(
        DelayedBindUtils::spawn_bind_scenarios(),
        tatami_test::standard_test_access_parameter_combinations(),
        ::testing::Values(
            std::make_pair(0, 0.6), 
            std::make_pair(0.25, 0.75), 
            std::make_pair(0.55, 1)
        )
    )
);

/****************************
 ****************************/

class DelayedBindIndexedAccessTest : public ::testing::TestWithParam<std::tuple<typename DelayedBindUtils::SimulationParameters, tatami_test::StandardTestAccessParameters, std::pair<double, int> > >, public DelayedBindUtils {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(DelayedBindIndexedAccessTest, Basic) {
    auto tparam = GetParam();
    auto params = tatami_test::convert_access_parameters(std::get<1>(tparam));

    auto interval_info = std::get<2>(tparam);
    size_t len = (params.use_row ? manual->ncol() : manual->nrow());
    size_t FIRST = interval_info.first * len, STEP = interval_info.second;

    tatami_test::test_indexed_access(params, bound_sparse.get(), manual.get(), FIRST, STEP);
    tatami_test::test_indexed_access(params, bound_dense.get(), manual.get(), FIRST, STEP);

    if (params.use_oracle) {
        tatami_test::test_indexed_access(params, forced_bound_sparse.get(), manual.get(), FIRST, STEP);
        tatami_test::test_indexed_access(params, forced_bound_dense.get(), manual.get(), FIRST, STEP);
    }
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindIndexedAccessTest,
    ::testing::Combine(
        DelayedBindUtils::spawn_bind_scenarios(),
        tatami_test::standard_test_access_parameter_combinations(),
        ::testing::Values(
            std::pair<double, int>(0, 5), 
            std::pair<double, int>(0.33, 3),
            std::pair<double, int>(0.5, 4)
        )
    )
);
