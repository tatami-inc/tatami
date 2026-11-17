#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <cstdint>
#include <random>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/bind/DelayedBind.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

class DelayedBindCore {
public:
    typedef std::tuple<std::vector<int>, bool> SimulationParameters;

    static auto spawn_bind_scenarios () {
        return ::testing::Combine(
            ::testing::Values(
                std::vector<int>{ 89 },
                std::vector<int>{ 78, 20 },
                std::vector<int>{ 54, 37 },
                std::vector<int>{ 52, 7, 23 },
                std::vector<int>{ 15, 41, 92 },
                std::vector<int>{ 25, 0, 77 }
            ),
            ::testing::Values(true, false) // bind by row or by column
        );
    }

protected:
    inline static int otherdim = 97;
    inline static std::shared_ptr<tatami::NumericMatrix> bound_dense, bound_sparse, manual;
    inline static std::shared_ptr<tatami::NumericMatrix> forced_bound_dense, forced_bound_sparse;
    inline static std::shared_ptr<tatami::NumericMatrix> uns_bound_dense, uns_bound_sparse;
    inline static SimulationParameters last_params;

    static void assemble(const SimulationParameters& sim_params) {
        if (bound_dense && sim_params == last_params) {
            return;
        }
        last_params = sim_params;

        const auto& lengths = std::get<0>(sim_params);
        bool row = std::get<1>(sim_params);

        std::vector<double> concat;
        int n_total = 0;
        std::vector<std::shared_ptr<tatami::NumericMatrix> > collected_dense, collected_sparse;
        std::vector<std::shared_ptr<tatami::NumericMatrix> > forced_collected_dense, forced_collected_sparse;
        std::vector<std::shared_ptr<tatami::NumericMatrix> > uns_collected_sparse;

        const std::size_t num_mats = lengths.size();
        for (std::size_t i = 0; i < num_mats; ++i) {
            auto to_add = tatami_test::simulate_vector<double>(lengths[i], otherdim, [&]{
                tatami_test::SimulateVectorOptions opt;
                opt.density = 0.2;
                opt.lower = -10;
                opt.upper = 10;
                opt.seed = tatami_test::create_seed(i, sim_params);
                return opt;
            }());

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
            uns_collected_sparse.emplace_back(std::make_shared<tatami_test::ReversedIndicesWrapper<double, int> >(collected_sparse.back()));
        }

        bound_dense.reset(new tatami::DelayedBind(std::move(collected_dense), row));
        bound_sparse.reset(new tatami::DelayedBind(std::move(collected_sparse), row));
        forced_bound_dense.reset(new tatami::DelayedBind(std::move(forced_collected_dense), row));
        forced_bound_sparse.reset(new tatami::DelayedBind(std::move(forced_collected_sparse), row));
        uns_bound_sparse.reset(new tatami::DelayedBind(std::move(uns_collected_sparse), row));

        // If there aren't any matrices to bind, DelayedBind defaults to a 0-by-0 matrix.
        const int actual_otherdim = (lengths.empty() ? 0 : otherdim);
        if (row) {
            manual.reset(new tatami::DenseRowMatrix<double, int>(n_total, actual_otherdim, std::move(concat)));
        } else {
            manual.reset(new tatami::DenseColumnMatrix<double, int>(actual_otherdim, n_total, std::move(concat)));
        }
    }
};

/****************************
 ****************************/

class DelayedBindBasicTest : public ::testing::Test, public DelayedBindCore {};

TEST_F(DelayedBindBasicTest, ByRow) {
    assemble(SimulationParameters({ 10, 20, 5 }, true));

    EXPECT_EQ(bound_dense->nrow(), 35);
    EXPECT_EQ(bound_dense->ncol(), otherdim);
    EXPECT_FALSE(bound_dense->is_sparse());
    EXPECT_EQ(bound_dense->is_sparse_proportion(), 0);
    EXPECT_TRUE(bound_dense->prefer_rows());
    EXPECT_EQ(bound_dense->prefer_rows_proportion(), 1);

    EXPECT_EQ(bound_sparse->nrow(), 35);
    EXPECT_EQ(bound_sparse->ncol(), otherdim);
    EXPECT_TRUE(bound_sparse->is_sparse());
    EXPECT_EQ(bound_sparse->is_sparse_proportion(), 1);
    EXPECT_FALSE(bound_sparse->prefer_rows());
    EXPECT_EQ(bound_sparse->prefer_rows_proportion(), 0);

    EXPECT_FALSE(bound_dense->uses_oracle(true));
    EXPECT_FALSE(bound_sparse->uses_oracle(true));
}

TEST_F(DelayedBindBasicTest, ByColumn) {
    assemble(SimulationParameters({ 10, 20, 5 }, false));

    EXPECT_EQ(bound_dense->nrow(), otherdim);
    EXPECT_EQ(bound_dense->ncol(), 35);
    EXPECT_FALSE(bound_dense->is_sparse());
    EXPECT_FALSE(bound_dense->prefer_rows());

    EXPECT_EQ(bound_sparse->nrow(), otherdim);
    EXPECT_EQ(bound_sparse->ncol(), 35);
    EXPECT_TRUE(bound_sparse->is_sparse());
    EXPECT_FALSE(bound_sparse->prefer_rows());

    EXPECT_FALSE(bound_dense->uses_oracle(false));
    EXPECT_FALSE(bound_sparse->uses_oracle(false));
}

TEST_F(DelayedBindBasicTest, InconsistentBinds) {
    assemble(SimulationParameters({ 10, 20, 5 }, true));

    // Bound_sparse is CSC, bound_dense is row-major.
    tatami::DelayedBind combined(std::vector<std::shared_ptr<tatami::NumericMatrix> >{ bound_sparse, bound_dense }, false);

    EXPECT_FLOAT_EQ(combined.is_sparse_proportion(), 0.5);
    EXPECT_FALSE(combined.is_sparse());

    EXPECT_FLOAT_EQ(combined.prefer_rows_proportion(), 0.5);
    EXPECT_FALSE(combined.prefer_rows());
}

TEST_F(DelayedBindBasicTest, ConstOverloads) {
    assemble(SimulationParameters({ 10, 50 }, true));

    std::vector<std::shared_ptr<const tatami::NumericMatrix> > const_collected({ bound_dense, bound_sparse });
    tatami::DelayedBind const_combined(std::move(const_collected), true);

    // Some cursory checks.
    EXPECT_EQ(const_combined.nrow(), 120); // i.e., (10 + 50) * 2 
    EXPECT_EQ(const_combined.ncol(), otherdim);
}

TEST_F(DelayedBindBasicTest, MakeHelpers) {
    assemble(SimulationParameters({ 10, 50 }, true));

    std::vector<std::shared_ptr<tatami::NumericMatrix> > collected({ bound_dense, bound_sparse });
    auto tdense1a = tatami::make_DelayedBind(collected, true);
    EXPECT_EQ(tdense1a->nrow(), 120);
    EXPECT_EQ(tdense1a->ncol(), otherdim);

    auto tdense1b = tatami::make_DelayedBind<0>(collected);
    EXPECT_EQ(tdense1b->nrow(), 120);
    EXPECT_EQ(tdense1b->ncol(), otherdim);

    std::vector<std::shared_ptr<const tatami::NumericMatrix> > const_collected({ bound_dense, bound_sparse });
    auto tdense2a = tatami::make_DelayedBind(const_collected, false);
    EXPECT_EQ(tdense2a->nrow(), 60);
    EXPECT_EQ(tdense2a->ncol(), otherdim * 2);

    auto tdense2b = tatami::make_DelayedBind<1>(const_collected);
    EXPECT_EQ(tdense2b->nrow(), 60);
    EXPECT_EQ(tdense2b->ncol(), otherdim * 2);
}

TEST(DelayedBind, ErrorCheck) {
    std::vector<std::shared_ptr<tatami::NumericMatrix> > collected;
    collected.emplace_back(new tatami::DenseRowMatrix<double, int>(10, 20, std::vector<double>(200)));
    collected.emplace_back(new tatami::DenseRowMatrix<double, int>(20, 10, std::vector<double>(200)));
    tatami_test::throws_error([&]() { tatami::DelayedBind(collected, true); }, "same number of columns");
    tatami_test::throws_error([&]() { tatami::DelayedBind(collected, false); }, "same number of rows");
}

TEST(DelayedBind, PartialOracleUsage) {
    auto p1 = std::make_shared<tatami::DenseRowMatrix<double, int> >(10, 20, std::vector<double>(200));
    auto p2 = std::make_shared<tatami::DenseRowMatrix<double, int> >(10, 20, std::vector<double>(200));

    {
        tatami::DelayedBind combined(std::vector<std::shared_ptr<tatami::Matrix<double, int> > >{ p1, p2 }, true); 
        EXPECT_FALSE(combined.uses_oracle(true));
        EXPECT_FALSE(combined.uses_oracle(false));
    }

    // Only one component uses an oracle.
    {
        tatami::DelayedBind combined(std::vector<std::shared_ptr<tatami::Matrix<double, int> > >{
            std::make_shared<tatami_test::ForcedOracleWrapper<double, int> >(p1),
            p2
        }, true);
        EXPECT_TRUE(combined.uses_oracle(true));
        EXPECT_TRUE(combined.uses_oracle(false));
    }

    // Both components use an oracle.
    {
        tatami::DelayedBind combined(std::vector<std::shared_ptr<tatami::Matrix<double, int> > >{
            std::make_shared<tatami_test::ForcedOracleWrapper<double, int> >(p1),
            std::make_shared<tatami_test::ForcedOracleWrapper<double, int> >(p2)
        }, true);
        EXPECT_TRUE(combined.uses_oracle(true));
        EXPECT_TRUE(combined.uses_oracle(false));
    }
}

/****************************
 ****************************/

class DelayedBindNoMatrixTest : 
    public ::testing::TestWithParam<std::tuple<bool, bool, tatami_test::TestAccessOptionsCombination> >,
    public DelayedBindCore
{
protected:
    void SetUp() {
        assemble({ std::vector<int>{}, std::get<0>(GetParam()) });
    }
};

TEST_P(DelayedBindNoMatrixTest, Full) {
    auto tparam = GetParam();
    const bool use_row = std::get<1>(tparam);
    auto options = tatami_test::convert_to_test_access_options(std::get<2>(tparam));

    tatami_test::test_access(*bound_dense, *manual, use_row, options);
    tatami_test::test_access(*bound_sparse, *manual, use_row, options);

    if (options.use_oracle) {
        tatami_test::test_access(*forced_bound_dense, *manual, use_row, options);
        tatami_test::test_access(*forced_bound_sparse, *manual, use_row, options);
    }

    // No need to test unsorted sparse extraction here.
    // We don't even run any extractions so there's no chance to even consider whether the indices are unsorted or not.
}

TEST_P(DelayedBindNoMatrixTest, Block) {
    auto tparam = GetParam();
    const bool use_row = std::get<1>(tparam);
    auto options = tatami_test::convert_to_test_access_options(std::get<2>(tparam));

    tatami_test::test_access(*bound_dense, *manual, use_row, 0, 0, options);
    tatami_test::test_access(*bound_sparse, *manual, use_row, 0, 0, options);

    if (options.use_oracle) {
        tatami_test::test_access(*forced_bound_dense, *manual, use_row, 0, 0, options);
        tatami_test::test_access(*forced_bound_sparse, *manual, use_row, 0, 0, options);
    }
}

TEST_P(DelayedBindNoMatrixTest, Index) {
    auto tparam = GetParam();
    const bool use_row = std::get<1>(tparam);
    auto options = tatami_test::convert_to_test_access_options(std::get<2>(tparam));

    auto empty = std::make_shared<const std::vector<int> >();
    tatami_test::test_access(*bound_dense, *manual, use_row, empty, options);
    tatami_test::test_access(*bound_sparse, *manual, use_row, empty, options);

    if (options.use_oracle) {
        tatami_test::test_access(*forced_bound_dense, *manual, use_row, empty, options);
        tatami_test::test_access(*forced_bound_sparse, *manual, use_row, empty, options);
    }
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindNoMatrixTest,
    ::testing::Combine(
        ::testing::Values(true, false), // bind by row
        ::testing::Values(true, false), // access by row
        tatami_test::combine_test_access_options()
    )
);

/****************************
 ****************************/

class DelayedBindZeroExtentTest : 
    public ::testing::TestWithParam<std::tuple<bool, bool, tatami_test::TestAccessOptionsCombination> >,
    public DelayedBindCore
{
protected:
    void SetUp() {
        assemble({ std::vector<int>{ 0, 0 }, std::get<0>(GetParam()) });
    }
};

TEST_P(DelayedBindZeroExtentTest, Full) {
    auto tparam = GetParam();
    const bool use_row = std::get<1>(tparam);
    auto options = tatami_test::convert_to_test_access_options(std::get<2>(tparam));

    tatami_test::test_access(*bound_dense, *manual, use_row, options);
    tatami_test::test_access(*bound_sparse, *manual, use_row, options);

    if (options.use_oracle) {
        tatami_test::test_access(*forced_bound_dense, *manual, use_row, options);
        tatami_test::test_access(*forced_bound_sparse, *manual, use_row, options);
    }

    tatami_test::test_unsorted_access(*uns_bound_sparse, use_row, options);
}

TEST_P(DelayedBindZeroExtentTest, Block) {
    auto tparam = GetParam();
    const bool use_row = std::get<1>(tparam);
    auto options = tatami_test::convert_to_test_access_options(std::get<2>(tparam));

    tatami_test::test_access(*bound_dense, *manual, use_row, 0, 0, options);
    tatami_test::test_access(*bound_sparse, *manual, use_row, 0, 0, options);

    if (options.use_oracle) {
        tatami_test::test_access(*forced_bound_dense, *manual, use_row, 0, 0, options);
        tatami_test::test_access(*forced_bound_sparse, *manual, use_row, 0, 0, options);
    }

    tatami_test::test_unsorted_access(*uns_bound_sparse, use_row, 0, 0, options);
}

TEST_P(DelayedBindZeroExtentTest, Index) {
    auto tparam = GetParam();
    const bool use_row = std::get<1>(tparam);
    auto options = tatami_test::convert_to_test_access_options(std::get<2>(tparam));

    auto empty = std::make_shared<const std::vector<int> >();
    tatami_test::test_access(*bound_dense, *manual, use_row, empty, options);
    tatami_test::test_access(*bound_sparse, *manual, use_row, empty, options);

    if (options.use_oracle) {
        tatami_test::test_access(*forced_bound_dense, *manual, use_row, empty, options);
        tatami_test::test_access(*forced_bound_sparse, *manual, use_row, empty, options);
    }

    tatami_test::test_unsorted_access(*uns_bound_sparse, use_row, empty, options);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindZeroExtentTest,
    ::testing::Combine(
        ::testing::Values(true, false), // bind by row or by column
        ::testing::Values(true, false), // access by row or column
        tatami_test::combine_test_access_options()
    )
);

/****************************
 ****************************/

class DelayedBindFullTest : 
    public ::testing::TestWithParam<std::tuple<typename DelayedBindCore::SimulationParameters, bool, tatami_test::TestAccessOptionsCombination> >,
    public DelayedBindCore 
{
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(DelayedBindFullTest, Basic) {
    auto tparam = GetParam();
    auto use_row = std::get<1>(tparam);
    auto options = tatami_test::convert_to_test_access_options(std::get<2>(tparam));

    tatami_test::test_access(*bound_dense, *manual, use_row, options);
    tatami_test::test_access(*bound_sparse, *manual, use_row, options);

    if (options.use_oracle) {
        tatami_test::test_access(*forced_bound_dense, *manual, use_row, options);
        tatami_test::test_access(*forced_bound_sparse, *manual, use_row, options);
    }

    tatami_test::test_unsorted_access(*uns_bound_sparse, use_row, options);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindFullTest,
    ::testing::Combine(
        DelayedBindCore::spawn_bind_scenarios(),
        ::testing::Values(false, true),
        tatami_test::combine_test_access_options()
    )
);

/****************************
 ****************************/

class DelayedBindBlockTest : 
    public ::testing::TestWithParam<
        std::tuple<
            typename DelayedBindCore::SimulationParameters,
            bool,
            std::pair<double, double>,
            tatami_test::TestAccessOptionsCombination
        >
    >,
    public DelayedBindCore
{
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(DelayedBindBlockTest, Basic) {
    auto tparam = GetParam();
    auto use_row = std::get<1>(tparam);
    auto raw_block = std::get<2>(tparam);
    auto block = tatami_test::define_non_target_block(*manual, use_row, raw_block.first, raw_block.second);
    auto options = tatami_test::convert_to_test_access_options(std::get<3>(tparam)); 

    tatami_test::test_access(*bound_dense, *manual, use_row, block.first, block.second, options);
    tatami_test::test_access(*bound_sparse, *manual, use_row, block.first, block.second, options);

    if (options.use_oracle) {
        tatami_test::test_access(*forced_bound_dense, *manual, use_row, block.first, block.second, options);
        tatami_test::test_access(*forced_bound_sparse, *manual, use_row, block.first, block.second, options);
    }

    tatami_test::test_unsorted_access(*uns_bound_sparse, use_row, block.first, block.second, options);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindBlockTest,
    ::testing::Combine(
        DelayedBindCore::spawn_bind_scenarios(),
        ::testing::Values(false, true),
        ::testing::Values(
            std::make_pair(0.0, 0.6),
            std::make_pair(0.25, 0.51),
            std::make_pair(0.55, 0.45)
        ),
        tatami_test::combine_test_access_options()
    )
);

/****************************
 ****************************/

class DelayedBindIndexTest :
    public ::testing::TestWithParam<
        std::tuple<
            typename DelayedBindCore::SimulationParameters,
            bool,
            std::tuple<double, double, double>,
            tatami_test::TestAccessOptionsCombination
        >
    >,
    public DelayedBindCore
{
protected:
    void SetUp() {
        assemble(SimulationParameters(std::get<0>(GetParam())));
    }
};

TEST_P(DelayedBindIndexTest, Basic) {
    auto tparam = GetParam();
    bool use_row = std::get<1>(tparam);
    auto raw_index = std::get<2>(tparam);
    auto indices_ptr = tatami_test::define_non_target_indices(*manual, use_row, std::get<0>(raw_index), std::get<1>(raw_index), std::get<2>(raw_index), {});
    auto options = tatami_test::convert_to_test_access_options(std::get<3>(tparam)); 

    tatami_test::test_access(*bound_dense, *manual, use_row, indices_ptr, options);
    tatami_test::test_access(*bound_sparse, *manual, use_row, indices_ptr, options);

    if (options.use_oracle) {
        tatami_test::test_access(*forced_bound_dense, *manual, use_row, indices_ptr, options);
        tatami_test::test_access(*forced_bound_sparse, *manual, use_row, indices_ptr, options);
    }

    tatami_test::test_unsorted_access(*uns_bound_sparse, use_row, indices_ptr, options);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindIndexTest,
    ::testing::Combine(
        DelayedBindCore::spawn_bind_scenarios(),
        ::testing::Values(false, true),
        ::testing::Values(
            std::make_tuple(0.0, 0.45, 0.25),
            std::make_tuple(0.25, 0.55, 0.2),
            std::make_tuple(0.6, 0.4, 0.3)
        ),
        tatami_test::combine_test_access_options()
    )
);

/****************************
 ****************************/

struct ParallelBlockByComponent {
    ParallelBlockByComponent(int start_mat_arg, double start_prop_arg, int end_mat_arg, double end_prop_arg) :
        start_mat(start_mat_arg), start_prop(start_prop_arg), end_mat(end_mat_arg), end_prop(end_prop_arg)
    {}
    int start_mat;
    double start_prop;
    int end_mat;
    double end_prop;
};

static std::pair<int, int> derive_parallel_block(const ParallelBlockByComponent& params, const std::vector<int>& sizes) {
    const int start_offset = std::accumulate(sizes.begin(), sizes.begin() + params.start_mat, 0); 
    const int end_offset = start_offset + std::accumulate(sizes.begin() + params.start_mat, sizes.begin() + params.end_mat, 0); 
    const auto block_start = start_offset + std::round(sizes[params.start_mat] * params.start_prop);
    const auto block_end = end_offset + std::round(sizes[params.end_mat] * params.end_prop);
    return std::pair<int, int>(block_start, block_end - block_start);
}

struct ParallelIndicesByComponent {
    ParallelIndicesByComponent(int start_mat_arg, double start_prop_arg, int end_mat_arg, double end_prop_arg, double probability_arg) :
        range(start_mat_arg, start_prop_arg, end_mat_arg, end_prop_arg), probability(probability_arg)
    {}
    ParallelBlockByComponent range;
    double probability;
};

static tatami::VectorPtr<int> derive_parallel_indices(const ParallelIndicesByComponent& params, const std::vector<int>& sizes) {
    auto block = derive_parallel_block(params.range, sizes);
    return tatami_test::define_non_target_indices(block.first, block.second, params.probability, {});
}

/****************************
 ****************************/

// Check parallel partial access in a DelayedBind with two component matrices.
// This checks the optimizations when only one component is used.

class DelayedBindTwoMatrixParallelBlockTest :
    public ::testing::TestWithParam<std::tuple<bool, ParallelBlockByComponent, tatami_test::TestAccessOptionsCombination> >,
    public DelayedBindCore
{
protected:
    void SetUp() {
        assemble(SimulationParameters(std::vector<int>{ 61, 47 }, std::get<0>(GetParam())));
    }
};

TEST_P(DelayedBindTwoMatrixParallelBlockTest, Basic) {
    auto tparam = GetParam();
    auto block = derive_parallel_block(std::get<1>(tparam), std::get<0>(last_params));
    auto options = tatami_test::convert_to_test_access_options(std::get<2>(tparam)); 

    // Access needs to be parallel to the bind dimension.
    auto use_row = !std::get<1>(last_params);

    tatami_test::test_access(*bound_dense, *manual, use_row, block.first, block.second, options);
    tatami_test::test_access(*bound_sparse, *manual, use_row, block.first, block.second, options);

    if (options.use_oracle) {
        tatami_test::test_access(*forced_bound_dense, *manual, use_row, block.first, block.second, options);
        tatami_test::test_access(*forced_bound_sparse, *manual, use_row, block.first, block.second, options);
    }

    tatami_test::test_unsorted_access(*uns_bound_sparse, use_row, block.first, block.second, options);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindTwoMatrixParallelBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false),
        ::testing::Values(
            // First matrix only.
            ParallelBlockByComponent(0, 0.0, 0, 0.4),
            ParallelBlockByComponent(0, 0.3, 0, 1.0),
            // Second matrix only.
            ParallelBlockByComponent(1, 0.0, 1, 0.5),
            ParallelBlockByComponent(1, 0.4, 1, 1.0),
            // Both matrices, as a control.
            ParallelBlockByComponent(0, 0.2, 1, 0.3),
            ParallelBlockByComponent(0, 0.7, 1, 0.8)
        ),
        tatami_test::combine_test_access_options()
    )
);

class DelayedBindTwoMatrixParallelIndexTest :
    public ::testing::TestWithParam<std::tuple<bool, ParallelIndicesByComponent, tatami_test::TestAccessOptionsCombination> >,
    public DelayedBindCore
{
protected:
    void SetUp() {
        assemble(SimulationParameters(std::vector<int>{ 38, 62 }, std::get<0>(GetParam())));
    }
};

TEST_P(DelayedBindTwoMatrixParallelIndexTest, Basic) {
    auto tparam = GetParam();
    auto indices_ptr = derive_parallel_indices(std::get<1>(tparam), std::get<0>(last_params));
    auto options = tatami_test::convert_to_test_access_options(std::get<2>(tparam));

    // Access needs to be parallel to the bind dimension.
    auto use_row = !std::get<1>(last_params);

    tatami_test::test_access(*bound_dense, *manual, use_row, indices_ptr, options);
    tatami_test::test_access(*bound_sparse, *manual, use_row, indices_ptr, options);

    if (options.use_oracle) {
        tatami_test::test_access(*forced_bound_dense, *manual, use_row, indices_ptr, options);
        tatami_test::test_access(*forced_bound_sparse, *manual, use_row, indices_ptr, options);
    }

    tatami_test::test_unsorted_access(*uns_bound_sparse, use_row, indices_ptr, options);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindTwoMatrixParallelIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false),
        ::testing::Values(
            // First matrix only.
            ParallelIndicesByComponent(0, 0.0, 0, 0.6, 0.3),
            ParallelIndicesByComponent(0, 0.5, 0, 1.0, 0.3),
            // Second matrix only.
            ParallelIndicesByComponent(1, 0.0, 1, 0.4, 0.25),
            ParallelIndicesByComponent(1, 0.6, 1, 1.0, 0.25),
            // Both matrices, as a control.
            ParallelIndicesByComponent(0, 0.2, 1, 0.8, 0.3),
            ParallelIndicesByComponent(0, 0.4, 1, 0.4, 0.3)
        ),
        tatami_test::combine_test_access_options()
    )
);

/****************************
 ****************************/

// Check parallel partial access in a DelayedBind with three component matrices.
// This checks the optimizations when only one component is used.

class DelayedBindThreeMatrixParallelBlockTest :
    public ::testing::TestWithParam<std::tuple<bool, ParallelBlockByComponent, tatami_test::TestAccessOptionsCombination> >,
    public DelayedBindCore
{
protected:
    void SetUp() {
        assemble(SimulationParameters(std::vector<int>{ 22, 47, 131 }, std::get<0>(GetParam())));
    }
};

TEST_P(DelayedBindThreeMatrixParallelBlockTest, Basic) {
    auto tparam = GetParam();
    auto block = derive_parallel_block(std::get<1>(tparam), std::get<0>(last_params));
    auto options = tatami_test::convert_to_test_access_options(std::get<2>(tparam)); 

    // Access needs to be parallel to the bind dimension.
    auto use_row = !std::get<1>(last_params);

    tatami_test::test_access(*bound_dense, *manual, use_row, block.first, block.second, options);
    tatami_test::test_access(*bound_sparse, *manual, use_row, block.first, block.second, options);

    if (options.use_oracle) {
        tatami_test::test_access(*forced_bound_dense, *manual, use_row, block.first, block.second, options);
        tatami_test::test_access(*forced_bound_sparse, *manual, use_row, block.first, block.second, options);
    }

    tatami_test::test_unsorted_access(*uns_bound_sparse, use_row, block.first, block.second, options);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindThreeMatrixParallelBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false),
        ::testing::Values(
            // First matrix only.
            ParallelBlockByComponent(0, 0.0, 0, 0.7),
            ParallelBlockByComponent(0, 0.7, 0, 1.0),
            // Second matrix only.
            ParallelBlockByComponent(1, 0.0, 1, 0.4),
            ParallelBlockByComponent(1, 0.5, 1, 1.0),
            // Third matrix only.
            ParallelBlockByComponent(2, 0.0, 2, 0.3),
            ParallelBlockByComponent(2, 0.5, 2, 1.0),
            // Multiple matrices, as a control.
            ParallelBlockByComponent(0, 0.5, 1, 0.4), // First two matrices.
            ParallelBlockByComponent(1, 0.2, 2, 0.8), // Last two matrices.
            ParallelBlockByComponent(0, 0.6, 2, 0.3) // All matrices.
        ),
        tatami_test::combine_test_access_options()
    )
);

class DelayedBindThreeMatrixParallelIndexTest :
    public ::testing::TestWithParam<std::tuple<bool, ParallelIndicesByComponent, tatami_test::TestAccessOptionsCombination> >,
    public DelayedBindCore
{
protected:
    void SetUp() {
        assemble(SimulationParameters(std::vector<int>{ 120, 60, 20 }, std::get<0>(GetParam())));
    }
};

TEST_P(DelayedBindThreeMatrixParallelIndexTest, Basic) {
    auto tparam = GetParam();
    auto indices_ptr = derive_parallel_indices(std::get<1>(tparam), std::get<0>(last_params));
    auto options = tatami_test::convert_to_test_access_options(std::get<2>(tparam));

    // Access needs to be parallel to the bind dimension.
    auto use_row = !std::get<1>(last_params);

    tatami_test::test_access(*bound_dense, *manual, use_row, indices_ptr, options);
    tatami_test::test_access(*bound_sparse, *manual, use_row, indices_ptr, options);

    if (options.use_oracle) {
        tatami_test::test_access(*forced_bound_dense, *manual, use_row, indices_ptr, options);
        tatami_test::test_access(*forced_bound_sparse, *manual, use_row, indices_ptr, options);
    }

    tatami_test::test_unsorted_access(*uns_bound_sparse, use_row, indices_ptr, options);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindThreeMatrixParallelIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false),
        ::testing::Values(
            // First matrix only.
            ParallelIndicesByComponent(0, 0.0, 0, 0.2, 0.5),
            ParallelIndicesByComponent(0, 0.4, 0, 1.0, 0.3),
            // Second matrix only.
            ParallelIndicesByComponent(1, 0.0, 1, 0.8, 0.25),
            ParallelIndicesByComponent(1, 0.7, 1, 1.0, 0.5),
            // Third matrix only.
            ParallelIndicesByComponent(2, 0.0, 2, 0.8, 0.3),
            ParallelIndicesByComponent(2, 0.1, 2, 0.9, 0.3),
            // Multiple matrices, as a control.
            ParallelIndicesByComponent(0, 0.8, 1, 0.5, 0.3), // first two matrices
            ParallelIndicesByComponent(1, 0.2, 2, 0.8, 0.3), // last two matrices
            ParallelIndicesByComponent(0, 0.7, 2, 0.5, 0.3) // All matrices
        ),
        tatami_test::combine_test_access_options()
    )
);

/****************************
 ****************************/

// Check for correct use of a single matrix's extractor when the predictions only cover a single matrix.
class DelayedBindLimitedPredictionsTest : 
    public ::testing::TestWithParam<std::vector<int> >,
    public DelayedBindCore
{};

TEST_P(DelayedBindLimitedPredictionsTest, Basic) {
    auto predictions = GetParam();
    assemble(std::make_tuple(std::vector<int>{ 100, 300, 200 }, true));

    auto observed_d = forced_bound_dense->dense(true, std::make_shared<tatami::FixedVectorOracle<int> >(predictions), {});
    auto observed_s = forced_bound_sparse->sparse(true, std::make_shared<tatami::FixedVectorOracle<int> >(predictions), {});
    auto expected = manual->dense(true, std::make_shared<tatami::FixedVectorOracle<int> >(predictions), {});

    std::vector<double> exp_buffer(otherdim), obs_buffer(otherdim), obs_buffer_d(otherdim);
    std::vector<int> obs_buffer_i(otherdim);
    for ([[maybe_unused]] auto p : predictions) {
        auto optr = observed_d->fetch(obs_buffer.data());
        tatami::copy_n(optr, otherdim, obs_buffer.data());

        auto eptr = expected->fetch(exp_buffer.data());
        tatami::copy_n(eptr, otherdim, exp_buffer.data());
        EXPECT_EQ(exp_buffer, obs_buffer);

        std::fill(obs_buffer.begin(), obs_buffer.end(), 0);
        auto range = observed_s->fetch(obs_buffer_d.data(), obs_buffer_i.data());
        for (int i = 0; i < range.number; ++i) {
            obs_buffer[range.index[i]] = range.value[i];
        }
        EXPECT_EQ(exp_buffer, obs_buffer);
    }
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindLimitedPredictionsTest,
    ::testing::Values(
        std::vector<int>{ 56, 27, 11, 0, 88, 21, 5, 64 }, // Only using the first matrix, so no need to wrap the oracle.
        std::vector<int>{ 115, 223, 210, 190, 381, 174, 165, 324 }, // Using the second matrix, wrap the oracle.
        std::vector<int>{} // Correct behavior with an empty input.
    )
);
