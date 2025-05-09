#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>
#include <random>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/subset/make_DelayedSubset.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

class SubsetCoreUtils {
protected:
    inline static size_t NR = 90, NC = 170;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, sparse;

    static void assemble() {
        if (dense) {
            return;
        }

        auto simulated = tatami_test::simulate_vector<double>(NR * NC, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.1;
            opt.seed = 120983712;
            return opt;
        }());

        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(NR, NC, std::move(simulated)));
        sparse = tatami::convert_to_compressed_sparse<false, double, int>(dense.get()); // column-major.
    }

protected:
    template<typename INDEX>
    static std::vector<INDEX> spawn_indices(size_t step, size_t max, bool duplicates, bool sorted) {
        std::vector<INDEX> output;
        for (size_t i = step; i < max; i += step) {
            output.push_back(i);
        }

        std::mt19937_64 rng(step + max + 10 * duplicates + sorted);

        if (duplicates) {
            for (size_t i = 0, end = output.size(); i < end; ++i) {
                output.insert(output.end(), rng() % 4, output[i]);
            }
            if (sorted) {
                std::sort(output.begin(), output.end());
            }
        }

        if (!sorted) {
            std::shuffle(output.begin(), output.end(), rng);
        }
        return output;
    }

    template<typename V>
    static std::shared_ptr<tatami::NumericMatrix> reference_on_rows(const tatami::NumericMatrix* src, const V& sub) {
        std::vector<double> reference(sub.size() * NC);
        auto ptr = reference.data();
        auto wrk = src->dense_row();

        for (auto r : sub) {
            auto src = wrk->fetch(r, ptr);
            tatami::copy_n(src, NC, ptr);
            ptr += NC;
        }

        return std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(sub.size(), NC, std::move(reference)));
    }

    template<typename V>
    static std::shared_ptr<tatami::NumericMatrix> reference_on_columns(const tatami::NumericMatrix* src, const V& sub) {
        std::vector<double> reference(sub.size() * NR);
        auto ptr = reference.data();
        std::vector<double> buffer(NC);
        auto wrk = src->dense_row();

        for (size_t r = 0; r < NR; ++r) {
            auto full = wrk->fetch(r, buffer.data());
            for (auto s : sub) {
                *ptr = full[s];
                ++ptr;
            }
        }

        return std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(NR, sub.size(), std::move(reference)));
    }
};

class SubsetUtils : public SubsetCoreUtils {
public:
    typedef std::tuple<bool, size_t, bool, bool> SimulationParameters;

    static auto simulation_parameter_combinations() {
        return ::testing::Combine(
            ::testing::Values(true, false), // whether to subset on the rows or columns.
            ::testing::Values(2, 5, 10), // step size.
            ::testing::Values(false, true), // whether to support duplicate indices.
            ::testing::Values(true, false) // whether to require sorted indices.
        );
    }

protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_subbed, sparse_subbed, uns_sparse_subbed, ref;
    inline static SimulationParameters last_params;
    inline static std::vector<size_t> sub;

    static void assemble(SimulationParameters sim_params) {
        if (ref && last_params == sim_params) {
            return;
        }
        last_params = sim_params;

        SubsetCoreUtils::assemble();

        auto byrow = std::get<0>(sim_params);
        auto step_size = std::get<1>(sim_params);
        auto duplicates = std::get<2>(sim_params);
        auto sorted = std::get<3>(sim_params);

        if (byrow) {
            sub = spawn_indices<size_t>(step_size, NR, duplicates, sorted);
            ref = SubsetCoreUtils::reference_on_rows(dense.get(), sub);
            dense_subbed = tatami::make_DelayedSubset<0>(dense, sub);
            sparse_subbed = tatami::make_DelayedSubset<0>(sparse, sub);
            uns_sparse_subbed = tatami::make_DelayedSubset<0, double, int>(std::make_shared<const tatami_test::ReversedIndicesWrapper<double, int> >(sparse), sub);

        } else {
            sub = spawn_indices<size_t>(step_size, NC, duplicates, sorted);
            ref = SubsetCoreUtils::reference_on_columns(dense.get(), sub);
            dense_subbed = tatami::make_DelayedSubset<1>(dense, sub);
            sparse_subbed = tatami::make_DelayedSubset<1>(sparse, sub);
            uns_sparse_subbed = tatami::make_DelayedSubset<1, double, int>(std::make_shared<const tatami_test::ReversedIndicesWrapper<double, int> >(sparse), sub);
        }
    }
};

class SubsetTest : public ::testing::TestWithParam<SubsetUtils::SimulationParameters>, public SubsetUtils {
protected:
    void SetUp() {
        assemble(GetParam());
    }
};

TEST_P(SubsetTest, Basic) {
    auto byrow = std::get<0>(last_params);
    if (byrow) {
        EXPECT_EQ(sub.size(), dense_subbed->nrow());
        EXPECT_EQ(dense->ncol(), dense_subbed->ncol());
    } else {
        EXPECT_EQ(dense->nrow(), dense_subbed->nrow());
        EXPECT_EQ(sub.size(), dense_subbed->ncol());
    }

    EXPECT_EQ(dense->is_sparse(), dense_subbed->is_sparse());
    EXPECT_EQ(dense->is_sparse_proportion(), dense_subbed->is_sparse_proportion());
    EXPECT_EQ(sparse->is_sparse(), sparse_subbed->is_sparse());
    EXPECT_EQ(sparse->is_sparse_proportion(), sparse_subbed->is_sparse_proportion());

    EXPECT_TRUE(dense_subbed->prefer_rows());
    EXPECT_EQ(dense->prefer_rows_proportion(), dense_subbed->prefer_rows_proportion());
    EXPECT_FALSE(sparse_subbed->prefer_rows());
    EXPECT_EQ(sparse->prefer_rows_proportion(), sparse_subbed->prefer_rows_proportion());

    EXPECT_FALSE(dense_subbed->uses_oracle(false));
}

INSTANTIATE_TEST_SUITE_P(
    DelayedSubset,
    SubsetTest,
    SubsetUtils::simulation_parameter_combinations()
);

/****************************************************
 ****************************************************/

class SubsetFullAccessTest : 
    public ::testing::TestWithParam<std::tuple<typename SubsetUtils::SimulationParameters, tatami_test::StandardTestAccessOptions> >, 
    public SubsetUtils {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(SubsetFullAccessTest, Basic) {
    auto tparam = GetParam();
    auto options = tatami_test::convert_test_access_options(std::get<1>(tparam));
    tatami_test::test_full_access(*dense_subbed, *ref, options);
    tatami_test::test_full_access(*sparse_subbed, *ref, options);
    tatami_test::test_unsorted_full_access(*uns_sparse_subbed, options);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedSubset,
    SubsetFullAccessTest,
    ::testing::Combine(
        SubsetUtils::simulation_parameter_combinations(),
        tatami_test::standard_test_access_options_combinations()
    )
);

/****************************************************
 ****************************************************/

class SubsetBlockAccessTest : 
    public ::testing::TestWithParam<std::tuple<typename SubsetUtils::SimulationParameters, tatami_test::StandardTestAccessOptions, std::pair<double, double> > >, 
    public SubsetUtils {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(SubsetBlockAccessTest, Block) {
    auto tparam = GetParam();
    auto options = tatami_test::convert_test_access_options(std::get<1>(tparam));
    auto interval_info = std::get<2>(tparam);
    tatami_test::test_block_access(*dense_subbed, *ref, interval_info.first, interval_info.second, options);
    tatami_test::test_block_access(*sparse_subbed, *ref, interval_info.first, interval_info.second, options);
    tatami_test::test_unsorted_block_access(*uns_sparse_subbed, interval_info.first, interval_info.second, options);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedSubset,
    SubsetBlockAccessTest,
    ::testing::Combine(
        SubsetUtils::simulation_parameter_combinations(),
        tatami_test::standard_test_access_options_combinations(),
        ::testing::Values(
            std::make_pair(0.0, 0.66), 
            std::make_pair(0.33, 0.457),
            std::make_pair(0.5, 0.5)
        )
    )
);

/****************************************************
 ****************************************************/

class SubsetIndexedAccessTest : 
    public ::testing::TestWithParam<std::tuple<typename SubsetUtils::SimulationParameters, tatami_test::StandardTestAccessOptions, std::pair<double, double> > >,
    public SubsetUtils {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(SubsetIndexedAccessTest, Indexed) {
    auto tparam = GetParam();
    auto options = tatami_test::convert_test_access_options(std::get<1>(tparam));
    auto interval_info = std::get<2>(tparam);
    tatami_test::test_indexed_access(*dense_subbed, *ref, interval_info.first, interval_info.second, options);
    tatami_test::test_indexed_access(*sparse_subbed, *ref, interval_info.first, interval_info.second, options);
    tatami_test::test_unsorted_indexed_access(*uns_sparse_subbed, interval_info.first, interval_info.second, options);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedSubset,
    SubsetIndexedAccessTest,
    ::testing::Combine(
        SubsetUtils::simulation_parameter_combinations(),
        tatami_test::standard_test_access_options_combinations(),
        ::testing::Values(
            std::make_pair(0.0, 0.2), 
            std::make_pair(0.33, 0.3),
            std::make_pair(0.5, 0.4)
        )
    )
);

/****************************************************
 ****************************************************/

class SubsetConstructorTest : 
    public ::testing::TestWithParam<std::tuple<bool, bool> >, 
    public SubsetCoreUtils {};

TEST_P(SubsetConstructorTest, SortedUnique) {
    auto param = GetParam();
    bool duplicate = std::get<0>(param);
    bool sorted = std::get<1>(param);
    auto sub = SubsetCoreUtils::spawn_indices<int>(true, 5, duplicate, sorted);

    if (sorted && !duplicate) {
        tatami::DelayedSubsetSortedUnique<double, int, decltype(sub)> manual(dense, sub, true);
        auto ref = SubsetCoreUtils::reference_on_rows(dense.get(), sub);
        tatami_test::test_simple_row_access(manual, *ref);
        tatami_test::test_simple_column_access(manual, *ref);
    } else {
        tatami_test::throws_error([&]() {
            tatami::DelayedSubsetSortedUnique<double, int, decltype(sub)> manual(dense, sub, true);
        }, "unique");
    }
}

TEST_P(SubsetConstructorTest, Sorted) {
    auto param = GetParam();
    bool duplicate = std::get<0>(param);
    bool sorted = std::get<1>(param);
    auto sub = SubsetCoreUtils::spawn_indices<int>(true, 5, duplicate, sorted);

    if (sorted) {
        tatami::DelayedSubsetSorted<double, int, decltype(sub)> manual(dense, sub, true);
        auto ref = SubsetCoreUtils::reference_on_rows(dense.get(), sub);
        tatami_test::test_simple_row_access(manual, *ref);
        tatami_test::test_simple_column_access(manual, *ref);
    } else {
        tatami_test::throws_error([&]() {
            tatami::DelayedSubsetSorted<double, int, decltype(sub)> manual(dense, sub, true); 
        }, "sorted");
    }
}

TEST_P(SubsetConstructorTest, Unique) {
    auto param = GetParam();
    bool duplicate = std::get<0>(param);
    bool sorted = std::get<1>(param);
    auto sub = SubsetCoreUtils::spawn_indices<int>(true, 5, duplicate, sorted);

    if (!duplicate) {
        tatami::DelayedSubsetUnique<double, int, decltype(sub)> manual(dense, sub, true);
        auto ref = SubsetCoreUtils::reference_on_rows(dense.get(), sub);
        tatami_test::test_simple_row_access(manual, *ref);
        tatami_test::test_simple_column_access(manual, *ref);
    } else {
        tatami_test::throws_error([&]() {
            tatami::DelayedSubsetUnique<double, int, decltype(sub)> manual(dense, sub, true);
        }, "unique");
    }
}

TEST_P(SubsetConstructorTest, Any) {
    auto param = GetParam();
    bool duplicate = std::get<0>(param);
    bool sorted = std::get<1>(param);
    auto sub = SubsetCoreUtils::spawn_indices<int>(true, 5, duplicate, sorted);

    tatami::DelayedSubset<double, int, decltype(sub)> manual(dense, sub, true);
    auto ref = reference_on_rows(dense.get(), sub);
    tatami_test::test_simple_row_access(manual, *ref);
    tatami_test::test_simple_column_access(manual, *ref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedSubset,
    SubsetConstructorTest,
    ::testing::Combine(
        ::testing::Values(false, true), // whether to support duplicate indices.
        ::testing::Values(true, false)  // whether to require sorted indices.
    )
);

/****************************************************
 ****************************************************/

TEST(DelayedSubset, ConstOverload) {
    int NR = 9, NC = 7;
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(NR, NC, std::vector<double>(NR * NC)));
    std::vector<int> subset{ 1, 3, 5 };

    auto sub = tatami::make_DelayedSubset<0>(dense, subset);
    EXPECT_EQ(sub->ncol(), NC);
    EXPECT_EQ(sub->nrow(), subset.size());
}

TEST(DelayedSubset, ArrayView) {
    int NR = 9, NC = 7;
    auto dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(NR, NC, std::vector<double>(NR * NC)));

    std::vector<int> subset{ 1, 3, 5 };
    tatami::ArrayView<int> aview(subset.data(), subset.size());

    auto sub = tatami::make_DelayedSubset<0>(dense, subset);
    EXPECT_EQ(sub->ncol(), NC);
    EXPECT_EQ(sub->nrow(), subset.size());
}
