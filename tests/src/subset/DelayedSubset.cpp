#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <vector>
#include <memory>
#include <tuple>
#include <random>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/subset/make_DelayedSubset.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

class SubsetTestCore {
protected:
    size_t NR= 90, NC = 170;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse;

    std::vector<size_t> sub;
    std::shared_ptr<tatami::NumericMatrix> dense_subbed, sparse_subbed, ref;

protected:
    void assemble(bool byrow, size_t step_size, bool duplicates, bool sorted) {
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1)));
        sparse = tatami::convert_to_compressed_sparse<false>(dense.get()); // column-major.

        if (byrow) {
            sub = spawn_indices<size_t>(step_size, NR, duplicates, sorted);
            dense_subbed = tatami::make_DelayedSubset<0>(dense, sub);
            sparse_subbed = tatami::make_DelayedSubset<0>(sparse, sub);
            ref = reference_on_rows(sub);
        } else {
            sub = spawn_indices<size_t>(step_size, NC, duplicates, sorted);
            dense_subbed = tatami::make_DelayedSubset<1>(dense, sub);
            sparse_subbed = tatami::make_DelayedSubset<1>(sparse, sub);
            ref = reference_on_columns(sub);
        }
    }

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
    std::shared_ptr<tatami::NumericMatrix> reference_on_rows(const V& sub) const {
        std::vector<double> reference(sub.size() * NC);
        auto ptr = reference.data();
        auto wrk = dense->dense_row();

        for (auto r : sub) {
            auto src = wrk->fetch(r, ptr);
            tatami::copy_n(src, NC, ptr);
            ptr += NC;
        }

        return std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sub.size(), NC, std::move(reference)));
    }

    template<typename V>
    std::shared_ptr<tatami::NumericMatrix> reference_on_columns(const V& sub) const {
        std::vector<double> reference(sub.size() * NR);
        auto ptr = reference.data();
        std::vector<double> buffer(NC);
        auto wrk = dense->dense_row();

        for (size_t r = 0; r < NR; ++r) {
            auto full = wrk->fetch(r, buffer.data());
            for (auto s : sub) {
                *ptr = full[s];
                ++ptr;
            }
        }

        return std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, sub.size(), std::move(reference)));
    }
};

/****************************************************
 ****************************************************/

class SubsetFullAccessTest : 
    public ::testing::TestWithParam<std::tuple<bool, size_t, bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t> >,
    public SubsetTestCore {};

TEST_P(SubsetFullAccessTest, Full) {
    auto tparam = GetParam();
    bool byrow = std::get<0>(tparam);
    assemble(byrow, std::get<1>(tparam), std::get<2>(tparam), std::get<3>(tparam));

    if (byrow) {
        EXPECT_EQ(sub.size(), dense_subbed->nrow());
        EXPECT_EQ(dense->ncol(), dense_subbed->ncol());
    } else {
        EXPECT_EQ(dense->nrow(), dense_subbed->nrow());
        EXPECT_EQ(sub.size(), dense_subbed->ncol());
    }

    EXPECT_EQ(dense->sparse(), dense_subbed->sparse());
    EXPECT_EQ(dense->sparse_proportion(), dense_subbed->sparse_proportion());
    EXPECT_EQ(sparse->sparse(), sparse_subbed->sparse());
    EXPECT_EQ(sparse->sparse_proportion(), sparse_subbed->sparse_proportion());

    EXPECT_TRUE(dense_subbed->prefer_rows());
    EXPECT_EQ(dense->prefer_rows_proportion(), dense_subbed->prefer_rows_proportion());
    EXPECT_FALSE(sparse_subbed->prefer_rows());
    EXPECT_EQ(sparse->prefer_rows_proportion(), sparse_subbed->prefer_rows_proportion());

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<4>(tparam);
    params.use_oracle = std::get<5>(tparam);
    params.order = std::get<6>(tparam);
    params.jump = std::get<7>(tparam);

    tatami_test::test_full_access(params, dense_subbed.get(), ref.get());
    tatami_test::test_full_access(params, sparse_subbed.get(), ref.get());
}

INSTANTIATE_TEST_SUITE_P(
    DelayedSubset,
    SubsetFullAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // whether to subset on the rows or columns.
        ::testing::Values(2, 5, 10), // step size.
        ::testing::Values(false, true), // whether to support duplicate indices.
        ::testing::Values(true, false), // whether to require sorted indices.
        ::testing::Values(true, false), // row or column access.
        ::testing::Values(true, false), // whether to use an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3) // jump, to test the workspace memory.
    )
);

/****************************************************
 ****************************************************/

class SubsetSlicedAccessTest : 
    public ::testing::TestWithParam<std::tuple<bool, size_t, bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > >,
    public SubsetTestCore {};

TEST_P(SubsetSlicedAccessTest, Sliced) {
    auto tparam = GetParam();
    bool byrow = std::get<0>(tparam);
    assemble(byrow, std::get<1>(tparam), std::get<2>(tparam), std::get<3>(tparam));

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<4>(tparam);
    params.use_oracle = std::get<5>(tparam);
    params.order = std::get<6>(tparam);
    params.jump = std::get<7>(tparam);

    auto interval_info = std::get<8>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;

    tatami_test::test_block_access(params, dense_subbed.get(), ref.get(), FIRST, LAST);
    tatami_test::test_block_access(params, sparse_subbed.get(), ref.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedSubset,
    SubsetSlicedAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // whether to subset on the rows or columns.
        ::testing::Values(2, 5, 10), // step size.
        ::testing::Values(false, true), // whether to support duplicate indices.
        ::testing::Values(true, false), // whether to require sorted indices.
        ::testing::Values(true, false), // row or column access.
        ::testing::Values(true, false), // whether to use an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace memory.
        ::testing::Values(
            std::make_pair(0.0, 0.66), 
            std::make_pair(0.33, 0.787),
            std::make_pair(0.5, 1.0)
        )
    )
);

/****************************************************
 ****************************************************/

class SubsetIndexedAccessTest : 
    public ::testing::TestWithParam<std::tuple<bool, size_t, bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, int> > >,
    public SubsetTestCore {};

TEST_P(SubsetIndexedAccessTest, Indexed) {
    auto tparam = GetParam();
    bool byrow = std::get<0>(tparam);
    assemble(byrow, std::get<1>(tparam), std::get<2>(tparam), std::get<3>(tparam));

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<4>(tparam);
    params.use_oracle = std::get<5>(tparam);
    params.order = std::get<6>(tparam);
    params.jump = std::get<7>(tparam);

    auto interval_info = std::get<8>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, STEP = interval_info.second;

    tatami_test::test_indexed_access(params, dense_subbed.get(), ref.get(), FIRST, STEP);
    tatami_test::test_indexed_access(params, sparse_subbed.get(), ref.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedSubset,
    SubsetIndexedAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // subset by row.
        ::testing::Values(2, 5, 10), // step size.
        ::testing::Values(false, true), // whether to support duplicate indices.
        ::testing::Values(true, false), // whether to require sorted indices.
        ::testing::Values(true, false), // row or column access.
        ::testing::Values(true, false), // whether to use an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace memory.
        ::testing::Values(
            std::pair<double, int>(0, 5), 
            std::pair<double, int>(0.33, 3),
            std::pair<double, int>(0.5, 2)
        )
    )
);

///****************************************************
// ****************************************************/
//
//class SubsetConstructorTest : public SubsetTest<std::tuple<bool, bool> > {};
//
//TEST_P(SubsetConstructorTest, SortedUnique) {
//    auto param = GetParam();
//    bool duplicate = std::get<0>(param);
//    bool sorted = std::get<1>(param);
//    auto sub = spawn_indices<size_t>(5, NR, duplicate, sorted);
//
//    if (sorted && !duplicate) {
//        tatami::DelayedSubsetSortedUnique<0, double, int, decltype(sub)> manual(dense, sub);
//        auto ref = reference_on_rows(sub);
//        tatami_test::test_simple_row_access(&manual, ref.get(), true, 1);
//        tatami_test::test_simple_column_access(&manual, ref.get(), true, 1);
//    } else {
//        try {
//            tatami::DelayedSubsetSortedUnique<0, double, int, decltype(sub)> manual(dense, sub);
//            FAIL() << "expected exception during construction";
//        } catch (std::exception& e) {
//            EXPECT_THAT(e.what(), ::testing::HasSubstr("unique"));
//        }
//    }
//}
//
//TEST_P(SubsetConstructorTest, Sorted) {
//    auto param = GetParam();
//    bool duplicate = std::get<0>(param);
//    bool sorted = std::get<1>(param);
//    auto sub = spawn_indices<size_t>(5, NR, duplicate, sorted);
//
//    if (sorted) {
//        tatami::DelayedSubsetSorted<0, double, int, decltype(sub)> manual(dense, sub);
//        auto ref = reference_on_rows(sub);
//        tatami_test::test_simple_row_access(&manual, ref.get(), true, 1);
//        tatami_test::test_simple_column_access(&manual, ref.get(), true, 1);
//    } else {
//        try {
//            tatami::DelayedSubsetSortedUnique<0, double, int, decltype(sub)> manual(dense, sub); // '<' breaks EXPECT_ANY_THROW macro.
//            FAIL() << "expected exception during construction";
//        } catch (std::exception& e) {
//            EXPECT_THAT(e.what(), ::testing::HasSubstr("sorted"));
//        }
//    }
//}
//
//TEST_P(SubsetConstructorTest, Unique) {
//    auto param = GetParam();
//    bool duplicate = std::get<0>(param);
//    bool sorted = std::get<1>(param);
//    auto sub = spawn_indices<size_t>(5, NR, duplicate, sorted);
//
//    if (!duplicate) {
//        tatami::DelayedSubsetUnique<0, double, int, decltype(sub)> manual(dense, sub);
//        auto ref = reference_on_rows(sub);
//        tatami_test::test_simple_row_access(&manual, ref.get(), true, 1);
//        tatami_test::test_simple_column_access(&manual, ref.get(), true, 1);
//    } else {
//        try {
//            tatami::DelayedSubsetSortedUnique<0, double, int, decltype(sub)> manual(dense, sub);
//            FAIL() << "expected exception during construction";
//        } catch (std::exception& e) {
//            EXPECT_THAT(e.what(), ::testing::HasSubstr("sorted"));
//        }
//    }
//}
//
//TEST_P(SubsetConstructorTest, Any) {
//    auto param = GetParam();
//    bool duplicate = std::get<0>(param);
//    bool sorted = std::get<1>(param);
//    auto sub = spawn_indices<size_t>(5, NR, duplicate, sorted);
//
//    tatami::DelayedSubset<0, double, int, decltype(sub)> manual(dense, sub);
//    auto ref = reference_on_rows(sub);
//    tatami_test::test_simple_row_access(&manual, ref.get(), true, 1);
//    tatami_test::test_simple_column_access(&manual, ref.get(), true, 1);
//}
//
//INSTANTIATE_TEST_SUITE_P(
//    DelayedSubset,
//    SubsetConstructorTest,
//    ::testing::Combine(
//        ::testing::Values(false, true), // whether to support duplicate indices.
//        ::testing::Values(true, false)  // whether to require sorted indices.
//    )
//);

///****************************************************
// ****************************************************/
//
//TEST(DelayedSubset, ConstOverload) {
//    int NR = 9, NC = 7;
//    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1)));
//    std::vector<int> subset{ 1, 3, 5 };
//
//    auto sub = tatami::make_DelayedSubset<0>(dense, subset);
//    EXPECT_EQ(sub->ncol(), NC);
//    EXPECT_EQ(sub->nrow(), subset.size());
//}
//
//TEST(DelayedSubset, ArrayView) {
//    int NR = 9, NC = 7;
//    auto dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1)));
//
//    std::vector<int> subset{ 1, 3, 5 };
//    tatami::ArrayView<int> aview(subset.data(), subset.size());
//
//    auto sub = tatami::make_DelayedSubset<0>(dense, subset);
//    EXPECT_EQ(sub->ncol(), NC);
//    EXPECT_EQ(sub->nrow(), subset.size());
//}
