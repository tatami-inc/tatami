#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <vector>
#include <memory>
#include <tuple>
#include <random>

#include "tatami/base/dense/DenseMatrix.hpp"
#include "tatami/base/subset/make_DelayedSubset.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../_tests/test_column_access.h"
#include "../_tests/test_row_access.h"
#include "../_tests/simulate_vector.h"

template<class PARAM> 
class SubsetTest : public ::testing::TestWithParam<PARAM> {
protected:
    size_t NR= 90, NC = 170;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse;

protected:
    void SetUp() {
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, simulate_sparse_vector<double>(NR * NC, 0.1)));
        sparse = tatami::convert_to_sparse<false>(dense.get()); // column-major.
        return;
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
                output.insert(output.end(), rng() % 3, output[i]);
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
            wrk->fetch_copy(r, ptr);
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

class SubsetFullAccessTest : public SubsetTest<std::tuple<size_t, bool, bool, bool, size_t> > {};

TEST_P(SubsetFullAccessTest, OnRow) {
    auto param = GetParam();
    auto sub = spawn_indices<size_t>(std::get<0>(param), NR, std::get<1>(param), std::get<2>(param));
    auto dense_subbed = tatami::make_DelayedSubset<0>(dense, sub);
    auto sparse_subbed = tatami::make_DelayedSubset<0>(sparse, sub);
    auto sparse_subbed2 = tatami::make_DelayedSubset<0>(sparse, std::vector<int>(sub.begin(), sub.end())); // check that eliding the copy works.
    auto ref = reference_on_rows(sub);

    EXPECT_EQ(sub.size(), dense_subbed->nrow());
    EXPECT_EQ(dense->ncol(), dense_subbed->ncol());
    EXPECT_EQ(dense->sparse(), dense_subbed->sparse());
    EXPECT_EQ(sparse->sparse(), sparse_subbed->sparse());
    EXPECT_TRUE(dense_subbed->prefer_rows());
    EXPECT_FALSE(sparse_subbed->prefer_rows());

    size_t FORWARD = std::get<3>(param);
    size_t JUMP = std::get<4>(param);
    test_simple_row_access(dense_subbed.get(), ref.get(), FORWARD, JUMP);
    test_simple_row_access(sparse_subbed.get(), ref.get(), FORWARD, JUMP);
    test_simple_row_access(sparse_subbed2.get(), ref.get(), FORWARD, JUMP);

    test_simple_column_access(dense_subbed.get(), ref.get(), FORWARD, JUMP);
    test_simple_column_access(sparse_subbed.get(), ref.get(), FORWARD, JUMP);
    test_simple_column_access(sparse_subbed2.get(), ref.get(), FORWARD, JUMP);
}

TEST_P(SubsetFullAccessTest, OnColumn) {
    auto param = GetParam();
    auto sub = spawn_indices<size_t>(std::get<0>(param), NC, std::get<1>(param), std::get<2>(param));
    auto dense_subbed = tatami::make_DelayedSubset<1>(dense, sub);
    auto sparse_subbed = tatami::make_DelayedSubset<1>(sparse, sub);
    auto sparse_subbed2 = tatami::make_DelayedSubset<1>(sparse, std::vector<int>(sub.begin(), sub.end())); // check that eliding the copy works.
    auto ref = reference_on_columns(sub);

    EXPECT_EQ(dense->nrow(), dense_subbed->nrow());
    EXPECT_EQ(sub.size(), dense_subbed->ncol());
    EXPECT_EQ(dense->sparse(), dense_subbed->sparse());
    EXPECT_EQ(sparse->sparse(), sparse_subbed->sparse());
    EXPECT_TRUE(dense_subbed->prefer_rows());
    EXPECT_FALSE(sparse_subbed->prefer_rows());

    size_t FORWARD = std::get<3>(param);
    size_t JUMP = std::get<4>(param);
    test_simple_row_access(dense_subbed.get(), ref.get(), FORWARD, JUMP);
    test_simple_row_access(sparse_subbed.get(), ref.get(), FORWARD, JUMP);
    test_simple_row_access(sparse_subbed2.get(), ref.get(), FORWARD, JUMP);

    test_simple_column_access(dense_subbed.get(), ref.get(), FORWARD, JUMP);
    test_simple_column_access(sparse_subbed.get(), ref.get(), FORWARD, JUMP);
    test_simple_column_access(sparse_subbed2.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubset,
    SubsetFullAccessTest,
    ::testing::Combine(
        ::testing::Values(1, 5, 10), // step size.
        ::testing::Values(false, true), // whether to support duplicate indices.
        ::testing::Values(true, false), // whether to require sorted indices.
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 3) // jump, to test the workspace memory.
    )
);

/****************************************************
 ****************************************************/

class SubsetSlicedAccessTest : public SubsetTest<std::tuple<size_t, bool, bool, size_t, std::vector<double> > > {};

TEST_P(SubsetSlicedAccessTest, OnRow) {
    auto param = GetParam();
    auto sub = spawn_indices<size_t>(std::get<0>(param), NR, std::get<1>(param), std::get<2>(param));

    auto dense_subbed = tatami::make_DelayedSubset<0>(dense, sub);
    auto sparse_subbed = tatami::make_DelayedSubset<0>(sparse, sub);
    auto sparse_subbed2 = tatami::make_DelayedSubset<0>(sparse, std::vector<int>(sub.begin(), sub.end())); // check that eliding the copy works.
    auto ref = reference_on_rows(sub);

    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t RFIRST = interval_info[0] * sub.size(), RLAST = interval_info[1] * sub.size();
    size_t CFIRST = interval_info[0] * NC, CLAST = interval_info[1] * NC;

    test_sliced_row_access(dense_subbed.get(), ref.get(), true, JUMP, CFIRST, CLAST);
    test_sliced_row_access(sparse_subbed.get(), ref.get(), true, JUMP, CFIRST, CLAST);
    test_sliced_row_access(sparse_subbed2.get(), ref.get(), true, JUMP, CFIRST, CLAST);

    test_sliced_column_access(dense_subbed.get(), ref.get(), true, JUMP, RFIRST, RLAST);
    test_sliced_column_access(sparse_subbed.get(), ref.get(), true, JUMP, RFIRST, RLAST);
    test_sliced_column_access(sparse_subbed2.get(), ref.get(), true, JUMP, RFIRST, RLAST);
}

TEST_P(SubsetSlicedAccessTest, OnColumn) {
    auto param = GetParam();
    auto sub = spawn_indices<size_t>(std::get<0>(param), NC, std::get<1>(param), std::get<2>(param));

    auto dense_subbed = tatami::make_DelayedSubset<1>(dense, sub);
    auto sparse_subbed = tatami::make_DelayedSubset<1>(sparse, sub);
    auto sparse_subbed2 = tatami::make_DelayedSubset<1>(sparse, std::vector<int>(sub.begin(), sub.end())); // check that eliding the copy works.
    auto ref = reference_on_columns(sub);

    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t RFIRST = interval_info[0] * NR, RLAST = interval_info[1] * NR;
    size_t CFIRST = interval_info[0] * sub.size(), CLAST = interval_info[1] * sub.size();

    test_sliced_row_access(dense_subbed.get(), ref.get(), true, JUMP, CFIRST, CLAST);
    test_sliced_row_access(sparse_subbed.get(), ref.get(), true, JUMP, CFIRST, CLAST);
    test_sliced_row_access(sparse_subbed2.get(), ref.get(), true, JUMP, CFIRST, CLAST);

    test_sliced_column_access(dense_subbed.get(), ref.get(), true, JUMP, RFIRST, RLAST);
    test_sliced_column_access(sparse_subbed.get(), ref.get(), true, JUMP, RFIRST, RLAST);
    test_sliced_column_access(sparse_subbed2.get(), ref.get(), true, JUMP, RFIRST, RLAST);
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubset,
    SubsetSlicedAccessTest,
    ::testing::Combine(
        ::testing::Values(1, 5, 10), // step size.
        ::testing::Values(false, true), // whether to support duplicate indices.
        ::testing::Values(true, false), // whether to require sorted indices.
        ::testing::Values(1, 3), // jump, to check the workspace memory
        ::testing::Values(
            std::vector<double>({ 0, 0.66 }), 
            std::vector<double>({ 0.33, 0.787 }),
            std::vector<double>({ 0.5, 1 })
        )
    )
);

/****************************************************
 ****************************************************/

class SubsetIndexedAccessTest : public SubsetTest<std::tuple<size_t, bool, bool, size_t, std::vector<double> > > {};

TEST_P(SubsetIndexedAccessTest, OnRow) {
    auto param = GetParam();
    auto sub = spawn_indices<size_t>(std::get<0>(param), NR, std::get<1>(param), std::get<2>(param));

    auto dense_subbed = tatami::make_DelayedSubset<0>(dense, sub);
    auto sparse_subbed = tatami::make_DelayedSubset<0>(sparse, sub);
    auto ref = reference_on_rows(sub);

    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t RFIRST = interval_info[0] * sub.size(), 
        CFIRST = interval_info[0] * NC, 
        STEP = interval_info[1];

    test_indexed_row_access(dense_subbed.get(), ref.get(), true, JUMP, CFIRST, STEP);
    test_indexed_row_access(sparse_subbed.get(), ref.get(), true, JUMP, CFIRST, STEP);

    test_indexed_column_access(dense_subbed.get(), ref.get(), true, JUMP, RFIRST, STEP);
    test_indexed_column_access(sparse_subbed.get(), ref.get(), true, JUMP, RFIRST, STEP);
}

TEST_P(SubsetIndexedAccessTest, OnColumn) {
    auto param = GetParam();
    auto sub = spawn_indices<size_t>(std::get<0>(param), NC, std::get<1>(param), std::get<2>(param));
    auto dense_subbed = tatami::make_DelayedSubset<1>(dense, sub);
    auto sparse_subbed = tatami::make_DelayedSubset<1>(sparse, sub);
    auto ref = reference_on_columns(sub);

    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t RFIRST = interval_info[0] * NR, 
        CFIRST = interval_info[0] * sub.size(), 
        STEP = interval_info[1];

    test_indexed_row_access(dense_subbed.get(), ref.get(), true, JUMP, CFIRST, STEP);
    test_indexed_row_access(sparse_subbed.get(), ref.get(), true, JUMP, CFIRST, STEP);

    test_indexed_column_access(dense_subbed.get(), ref.get(), true, JUMP, RFIRST, STEP);
    test_indexed_column_access(sparse_subbed.get(), ref.get(), true, JUMP, RFIRST, STEP);
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubset,
    SubsetIndexedAccessTest,
    ::testing::Combine(
        ::testing::Values(1, 5, 10), // step size.
        ::testing::Values(false, true), // whether to support duplicate indices.
        ::testing::Values(true, false), // whether to require sorted indices.
        ::testing::Values(1, 3), // jump, to check the workspace memory
        ::testing::Values(
            std::vector<double>({ 0, 5 }), 
            std::vector<double>({ 0.33, 3 }),
            std::vector<double>({ 0.5, 2 })
        )
    )
);

/****************************************************
 ****************************************************/

class SubsetConstructorTest : public SubsetTest<std::tuple<bool, bool> > {};

TEST_P(SubsetConstructorTest, SortedUnique) {
    auto param = GetParam();
    bool duplicate = std::get<0>(param);
    bool sorted = std::get<1>(param);
    auto sub = spawn_indices<size_t>(5, NR, duplicate, sorted);

    if (sorted && !duplicate) {
        tatami::DelayedSubsetSortedUnique<0, double, int, decltype(sub)> manual(dense, sub);
        auto ref = reference_on_rows(sub);
        test_simple_row_access(&manual, ref.get(), true, 1);
        test_simple_column_access(&manual, ref.get(), true, 1);
    } else {
        try {
            tatami::DelayedSubsetSortedUnique<0, double, int, decltype(sub)> manual(dense, sub);
            FAIL() << "expected exception during construction";
        } catch (std::exception& e) {
            EXPECT_THAT(e.what(), ::testing::HasSubstr("unique"));
        }
    }
}

TEST_P(SubsetConstructorTest, Sorted) {
    auto param = GetParam();
    bool duplicate = std::get<0>(param);
    bool sorted = std::get<1>(param);
    auto sub = spawn_indices<size_t>(5, NR, duplicate, sorted);

    if (sorted) {
        tatami::DelayedSubsetSorted<0, double, int, decltype(sub)> manual(dense, sub);
        auto ref = reference_on_rows(sub);
        test_simple_row_access(&manual, ref.get(), true, 1);
        test_simple_column_access(&manual, ref.get(), true, 1);
    } else {
        try {
            tatami::DelayedSubsetSortedUnique<0, double, int, decltype(sub)> manual(dense, sub); // '<' breaks EXPECT_ANY_THROW macro.
            FAIL() << "expected exception during construction";
        } catch (std::exception& e) {
            EXPECT_THAT(e.what(), ::testing::HasSubstr("sorted"));
        }
    }
}

TEST_P(SubsetConstructorTest, Unique) {
    auto param = GetParam();
    bool duplicate = std::get<0>(param);
    bool sorted = std::get<1>(param);
    auto sub = spawn_indices<size_t>(5, NR, duplicate, sorted);

    if (!duplicate) {
        tatami::DelayedSubsetUnique<0, double, int, decltype(sub)> manual(dense, sub);
        auto ref = reference_on_rows(sub);
        test_simple_row_access(&manual, ref.get(), true, 1);
        test_simple_column_access(&manual, ref.get(), true, 1);
    } else {
        try {
            tatami::DelayedSubsetSortedUnique<0, double, int, decltype(sub)> manual(dense, sub);
            FAIL() << "expected exception during construction";
        } catch (std::exception& e) {
            EXPECT_THAT(e.what(), ::testing::HasSubstr("sorted"));
        }
    }
}

TEST_P(SubsetConstructorTest, Any) {
    auto param = GetParam();
    bool duplicate = std::get<0>(param);
    bool sorted = std::get<1>(param);
    auto sub = spawn_indices<size_t>(5, NR, duplicate, sorted);

    tatami::DelayedSubset<0, double, int, decltype(sub)> manual(dense, sub);
    auto ref = reference_on_rows(sub);
    test_simple_row_access(&manual, ref.get(), true, 1);
    test_simple_column_access(&manual, ref.get(), true, 1);
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubset,
    SubsetConstructorTest,
    ::testing::Combine(
        ::testing::Values(false, true), // whether to support duplicate indices.
        ::testing::Values(true, false)  // whether to require sorted indices.
    )
);
