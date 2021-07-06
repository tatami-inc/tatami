#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedSubset.hpp"
#include "tatami/base/DelayedSubsetBlock.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../data/data.h"
#include "TestCore.h"

template<class PARAM> 
class SubsetBlockTest : public TestCore<::testing::TestWithParam<PARAM> > {
protected:
    std::shared_ptr<tatami::numeric_matrix> dense, sparse, ref, dense_block, sparse_block;

    size_t first, last;
    std::vector<size_t> sub;
protected:
    void SetUp() {
        dense = std::shared_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
        sparse = tatami::convert_to_sparse(dense.get(), false); // column-major.
        return;
    }

    void extra_assemble(const PARAM& param) {
        size_t full;
        if (std::get<0>(param)) {
            full = dense->nrow(); 
        } else {
            full = dense->ncol();
        }

        first = static_cast<double>(full) * std::get<1>(param).first;
        last = static_cast<double>(full) * std::get<1>(param).second;
        sub.resize(last - first);
        for (size_t i = 0; i < sub.size(); ++i) { sub[i] = i + first; }

        if (std::get<0>(param)) {
            ref = tatami::make_DelayedSubset<0>(dense, sub);
            dense_block = tatami::make_DelayedSubsetBlock<0>(dense, first, last);
            sparse_block = tatami::make_DelayedSubsetBlock<0>(sparse, first, last);
        } else {
            ref = tatami::make_DelayedSubset<1>(dense, sub);
            dense_block = tatami::make_DelayedSubsetBlock<1>(dense, first, last);
            sparse_block = tatami::make_DelayedSubsetBlock<1>(sparse, first, last);
        }
    }
};

/*****************************
 *****************************/

using SubsetBlockFullAccessTest = SubsetBlockTest<std::tuple<bool, std::pair<double, double>, size_t> >;

TEST_P(SubsetBlockFullAccessTest, RowAccess) {
    auto param = GetParam();
    extra_assemble(param);

    EXPECT_EQ(ref->nrow(), dense_block->nrow());
    EXPECT_EQ(ref->ncol(), dense_block->ncol());
    if (std::get<0>(param)) {
        EXPECT_EQ(sub.size(), dense_block->nrow());
        EXPECT_EQ(dense->ncol(), dense_block->ncol());
    } else {
        EXPECT_EQ(dense->nrow(), dense_block->nrow());
        EXPECT_EQ(sub.size(), dense_block->ncol());
    }

    EXPECT_TRUE(dense_block->prefer_rows());
    EXPECT_FALSE(sparse_block->prefer_rows());
    EXPECT_FALSE(dense_block->sparse());
    EXPECT_TRUE(sparse_block->sparse());

    auto work_dense = dense_block->new_workspace(true);
    auto work_sparse = dense_block->new_workspace(true);

    for (size_t i = 0; i < dense_block->nrow(); i += std::get<2>(param)) {
        auto expected = extract_dense<true>(ref.get(), i);
        auto output = extract_dense<true>(dense_block.get(), i);
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        auto output2 = extract_dense<true>(sparse_block.get(), i);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<true>(sparse_block.get(), i);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<true>(dense_block.get(), i, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputSW = extract_sparse<true>(sparse_block.get(), i, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

TEST_P(SubsetBlockFullAccessTest, ColumnAccess) {
    auto param = GetParam();
    extra_assemble(param);

    auto work_dense = dense_block->new_workspace(false);
    auto work_sparse = dense_block->new_workspace(false);

    for (size_t i = 0; i < dense_block->ncol(); i += std::get<2>(param)) {
        auto expected = extract_dense<false>(ref.get(), i);
        auto output = extract_dense<false>(dense_block.get(), i);
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        auto output2 = extract_dense<false>(sparse_block.get(), i);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<false>(sparse_block.get(), i);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<false>(dense_block.get(), i, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputSW = extract_sparse<false>(sparse_block.get(), i, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubsetBlock,
    SubsetBlockFullAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row or column subsetting, respectively.
        ::testing::Values(
            std::make_pair(0.0, 0.5),
            std::make_pair(0.25, 0.8),
            std::make_pair(0.4, 1)
        ),
        ::testing::Values(1, 3) // jump, to check the workspace memory
    )
);

/*****************************
 *****************************/

using SubsetBlockSlicedAccessTest = SubsetBlockTest<std::tuple<bool, std::pair<double, double>, size_t, std::vector<size_t> > >;

TEST_P(SubsetBlockSlicedAccessTest, RowAccess) {
    auto param = GetParam();
    extra_assemble(param);

    auto slice = std::get<3>(GetParam());
    size_t FIRST = slice[0], LEN = slice[1], SHIFT = slice[2];

    auto work_dense = dense_block->new_workspace(true);
    auto work_sparse = dense_block->new_workspace(true);

    for (size_t i = 0; i < dense_block->nrow(); i += std::get<2>(param), FIRST += SHIFT) {
        auto interval = wrap_intervals(FIRST, FIRST + SHIFT, dense_block->ncol());
        size_t first = interval.first, last = interval.second;

        auto expected = extract_dense<true>(ref.get(), i, first, last);
        auto output = extract_dense<true>(dense_block.get(), i, first, last);
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        auto output2 = extract_dense<true>(sparse_block.get(), i, first, last);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<true>(sparse_block.get(), i, first, last);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<true>(dense_block.get(), i, first, last, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputSW = extract_sparse<true>(sparse_block.get(), i, first, last, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

TEST_P(SubsetBlockSlicedAccessTest, ColumnAccess) {
    auto param = GetParam();
    extra_assemble(param);

    auto slice = std::get<3>(GetParam());
    size_t FIRST = slice[0], LEN = slice[1], SHIFT = slice[2];

    auto work_dense = dense_block->new_workspace(false);
    auto work_sparse = dense_block->new_workspace(false);

    for (size_t i = 0; i < dense_block->ncol(); i += std::get<2>(param), FIRST += SHIFT) {
        auto interval = wrap_intervals(FIRST, FIRST + SHIFT, dense_block->nrow());
        size_t first = interval.first, last = interval.second;

        auto expected = extract_dense<false>(ref.get(), i, first, last);
        auto output = extract_dense<false>(dense_block.get(), i, first, last);
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        auto output2 = extract_dense<false>(sparse_block.get(), i, first, last);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<false>(sparse_block.get(), i, first, last);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<false>(dense_block.get(), i, first, last, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputSW = extract_sparse<false>(sparse_block.get(), i, first, last, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubsetBlock,
    SubsetBlockSlicedAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row or column subsetting, respectively.
        ::testing::Values(
            std::make_pair(0.0, 0.5),
            std::make_pair(0.25, 0.8),
            std::make_pair(0.4, 1)
        ),
        ::testing::Values(1, 3), // jump, to check the workspace memory
        ::testing::Values(
            std::vector<size_t>({ 0, 10, 1 }), // start, length, shift
            std::vector<size_t>({ 3, 5, 0 })
        )        
    )
);
