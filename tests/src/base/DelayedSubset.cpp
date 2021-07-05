#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedSubset.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../data/data.h"
#include "TestCore.h"

template<class PARAM> 
class SubsetTest : public TestCore0<::testing::TestWithParam<PARAM> > {
protected:
    std::shared_ptr<tatami::numeric_matrix> dense;
    std::shared_ptr<tatami::typed_matrix<double, int> > sparse;
    std::shared_ptr<tatami::workspace> work_dense, work_sparse;

protected:
    void assemble(size_t nr, size_t nc, const std::vector<double>& source) {
        dense = std::shared_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(nr, nc, source));
        sparse = tatami::convert_to_sparse(dense.get(), false); // column-major.
        return;
    }

    void SetUp() {
        assemble(sparse_nrow, sparse_ncol, sparse_matrix);
    }
};

/****************************************************
 ****************************************************/

class SubsetRowFullColumnAccessTest : public SubsetTest<std::vector<size_t> > {};

TEST_P(SubsetRowFullColumnAccessTest, Access) {
    std::vector<size_t> sub = GetParam();
    std::vector<double> buffer_full(dense->nrow());

    auto dense_subbed = tatami::make_DelayedSubset<0>(dense, sub);
    auto sparse_subbed = tatami::make_DelayedSubset<0>(sparse, sub);

    set_sizes(0, sub.size());
    EXPECT_EQ(sub.size(), dense_subbed->nrow());
    EXPECT_EQ(dense->ncol(), dense_subbed->ncol());

    EXPECT_TRUE(dense_subbed->prefer_rows());
    EXPECT_FALSE(sparse_subbed->prefer_rows());

    work_dense = dense_subbed->new_workspace(false);
    work_sparse = sparse_subbed->new_workspace(false);

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_output();
        fill_output(sparse_subbed->column(i, output.data()));

        // Reference extraction, the simple way.
        wipe_expected();
        auto X = sparse->column(i, buffer_full.data());
        auto eIt = expected.begin();
        for (auto s : sub) {
            (*eIt) = X[s];
            ++eIt;
        }
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed->column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed->sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed->sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_subbed->column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubset,
    SubsetRowFullColumnAccessTest,
    ::testing::Values(
        std::vector<size_t>({ 0, 3, 3, 13, 5, 2, 19, 4, 6, 11, 19, 8 }), // with duplicates
        std::vector<size_t>({ 1, 2, 3, 5, 9, 13, 17 }), // ordered, no duplicates
        std::vector<size_t>({ 8, 9, 10, 11}) // consecutive
    )
);

/****************************************************
 ****************************************************/

class SubsetRowSlicedColumnAccessTest : public SubsetTest<std::tuple<std::vector<size_t>, std::vector<size_t> > > {};

TEST_P(SubsetRowSlicedColumnAccessTest, Access) {
    std::vector<size_t> sub = std::get<0>(GetParam()); 
    std::vector<double> buffer_full(dense->nrow());

    auto dense_subbed = tatami::make_DelayedSubset<0>(dense, sub);
    auto sparse_subbed = tatami::make_DelayedSubset<0>(sparse, sub);

    auto slice = std::get<1>(GetParam());
    size_t FIRST = slice[0], LEN = slice[1], SHIFT = slice[2];

    work_dense = dense_subbed->new_workspace(false);
    work_sparse = sparse_subbed->new_workspace(false);

    for (size_t i = 0; i < dense->ncol(); ++i) {
        set_sizes(FIRST, std::min(FIRST + LEN, sub.size()));

        wipe_output();
        fill_output(sparse_subbed->column(i, output.data(), first, last));

        // Reference extraction, the simple way.
        wipe_expected();
        auto X = sparse->column(i, buffer_full.data());
        auto eIt = expected.begin();
        for (auto s = first; s < last; ++s) {
            (*eIt) = X[sub[s]];
            ++eIt;
        }
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed->column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed->sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);        

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed->sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        FIRST += SHIFT;
        FIRST %= sub.size();
    }
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubset,
    SubsetRowSlicedColumnAccessTest,
    ::testing::Combine(
        ::testing::Values(
            std::vector<size_t>({ 17, 18, 11, 18, 15, 17, 13, 18, 11, 9, 6, 3, 6, 18, 1 }), // with duplicates
            std::vector<size_t>({ 2, 3, 5, 7, 9, 12, 13 }), // ordered, no duplicates
            std::vector<size_t>({ 4, 5, 6, 7, 8, 9, 10 }) // consecutive
        ),
        ::testing::Values(
            std::vector<size_t>({ 0, 6, 13 }), // start, length, shift
            std::vector<size_t>({ 1, 7, 3 }),
            std::vector<size_t>({ 3, 18, 0 })
        )
    )
);

/****************************************************
 ****************************************************/

class SubsetRowFullRowAccessTest : public SubsetTest<std::vector<size_t> > {};

TEST_P(SubsetRowFullRowAccessTest, Access) {
    std::vector<size_t> sub = GetParam();
    std::vector<double> buffer_full(dense->nrow());

    auto dense_subbed = tatami::make_DelayedSubset<0>(dense, sub);
    auto sparse_subbed = tatami::make_DelayedSubset<0>(sparse, sub);

    set_sizes(0, dense->ncol());

    work_dense = dense_subbed->new_workspace(true);
    work_sparse = sparse_subbed->new_workspace(true);

    for (size_t i = 0; i < sub.size(); ++i) {
        wipe_output();
        fill_output(sparse_subbed->row(i, output.data()));

        wipe_expected();
        fill_expected(sparse->row(sub[i], expected.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed->row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed->sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed->sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_subbed->row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubset,
    SubsetRowFullRowAccessTest,
    ::testing::Values(
        std::vector<size_t>({ 13, 4, 17, 0, 17, 1, 19, 6, 1 }),
        std::vector<size_t>({ 0, 5, 10, 12, 16, 18, 19 }),
        std::vector<size_t>({ 15, 16, 17, 18, 19 }) // consecutive
    )
);

/****************************************************
 ****************************************************/

class SubsetColumnFullRowAccessTest : public SubsetTest<std::vector<size_t> > {};

TEST_P(SubsetColumnFullRowAccessTest, Access) {
    std::vector<size_t> sub = GetParam();
    std::vector<double> buffer_full(dense->ncol());

    auto dense_subbed = tatami::make_DelayedSubset<1>(dense, sub);
    auto sparse_subbed = tatami::make_DelayedSubset<1>(sparse, sub);

    set_sizes(0, sub.size());
    EXPECT_EQ(sub.size(), dense_subbed->ncol());
    EXPECT_EQ(dense->nrow(), dense_subbed->nrow());

    work_dense = dense_subbed->new_workspace(true);
    work_sparse = sparse_subbed->new_workspace(true);

    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_output();
        fill_output(sparse_subbed->row(i, output.data()));

        // Reference extraction, the simple way.
        wipe_expected();
        auto X = sparse->row(i, buffer_full.data());
        auto eIt = expected.begin();
        for (auto s : sub) {
            (*eIt) = X[s];
            ++eIt;
        }
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed->row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed->sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed->sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_subbed->row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubset,
    SubsetColumnFullRowAccessTest,
    ::testing::Values(
        std::vector<size_t>({ 3, 9, 1, 0, 9, 5, 8, 3, 1, 8, 7 }),
        std::vector<size_t>({ 0, 1, 2, 3, 5, 8 }),
        std::vector<size_t>({ 2, 3, 4, 5 })
    )
);

/****************************************************
 ****************************************************/

class SubsetColumnSlicedRowAccessTest : public SubsetTest<std::tuple<std::vector<size_t>, std::vector<size_t> > > {};

TEST_P(SubsetColumnSlicedRowAccessTest, Access) {
    std::vector<size_t> sub = std::get<0>(GetParam());
    std::vector<double> buffer_full(dense->ncol());

    auto dense_subbed = tatami::make_DelayedSubset<1>(dense, sub);
    auto sparse_subbed = tatami::make_DelayedSubset<1>(sparse, sub);

    work_dense = dense_subbed->new_workspace(true);
    work_sparse = sparse_subbed->new_workspace(true);

    auto slice = std::get<1>(GetParam());
    size_t FIRST = slice[0], LEN = slice[1], SHIFT = slice[2];

    for (size_t i = 0; i < dense->nrow(); ++i) {
        set_sizes(FIRST, std::min(FIRST + LEN, sub.size()));

        wipe_output();
        fill_output(sparse_subbed->row(i, output.data(), first, last));

        // Reference extraction, the simple way.
        wipe_expected();
        auto X = sparse->row(i, buffer_full.data());
        auto eIt = expected.begin();
        for (auto s = first; s < last && s < sub.size(); ++s) {
            (*eIt) = X[sub[s]];
            ++eIt;
        }
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed->row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed->sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);        

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed->sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        first += SHIFT;
        first %= sub.size();
    }
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubset,
    SubsetColumnSlicedRowAccessTest,
    ::testing::Combine(
        ::testing::Values(
            std::vector<size_t>({ 2, 2, 4, 8, 0, 7, 3, 1, 1, 2, 7, 8, 9, 9, 4, 5, 8, 5, 6, 2, 0 }),
            std::vector<size_t>({ 2, 3, 5, 7, 9 }), // ordered, no duplicates
            std::vector<size_t>({ 3, 4, 5, 6, 7, 8, 9 }) // consecutive
        ),
        ::testing::Values(
            std::vector<size_t>({ 0, 6, 1 }), // start, length, shift
            std::vector<size_t>({ 5, 5, 2 }),
            std::vector<size_t>({ 3, 7, 0 })
        )
    )
);

/****************************************************
 ****************************************************/

class SubsetColumnFullColumnAccessTest : public SubsetTest<std::vector<size_t> > {};

TEST_P(SubsetColumnFullColumnAccessTest, Access) {
    auto sub = GetParam();
    std::vector<double> buffer_full(dense->ncol());

    auto dense_subbed = tatami::make_DelayedSubset<1>(dense, sub);
    auto sparse_subbed = tatami::make_DelayedSubset<1>(sparse, sub);

    set_sizes(0, sparse->nrow());

    work_dense = dense_subbed->new_workspace(false);
    work_sparse = sparse_subbed->new_workspace(false);

    for (size_t i = 0; i < sub.size(); ++i) {
        wipe_output();
        fill_output(sparse_subbed->column(i, output.data()));

        wipe_expected();
        fill_expected(sparse->column(sub[i], expected.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed->column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed->sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed->sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_subbed->column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubset,
    SubsetColumnFullColumnAccessTest,
    ::testing::Values(
        std::vector<size_t>({ 7, 8, 0, 5, 1, 4, 1 }),
        std::vector<size_t>({ 1, 3, 5, 7 }),
        std::vector<size_t>({ 3, 4, 5 })
    )
);
