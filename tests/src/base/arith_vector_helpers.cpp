#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedIsometricOp.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../data/data.h"
#include "TestCore.h"

template<class PARAM>
class ArithVectorTest : public TestCore<::testing::TestWithParam<PARAM> > {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense, sparse;
protected:
    void SetUp() {
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
        sparse = tatami::convert_to_sparse(dense.get(), false);
        return;
    }
};

std::vector<double> create_vector(size_t n, double starter = 0, double jump = 1){ 
    std::vector<double> output(n, starter);
    for (size_t i = 1; i < n; ++i) {
        output[i] = output[i-1] + jump;
    }
    return output;
}

/****************************
 ********* ADDITION *********
 ****************************/

template<class PARAM>
class ArithVectorAdditionTest : public ArithVectorTest<PARAM> {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;
    std::vector<double> vec;

    void extra_assemble(const PARAM& param) {
        bool row = std::get<2>(param);
        vec = create_vector(
            row ? (this->dense)->nrow() : (this->dense)->ncol(),
            std::get<0>(param), 
            std::get<1>(param)
        );
        if (row) {
            auto op = tatami::make_DelayedAddVectorHelper<0>(vec);
            dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
        } else {
            auto op = tatami::make_DelayedAddVectorHelper<1>(vec);
            dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
        }
    }

    void add_value(double& ref, const PARAM& param, size_t r, size_t c) {
        ref += vec[std::get<2>(param) ? r : c];
        return;
    }
};

using ArithVectorAdditionFullTest = ArithVectorAdditionTest<std::tuple<double, double, bool, size_t> >;

TEST_P(ArithVectorAdditionFullTest, ColumnAccess) {
    auto param = GetParam();
    extra_assemble(param);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense->prefer_rows());
    EXPECT_FALSE(sparse->prefer_rows());

    auto work_dense = dense_mod->new_workspace(false);
    auto work_sparse = sparse_mod->new_workspace(false);

    for (size_t c = 0; c < dense->ncol(); c += std::get<3>(param)) {
        auto expected = extract_dense<false>(dense.get(), c);
        for (size_t r = 0; r < dense->nrow(); ++r) {
            add_value(expected[r], param, r, c);
        }

        auto output = extract_dense<false>(dense_mod.get(), c);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<false>(sparse_mod.get(), c);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<false>(sparse_mod.get(), c);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<false>(dense_mod.get(), c, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputSW = extract_sparse<false>(sparse_mod.get(), c, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

TEST_P(ArithVectorAdditionFullTest, RowAccess) {
    auto param = GetParam();
    extra_assemble(param);
    auto work_dense = dense_mod->new_workspace(true);
    auto work_sparse = sparse_mod->new_workspace(true);

    for (size_t r = 0; r < dense->nrow(); r += std::get<3>(param)) {
        auto expected = extract_dense<true>(dense.get(), r);
        for (size_t c = 0; c < dense->ncol(); ++c) {
            add_value(expected[c], param, r, c);
        }

        auto output = extract_dense<true>(dense_mod.get(), r);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<true>(sparse_mod.get(), r);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<true>(sparse_mod.get(), r);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<true>(dense_mod.get(), r, work_dense.get());
        EXPECT_EQ(output, expected);

        auto outputSW = extract_sparse<true>(sparse_mod.get(), r, work_sparse.get());
        EXPECT_EQ(outputDW, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorAdditionFullTest,
    ::testing::Combine(
        ::testing::Values(5),
        ::testing::Values(0.5),
        ::testing::Values(true, false), // add by row, or by column
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorAdditionSubsetTest = ArithVectorAdditionTest<std::tuple<double, double, bool, size_t, std::vector<size_t> > >;

TEST_P(ArithVectorAdditionSubsetTest, ColumnAccess) {
    auto param = GetParam();
    extra_assemble(param);
    auto work_dense = dense_mod->new_workspace(false);
    auto work_sparse = sparse_mod->new_workspace(false);

    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    for (size_t c = 0; c < dense->ncol(); c += JUMP, FIRST += SHIFT) {
        auto interval = wrap_intervals(FIRST, FIRST + LEN, dense->nrow());
        size_t start = interval.first, end = interval.second;

        auto expected = extract_dense<false>(dense.get(), c, start, end);
        for (size_t r = start; r < end; ++r) {
            add_value(expected[r - start], param, r, c);
        }

        auto output = extract_dense<false>(dense_mod.get(), c, start, end);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<false>(sparse_mod.get(), c, start, end);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<false>(sparse_mod.get(), c, start, end);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<false>(dense_mod.get(), c, start, end, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputSW = extract_sparse<false>(sparse_mod.get(), c, start, end, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

TEST_P(ArithVectorAdditionSubsetTest, RowAccess) {
    auto param = GetParam();
    extra_assemble(param);
    auto work_dense = dense_mod->new_workspace(true);
    auto work_sparse = sparse_mod->new_workspace(true);

    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    for (size_t r = 0; r < dense->nrow(); r += JUMP, FIRST += SHIFT) {
        auto interval = wrap_intervals(FIRST, FIRST + LEN, dense->ncol());
        size_t start = interval.first, end = interval.second;

        auto expected = extract_dense<true>(dense.get(), r, start, end);
        for (size_t c = start; c < end; ++c) {
            add_value(expected[c - start], param, r, c);
        }

        auto output = extract_dense<true>(dense_mod.get(), r, start, end);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<true>(sparse_mod.get(), r, start, end);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<true>(sparse_mod.get(), r, start, end);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<true>(dense_mod.get(), r, start, end, work_dense.get());
        EXPECT_EQ(output, expected);

        auto outputSW = extract_sparse<true>(sparse_mod.get(), r, start, end, work_sparse.get());
        EXPECT_EQ(outputDW, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorAdditionSubsetTest,
    ::testing::Combine(
        ::testing::Values(5),
        ::testing::Values(0.5),
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<size_t>({ 0, 8, 3 }), // overlapping shifts
            std::vector<size_t>({ 1, 4, 4 }), // non-overlapping shifts
            std::vector<size_t>({ 3, 10, 0 })
        )
    )
);

/*******************************
 ********* SUBTRACTION *********
 *******************************/

template<class PARAM>
class ArithVectorSubtractionTest : public ArithVectorTest<PARAM> {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;
    std::vector<double> vec;

    void extra_assemble(const PARAM& param) {
        bool row = std::get<2>(param);
        vec = create_vector(
            row ? (this->dense)->nrow() : (this->dense)->ncol(),
            std::get<0>(param), 
            std::get<1>(param)
        );

        bool right = std::get<3>(param);
        if (row) {
            if (right) {
                auto op = tatami::make_DelayedSubtractVectorHelper<true, 0>(vec);
                dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
            } else {
                auto op = tatami::make_DelayedSubtractVectorHelper<false, 0>(vec);
                dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
            }
        } else {
            if (right) {
                auto op = tatami::make_DelayedSubtractVectorHelper<true, 1>(vec);
                dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
            } else {
                auto op = tatami::make_DelayedSubtractVectorHelper<false, 1>(vec);
                dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
            }
        }
    }

    void subtract_value(double& ref, const PARAM& param, size_t r, size_t c) {
        const double& val = vec[std::get<2>(param) ? r : c];
        if (std::get<3>(param)) {
            ref -= val;
        } else {
            ref = val - ref;
        }
        return;
    }
};

using ArithVectorSubtractionFullTest = ArithVectorSubtractionTest<std::tuple<double, double, bool, bool, size_t> >;

TEST_P(ArithVectorSubtractionFullTest, ColumnAccess) {
    auto param = GetParam();
    extra_assemble(param);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense->prefer_rows());
    EXPECT_FALSE(sparse->prefer_rows());

    auto work_dense = dense_mod->new_workspace(false);
    auto work_sparse = sparse_mod->new_workspace(false);

    for (size_t c = 0; c < dense->ncol(); c += std::get<4>(param)) {
        auto expected = extract_dense<false>(dense.get(), c);
        for (size_t r = 0; r < dense->nrow(); ++r) {
            subtract_value(expected[r], param, r, c);
        }

        auto output = extract_dense<false>(dense_mod.get(), c);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<false>(sparse_mod.get(), c);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<false>(sparse_mod.get(), c);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<false>(dense_mod.get(), c, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputSW = extract_sparse<false>(sparse_mod.get(), c, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

TEST_P(ArithVectorSubtractionFullTest, RowAccess) {
    auto param = GetParam();
    extra_assemble(param);
    auto work_dense = dense_mod->new_workspace(true);
    auto work_sparse = sparse_mod->new_workspace(true);

    for (size_t r = 0; r < dense->nrow(); r += std::get<4>(param)) {
        auto expected = extract_dense<true>(dense.get(), r);
        for (size_t c = 0; c < dense->ncol(); ++c) {
            subtract_value(expected[c], param, r, c);
        }

        auto output = extract_dense<true>(dense_mod.get(), r);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<true>(sparse_mod.get(), r);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<true>(sparse_mod.get(), r);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<true>(dense_mod.get(), r, work_dense.get());
        EXPECT_EQ(output, expected);

        auto outputSW = extract_sparse<true>(sparse_mod.get(), r, work_sparse.get());
        EXPECT_EQ(outputDW, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorSubtractionFullTest,
    ::testing::Combine(
        ::testing::Values(1),
        ::testing::Values(0.33),
        ::testing::Values(true, false), // by row or by column
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorSubtractionSubsetTest = ArithVectorSubtractionTest<std::tuple<double, double, bool, bool, size_t, std::vector<size_t> > >;

TEST_P(ArithVectorSubtractionSubsetTest, ColumnAccess) {
    auto param = GetParam();
    extra_assemble(param);
    auto work_dense = dense_mod->new_workspace(false);
    auto work_sparse = sparse_mod->new_workspace(false);

    size_t JUMP = std::get<4>(param);
    auto interval_info = std::get<5>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    for (size_t c = 0; c < dense->ncol(); c += JUMP, FIRST += SHIFT) {
        auto interval = wrap_intervals(FIRST, FIRST + LEN, dense->nrow());
        size_t start = interval.first, end = interval.second;

        auto expected = extract_dense<false>(dense.get(), c, start, end);
        for (size_t r = start; r < end; ++r) {
            subtract_value(expected[r - start], param, r, c);
        }

        auto output = extract_dense<false>(dense_mod.get(), c, start, end);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<false>(sparse_mod.get(), c, start, end);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<false>(sparse_mod.get(), c, start, end);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<false>(dense_mod.get(), c, start, end, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputSW = extract_sparse<false>(sparse_mod.get(), c, start, end, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

TEST_P(ArithVectorSubtractionSubsetTest, RowAccess) {
    auto param = GetParam();
    extra_assemble(param);
    auto work_dense = dense_mod->new_workspace(true);
    auto work_sparse = sparse_mod->new_workspace(true);

    size_t JUMP = std::get<4>(param);
    auto interval_info = std::get<5>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    for (size_t r = 0; r < dense->nrow(); r += JUMP, FIRST += SHIFT) {
        auto interval = wrap_intervals(FIRST, FIRST + LEN, dense->ncol());
        size_t start = interval.first, end = interval.second;

        auto expected = extract_dense<true>(dense.get(), r, start, end);
        for (size_t c = start; c < end; ++c) {
            subtract_value(expected[c - start], param, r, c);
        }

        auto output = extract_dense<true>(dense_mod.get(), r, start, end);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<true>(sparse_mod.get(), r, start, end);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<true>(sparse_mod.get(), r, start, end);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<true>(dense_mod.get(), r, start, end, work_dense.get());
        EXPECT_EQ(output, expected);

        auto outputSW = extract_sparse<true>(sparse_mod.get(), r, start, end, work_sparse.get());
        EXPECT_EQ(outputDW, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorSubtractionSubsetTest,
    ::testing::Combine(
        ::testing::Values(5),
        ::testing::Values(0.5),
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<size_t>({ 0, 8, 3 }), // overlapping shifts
            std::vector<size_t>({ 1, 4, 4 }), // non-overlapping shifts
            std::vector<size_t>({ 3, 10, 0 })
        )
    )
);

/**********************************
 ********* MULTIPLICATION *********
 **********************************/

template<class PARAM>
class ArithVectorMultiplicationTest : public ArithVectorTest<PARAM> {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;
    std::vector<double> vec;

    void extra_assemble(const PARAM& param) {
        bool row = std::get<2>(param);
        vec = create_vector(
            row ? (this->dense)->nrow() : (this->dense)->ncol(),
            std::get<0>(param), 
            std::get<1>(param)
        );
        if (row) {
            auto op = tatami::make_DelayedMultiplyVectorHelper<0>(vec);
            dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
        } else {
            auto op = tatami::make_DelayedMultiplyVectorHelper<1>(vec);
            dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
        }
    }

    void multiply_value(double& ref, const PARAM& param, size_t r, size_t c) {
        ref *= vec[std::get<2>(param) ? r : c];
        return;
    }
};

using ArithVectorMultiplicationFullTest = ArithVectorMultiplicationTest<std::tuple<double, double, bool, size_t> >;

TEST_P(ArithVectorMultiplicationFullTest, ColumnAccess) {
    auto param = GetParam();
    extra_assemble(param);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense->prefer_rows());
    EXPECT_FALSE(sparse->prefer_rows());

    auto work_dense = dense_mod->new_workspace(false);
    auto work_sparse = sparse_mod->new_workspace(false);

    for (size_t c = 0; c < dense->ncol(); c += std::get<3>(param)) {
        auto expected = extract_dense<false>(dense.get(), c);
        for (size_t r = 0; r < dense->nrow(); ++r) {
            multiply_value(expected[r], param, r, c);
        }

        auto output = extract_dense<false>(dense_mod.get(), c);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<false>(sparse_mod.get(), c);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<false>(sparse_mod.get(), c);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<false>(dense_mod.get(), c, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputSW = extract_sparse<false>(sparse_mod.get(), c, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

TEST_P(ArithVectorMultiplicationFullTest, RowAccess) {
    auto param = GetParam();
    extra_assemble(param);
    auto work_dense = dense_mod->new_workspace(true);
    auto work_sparse = sparse_mod->new_workspace(true);

    for (size_t r = 0; r < dense->nrow(); r += std::get<3>(param)) {
        auto expected = extract_dense<true>(dense.get(), r);
        for (size_t c = 0; c < dense->ncol(); ++c) {
            multiply_value(expected[c], param, r, c);
        }

        auto output = extract_dense<true>(dense_mod.get(), r);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<true>(sparse_mod.get(), r);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<true>(sparse_mod.get(), r);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<true>(dense_mod.get(), r, work_dense.get());
        EXPECT_EQ(output, expected);

        auto outputSW = extract_sparse<true>(sparse_mod.get(), r, work_sparse.get());
        EXPECT_EQ(outputDW, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorMultiplicationFullTest,
    ::testing::Combine(
        ::testing::Values(5),
        ::testing::Values(0.5),
        ::testing::Values(true, false),
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorMultiplicationSubsetTest = ArithVectorMultiplicationTest<std::tuple<double, double, bool, size_t, std::vector<size_t> > >;

TEST_P(ArithVectorMultiplicationSubsetTest, ColumnAccess) {
    auto param = GetParam();
    extra_assemble(param);
    auto work_dense = dense_mod->new_workspace(false);
    auto work_sparse = sparse_mod->new_workspace(false);

    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    for (size_t c = 0; c < dense->ncol(); c += JUMP, FIRST += SHIFT) {
        auto interval = wrap_intervals(FIRST, FIRST + LEN, dense->nrow());
        size_t start = interval.first, end = interval.second;

        auto expected = extract_dense<false>(dense.get(), c, start, end);
        for (size_t r = start; r < end; ++r) {
            multiply_value(expected[r - start], param, r, c);
        }

        auto output = extract_dense<false>(dense_mod.get(), c, start, end);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<false>(sparse_mod.get(), c, start, end);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<false>(sparse_mod.get(), c, start, end);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<false>(dense_mod.get(), c, start, end, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputSW = extract_sparse<false>(sparse_mod.get(), c, start, end, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

TEST_P(ArithVectorMultiplicationSubsetTest, RowAccess) {
    auto param = GetParam();
    extra_assemble(param);
    auto work_dense = dense_mod->new_workspace(true);
    auto work_sparse = sparse_mod->new_workspace(true);

    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    for (size_t r = 0; r < dense->nrow(); r += JUMP, FIRST += SHIFT) {
        auto interval = wrap_intervals(FIRST, FIRST + LEN, dense->ncol());
        size_t start = interval.first, end = interval.second;

        auto expected = extract_dense<true>(dense.get(), r, start, end);
        for (size_t c = start; c < end; ++c) {
            multiply_value(expected[c - start], param, r, c);
        }

        auto output = extract_dense<true>(dense_mod.get(), r, start, end);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<true>(sparse_mod.get(), r, start, end);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<true>(sparse_mod.get(), r, start, end);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<true>(dense_mod.get(), r, start, end, work_dense.get());
        EXPECT_EQ(output, expected);

        auto outputSW = extract_sparse<true>(sparse_mod.get(), r, start, end, work_sparse.get());
        EXPECT_EQ(outputDW, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorMultiplicationSubsetTest,
    ::testing::Combine(
        ::testing::Values(5),
        ::testing::Values(0.5),
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<size_t>({ 0, 8, 3 }), // overlapping shifts
            std::vector<size_t>({ 1, 4, 4 }), // non-overlapping shifts
            std::vector<size_t>({ 3, 10, 0 })
        )
    )
);

/****************************
 ********* DIVISION *********
 ****************************/

template<class PARAM>
class ArithVectorDivisionTest : public ArithVectorTest<PARAM> {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;
    std::vector<double> vec;

    void extra_assemble(const PARAM& param) {
        bool row = std::get<2>(param);
        vec = create_vector(
            row ? (this->dense)->nrow() : (this->dense)->ncol(),
            std::get<0>(param), 
            std::get<1>(param)
        );

        bool right = std::get<3>(param);
        if (!right) {
            // Some shenanigans to avoid division by zero.
            tatami::DelayedExpHelper<double> op2a;
            this->dense = tatami::make_DelayedIsometricOp(this->dense, op2a);
            this->sparse = tatami::make_DelayedIsometricOp(this->sparse, op2a);
        }

        if (row) {
            if (right) {
                auto op = tatami::make_DelayedDivideVectorHelper<true, 0>(vec);
                dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
            } else {
                auto op = tatami::make_DelayedDivideVectorHelper<false, 0>(vec);
                dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
            }
        } else {
            if (right) {
                auto op = tatami::make_DelayedDivideVectorHelper<true, 1>(vec);
                dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
            } else {
                auto op = tatami::make_DelayedDivideVectorHelper<false, 1>(vec);
                dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
            }
        }
    }

    void divide_value(double& ref, const PARAM& param, size_t r, size_t c) {
        const double& val = vec[std::get<2>(param) ? r : c];
        if (std::get<3>(param)) {
            ref /= val;
        } else {
            ref = val / ref;
        }
        return;
    }
};

using ArithVectorDivisionFullTest = ArithVectorDivisionTest<std::tuple<double, double, bool, bool, size_t> >;

TEST_P(ArithVectorDivisionFullTest, ColumnAccess) {
    auto param = GetParam();
    extra_assemble(param);

    EXPECT_FALSE(dense_mod->sparse());
    if (std::get<3>(param)) {
        // only true if we're dividing on the right; 
        // divisions on the left need to undo the sparsity
        // to avoid divide by zero errors.
        EXPECT_TRUE(sparse_mod->sparse());
    }
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense->prefer_rows());
    EXPECT_FALSE(sparse->prefer_rows());

    auto work_dense = dense_mod->new_workspace(false);
    auto work_sparse = sparse_mod->new_workspace(false);

    for (size_t c = 0; c < dense->ncol(); c += std::get<4>(param)) {
        auto expected = extract_dense<false>(dense.get(), c);
        for (size_t r = 0; r < dense->nrow(); ++r) {
            divide_value(expected[r], param, r, c);
        }

        auto output = extract_dense<false>(dense_mod.get(), c);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<false>(sparse_mod.get(), c);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<false>(sparse_mod.get(), c);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<false>(dense_mod.get(), c, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputSW = extract_sparse<false>(sparse_mod.get(), c, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

TEST_P(ArithVectorDivisionFullTest, RowAccess) {
    auto param = GetParam();
    extra_assemble(param);
    auto work_dense = dense_mod->new_workspace(true);
    auto work_sparse = sparse_mod->new_workspace(true);

    for (size_t r = 0; r < dense->nrow(); r += std::get<4>(param)) {
        auto expected = extract_dense<true>(dense.get(), r);
        for (size_t c = 0; c < dense->ncol(); ++c) {
            divide_value(expected[c], param, r, c);
        }

        auto output = extract_dense<true>(dense_mod.get(), r);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<true>(sparse_mod.get(), r);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<true>(sparse_mod.get(), r);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<true>(dense_mod.get(), r, work_dense.get());
        EXPECT_EQ(output, expected);

        auto outputSW = extract_sparse<true>(sparse_mod.get(), r, work_sparse.get());
        EXPECT_EQ(outputDW, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorDivisionFullTest,
    ::testing::Combine(
        ::testing::Values(1),
        ::testing::Values(0.33),
        ::testing::Values(true, false), // by row or by column
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorDivisionSubsetTest = ArithVectorDivisionTest<std::tuple<double, double, bool, bool, size_t, std::vector<size_t> > >;

TEST_P(ArithVectorDivisionSubsetTest, ColumnAccess) {
    auto param = GetParam();
    extra_assemble(param);
    auto work_dense = dense_mod->new_workspace(false);
    auto work_sparse = sparse_mod->new_workspace(false);

    size_t JUMP = std::get<4>(param);
    auto interval_info = std::get<5>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    for (size_t c = 0; c < dense->ncol(); c += JUMP, FIRST += SHIFT) {
        auto interval = wrap_intervals(FIRST, FIRST + LEN, dense->nrow());
        size_t start = interval.first, end = interval.second;

        auto expected = extract_dense<false>(dense.get(), c, start, end);
        for (size_t r = start; r < end; ++r) {
            divide_value(expected[r - start], param, r, c);
        }

        auto output = extract_dense<false>(dense_mod.get(), c, start, end);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<false>(sparse_mod.get(), c, start, end);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<false>(sparse_mod.get(), c, start, end);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<false>(dense_mod.get(), c, start, end, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputSW = extract_sparse<false>(sparse_mod.get(), c, start, end, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

TEST_P(ArithVectorDivisionSubsetTest, RowAccess) {
    auto param = GetParam();
    extra_assemble(param);
    auto work_dense = dense_mod->new_workspace(true);
    auto work_sparse = sparse_mod->new_workspace(true);

    size_t JUMP = std::get<4>(param);
    auto interval_info = std::get<5>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    for (size_t r = 0; r < dense->nrow(); r += JUMP, FIRST += SHIFT) {
        auto interval = wrap_intervals(FIRST, FIRST + LEN, dense->ncol());
        size_t start = interval.first, end = interval.second;

        auto expected = extract_dense<true>(dense.get(), r, start, end);
        for (size_t c = start; c < end; ++c) {
            divide_value(expected[c - start], param, r, c);
        }

        auto output = extract_dense<true>(dense_mod.get(), r, start, end);
        EXPECT_EQ(output, expected);

        // Works with different backends.
        auto output2 = extract_dense<true>(sparse_mod.get(), r, start, end);
        EXPECT_EQ(output2, expected);

        // Works in sparse mode as well.
        auto outputS = extract_sparse<true>(sparse_mod.get(), r, start, end);
        EXPECT_EQ(outputS, expected);

        // Passes along the workspace.
        auto outputDW = extract_dense<true>(dense_mod.get(), r, start, end, work_dense.get());
        EXPECT_EQ(output, expected);

        auto outputSW = extract_sparse<true>(sparse_mod.get(), r, start, end, work_sparse.get());
        EXPECT_EQ(outputDW, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorDivisionSubsetTest,
    ::testing::Combine(
        ::testing::Values(5),
        ::testing::Values(0.5),
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<size_t>({ 0, 8, 3 }), // overlapping shifts
            std::vector<size_t>({ 1, 4, 4 }), // non-overlapping shifts
            std::vector<size_t>({ 3, 10, 0 })
        )
    )
);


