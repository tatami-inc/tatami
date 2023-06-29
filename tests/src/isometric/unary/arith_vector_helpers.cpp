#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOp.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class ArithVectorUtils {
protected:
    size_t nrow = 291, ncol = 188;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    std::vector<double> simulated;
protected:
    void assemble() {
        simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1);
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated));
        sparse = tatami::convert_to_sparse<false>(dense.get()); // column major.
        return;
    }

    static std::vector<double> create_vector(size_t n, double starter = 0, double jump = 1) { 
        std::vector<double> output(n, starter);
        for (size_t i = 1; i < n; ++i) {
            output[i] = output[i-1] + jump;
        }
        return output;
    }
};

/****************************
 ********* ADDITION *********
 ****************************/

template<class PARAM>
class ArithVectorAdditionTest : public ::testing::TestWithParam<PARAM>, public ArithVectorUtils {
protected:
    void SetUp() {
        assemble();
    }

    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    std::vector<double> vec;

    void extra_assemble(bool row) {
        vec = this->create_vector(row ? this->nrow : this->ncol, 5, 0.5);

        if (row) {
            auto op = tatami::make_DelayedAddVectorHelper<0>(vec);
            dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
        } else {
            auto op = tatami::make_DelayedAddVectorHelper<1>(vec);
            dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
        }

        auto refvec = this->simulated;
        for (size_t r = 0; r < this->nrow; ++r) {
            for (size_t c = 0; c < this->ncol; ++c) {
                refvec[r * this->ncol + c] += vec[row ? r : c];
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));
    }
};

using ArithVectorAdditionFullTest = ArithVectorAdditionTest<std::tuple<bool, bool, size_t> >;

TEST_P(ArithVectorAdditionFullTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));
    bool FORWARD = std::get<1>(param);
    int JUMP = std::get<2>(param);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense_mod->sparse_proportion(), 0);
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(sparse_mod->sparse_proportion(), 0);
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_EQ(dense_mod->prefer_rows_proportion(), 1);
    EXPECT_FALSE(sparse_mod->prefer_rows());
    EXPECT_EQ(sparse_mod->prefer_rows_proportion(), 0);

    tatami_test::test_simple_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);

    tatami_test::test_simple_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorAdditionFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or by column
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorAdditionBlockTest = ArithVectorAdditionTest<std::tuple<bool, bool, size_t, std::vector<double> > >;

TEST_P(ArithVectorAdditionBlockTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));

    bool FORWARD = std::get<1>(param);
    size_t JUMP = std::get<2>(param);
    auto interval_info = std::get<3>(param);

    {
        size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;
        tatami_test::test_sliced_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }

    {
        size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;
        tatami_test::test_sliced_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorAdditionBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.35 }),
            std::vector<double>({ 0.38, 0.61 }), 
            std::vector<double>({ 0.777, 1 })

        )
    )
);

using ArithVectorAdditionIndexTest = ArithVectorAdditionTest<std::tuple<bool, bool, size_t, std::vector<double> > >;

TEST_P(ArithVectorAdditionIndexTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));

    bool FORWARD = std::get<1>(param);
    size_t JUMP = std::get<2>(param);
    auto interval_info = std::get<3>(param);

    {
        size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;
        tatami_test::test_indexed_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }

    {
        size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;
        tatami_test::test_indexed_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorAdditionIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.077 }),
            std::vector<double>({ 0.21, 0.09 }), 
            std::vector<double>({ 0.56, 0.01 })
        )
    )
);

/*******************************
 ********* SUBTRACTION *********
 *******************************/

template<class PARAM>
class ArithVectorSubtractionTest : public ::testing::TestWithParam<PARAM>, public ArithVectorUtils {
protected:
    void SetUp() {
        assemble();
    }

    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    std::vector<double> vec;

    void extra_assemble(bool row, bool right) {
        vec = this->create_vector(row ? this->nrow : this->ncol, -10, 2.222);

        if (row) {
            if (right) {
                auto op = tatami::make_DelayedSubtractVectorHelper<true, 0>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            } else {
                auto op = tatami::make_DelayedSubtractVectorHelper<false, 0>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            }
        } else {
            if (right) {
                auto op = tatami::make_DelayedSubtractVectorHelper<true, 1>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            } else {
                auto op = tatami::make_DelayedSubtractVectorHelper<false, 1>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            }
        }

        auto refvec = this->simulated;
        for (size_t r = 0; r < this->nrow; ++r) {
            for (size_t c = 0; c < this->ncol; ++c) {
                auto& x = refvec[r * this->ncol + c];
                if (right) {
                    x -= vec[row ? r : c];
                } else {
                    x = vec[row ? r : c] - x;
                }
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));
    }
};

using ArithVectorSubtractionFullTest = ArithVectorSubtractionTest<std::tuple<bool, bool, bool, size_t> >;

TEST_P(ArithVectorSubtractionFullTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));
    bool FORWARD = std::get<2>(param);
    int JUMP = std::get<3>(param);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());

    tatami_test::test_simple_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);

    tatami_test::test_simple_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorSubtractionFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // by row or by column
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorSubtractionBlockTest = ArithVectorSubtractionTest<std::tuple<bool, bool, bool, size_t, std::vector<double> > >;

TEST_P(ArithVectorSubtractionBlockTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);

    {
        size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;
        tatami_test::test_sliced_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }

    {
        size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;
        tatami_test::test_sliced_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorSubtractionBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.541 }), 
            std::vector<double>({ 0.111, 0.999 }), 
            std::vector<double>({ 0.42, 1 })
        )
    )
);

using ArithVectorSubtractionIndexTest = ArithVectorSubtractionTest<std::tuple<bool, bool, bool, size_t, std::vector<double> > >;

TEST_P(ArithVectorSubtractionIndexTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);

    {
        size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;
        tatami_test::test_indexed_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }

    {
        size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;
        tatami_test::test_indexed_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorSubtractionIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.077 }),
            std::vector<double>({ 0.21, 0.09 }), 
            std::vector<double>({ 0.56, 0.01 })
        )
    )
);

/**********************************
 ********* MULTIPLICATION *********
 **********************************/

template<class PARAM>
class ArithVectorMultiplicationTest : public ::testing::TestWithParam<PARAM>, public ArithVectorUtils {
protected:
    void SetUp() {
        assemble();
    }

    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    std::vector<double> vec;

    void extra_assemble(bool row) {
        vec = this->create_vector(row ? this->nrow : this->ncol, 99, -2.5);

        if (row) {
            auto op = tatami::make_DelayedMultiplyVectorHelper<0>(vec);
            dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
        } else {
            auto op = tatami::make_DelayedMultiplyVectorHelper<1>(vec);
            dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
        }

        auto refvec = this->simulated;
        for (size_t r = 0; r < this->nrow; ++r) {
            for (size_t c = 0; c < this->ncol; ++c) {
                refvec[r * this->ncol + c] *= vec[row ? r : c];
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));
    }
};

using ArithVectorMultiplicationFullTest = ArithVectorMultiplicationTest<std::tuple<bool, bool, size_t> >;

TEST_P(ArithVectorMultiplicationFullTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));
    bool FORWARD = std::get<1>(param);
    int JUMP = std::get<2>(param);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense_mod->sparse_proportion(), 0);
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(sparse_mod->sparse_proportion(), 1);
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());

    tatami_test::test_simple_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);

    tatami_test::test_simple_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorMultiplicationFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // by row or by column
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorMultiplicationBlockTest = ArithVectorMultiplicationTest<std::tuple<bool, bool, size_t, std::vector<double> > >;

TEST_P(ArithVectorMultiplicationBlockTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));

    bool FORWARD = std::get<1>(param);
    size_t JUMP = std::get<2>(param);
    auto interval_info = std::get<3>(param);

    {
        size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;
        tatami_test::test_sliced_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }

    {
        size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;
        tatami_test::test_sliced_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorMultiplicationBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.35 }),
            std::vector<double>({ 0.38, 0.61 }), 
            std::vector<double>({ 0.777, 1 })

        )
    )
);

using ArithVectorMultiplicationIndexTest = ArithVectorMultiplicationTest<std::tuple<bool, bool, size_t, std::vector<double> > >;

TEST_P(ArithVectorMultiplicationIndexTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));

    bool FORWARD = std::get<1>(param);
    size_t JUMP = std::get<2>(param);
    auto interval_info = std::get<3>(param);

    {
        size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;
        tatami_test::test_indexed_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }

    {
        size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;
        tatami_test::test_indexed_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorMultiplicationIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.077 }),
            std::vector<double>({ 0.21, 0.09 }), 
            std::vector<double>({ 0.56, 0.01 })
        )
    )
);

/****************************
 ********* DIVISION *********
 ****************************/

template<class PARAM>
class ArithVectorDivisionTest : public ::testing::TestWithParam<PARAM>, public ArithVectorUtils {
protected:
    void SetUp() {
        assemble();
    }

    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    std::vector<double> vec;

    void extra_assemble(bool row, bool right) {
        vec = this->create_vector(row ? this->nrow : this->ncol, -19, 2.11);

        if (row) {
            if (right) {
                auto op = tatami::make_DelayedDivideVectorHelper<true, 0>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            } else {
                auto op = tatami::make_DelayedDivideVectorHelper<false, 0>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            }
        } else {
            if (right) {
                auto op = tatami::make_DelayedDivideVectorHelper<true, 1>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            } else {
                auto op = tatami::make_DelayedDivideVectorHelper<false, 1>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            }
        }

        auto refvec = this->simulated;
        for (size_t r = 0; r < this->nrow; ++r) {
            for (size_t c = 0; c < this->ncol; ++c) {
                auto& x = refvec[r * this->ncol + c];
                auto val = vec[row ? r : c];
                if (right) {
                    x /= val;
                } else {
                    if (x) {
                        x = val / x;
                    } else if (val > 0) {
                        x = std::numeric_limits<double>::infinity();
                    } else {
                        x = -std::numeric_limits<double>::infinity();
                    }
                }
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));
    }
};

using ArithVectorDivisionFullTest = ArithVectorDivisionTest<std::tuple<bool, bool, bool, size_t> >;

TEST_P(ArithVectorDivisionFullTest, Basic) {
    auto param = GetParam();
    auto RIGHT = std::get<1>(param);
    extra_assemble(std::get<0>(param), RIGHT);
    bool FORWARD = std::get<2>(param);
    int JUMP = std::get<3>(param);

    EXPECT_FALSE(dense_mod->sparse());
    if (RIGHT) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());

    tatami_test::test_simple_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);

    tatami_test::test_simple_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorDivisionFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // by row or by column
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorDivisionBlockTest = ArithVectorDivisionTest<std::tuple<bool, bool, bool, size_t, std::vector<double> > >;

TEST_P(ArithVectorDivisionBlockTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);

    {
        size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;
        tatami_test::test_sliced_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }

    {
        size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;
        tatami_test::test_sliced_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorDivisionBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.541 }), 
            std::vector<double>({ 0.111, 0.999 }), 
            std::vector<double>({ 0.42, 1 })
        )
    )
);

using ArithVectorDivisionIndexTest = ArithVectorDivisionTest<std::tuple<bool, bool, bool, size_t, std::vector<double> > >;

TEST_P(ArithVectorDivisionIndexTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);

    {
        size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;
        tatami_test::test_indexed_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }

    {
        size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;
        tatami_test::test_indexed_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorDivisionIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.077 }),
            std::vector<double>({ 0.21, 0.09 }), 
            std::vector<double>({ 0.56, 0.01 })
        )
    )
);

/****************************
 ********** POWER ***********
 ****************************/

template<class PARAM>
class ArithVectorPowerTest : public ::testing::TestWithParam<PARAM>, public ArithVectorUtils {
protected:
    void SetUp() {
        assemble();
    }

    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    std::vector<double> vec;

    void extra_assemble(bool row, bool right) {
        vec = this->create_vector(row ? this->nrow : this->ncol, 0.01, 2.11);

        tatami::DelayedAbsHelper op0;
        dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op0);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op0);
        if (row) {
            if (right) {
                auto op = tatami::make_DelayedPowerVectorHelper<true, 0>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense_mod, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse_mod, op);
            } else {
                auto op = tatami::make_DelayedPowerVectorHelper<false, 0>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense_mod, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse_mod, op);
            }
        } else {
            if (right) {
                auto op = tatami::make_DelayedPowerVectorHelper<true, 1>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense_mod, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse_mod, op);
            } else {
                auto op = tatami::make_DelayedPowerVectorHelper<false, 1>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense_mod, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse_mod, op);
            }
        }

        auto refvec = this->simulated;
        for (size_t r = 0; r < this->nrow; ++r) {
            for (size_t c = 0; c < this->ncol; ++c) {
                auto& x = refvec[r * this->ncol + c];
                auto val = vec[row ? r : c];
                if (right) {
                    x = std::pow(std::abs(x), val);
                } else {
                    x = std::pow(val, std::abs(x));
                }
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));
    }
};

using ArithVectorPowerFullTest = ArithVectorPowerTest<std::tuple<bool, bool, bool, size_t> >;

TEST_P(ArithVectorPowerFullTest, Basic) {
    auto param = GetParam();
    auto RIGHT = std::get<1>(param);
    extra_assemble(std::get<0>(param), RIGHT);
    bool FORWARD = std::get<2>(param);
    int JUMP = std::get<3>(param);

    EXPECT_FALSE(dense_mod->sparse());
    if (RIGHT) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());

    tatami_test::test_simple_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);

    tatami_test::test_simple_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorPowerFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // by row or by column
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorPowerBlockTest = ArithVectorPowerTest<std::tuple<bool, bool, bool, size_t, std::vector<double> > >;

TEST_P(ArithVectorPowerBlockTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);

    {
        size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;
        tatami_test::test_sliced_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }

    {
        size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;
        tatami_test::test_sliced_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorPowerBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.541 }),
            std::vector<double>({ 0.111, 0.999 }),
            std::vector<double>({ 0.42, 1 })
        )
    )
);

using ArithVectorPowerIndexTest = ArithVectorPowerTest<std::tuple<bool, bool, bool, size_t, std::vector<double> > >;

TEST_P(ArithVectorPowerIndexTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);

    {
        size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;
        tatami_test::test_indexed_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }

    {
        size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;
        tatami_test::test_indexed_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorPowerIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.077 }),
            std::vector<double>({ 0.21, 0.09 }),
            std::vector<double>({ 0.56, 0.01 })
        )
    )
);

/****************************
 ********** MODULO **********
 ****************************/

template<class PARAM>
class ArithVectorModuloTest : public ::testing::TestWithParam<PARAM>, public ArithVectorUtils {
protected:
    void SetUp() {
        assemble();
    }

    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    std::vector<double> vec;

    void extra_assemble(bool row, bool right) {
        vec = this->create_vector(row ? this->nrow : this->ncol, -19, 2.11);

        if (row) {
            if (right) {
                auto op = tatami::make_DelayedModuloVectorHelper<true, 0>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            } else {
                auto op = tatami::make_DelayedModuloVectorHelper<false, 0>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            }
        } else {
            if (right) {
                auto op = tatami::make_DelayedModuloVectorHelper<true, 1>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            } else {
                auto op = tatami::make_DelayedModuloVectorHelper<false, 1>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            }
        }

        auto refvec = this->simulated;
        for (size_t r = 0; r < this->nrow; ++r) {
            for (size_t c = 0; c < this->ncol; ++c) {
                auto& x = refvec[r * this->ncol + c];
                auto val = vec[row ? r : c];
                if (right) {
                    x = std::fmod(x, val);
                } else {
                    x = std::fmod(val, x);
                }
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));
    }
};

using ArithVectorModuloFullTest = ArithVectorModuloTest<std::tuple<bool, bool, bool, size_t> >;

TEST_P(ArithVectorModuloFullTest, Basic) {
    auto param = GetParam();
    auto RIGHT = std::get<1>(param);
    extra_assemble(std::get<0>(param), RIGHT);
    bool FORWARD = std::get<2>(param);
    int JUMP = std::get<3>(param);

    EXPECT_FALSE(dense_mod->sparse());
    if (RIGHT) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());

    tatami_test::test_simple_column_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP);

    tatami_test::test_simple_row_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorModuloFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // by row or by column
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorModuloBlockTest = ArithVectorModuloTest<std::tuple<bool, bool, bool, size_t, std::vector<double> > >;

TEST_P(ArithVectorModuloBlockTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);

    {
        size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;
        tatami_test::test_sliced_column_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }

    {
        size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;
        tatami_test::test_sliced_row_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorModuloBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.541 }),
            std::vector<double>({ 0.111, 0.999 }),
            std::vector<double>({ 0.42, 1 })
        )
    )
);

using ArithVectorModuloIndexTest = ArithVectorModuloTest<std::tuple<bool, bool, bool, size_t, std::vector<double> > >;

TEST_P(ArithVectorModuloIndexTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);

    {
        size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;
        tatami_test::test_indexed_column_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }

    {
        size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;
        tatami_test::test_indexed_row_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorModuloIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.077 }),
            std::vector<double>({ 0.21, 0.09 }),
            std::vector<double>({ 0.56, 0.01 })
        )
    )
);

/****************************
 ***** INTEGER DIVISION *****
 ****************************/

template<class PARAM>
class ArithVectorIntegerDivisionTest : public ::testing::TestWithParam<PARAM>, public ArithVectorUtils {
protected:
    void SetUp() {
        assemble();
    }

    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    std::vector<double> vec;

    void extra_assemble(bool row, bool right) {
        vec = this->create_vector(row ? this->nrow : this->ncol, -19, 2.11);

        if (row) {
            if (right) {
                auto op = tatami::make_DelayedIntegerDivideVectorHelper<true, 0>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            } else {
                auto op = tatami::make_DelayedIntegerDivideVectorHelper<false, 0>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            }
        } else {
            if (right) {
                auto op = tatami::make_DelayedIntegerDivideVectorHelper<true, 1>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            } else {
                auto op = tatami::make_DelayedIntegerDivideVectorHelper<false, 1>(vec);
                dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
                sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            }
        }

        auto refvec = this->simulated;
        for (size_t r = 0; r < this->nrow; ++r) {
            for (size_t c = 0; c < this->ncol; ++c) {
                auto& x = refvec[r * this->ncol + c];
                auto val = vec[row ? r : c];
                // x == (x %% y) + y * (x %/% y)
                if (right) {
                    x = std::floor(careful_division(x, val));
                } else {
                    x = std::floor(careful_division(val, x));
                }
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));
    }
};

using ArithVectorIntegerDivisionFullTest = ArithVectorIntegerDivisionTest<std::tuple<bool, bool, bool, size_t> >;

TEST_P(ArithVectorIntegerDivisionFullTest, Basic) {
    auto param = GetParam();
    auto RIGHT = std::get<1>(param);
    extra_assemble(std::get<0>(param), RIGHT);
    bool FORWARD = std::get<2>(param);
    int JUMP = std::get<3>(param);

    EXPECT_FALSE(dense_mod->sparse());
    if (RIGHT) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());

    tatami_test::test_simple_column_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP);

    tatami_test::test_simple_row_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorIntegerDivisionFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // by row or by column
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorIntegerDivisionBlockTest = ArithVectorIntegerDivisionTest<std::tuple<bool, bool, bool, size_t, std::vector<double> > >;

TEST_P(ArithVectorIntegerDivisionBlockTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);

    {
        size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;
        tatami_test::test_sliced_column_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_column_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }

    {
        size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;
        tatami_test::test_sliced_row_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
        tatami_test::test_sliced_row_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorIntegerDivisionBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.541 }),
            std::vector<double>({ 0.111, 0.999 }),
            std::vector<double>({ 0.42, 1 })
        )
    )
);

using ArithVectorIntegerDivisionIndexTest = ArithVectorIntegerDivisionTest<std::tuple<bool, bool, bool, size_t, std::vector<double> > >;

TEST_P(ArithVectorIntegerDivisionIndexTest, Basic) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);

    {
        size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;
        tatami_test::test_indexed_column_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_column_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }

    {
        size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;
        tatami_test::test_indexed_row_access<true>(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
        tatami_test::test_indexed_row_access<true>(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorIntegerDivisionIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.077 }),
            std::vector<double>({ 0.21, 0.09 }),
            std::vector<double>({ 0.56, 0.01 })
        )
    )
);

/**************************
 ********* ORACLE *********
 **************************/

class ArithVectorOracleTest : public ::testing::TestWithParam<std::tuple<bool, bool> >, public ArithVectorUtils {
protected:
    void SetUp() {
        assemble();
    }

    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, wrapped_dense_mod, wrapped_sparse_mod;
    std::vector<double> vec;

    void extra_assemble(bool row) {
        vec = this->create_vector(row ? this->nrow : this->ncol, 5, 0.5);

        if (row) {
            auto op = tatami::make_DelayedAddVectorHelper<0>(vec);
            dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            wrapped_dense_mod = tatami::make_DelayedUnaryIsometricOp(tatami_test::make_CrankyMatrix(this->dense), op);
            wrapped_sparse_mod = tatami::make_DelayedUnaryIsometricOp(tatami_test::make_CrankyMatrix(this->sparse), op);
        } else {
            auto op = tatami::make_DelayedAddVectorHelper<1>(vec);
            dense_mod = tatami::make_DelayedUnaryIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedUnaryIsometricOp(this->sparse, op);
            wrapped_dense_mod = tatami::make_DelayedUnaryIsometricOp(tatami_test::make_CrankyMatrix(this->dense), op);
            wrapped_sparse_mod = tatami::make_DelayedUnaryIsometricOp(tatami_test::make_CrankyMatrix(this->sparse), op);
        }
    }
};

TEST_P(ArithVectorOracleTest, Validate) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));
    EXPECT_FALSE(dense_mod->uses_oracle(true));
    EXPECT_TRUE(wrapped_dense_mod->uses_oracle(true));

    tatami_test::test_oracle_column_access(wrapped_dense_mod.get(), dense_mod.get(), std::get<1>(param));
    tatami_test::test_oracle_column_access(wrapped_sparse_mod.get(), sparse_mod.get(), std::get<1>(param));

    tatami_test::test_oracle_row_access(wrapped_dense_mod.get(), dense_mod.get(), std::get<1>(param));
    tatami_test::test_oracle_row_access(wrapped_sparse_mod.get(), sparse_mod.get(), std::get<1>(param));
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorOracleTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false)  // use random or consecutive oracle.
    )
);

/***********************************
 ********* CONST OVERLOADS *********
 ***********************************/

TEST(ArithVector, ConstOverload) {
    int nrow = 23, ncol = 42;
    auto simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1);
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated));

    auto vec = std::vector<double>(nrow);
    auto op = tatami::make_DelayedAddVectorHelper<0>(vec);
    auto mat = tatami::make_DelayedUnaryIsometricOp(dense, std::move(op));

    // cursory checks.
    EXPECT_EQ(mat->nrow(), dense->nrow());
    EXPECT_EQ(mat->ncol(), dense->ncol());
}

/**************************
 ********* ZEROED *********
 **************************/

class ArithVectorZeroedTest : public ::testing::TestWithParam<bool>, public ArithVectorUtils {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;

    void SetUp() {
        assemble();
    }
};

TEST_P(ArithVectorZeroedTest, Addition) {
    std::vector<double> zeroed(GetParam() ? nrow : ncol);

    if (GetParam()) {
        auto op = tatami::make_DelayedAddVectorHelper<0>(std::move(zeroed));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    } else {
        auto op = tatami::make_DelayedAddVectorHelper<1>(std::move(zeroed));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    }

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());

    tatami_test::test_simple_column_access(dense_mod.get(), dense.get(), true, 1);
    tatami_test::test_simple_column_access(sparse_mod.get(), sparse.get(), true, 1); 

    tatami_test::test_simple_row_access(dense_mod.get(), dense.get(), true, 1);
    tatami_test::test_simple_row_access(sparse_mod.get(), sparse.get(), true, 1);
}

TEST_P(ArithVectorZeroedTest, Subtraction) {
    std::vector<double> zeroed(GetParam() ? nrow : ncol);

    {
        if (GetParam()) {
            auto op = tatami::make_DelayedSubtractVectorHelper<true, 0>(zeroed);
            dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
            sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
        } else {
            auto op = tatami::make_DelayedSubtractVectorHelper<true, 1>(zeroed);
            dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
            sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
        }

        EXPECT_FALSE(dense_mod->sparse());
        EXPECT_TRUE(sparse_mod->sparse());

        tatami_test::test_simple_column_access(dense_mod.get(), dense.get(), true, 1);
        tatami_test::test_simple_column_access(sparse_mod.get(), sparse.get(), true, 1); 

        tatami_test::test_simple_row_access(dense_mod.get(), dense.get(), true, 1);
        tatami_test::test_simple_row_access(sparse_mod.get(), sparse.get(), true, 1);
    }

    {
        if (GetParam()) {
            auto op = tatami::make_DelayedSubtractVectorHelper<false, 0>(zeroed);
            dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
            sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
        } else {
            auto op = tatami::make_DelayedSubtractVectorHelper<false, 1>(zeroed);
            dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
            sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
        }

        auto copy = simulated;
        for (auto& x : copy) {
            x *= -1;
        }
        tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

        EXPECT_FALSE(dense_mod->sparse());
        EXPECT_TRUE(sparse_mod->sparse());

        tatami_test::test_simple_column_access(dense_mod.get(), &ref, true, 1);
        tatami_test::test_simple_column_access(sparse_mod.get(), &ref, true, 1); 

        tatami_test::test_simple_row_access(dense_mod.get(), &ref, true, 1);
        tatami_test::test_simple_row_access(sparse_mod.get(), &ref, true, 1);
    }
}

TEST_P(ArithVectorZeroedTest, Multiplication) {
    std::vector<double> zeroed(GetParam() ? nrow : ncol);

    if (GetParam()) {
        auto op = tatami::make_DelayedMultiplyVectorHelper<0>(std::move(zeroed));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    } else {
        auto op = tatami::make_DelayedMultiplyVectorHelper<1>(std::move(zeroed));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    }

    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::vector<double>(nrow * ncol));

    tatami_test::test_simple_column_access(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_column_access(sparse_mod.get(), &ref, true, 1); 

    tatami_test::test_simple_row_access(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_row_access(sparse_mod.get(), &ref, true, 1);
}

TEST_P(ArithVectorZeroedTest, DivisionAllZero) {
    std::vector<double> zeroed(GetParam() ? nrow : ncol);

    if (GetParam()) {
        auto op = tatami::make_DelayedDivideVectorHelper<true, 0>(std::move(zeroed));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    } else {
        auto op = tatami::make_DelayedDivideVectorHelper<true, 1>(std::move(zeroed));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    }

    auto copy = simulated;
    for (auto& x : copy) {
        x = careful_division(x, 0.0);
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());

    // Turning on NaN protection.
    tatami_test::test_simple_column_access<true>(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_column_access<true>(sparse_mod.get(), &ref, true, 1); 

    tatami_test::test_simple_row_access<true>(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_row_access<true>(sparse_mod.get(), &ref, true, 1);
}

TEST_P(ArithVectorZeroedTest, DivisionOneZero) {
    // But actually, even 1 zero is enough to break sparsity.
    std::vector<double> solo_zero(GetParam() ? nrow : ncol, 1);
    solo_zero[0] = 0;
    auto copy = simulated;

    if (GetParam()) {
        auto op = tatami::make_DelayedDivideVectorHelper<true, 0>(std::move(solo_zero));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

        for (int c = 0; c < ncol; ++c) {
            copy[c] = careful_division(copy[c], 0.0);
        }
    } else {
        auto op = tatami::make_DelayedDivideVectorHelper<true, 1>(std::move(solo_zero));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

        for (int r = 0; r < nrow; ++r) {
            copy[r * ncol] = careful_division(copy[r * ncol], 0.0);
        }
    }

    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());

    // Turning on NaN protection.
    tatami_test::test_simple_column_access<true>(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_column_access<true>(sparse_mod.get(), &ref, true, 1); 

    tatami_test::test_simple_row_access<true>(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_row_access<true>(sparse_mod.get(), &ref, true, 1);
}

TEST_P(ArithVectorZeroedTest, PowerAllZero) {
    std::vector<double> zeroed(GetParam() ? nrow : ncol);

    if (GetParam()) {
        auto op = tatami::make_DelayedPowerVectorHelper<true, 0>(std::move(zeroed));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    } else {
        auto op = tatami::make_DelayedPowerVectorHelper<true, 1>(std::move(zeroed));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    }

    auto copy = simulated;
    for (auto& x : copy) {
        x = std::pow(x, 0.0);
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());

    // Turning on NaN protection.
    tatami_test::test_simple_column_access<true>(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_column_access<true>(sparse_mod.get(), &ref, true, 1);

    tatami_test::test_simple_row_access<true>(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_row_access<true>(sparse_mod.get(), &ref, true, 1);
}

TEST_P(ArithVectorZeroedTest, PowerOneZero) {
    // But actually, even 1 zero is enough to break sparsity.
    std::vector<double> solo_zero(GetParam() ? nrow : ncol, 1);
    solo_zero[0] = 0;
    auto copy = simulated;

    if (GetParam()) {
        auto op = tatami::make_DelayedPowerVectorHelper<true, 0>(std::move(solo_zero));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

        for (int c = 0; c < ncol; ++c) {
            copy[c] = std::pow(copy[c], 0.0);
        }
    } else {
        auto op = tatami::make_DelayedPowerVectorHelper<true, 1>(std::move(solo_zero));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

        for (int r = 0; r < nrow; ++r) {
            copy[r * ncol] = std::pow(copy[r * ncol], 0.0);
        }
    }

    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());

    // Turning on NaN protection.
    tatami_test::test_simple_column_access<true>(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_column_access<true>(sparse_mod.get(), &ref, true, 1);

    tatami_test::test_simple_row_access<true>(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_row_access<true>(sparse_mod.get(), &ref, true, 1);
}

TEST_P(ArithVectorZeroedTest, Modulo) {
    std::vector<double> zeroed(GetParam() ? nrow : ncol);

    if (GetParam()) {
        auto op = tatami::make_DelayedModuloVectorHelper<true, 0>(std::move(zeroed));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    } else {
        auto op = tatami::make_DelayedModuloVectorHelper<true, 1>(std::move(zeroed));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    }

    auto copy = simulated;
    for (auto& x : copy) {
        x = std::fmod(x, 0.0);
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());

    // Turning on NaN protection.
    tatami_test::test_simple_column_access<true>(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_column_access<true>(sparse_mod.get(), &ref, true, 1);

    tatami_test::test_simple_row_access<true>(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_row_access<true>(sparse_mod.get(), &ref, true, 1);
}

TEST_P(ArithVectorZeroedTest, IntegerDivision) {
    std::vector<double> zeroed(GetParam() ? nrow : ncol);

    if (GetParam()) {
        auto op = tatami::make_DelayedIntegerDivideVectorHelper<true, 0>(std::move(zeroed));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    } else {
        auto op = tatami::make_DelayedIntegerDivideVectorHelper<true, 1>(std::move(zeroed));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    }

    auto copy = simulated;
    for (auto& x : copy) {
        // x == (x %% y) + y * (x %/% y)
        x = std::floor(careful_division(x, 0.0));
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());

    // Turning on NaN protection.
    tatami_test::test_simple_column_access<true>(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_column_access<true>(sparse_mod.get(), &ref, true, 1);

    tatami_test::test_simple_row_access<true>(dense_mod.get(), &ref, true, 1);
    tatami_test::test_simple_row_access<true>(sparse_mod.get(), &ref, true, 1);
}

INSTANTIATE_TEST_CASE_P(
    ArithVector,
    ArithVectorZeroedTest,
    ::testing::Values(true, false) // add by row, or by column
);
