#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/dense/DenseMatrix.hpp"
#include "tatami/base/isometric/DelayedIsometricOp.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../_tests/test_column_access.h"
#include "../_tests/test_row_access.h"
#include "../_tests/test_oracle_access.h"
#include "../_tests/simulate_vector.h"

template<class PARAM>
class ArithVectorTest : public ::testing::TestWithParam<PARAM> {
protected:
    size_t nrow = 291, ncol = 188;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    std::vector<double> simulated;
protected:
    void SetUp() {
        simulated = simulate_sparse_vector<double>(nrow * ncol, 0.1);
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
class ArithVectorAdditionTest : public ArithVectorTest<PARAM> {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    std::vector<double> vec;

    void extra_assemble(bool row) {
        vec = this->create_vector(row ? this->nrow : this->ncol, 5, 0.5);

        if (row) {
            auto op = tatami::make_DelayedAddVectorHelper<0>(vec);
            dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
        } else {
            auto op = tatami::make_DelayedAddVectorHelper<1>(vec);
            dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
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

TEST_P(ArithVectorAdditionFullTest, Column) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));
    bool FORWARD = std::get<1>(param);
    int JUMP = std::get<2>(param);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());

    test_simple_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    test_simple_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
}

TEST_P(ArithVectorAdditionFullTest, Row) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));
    bool FORWARD = std::get<1>(param);
    int JUMP = std::get<2>(param);

    test_simple_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    test_simple_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
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

TEST_P(ArithVectorAdditionBlockTest, Column) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));

    bool FORWARD = std::get<1>(param);
    size_t JUMP = std::get<2>(param);
    auto interval_info = std::get<3>(param);
    size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;

    test_sliced_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    test_sliced_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
}

TEST_P(ArithVectorAdditionBlockTest, Row) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));

    bool FORWARD = std::get<1>(param);
    size_t JUMP = std::get<2>(param);
    auto interval_info = std::get<3>(param);
    size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;

    test_sliced_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    test_sliced_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
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

TEST_P(ArithVectorAdditionIndexTest, Column) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));

    bool FORWARD = std::get<1>(param);
    size_t JUMP = std::get<2>(param);
    auto interval_info = std::get<3>(param);
    size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;

    test_indexed_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    test_indexed_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
}

TEST_P(ArithVectorAdditionIndexTest, Row) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));

    bool FORWARD = std::get<1>(param);
    size_t JUMP = std::get<2>(param);
    auto interval_info = std::get<3>(param);
    size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;

    test_indexed_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    test_indexed_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
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
class ArithVectorSubtractionTest : public ArithVectorTest<PARAM> {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    std::vector<double> vec;

    void extra_assemble(bool row, bool right) {
        vec = this->create_vector(row ? this->nrow : this->ncol, -10, 2.222);

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

TEST_P(ArithVectorSubtractionFullTest, Column) {
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

    test_simple_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    test_simple_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
}

TEST_P(ArithVectorSubtractionFullTest, Row) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));
    bool FORWARD = std::get<2>(param);
    int JUMP = std::get<3>(param);

    test_simple_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    test_simple_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
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

TEST_P(ArithVectorSubtractionBlockTest, Column) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;

    test_sliced_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    test_sliced_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
}

TEST_P(ArithVectorSubtractionBlockTest, Row) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;

    test_sliced_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    test_sliced_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
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

TEST_P(ArithVectorSubtractionIndexTest, Column) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;

    test_indexed_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    test_indexed_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
}

TEST_P(ArithVectorSubtractionIndexTest, Row) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;

    test_indexed_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    test_indexed_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
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
class ArithVectorMultiplicationTest : public ArithVectorTest<PARAM> {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    std::vector<double> vec;

    void extra_assemble(bool row) {
        vec = this->create_vector(row ? this->nrow : this->ncol, 99, -2.5);

        if (row) {
            auto op = tatami::make_DelayedMultiplyVectorHelper<0>(vec);
            dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
        } else {
            auto op = tatami::make_DelayedMultiplyVectorHelper<1>(vec);
            dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
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

TEST_P(ArithVectorMultiplicationFullTest, Column) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));
    bool FORWARD = std::get<1>(param);
    int JUMP = std::get<2>(param);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());

    test_simple_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    test_simple_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
}

TEST_P(ArithVectorMultiplicationFullTest, Row) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));
    bool FORWARD = std::get<1>(param);
    int JUMP = std::get<2>(param);

    test_simple_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    test_simple_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
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

TEST_P(ArithVectorMultiplicationBlockTest, Column) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));

    bool FORWARD = std::get<1>(param);
    size_t JUMP = std::get<2>(param);
    auto interval_info = std::get<3>(param);
    size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;

    test_sliced_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    test_sliced_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
}

TEST_P(ArithVectorMultiplicationBlockTest, Row) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));

    bool FORWARD = std::get<1>(param);
    size_t JUMP = std::get<2>(param);
    auto interval_info = std::get<3>(param);
    size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;

    test_sliced_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    test_sliced_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
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

TEST_P(ArithVectorMultiplicationIndexTest, Column) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));

    bool FORWARD = std::get<1>(param);
    size_t JUMP = std::get<2>(param);
    auto interval_info = std::get<3>(param);
    size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;

    test_indexed_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    test_indexed_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
}

TEST_P(ArithVectorMultiplicationIndexTest, Row) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));

    bool FORWARD = std::get<1>(param);
    size_t JUMP = std::get<2>(param);
    auto interval_info = std::get<3>(param);
    size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;

    test_indexed_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    test_indexed_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
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
class ArithVectorDivisionTest : public ArithVectorTest<PARAM> {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;
    std::vector<double> vec;

    void extra_assemble(bool row, bool right) {
        vec = this->create_vector(row ? this->nrow : this->ncol, -19, 2.11);

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

        auto refvec = this->simulated;
        for (size_t r = 0; r < this->nrow; ++r) {
            for (size_t c = 0; c < this->ncol; ++c) {
                auto& x = refvec[r * this->ncol + c];
                if (right) {
                    x /= vec[row ? r : c];
                } else {
                    if (x) {
                        x = vec[row ? r : c] / x;
                    } else {
                        x = std::numeric_limits<double>::infinity();
                    }
                }
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(this->nrow, this->ncol, std::move(refvec)));
    }
};

using ArithVectorDivisionFullTest = ArithVectorDivisionTest<std::tuple<bool, bool, bool, size_t> >;

TEST_P(ArithVectorDivisionFullTest, Column) {
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

    test_simple_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    test_simple_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
}

TEST_P(ArithVectorDivisionFullTest, Row) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));
    bool FORWARD = std::get<2>(param);
    int JUMP = std::get<3>(param);

    test_simple_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP);
    test_simple_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP);
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

TEST_P(ArithVectorDivisionBlockTest, Column) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;

    test_sliced_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    test_sliced_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
}

TEST_P(ArithVectorDivisionBlockTest, Row) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;

    test_sliced_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    test_sliced_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
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

TEST_P(ArithVectorDivisionIndexTest, Column) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;

    test_indexed_column_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    test_indexed_column_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
}

TEST_P(ArithVectorDivisionIndexTest, Row) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param), std::get<1>(param));

    bool FORWARD = std::get<2>(param);
    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;

    test_indexed_row_access(dense_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    test_indexed_row_access(sparse_mod.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
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

/**************************
 ********* ORACLE *********
 **************************/

class ArithVectorOracleTest : public ArithVectorTest<std::tuple<bool, bool> > {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, wrapped_dense_mod, wrapped_sparse_mod;
    std::vector<double> vec;

    void extra_assemble(bool row) {
        vec = this->create_vector(row ? this->nrow : this->ncol, 5, 0.5);

        if (row) {
            auto op = tatami::make_DelayedAddVectorHelper<0>(vec);
            dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
            wrapped_dense_mod = tatami::make_DelayedIsometricOp(make_CrankyMatrix(this->dense), op);
            wrapped_sparse_mod = tatami::make_DelayedIsometricOp(make_CrankyMatrix(this->sparse), op);
        } else {
            auto op = tatami::make_DelayedAddVectorHelper<1>(vec);
            dense_mod = tatami::make_DelayedIsometricOp(this->dense, op);
            sparse_mod = tatami::make_DelayedIsometricOp(this->sparse, op);
            wrapped_dense_mod = tatami::make_DelayedIsometricOp(make_CrankyMatrix(this->dense), op);
            wrapped_sparse_mod = tatami::make_DelayedIsometricOp(make_CrankyMatrix(this->sparse), op);
        }
    }
};

TEST_P(ArithVectorOracleTest, Validate) {
    auto param = GetParam();
    extra_assemble(std::get<0>(param));
    EXPECT_FALSE(dense_mod->uses_oracle(true));
    EXPECT_TRUE(wrapped_dense_mod->uses_oracle(true));

    test_oracle_column_access(wrapped_dense_mod.get(), dense_mod.get(), std::get<1>(param));
    test_oracle_column_access(wrapped_sparse_mod.get(), sparse_mod.get(), std::get<1>(param));

    test_oracle_row_access(wrapped_dense_mod.get(), dense_mod.get(), std::get<1>(param));
    test_oracle_row_access(wrapped_sparse_mod.get(), sparse_mod.get(), std::get<1>(param));
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
    auto simulated = simulate_sparse_vector<double>(nrow * ncol, 0.1);
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated));

    auto vec = std::vector<double>(nrow);
    auto op = tatami::make_DelayedAddVectorHelper<0>(vec);
    auto mat = tatami::make_DelayedIsometricOp(dense, std::move(op));

    // cursory checks.
    EXPECT_EQ(mat->nrow(), dense->nrow());
    EXPECT_EQ(mat->ncol(), dense->ncol());
}
