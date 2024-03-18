#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOp.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

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
        sparse = tatami::convert_to_compressed_sparse<false>(dense.get()); // column major.
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

using ArithVectorAdditionFullTest = ArithVectorAdditionTest<std::tuple<bool, bool, bool, tatami_test::TestAccessOrder, size_t> >;

TEST_P(ArithVectorAdditionFullTest, Basic) {
    auto tparam = GetParam();
    extra_assemble(std::get<0>(tparam));

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

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<1>(tparam);
    params.use_oracle = std::get<2>(tparam);
    params.order = std::get<3>(tparam);
    params.jump = std::get<4>(tparam);

    tatami_test::test_full_access(params, dense_mod.get(), ref.get());
    tatami_test::test_full_access(params, sparse_mod.get(), ref.get());
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorAdditionFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or by column
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorAdditionBlockTest = ArithVectorAdditionTest<std::tuple<bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > >;

TEST_P(ArithVectorAdditionBlockTest, Basic) {
    auto tparam = GetParam();
    extra_assemble(std::get<0>(tparam));

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<1>(tparam);
    params.use_oracle = std::get<2>(tparam);
    params.order = std::get<3>(tparam);
    params.jump = std::get<4>(tparam);

    auto interval_info = std::get<5>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;

    tatami_test::test_block_access(params, dense_mod.get(), ref.get(), FIRST, LAST);
    tatami_test::test_block_access(params, sparse_mod.get(), ref.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorAdditionBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair(0.0, 0.35),
            std::make_pair(0.38, 0.61), 
            std::make_pair(0.777, 1.0)
        )
    )
);

using ArithVectorAdditionIndexTest = ArithVectorAdditionTest<std::tuple<bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > >;

TEST_P(ArithVectorAdditionIndexTest, Basic) {
    auto tparam = GetParam();
    extra_assemble(std::get<0>(tparam));

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<1>(tparam);
    params.use_oracle = std::get<2>(tparam);
    params.order = std::get<3>(tparam);
    params.jump = std::get<4>(tparam);

    auto interval_info = std::get<5>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, STEP = interval_info.second * len;

    tatami_test::test_indexed_access(params, dense_mod.get(), ref.get(), FIRST, STEP);
    tatami_test::test_indexed_access(params, sparse_mod.get(), ref.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorAdditionIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair(0.0, 0.077),
            std::make_pair(0.21, 0.09), 
            std::make_pair(0.56, 0.01)
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

using ArithVectorSubtractionFullTest = ArithVectorSubtractionTest<std::tuple<bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t> >;

TEST_P(ArithVectorSubtractionFullTest, Basic) {
    auto tparam = GetParam();
    extra_assemble(std::get<0>(tparam), std::get<1>(tparam));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);

    tatami_test::test_full_access(params, dense_mod.get(), ref.get());
    tatami_test::test_full_access(params, sparse_mod.get(), ref.get());
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorSubtractionFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // by row or by column
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorSubtractionBlockTest = ArithVectorSubtractionTest<std::tuple<bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > >;

TEST_P(ArithVectorSubtractionBlockTest, Basic) {
    auto tparam = GetParam();
    extra_assemble(std::get<0>(tparam), std::get<1>(tparam));

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);

    auto interval_info = std::get<6>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;
    tatami_test::test_block_access(params, dense_mod.get(), ref.get(), FIRST, LAST);
    tatami_test::test_block_access(params, sparse_mod.get(), ref.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorSubtractionBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // iterate forward or back
        ::testing::Values(true, false), // with or without an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair(0.0, 0.541), 
            std::make_pair(0.111, 0.999), 
            std::make_pair(0.42, 1.0)
        )
    )
);

using ArithVectorSubtractionIndexTest = ArithVectorSubtractionTest<std::tuple<bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > >;

TEST_P(ArithVectorSubtractionIndexTest, Basic) {
    auto tparam = GetParam();
    extra_assemble(std::get<0>(tparam), std::get<1>(tparam));

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);

    auto interval_info = std::get<6>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, STEP = interval_info.second * len;
    tatami_test::test_indexed_access(params, dense_mod.get(), ref.get(), FIRST, STEP);
    tatami_test::test_indexed_access(params, sparse_mod.get(), ref.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorSubtractionIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair(0.0, 0.077),
            std::make_pair(0.21, 0.09), 
            std::make_pair(0.56, 0.01)
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

using ArithVectorMultiplicationFullTest = ArithVectorMultiplicationTest<std::tuple<bool, bool, bool, tatami_test::TestAccessOrder, size_t> >;

TEST_P(ArithVectorMultiplicationFullTest, Basic) {
    auto tparam = GetParam();
    extra_assemble(std::get<0>(tparam));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense_mod->sparse_proportion(), 0);
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(sparse_mod->sparse_proportion(), 1);
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    EXPECT_TRUE(dense_mod->prefer_rows());
    EXPECT_FALSE(sparse_mod->prefer_rows());

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<1>(tparam);
    params.use_oracle = std::get<2>(tparam);
    params.order = std::get<3>(tparam);
    params.jump = std::get<4>(tparam);

    tatami_test::test_full_access(params, dense_mod.get(), ref.get());
    tatami_test::test_full_access(params, sparse_mod.get(), ref.get());

    tatami_test::test_simple_row_access(dense_mod.get(), ref.get());
    tatami_test::test_simple_row_access(sparse_mod.get(), ref.get());
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorMultiplicationFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // by row or by column
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorMultiplicationBlockTest = ArithVectorMultiplicationTest<std::tuple<bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > >;

TEST_P(ArithVectorMultiplicationBlockTest, Basic) {
    auto tparam = GetParam();
    extra_assemble(std::get<0>(tparam));

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<1>(tparam);
    params.use_oracle = std::get<2>(tparam);
    params.order = std::get<3>(tparam);
    params.jump = std::get<4>(tparam);

    auto interval_info = std::get<5>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;
    tatami_test::test_block_access(params, dense_mod.get(), ref.get(), FIRST, LAST);
    tatami_test::test_block_access(params, sparse_mod.get(), ref.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorMultiplicationBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair(0.0, 0.35),
            std::make_pair(0.38, 0.61), 
            std::make_pair(0.777, 1.0)
        )
    )
);

using ArithVectorMultiplicationIndexTest = ArithVectorMultiplicationTest<std::tuple<bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > >;

TEST_P(ArithVectorMultiplicationIndexTest, Basic) {
    auto tparam = GetParam();
    extra_assemble(std::get<0>(tparam));

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<1>(tparam);
    params.use_oracle = std::get<2>(tparam);
    params.order = std::get<3>(tparam);
    params.jump = std::get<4>(tparam);

    auto interval_info = std::get<5>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, STEP = interval_info.second * len;
    tatami_test::test_indexed_access(params, dense_mod.get(), ref.get(), FIRST, STEP);
    tatami_test::test_indexed_access(params, sparse_mod.get(), ref.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorMultiplicationIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair(0.0, 0.077),
            std::make_pair(0.21, 0.09), 
            std::make_pair(0.56, 0.01)
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

using ArithVectorDivisionFullTest = ArithVectorDivisionTest<std::tuple<bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t> >;

TEST_P(ArithVectorDivisionFullTest, Basic) {
    auto tparam = GetParam();
    auto RIGHT = std::get<1>(tparam);
    extra_assemble(std::get<0>(tparam), RIGHT);

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

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);

    tatami_test::test_full_access(params, dense_mod.get(), ref.get());
    tatami_test::test_full_access(params, sparse_mod.get(), ref.get());
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorDivisionFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // by row or by column
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorDivisionBlockTest = ArithVectorDivisionTest<std::tuple<bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > >;

TEST_P(ArithVectorDivisionBlockTest, Basic) {
    auto tparam = GetParam();
    extra_assemble(std::get<0>(tparam), std::get<1>(tparam));

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);

    auto interval_info = std::get<6>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;
    tatami_test::test_block_access(params, dense_mod.get(), ref.get(), FIRST, LAST);
    tatami_test::test_block_access(params, sparse_mod.get(), ref.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorDivisionBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair( 0, 0.541 ), 
            std::make_pair( 0.111, 0.999 ), 
            std::make_pair( 0.42, 1 )
        )
    )
);

using ArithVectorDivisionIndexTest = ArithVectorDivisionTest<std::tuple<bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > >;

TEST_P(ArithVectorDivisionIndexTest, Basic) {
    auto tparam = GetParam();
    extra_assemble(std::get<0>(tparam), std::get<1>(tparam));

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);

    auto interval_info = std::get<6>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, STEP = interval_info.second * len;

    tatami_test::test_indexed_access(params, dense_mod.get(), ref.get(), FIRST, STEP);
    tatami_test::test_indexed_access(params, sparse_mod.get(), ref.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorDivisionIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair( 0, 0.077 ),
            std::make_pair( 0.21, 0.09 ), 
            std::make_pair( 0.56, 0.01 )
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

using ArithVectorPowerFullTest = ArithVectorPowerTest<std::tuple<bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t> >;

TEST_P(ArithVectorPowerFullTest, Basic) {
    auto tparam = GetParam();
    auto RIGHT = std::get<1>(tparam);
    extra_assemble(std::get<0>(tparam), RIGHT);

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

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);

    tatami_test::test_full_access(params, dense_mod.get(), ref.get());
    tatami_test::test_full_access(params, sparse_mod.get(), ref.get());
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorPowerFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // by row or by column
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorPowerBlockTest = ArithVectorPowerTest<std::tuple<bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > >;

TEST_P(ArithVectorPowerBlockTest, Basic) {
    auto tparam = GetParam();
    extra_assemble(std::get<0>(tparam), std::get<1>(tparam));

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);

    auto interval_info = std::get<6>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;
    tatami_test::test_block_access(params, dense_mod.get(), ref.get(), FIRST, LAST);
    tatami_test::test_block_access(params, sparse_mod.get(), ref.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorPowerBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair( 0, 0.541 ),
            std::make_pair( 0.111, 0.999 ),
            std::make_pair( 0.42, 1 )
        )
    )
);

using ArithVectorPowerIndexTest = ArithVectorPowerTest<std::tuple<bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > >;

TEST_P(ArithVectorPowerIndexTest, Basic) {
    auto tparam = GetParam();
    extra_assemble(std::get<0>(tparam), std::get<1>(tparam));

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);

    auto interval_info = std::get<6>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, STEP = interval_info.second * len;

    tatami_test::test_indexed_access(params, dense_mod.get(), ref.get(), FIRST, STEP);
    tatami_test::test_indexed_access(params, sparse_mod.get(), ref.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorPowerIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair( 0, 0.077 ),
            std::make_pair( 0.21, 0.09 ),
            std::make_pair( 0.56, 0.01 )
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

using ArithVectorModuloFullTest = ArithVectorModuloTest<std::tuple<bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t> >;

TEST_P(ArithVectorModuloFullTest, Basic) {
    auto tparam = GetParam();
    auto RIGHT = std::get<1>(tparam);
    extra_assemble(std::get<0>(tparam), RIGHT);

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

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);
    params.has_nan = !RIGHT;

    tatami_test::test_full_access(params, dense_mod.get(), ref.get());
    tatami_test::test_full_access(params, sparse_mod.get(), ref.get());
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorModuloFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // by row or by column
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorModuloBlockTest = ArithVectorModuloTest<std::tuple<bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > >;

TEST_P(ArithVectorModuloBlockTest, Basic) {
    auto tparam = GetParam();
    auto RIGHT = std::get<1>(tparam);
    extra_assemble(std::get<0>(tparam), RIGHT);

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);
    params.has_nan = !RIGHT;

    auto interval_info = std::get<6>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len , LAST = interval_info.second * len;

    tatami_test::test_block_access(params, dense_mod.get(), ref.get(), FIRST, LAST);
    tatami_test::test_block_access(params, sparse_mod.get(), ref.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorModuloBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair( 0, 0.541 ),
            std::make_pair( 0.111, 0.999 ),
            std::make_pair( 0.42, 1 )
        )
    )
);

using ArithVectorModuloIndexTest = ArithVectorModuloTest<std::tuple<bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > >;

TEST_P(ArithVectorModuloIndexTest, Basic) {
    auto tparam = GetParam();
    auto RIGHT = std::get<1>(tparam);
    extra_assemble(std::get<0>(tparam), RIGHT);

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);
    params.has_nan = !RIGHT;

    auto interval_info = std::get<6>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, STEP = interval_info.second * len;

    tatami_test::test_indexed_access(params, dense_mod.get(), ref.get(), FIRST, STEP);
    tatami_test::test_indexed_access(params, sparse_mod.get(), ref.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorModuloIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair( 0, 0.077 ),
            std::make_pair( 0.21, 0.09 ),
            std::make_pair( 0.56, 0.01 )
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

using ArithVectorIntegerDivisionFullTest = ArithVectorIntegerDivisionTest<std::tuple<bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t> >;

TEST_P(ArithVectorIntegerDivisionFullTest, Basic) {
    auto tparam = GetParam();
    auto RIGHT = std::get<1>(tparam);
    extra_assemble(std::get<0>(tparam), RIGHT);

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

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);

    tatami_test::test_full_access(params, dense_mod.get(), ref.get());
    tatami_test::test_full_access(params, sparse_mod.get(), ref.get());
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorIntegerDivisionFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // by row or by column
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

using ArithVectorIntegerDivisionBlockTest = ArithVectorIntegerDivisionTest<std::tuple<bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > >;

TEST_P(ArithVectorIntegerDivisionBlockTest, Basic) {
    auto tparam = GetParam();
    auto RIGHT = std::get<1>(tparam);
    extra_assemble(std::get<0>(tparam), RIGHT);

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);

    auto interval_info = std::get<6>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;

    tatami_test::test_block_access(params, dense_mod.get(), ref.get(), FIRST, LAST);
    tatami_test::test_block_access(params, sparse_mod.get(), ref.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorIntegerDivisionBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair( 0, 0.541 ),
            std::make_pair( 0.111, 0.999 ),
            std::make_pair( 0.42, 1 )
        )
    )
);

using ArithVectorIntegerDivisionIndexTest = ArithVectorIntegerDivisionTest<std::tuple<bool, bool, bool, bool, tatami_test::TestAccessOrder, size_t, std::pair<double, double> > >;

TEST_P(ArithVectorIntegerDivisionIndexTest, Basic) {
    auto tparam = GetParam();
    auto RIGHT = std::get<1>(tparam);
    extra_assemble(std::get<0>(tparam), RIGHT);

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<2>(tparam);
    params.use_oracle = std::get<3>(tparam);
    params.order = std::get<4>(tparam);
    params.jump = std::get<5>(tparam);

    auto interval_info = std::get<6>(tparam);
    auto len = (params.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, STEP = interval_info.second * len;

    tatami_test::test_indexed_access(params, dense_mod.get(), ref.get(), FIRST, STEP);
    tatami_test::test_indexed_access(params, sparse_mod.get(), ref.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorIntegerDivisionIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // add by row, or add by column.
        ::testing::Values(true, false), // on the right or left
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair( 0, 0.077 ),
            std::make_pair( 0.21, 0.09 ),
            std::make_pair( 0.56, 0.01 )
        )
    )
);

/**************************
 ********* ZEROED *********
 **************************/

class ArithVectorZeroedTest : public ::testing::TestWithParam<bool>, public ArithVectorUtils {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod, ref;

    void SetUp() {
        assemble();
    }

    void test_simple_row_access_wt_nan(const tatami::NumericMatrix* test, const tatami::NumericMatrix* ref) {
        tatami_test::TestAccessParameters params;
        params.has_nan = true;
        params.use_row = true;
        tatami_test::test_full_access(params, test, ref);
    }

    void test_simple_column_access_wt_nan(const tatami::NumericMatrix* test, const tatami::NumericMatrix* ref) {
        tatami_test::TestAccessParameters params;
        params.has_nan = true;
        params.use_row = false;
        tatami_test::test_full_access(params, test, ref);
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

    tatami_test::test_simple_column_access(dense_mod.get(), dense.get());
    tatami_test::test_simple_column_access(sparse_mod.get(), sparse.get()); 

    tatami_test::test_simple_row_access(dense_mod.get(), dense.get());
    tatami_test::test_simple_row_access(sparse_mod.get(), sparse.get());
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

        tatami_test::test_simple_column_access(dense_mod.get(), dense.get());
        tatami_test::test_simple_column_access(sparse_mod.get(), sparse.get()); 

        tatami_test::test_simple_row_access(dense_mod.get(), dense.get());
        tatami_test::test_simple_row_access(sparse_mod.get(), sparse.get());
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

        tatami_test::test_simple_column_access(dense_mod.get(), &ref);
        tatami_test::test_simple_column_access(sparse_mod.get(), &ref);

        tatami_test::test_simple_row_access(dense_mod.get(), &ref);
        tatami_test::test_simple_row_access(sparse_mod.get(), &ref);
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

    tatami_test::test_simple_column_access(dense_mod.get(), &ref);
    tatami_test::test_simple_column_access(sparse_mod.get(), &ref);

    tatami_test::test_simple_row_access(dense_mod.get(), &ref);
    tatami_test::test_simple_row_access(sparse_mod.get(), &ref);
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
    test_simple_column_access_wt_nan(dense_mod.get(), &ref);
    test_simple_column_access_wt_nan(sparse_mod.get(), &ref);

    test_simple_row_access_wt_nan(dense_mod.get(), &ref);
    test_simple_row_access_wt_nan(sparse_mod.get(), &ref);
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

        for (size_t c = 0; c < ncol; ++c) {
            copy[c] = careful_division(copy[c], 0.0);
        }
    } else {
        auto op = tatami::make_DelayedDivideVectorHelper<true, 1>(std::move(solo_zero));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

        for (size_t r = 0; r < nrow; ++r) {
            copy[r * ncol] = careful_division(copy[r * ncol], 0.0);
        }
    }

    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());

    // Turning on NaN protection.
    test_simple_column_access_wt_nan(dense_mod.get(), &ref);
    test_simple_column_access_wt_nan(sparse_mod.get(), &ref);

    test_simple_row_access_wt_nan(dense_mod.get(), &ref);
    test_simple_row_access_wt_nan(sparse_mod.get(), &ref);
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
    test_simple_column_access_wt_nan(dense_mod.get(), &ref);
    test_simple_column_access_wt_nan(sparse_mod.get(), &ref);

    test_simple_row_access_wt_nan(dense_mod.get(), &ref);
    test_simple_row_access_wt_nan(sparse_mod.get(), &ref);
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

        for (size_t c = 0; c < ncol; ++c) {
            copy[c] = std::pow(copy[c], 0.0);
        }
    } else {
        auto op = tatami::make_DelayedPowerVectorHelper<true, 1>(std::move(solo_zero));
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

        for (size_t r = 0; r < nrow; ++r) {
            copy[r * ncol] = std::pow(copy[r * ncol], 0.0);
        }
    }

    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(copy));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());

    // Turning on NaN protection.
    test_simple_column_access_wt_nan(dense_mod.get(), &ref);
    test_simple_column_access_wt_nan(sparse_mod.get(), &ref);

    test_simple_row_access_wt_nan(dense_mod.get(), &ref);
    test_simple_row_access_wt_nan(sparse_mod.get(), &ref);
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
    test_simple_column_access_wt_nan(dense_mod.get(), &ref);
    test_simple_column_access_wt_nan(sparse_mod.get(), &ref);

    test_simple_row_access_wt_nan(dense_mod.get(), &ref);
    test_simple_row_access_wt_nan(sparse_mod.get(), &ref);
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
    test_simple_column_access_wt_nan(dense_mod.get(), &ref);
    test_simple_column_access_wt_nan(sparse_mod.get(), &ref);

    test_simple_row_access_wt_nan(dense_mod.get(), &ref);
    test_simple_row_access_wt_nan(sparse_mod.get(), &ref);
}

INSTANTIATE_TEST_SUITE_P(
    ArithVector,
    ArithVectorZeroedTest,
    ::testing::Values(true, false) // add by row, or by column
);
