#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/other/DelayedTranspose.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

class TransposeUtils {
protected:
    inline static size_t nrow = 199, ncol = 201;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, sparse, tdense, tsparse, ref;

    static void assemble() {
        if (ref) {
            return;
        }

        auto simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.05);
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulated));
        sparse = tatami::convert_to_compressed_sparse<false, double, int>(dense.get()); // column-major.
        tdense = tatami::make_DelayedTranspose(dense);
        tsparse = tatami::make_DelayedTranspose(sparse);

        std::vector<double> refvec(nrow * ncol);
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                refvec[c * nrow + r] = simulated[r * ncol + c];
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double, int>(ncol, nrow, refvec));
    }
};

class TransposeTest : public ::testing::Test, public TransposeUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_F(TransposeTest, Basic) {
    EXPECT_EQ(tdense->ncol(), nrow);
    EXPECT_EQ(tdense->nrow(), ncol);
    EXPECT_FALSE(tdense->is_sparse());
    EXPECT_EQ(tdense->is_sparse_proportion(), 0);
    EXPECT_FALSE(tdense->prefer_rows());
    EXPECT_EQ(tdense->prefer_rows_proportion(), 0);

    EXPECT_EQ(tsparse->ncol(), nrow);
    EXPECT_EQ(tsparse->nrow(), ncol);
    EXPECT_TRUE(tsparse->is_sparse());
    EXPECT_EQ(tsparse->is_sparse_proportion(), 1);
    EXPECT_TRUE(tsparse->prefer_rows());
    EXPECT_EQ(tsparse->prefer_rows_proportion(), 1);

    EXPECT_FALSE(tsparse->uses_oracle(false));
}

TEST_F(TransposeTest, ConstOverload) {
    std::shared_ptr<const tatami::NumericMatrix> const_dense(dense);
    auto tdense = tatami::make_DelayedTranspose(const_dense);

    // Cursory checks.
    EXPECT_EQ(dense->nrow(), tdense->ncol());
    EXPECT_EQ(dense->ncol(), tdense->nrow());
}

class TransposeFullTest : public ::testing::TestWithParam<tatami_test::StandardTestAccessParameters>, public TransposeUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(TransposeFullTest, Basic) {
    auto tparam = GetParam(); 
    auto params = tatami_test::convert_access_parameters(tparam);
    tatami_test::test_full_access(params, tdense.get(), ref.get());
    tatami_test::test_full_access(params, tsparse.get(), ref.get());
}

INSTANTIATE_TEST_SUITE_P(
    TransposeTest,
    TransposeFullTest,
    tatami_test::standard_test_access_parameter_combinations()
);

class TransposeBlockTest : public ::testing::TestWithParam<std::tuple<tatami_test::StandardTestAccessParameters, std::pair<double, double> > >, public TransposeUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(TransposeBlockTest, Basic) {
    auto tparam = GetParam(); 
    auto params = tatami_test::convert_access_parameters(std::get<0>(tparam));

    auto interval_info = std::get<1>(tparam);
    auto len = (params.use_row ? tdense->ncol() : tdense->nrow());
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;

    tatami_test::test_block_access(params, tdense.get(), ref.get(), FIRST, LAST);
    tatami_test::test_block_access(params, tsparse.get(), ref.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    TransposeTest,
    TransposeBlockTest,
    ::testing::Combine(
        tatami_test::standard_test_access_parameter_combinations(),
        ::testing::Values(
            std::make_pair(0, 0.44),
            std::make_pair(0.21, 0.89), 
            std::make_pair(0.33, 1)
        )
    )
);

class TransposeIndexTest : public ::testing::TestWithParam<std::tuple<tatami_test::StandardTestAccessParameters, std::pair<double, int> > >, public TransposeUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(TransposeIndexTest, Basic) {
    auto tparam = GetParam(); 
    auto params = tatami_test::convert_access_parameters(std::get<0>(tparam));

    auto interval_info = std::get<1>(tparam);
    auto len = (params.use_row ? tdense->ncol() : tdense->nrow());
    size_t FIRST = interval_info.first * len, STEP = interval_info.second;

    tatami_test::test_indexed_access(params, tdense.get(), ref.get(), FIRST, STEP);
    tatami_test::test_indexed_access(params, tsparse.get(), ref.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    TransposeTest,
    TransposeIndexTest,
    ::testing::Combine(
        tatami_test::standard_test_access_parameter_combinations(),
        ::testing::Values(
            std::make_pair(0, 5),
            std::make_pair(0.2, 10), 
            std::make_pair(0.7, 3)
        )
    )
);
