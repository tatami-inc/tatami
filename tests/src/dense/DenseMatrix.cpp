#include <gtest/gtest.h>

#include <vector>
#include <deque>
#include <numeric>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami_test/tatami_test.hpp"

TEST(DenseMatrix, Basic) {
    std::vector<double> contents(200);
    double counter = -105;
    for (auto& i : contents) { i = counter++; }

    tatami::DenseColumnMatrix<double> mat(10, 20, contents);

    EXPECT_FALSE(mat.sparse());
    EXPECT_EQ(mat.sparse_proportion(), 0);
    EXPECT_FALSE(mat.prefer_rows());
    EXPECT_EQ(mat.prefer_rows_proportion(), 0);

    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);

    {
        auto wrk = mat.dense_column();
        for (size_t i = 0, end = mat.ncol(); i < end; ++i) {
            auto start = contents.begin() + i * mat.nrow();
            std::vector<double> expected(start, start + mat.nrow());
            auto observed = tatami_test::fetch<double, int>(wrk.get(), i, mat.nrow());
            EXPECT_EQ(observed, expected);
        }
    }

    {
        auto wrk = mat.dense_row();
        for (size_t i = 0, end = mat.nrow(); i < end; ++i) {
            std::vector<double> expected(mat.ncol());
            for (size_t j = 0, jend = mat.ncol(); j < jend; ++j) {
                expected[j] = contents[j * mat.nrow() + i];
            }
            auto observed = tatami_test::fetch<double, int>(wrk.get(), i, mat.ncol());
            EXPECT_EQ(observed, expected);
        }
    }

    EXPECT_FALSE(mat.uses_oracle(true));
}

TEST(DenseMatrix, OddsAndEnds) {
    std::vector<double> contents(200);
    double counter = -105;
    for (auto& i : contents) { i = counter++; }

    // Checks run properly.
    contents.clear();
    EXPECT_ANY_THROW({
        tatami::DenseColumnMatrix<double> mat(10, 20, contents);
    });
    EXPECT_ANY_THROW({
        tatami::DenseColumnMatrix<double> mat(10, 20, std::move(contents));
    });

    std::deque<double> more_contents(200);
    std::iota(more_contents.begin(), more_contents.end(), 1);
    tatami::DenseColumnMatrix<double, int, std::deque<double> > mat2(10, 20, more_contents);
    EXPECT_EQ(more_contents.size(), 200);
}

/*************************************
 *************************************/

class DenseTestMethods {
protected:
    inline static size_t nrow = 200, ncol = 100;
    inline static std::shared_ptr<tatami::NumericMatrix> dense_row, dense_column;

    static void assemble() {
        if (dense_row) {
            return;
        }

        auto simulated = tatami_test::simulate_dense_vector<double>(nrow * ncol, 0.05);
        auto transposed = std::vector<double>(nrow * ncol);
        for (size_t c = 0; c < ncol; ++c) {
            for (size_t r = 0; r < nrow; ++r) {
                transposed[c * nrow + r] = simulated[r * ncol + c];
            }
        }

        dense_row.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(simulated)));
        dense_column.reset(new tatami::DenseColumnMatrix<double, int>(nrow, ncol, std::move(transposed)));
        return;
    }
};

class DenseUtilsTest : public ::testing::Test, public DenseTestMethods {
protected:
    void SetUp() {
        assemble();
    }
};

TEST_F(DenseUtilsTest, Basic) {
    size_t NC = dense_column->ncol(), NR = dense_column->nrow();
    EXPECT_EQ(NC, ncol);
    EXPECT_EQ(NR, nrow);

    EXPECT_FALSE(dense_column->prefer_rows());
    EXPECT_TRUE(dense_row->prefer_rows());

    EXPECT_FALSE(dense_column->sparse());
    EXPECT_FALSE(dense_row->sparse());
}

/*************************************
 *************************************/

class DenseFullAccessTest : 
    public ::testing::TestWithParam<typename tatami_test::StandardTestAccessParameters>,
    public DenseTestMethods {
protected:
    void SetUp() {
        assemble();
    }
};

TEST_P(DenseFullAccessTest, Full) {
    auto tparam = GetParam(); 
    auto params = tatami_test::convert_access_parameters(tparam);
    tatami_test::test_full_access(params, dense_row.get(), dense_column.get());
    tatami_test::test_full_access(params, dense_column.get(), dense_row.get());
}

INSTANTIATE_TEST_SUITE_P(
    DenseMatrix,
    DenseFullAccessTest,
    TATAMI_TEST_STANDARD_ACCESS_PARAMETER_COMBINATIONS
);

/*************************************
 *************************************/

class DenseSlicedAccessTest : 
    public ::testing::TestWithParam<std::tuple<typename tatami_test::StandardTestAccessParameters, std::pair<double, double> > >,
    public DenseTestMethods {
protected:
    void SetUp() {
        assemble();
    }
};

TEST_P(DenseSlicedAccessTest, Sliced) {
    auto tparam = GetParam(); 
    auto params = tatami_test::convert_access_parameters(std::get<0>(tparam));

    auto interval_info = std::get<1>(tparam);
    auto len = (params.use_row ? ncol : nrow);
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;

    tatami_test::test_block_access(params, dense_column.get(), dense_row.get(), FIRST, LAST);
    tatami_test::test_block_access(params, dense_row.get(), dense_column.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    DenseMatrix,
    DenseSlicedAccessTest,
    ::testing::Combine(
        TATAMI_TEST_STANDARD_ACCESS_PARAMETER_COMBINATIONS,
        ::testing::Values(
            std::make_pair(0.0, 0.45),
            std::make_pair(0.2, 0.8), 
            std::make_pair(0.7, 1.0)
        )
    )
);

/*************************************
 *************************************/

class DenseIndexedAccessTest :
    public ::testing::TestWithParam<std::tuple<typename tatami_test::StandardTestAccessParameters, std::pair<double, int> > >,
    public DenseTestMethods {
protected:
    void SetUp() {
        assemble();
    }
};

TEST_P(DenseIndexedAccessTest, Indexed) {
    auto tparam = GetParam(); 
    auto params = tatami_test::convert_access_parameters(std::get<0>(tparam));

    auto interval_info = std::get<1>(tparam);
    auto len = (params.use_row ? ncol : nrow);
    size_t FIRST = interval_info.first * len, STEP = interval_info.second;

    tatami_test::test_indexed_access(params, dense_column.get(), dense_row.get(), FIRST, STEP);
    tatami_test::test_indexed_access(params, dense_row.get(), dense_column.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    DenseMatrix,
    DenseIndexedAccessTest,
    ::testing::Combine(
        TATAMI_TEST_STANDARD_ACCESS_PARAMETER_COMBINATIONS,
        ::testing::Values(
            std::make_pair(0.0, 5),
            std::make_pair(0.2, 10), 
            std::make_pair(0.7, 3)
        )
    )
);

/*************************************
 *************************************/

TEST(DenseMatrix, TypeOverflow) {
    // Check for correct pointer arithmetic when indices are supplied in a smaller integer type;
    // the class should convert them into size_t before doing any work.
    std::vector<double> contents(20000);
    double counter = -105;
    for (auto& i : contents) { i = counter++; }

    tatami::DenseColumnMatrix<double> ref(100, 200, contents);
    tatami::DenseColumnMatrix<double, unsigned char> limited(100, 200, contents);

    EXPECT_EQ(limited.nrow(), 100);
    EXPECT_EQ(limited.ncol(), 200);

    {
        auto rwrk = ref.dense_column();
        auto lwrk = limited.dense_column();
        for (int i = 0; i < ref.ncol(); ++i) {
            auto expected = tatami_test::fetch<double, int>(rwrk.get(), i, ref.nrow());
            auto observed = tatami_test::fetch<double, unsigned char>(lwrk.get(), i, limited.nrow());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        auto rwrk = ref.dense_row();
        auto lwrk = limited.dense_row();
        for (int i = 0; i < ref.nrow(); ++i) {
            auto expected = tatami_test::fetch<double, int>(rwrk.get(), i, ref.ncol());
            auto observed = tatami_test::fetch<double, unsigned char>(lwrk.get(), i, limited.ncol());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        auto rwrk = ref.dense_column(59, 189);
        auto lwrk = limited.dense_column(59, 189);
        for (int i = 0; i < ref.ncol(); ++i) {
            auto expected = tatami_test::fetch<double, int>(rwrk.get(), i, ref.nrow());
            auto observed = tatami_test::fetch<double, unsigned char>(lwrk.get(), i, limited.nrow());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        auto rwrk = ref.dense_row(59, 89);
        auto lwrk = limited.dense_row(59, 89);
        for (int i = 0; i < ref.nrow(); ++i) {
            auto expected = tatami_test::fetch<double, int>(rwrk.get(), i, ref.ncol());
            auto observed = tatami_test::fetch<double, unsigned char>(lwrk.get(), i, limited.ncol());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        std::vector<int> indices{ 11, 33, 55, 99, 111, 122, 155, 177, 199 };
        std::vector<unsigned char> uindices(indices.begin(), indices.end());

        auto rwrk = ref.dense_column(std::move(indices));
        auto lwrk = limited.dense_column(std::move(uindices));
        for (int i = 0; i < ref.ncol(); ++i) {
            auto expected = tatami_test::fetch<double, int>(rwrk.get(), i, ref.nrow());
            auto observed = tatami_test::fetch<double, unsigned char>(lwrk.get(), i, limited.nrow());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        std::vector<int> indices{ 10, 20, 40, 60, 80, 99 };
        std::vector<unsigned char> uindices(indices.begin(), indices.end());

        auto rwrk = ref.dense_row(std::move(indices));
        auto lwrk = limited.dense_row(std::move(uindices));
        for (int i = 0; i < ref.nrow(); ++i) {
            auto expected = tatami_test::fetch<double, int>(rwrk.get(), i, ref.ncol());
            auto observed = tatami_test::fetch<double, unsigned char>(lwrk.get(), i, limited.ncol());
            EXPECT_EQ(expected, observed);
        }
    }
}
