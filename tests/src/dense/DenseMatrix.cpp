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

    tatami::DenseColumnMatrix<double, int> mat(10, 20, contents);

    EXPECT_FALSE(mat.is_sparse());
    EXPECT_EQ(mat.is_sparse_proportion(), 0);
    EXPECT_FALSE(mat.prefer_rows());
    EXPECT_EQ(mat.prefer_rows_proportion(), 0);

    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);

    {
        auto wrk = mat.dense_column();
        for (size_t i = 0, end = mat.ncol(); i < end; ++i) {
            auto start = contents.begin() + i * mat.nrow();
            std::vector<double> expected(start, start + mat.nrow());
            auto observed = tatami_test::fetch<double, int>(*wrk, i, mat.nrow());
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
            auto observed = tatami_test::fetch<double, int>(*wrk, i, mat.ncol());
            EXPECT_EQ(observed, expected);
        }
    }

    EXPECT_FALSE(mat.uses_oracle(true));
}

TEST(DenseMatrix, Empty) {
    std::vector<double> contents;
    tatami::DenseColumnMatrix<double, int> mat(0, 20, contents);
    EXPECT_EQ(mat.nrow(), 0);
    EXPECT_EQ(mat.ncol(), 20);

    tatami::DenseColumnMatrix<double, int> mat2(20, 0, contents);
    EXPECT_EQ(mat2.nrow(), 20);
    EXPECT_EQ(mat2.ncol(), 0);
}

TEST(DenseMatrix, Errors) {
    std::vector<double> contents;
    tatami_test::throws_error([&]() {
        tatami::DenseColumnMatrix<double, int> mat(10, 20, contents);
    }, "length of 'values' should be equal");

    contents.push_back(1);
    tatami_test::throws_error([&]() {
        tatami::DenseColumnMatrix<double, int> mat(0, 20, contents);
    }, "length of 'values' should be equal");
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

        auto simulated = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.seed = 438579088;
            return opt;
        }());

        auto transposed = std::vector<double>(nrow * ncol);
        for (size_t c = 0; c < ncol; ++c) {
            for (size_t r = 0; r < nrow; ++r) {
                transposed[c * nrow + r] = simulated[r * ncol + c];
            }
        }

        dense_row.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(simulated)));
        dense_column.reset(new tatami::DenseColumnMatrix<double, int>(nrow, ncol, std::move(transposed)));
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

    EXPECT_FALSE(dense_column->is_sparse());
    EXPECT_FALSE(dense_row->is_sparse());
}

/*************************************
 *************************************/

class DenseFullAccessTest : 
    public ::testing::TestWithParam<tatami_test::StandardTestAccessOptions>,
    public DenseTestMethods {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DenseFullAccessTest, Full) {
    auto tparam = GetParam(); 
    auto opt = tatami_test::convert_test_access_options(tparam);
    tatami_test::test_full_access(*dense_row, *dense_column, opt);
    tatami_test::test_full_access(*dense_column, *dense_row, opt);
}

INSTANTIATE_TEST_SUITE_P(
    DenseMatrix,
    DenseFullAccessTest,
    tatami_test::standard_test_access_options_combinations()
);

/*************************************
 *************************************/

class DenseBlockAccessTest : 
    public ::testing::TestWithParam<std::tuple<tatami_test::StandardTestAccessOptions, std::pair<double, double> > >,
    public DenseTestMethods {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DenseBlockAccessTest, Block) {
    auto tparam = GetParam(); 
    auto opts = tatami_test::convert_test_access_options(std::get<0>(tparam));
    auto interval_info = std::get<1>(tparam);
    tatami_test::test_block_access(*dense_column, *dense_row, interval_info.first, interval_info.second, opts);
    tatami_test::test_block_access(*dense_row, *dense_column, interval_info.first, interval_info.second, opts);
}

INSTANTIATE_TEST_SUITE_P(
    DenseMatrix,
    DenseBlockAccessTest,
    ::testing::Combine(
        tatami_test::standard_test_access_options_combinations(),
        ::testing::Values(
            std::make_pair(0.0, 0.45),
            std::make_pair(0.2, 0.6), 
            std::make_pair(0.7, 0.3)
        )
    )
);

/*************************************
 *************************************/

class DenseIndexedAccessTest :
    public ::testing::TestWithParam<std::tuple<typename tatami_test::StandardTestAccessOptions, std::pair<double, double> > >,
    public DenseTestMethods {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DenseIndexedAccessTest, Indexed) {
    auto tparam = GetParam(); 
    auto opts = tatami_test::convert_test_access_options(std::get<0>(tparam));
    auto interval_info = std::get<1>(tparam);
    tatami_test::test_indexed_access(*dense_column, *dense_row, interval_info.first, interval_info.second, opts);
    tatami_test::test_indexed_access(*dense_row, *dense_column, interval_info.first, interval_info.second, opts);
}

INSTANTIATE_TEST_SUITE_P(
    DenseMatrix,
    DenseIndexedAccessTest,
    ::testing::Combine(
        tatami_test::standard_test_access_options_combinations(),
        ::testing::Values(
            std::make_pair(0.0, 0.15),
            std::make_pair(0.2, 0.25), 
            std::make_pair(0.7, 0.3)
        )
    )
);

/*************************************
 *************************************/

TEST(DenseMatrix, IndexTypeOverflow) {
    // Check for correct pointer arithmetic when indices are supplied in a smaller integer type;
    // the class should convert them into size_t before doing any work.
    std::vector<double> contents(20000);
    double counter = -105;
    for (auto& i : contents) { i = counter++; }

    tatami::DenseColumnMatrix<double, int> ref(100, 200, contents);
    tatami::DenseColumnMatrix<double, unsigned char> limited(100, 200, contents);

    EXPECT_EQ(limited.nrow(), 100);
    EXPECT_EQ(limited.ncol(), 200);

    {
        auto rwrk = ref.dense_column();
        auto lwrk = limited.dense_column();
        for (int i = 0; i < ref.ncol(); ++i) {
            auto expected = tatami_test::fetch<double, int>(*rwrk, i, ref.nrow());
            auto observed = tatami_test::fetch<double, unsigned char>(*lwrk, i, limited.nrow());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        auto rwrk = ref.dense_row();
        auto lwrk = limited.dense_row();
        for (int i = 0; i < ref.nrow(); ++i) {
            auto expected = tatami_test::fetch<double, int>(*rwrk, i, ref.ncol());
            auto observed = tatami_test::fetch<double, unsigned char>(*lwrk, i, limited.ncol());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        auto rwrk = ref.dense_column(59, 189);
        auto lwrk = limited.dense_column(59, 189);
        for (int i = 0; i < ref.ncol(); ++i) {
            auto expected = tatami_test::fetch<double, int>(*rwrk, i, ref.nrow());
            auto observed = tatami_test::fetch<double, unsigned char>(*lwrk, i, limited.nrow());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        auto rwrk = ref.dense_row(59, 89);
        auto lwrk = limited.dense_row(59, 89);
        for (int i = 0; i < ref.nrow(); ++i) {
            auto expected = tatami_test::fetch<double, int>(*rwrk, i, ref.ncol());
            auto observed = tatami_test::fetch<double, unsigned char>(*lwrk, i, limited.ncol());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        std::vector<int> indices{ 11, 33, 55, 99, 111, 122, 155, 177, 199 };
        std::vector<unsigned char> uindices(indices.begin(), indices.end());

        auto rwrk = ref.dense_column(std::move(indices));
        auto lwrk = limited.dense_column(std::move(uindices));
        for (int i = 0; i < ref.ncol(); ++i) {
            auto expected = tatami_test::fetch<double, int>(*rwrk, i, ref.nrow());
            auto observed = tatami_test::fetch<double, unsigned char>(*lwrk, i, limited.nrow());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        std::vector<int> indices{ 10, 20, 40, 60, 80, 99 };
        std::vector<unsigned char> uindices(indices.begin(), indices.end());

        auto rwrk = ref.dense_row(std::move(indices));
        auto lwrk = limited.dense_row(std::move(uindices));
        for (int i = 0; i < ref.nrow(); ++i) {
            auto expected = tatami_test::fetch<double, int>(*rwrk, i, ref.ncol());
            auto observed = tatami_test::fetch<double, unsigned char>(*lwrk, i, limited.ncol());
            EXPECT_EQ(expected, observed);
        }
    }
}

TEST(DenseMatrix, DifferentValueType) {
    // Checking that everything works when there is no .data() method of the right type.
    std::deque<int> contents(20000);
    double counter = -105;
    for (auto& i : contents) { i = counter++; }

    tatami::DenseColumnMatrix<double, int> ref(100, 200, std::vector<double>(contents.begin(), contents.end()));
    tatami::DenseColumnMatrix<double, int, decltype(contents)> vstore(100, 200, contents);

    EXPECT_EQ(vstore.nrow(), 100);
    EXPECT_EQ(vstore.ncol(), 200);

    {
        auto rwrk = ref.dense_column();
        auto lwrk = vstore.dense_column();
        for (int i = 0; i < ref.ncol(); ++i) {
            auto expected = tatami_test::fetch<double, int>(*rwrk, i, ref.nrow());
            auto observed = tatami_test::fetch<double, int>(*lwrk, i, vstore.nrow());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        auto rwrk = ref.dense_row();
        auto lwrk = vstore.dense_row();
        for (int i = 0; i < ref.nrow(); ++i) {
            auto expected = tatami_test::fetch<double, int>(*rwrk, i, ref.ncol());
            auto observed = tatami_test::fetch<double, int>(*lwrk, i, vstore.ncol());
            EXPECT_EQ(expected, observed);
        }
    }

    // Same for blocks, just to make sure we get test coverage.
    {
        auto rwrk = ref.dense_column(55, 45);
        auto lwrk = vstore.dense_column(55, 45);
        for (int i = 0; i < ref.ncol(); ++i) {
            auto expected = tatami_test::fetch<double, int>(*rwrk, i, 45);
            auto observed = tatami_test::fetch<double, int>(*lwrk, i, 45);
            EXPECT_EQ(expected, observed);
        }
    }

    {
        auto rwrk = ref.dense_row(10, 70);
        auto lwrk = vstore.dense_row(10, 70);
        for (int i = 0; i < ref.nrow(); ++i) {
            auto expected = tatami_test::fetch<double, int>(*rwrk, i, 70);
            auto observed = tatami_test::fetch<double, int>(*lwrk, i, 70);
            EXPECT_EQ(expected, observed);
        }
    }
}
