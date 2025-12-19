#include <gtest/gtest.h>

#include <vector>
#include <deque>
#include <numeric>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/utils/copy.hpp"
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

    const int NR = mat.nrow();
    EXPECT_EQ(NR, 10);
    const int NC = mat.ncol();
    EXPECT_EQ(NC, 20);

    {
        auto wrk = mat.dense_column();
        std::vector<double> expected;
        auto buffer = sanisizer::create<std::vector<double> >(NR);

        for (int c = 0; c < NC; ++c) {
            auto start = contents.begin() + sanisizer::product_unsafe<std::size_t>(c, NR);
            expected.clear();
            expected.insert(expected.end(), start, start + NR);

            auto observed = wrk->fetch(c, buffer.data());
            tatami::copy_n(observed, NR, buffer.data());
            EXPECT_EQ(expected, buffer);
        }
    }

    {
        auto wrk = mat.dense_row();
        sanisizer::as_size_type<std::vector<double> >(NC);
        std::vector<double> expected(NC), buffer(NC);

        for (int r = 0; r < NR; ++r) {
            for (int c = 0; c < NC; ++c) {
                expected[c] = contents[sanisizer::nd_offset<std::size_t>(r, NR, c)];
            }

            auto observed = wrk->fetch(r, buffer.data());
            tatami::copy_n(observed, NC, buffer.data());
            EXPECT_EQ(expected, buffer);
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
    inline static int nrow = 200, ncol = 100;
    inline static std::shared_ptr<tatami::NumericMatrix> dense_row, dense_column;

    static void assemble() {
        if (dense_row) {
            return;
        }

        const auto full_size = sanisizer::product<std::size_t>(nrow, ncol);
        auto simulated = tatami_test::simulate_vector<double>(full_size, []{
            tatami_test::SimulateVectorOptions opt;
            opt.seed = sanisizer::cap<tatami_test::SeedType>(438579088);
            return opt;
        }());

        std::vector<double> transposed(simulated.size());
        for (int c = 0; c < ncol; ++c) {
            for (int r = 0; r < nrow; ++r) {
                transposed[sanisizer::nd_offset<std::size_t>(r, nrow, c)] = simulated[sanisizer::nd_offset<std::size_t>(c, ncol, r)];
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
    const auto NC = dense_column->ncol(), NR = dense_column->nrow();
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

    const int NR = 100, NC = 200;
    tatami::DenseColumnMatrix<double, int> ref(NR, NC, contents);
    tatami::DenseColumnMatrix<double, unsigned char> limited(static_cast<unsigned char>(NR), static_cast<unsigned char>(NC), contents);

    EXPECT_EQ(limited.nrow(), NR);
    EXPECT_EQ(limited.ncol(), NC);

    {
        auto rwrk = ref.dense_column();
        auto lwrk = limited.dense_column();
        sanisizer::as_size_type<std::vector<double> >(NR);
        std::vector<double> expected(NR), observed(NR);

        for (int c = 0; c < NC; ++c) {
            auto eptr = rwrk->fetch(c, expected.data());
            tatami::copy_n(eptr, NR, expected.data());
            auto optr = lwrk->fetch(static_cast<unsigned char>(c), observed.data());
            tatami::copy_n(optr, NR, observed.data());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        auto rwrk = ref.dense_row();
        auto lwrk = limited.dense_row();
        sanisizer::as_size_type<std::vector<double> >(NC);
        std::vector<double> expected(NC), observed(NC);

        for (int r = 0; r < NR; ++r) {
            auto eptr = rwrk->fetch(r, expected.data());
            tatami::copy_n(eptr, NC, expected.data());
            auto optr = lwrk->fetch(static_cast<unsigned char>(r), observed.data());
            tatami::copy_n(optr, NC, observed.data());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        const int start = 59, len = 130;
        auto rwrk = ref.dense_column(start, len);
        auto lwrk = limited.dense_column(static_cast<unsigned char>(start), static_cast<unsigned char>(len));
        sanisizer::as_size_type<std::vector<double> >(len);
        std::vector<double> expected(len), observed(len);

        for (int c = 0; c < NC; ++c) {
            auto eptr = rwrk->fetch(c, expected.data());
            tatami::copy_n(eptr, len, expected.data());
            auto optr = lwrk->fetch(static_cast<unsigned char>(c), observed.data());
            tatami::copy_n(optr, len, observed.data());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        const int start = 59, len = 20;
        auto rwrk = ref.dense_row(start, len);
        auto lwrk = limited.dense_row(static_cast<unsigned char>(start), static_cast<unsigned char>(len));
        sanisizer::as_size_type<std::vector<double> >(len);
        std::vector<double> expected(len), observed(len);

        for (int r = 0; r < NR; ++r) {
            auto eptr = rwrk->fetch(r, expected.data());
            tatami::copy_n(eptr, len, expected.data());
            auto optr = lwrk->fetch(static_cast<unsigned char>(r), observed.data());
            tatami::copy_n(optr, len, observed.data());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        std::vector<int> indices{ 11, 33, 55, 99, 111, 122, 155, 177, 199 };
        std::vector<unsigned char> uindices(indices.begin(), indices.end());
        const auto num_i = indices.size();

        auto rwrk = ref.dense_column(std::move(indices));
        auto lwrk = limited.dense_column(std::move(uindices));
        sanisizer::as_size_type<std::vector<double> >(num_i);
        std::vector<double> expected(num_i), observed(num_i);

        for (int c = 0; c < NC; ++c) {
            auto eptr = rwrk->fetch(c, expected.data());
            tatami::copy_n(eptr, num_i, expected.data());
            auto optr = lwrk->fetch(static_cast<unsigned char>(c), observed.data());
            tatami::copy_n(optr, num_i, observed.data());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        std::vector<int> indices{ 10, 20, 40, 60, 80, 99 };
        std::vector<unsigned char> uindices(indices.begin(), indices.end());
        const auto num_i = indices.size();

        auto rwrk = ref.dense_row(std::move(indices));
        auto lwrk = limited.dense_row(std::move(uindices));
        sanisizer::as_size_type<std::vector<double> >(num_i);
        std::vector<double> expected(num_i), observed(num_i);

        for (int r = 0; r < NR; ++r) {
            auto eptr = rwrk->fetch(r, expected.data());
            tatami::copy_n(eptr, num_i, expected.data());
            auto optr = lwrk->fetch(static_cast<unsigned char>(r), observed.data());
            tatami::copy_n(optr, num_i, observed.data());
            EXPECT_EQ(expected, observed);
        }
    }
}

TEST(DenseMatrix, DifferentValueType) {
    // Checking that everything works when there is no .data() method of the right type.
    std::deque<int> contents(20000);
    double counter = -105;
    for (auto& i : contents) { i = counter++; }

    const int NR = 100, NC = 200;
    tatami::DenseColumnMatrix<double, int> ref(NR, NC, std::vector<double>(contents.begin(), contents.end()));
    tatami::DenseColumnMatrix<double, int, decltype(contents)> vstore(NR, NC, contents);

    EXPECT_EQ(vstore.nrow(), NR);
    EXPECT_EQ(vstore.ncol(), NC);

    {
        auto rwrk = ref.dense_column();
        auto lwrk = vstore.dense_column();
        sanisizer::as_size_type<std::vector<double> >(NR);
        std::vector<double> expected(NR), observed(NR);

        for (int c = 0; c < NC; ++c) {
            auto eptr = rwrk->fetch(c, expected.data());
            tatami::copy_n(eptr, NR, expected.data());
            auto optr = lwrk->fetch(c, observed.data());
            tatami::copy_n(optr, NR, observed.data());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        auto rwrk = ref.dense_row();
        auto lwrk = vstore.dense_row();
        sanisizer::as_size_type<std::vector<double> >(NC);
        std::vector<double> expected(NC), observed(NC);

        for (int r = 0; r < NR; ++r) {
            auto eptr = rwrk->fetch(r, expected.data());
            tatami::copy_n(eptr, NC, expected.data());
            auto optr = lwrk->fetch(r, observed.data());
            tatami::copy_n(optr, NC, observed.data());
            EXPECT_EQ(expected, observed);
        }
    }

    // Same for blocks, just to get some coverage.
    {
        const int start = 55, len = 45;
        auto rwrk = ref.dense_column(start, len);
        auto lwrk = vstore.dense_column(start, len);
        sanisizer::as_size_type<std::vector<double> >(len);
        std::vector<double> expected(len), observed(len);

        for (int c = 0; c < NC; ++c) {
            auto eptr = rwrk->fetch(c, expected.data());
            tatami::copy_n(eptr, len, expected.data());
            auto optr = lwrk->fetch(c, observed.data());
            tatami::copy_n(optr, len, observed.data());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        const int start = 10, len = 70;
        auto rwrk = ref.dense_row(start, len);
        auto lwrk = vstore.dense_row(start, len);
        sanisizer::as_size_type<std::vector<double> >(len);
        std::vector<double> expected(len), observed(len);

        for (int r = 0; r < NR; ++r) {
            auto eptr = rwrk->fetch(r, expected.data());
            tatami::copy_n(eptr, len, expected.data());
            auto optr = lwrk->fetch(r, observed.data());
            tatami::copy_n(optr, len, observed.data());
            EXPECT_EQ(expected, observed);
        }
    }

    // Same for indices, just to get some coverage.
    {
        std::vector<int> indices{ 12, 24, 25, 48, 50, 75, 96 };
        const auto num_i = indices.size();

        auto rwrk = ref.dense_column(indices);
        auto lwrk = vstore.dense_column(std::move(indices));
        sanisizer::as_size_type<std::vector<double> >(num_i);
        std::vector<double> expected(num_i), observed(num_i);

        for (int c = 0; c < NC; ++c) {
            auto eptr = rwrk->fetch(c, expected.data());
            tatami::copy_n(eptr, num_i, expected.data());
            auto optr = lwrk->fetch(c, observed.data());
            tatami::copy_n(optr, num_i, observed.data());
            EXPECT_EQ(expected, observed);
        }
    }

    {
        std::vector<int> indices{ 10, 15, 20, 25, 30, 35, 40, 45, 50, 96, 97, 98, 99 };
        const auto num_i = indices.size();

        auto rwrk = ref.dense_row(indices);
        auto lwrk = vstore.dense_row(std::move(indices));
        sanisizer::as_size_type<std::vector<double> >(num_i);
        std::vector<double> expected(num_i), observed(num_i);

        for (int r = 0; r < NR; ++r) {
            auto eptr = rwrk->fetch(r, expected.data());
            tatami::copy_n(eptr, num_i, expected.data());
            auto optr = lwrk->fetch(r, observed.data());
            tatami::copy_n(optr, num_i, observed.data());
            EXPECT_EQ(expected, observed);
        }
    }
}
