#include <gtest/gtest.h>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/CompressedSparseMatrix.hpp"
#include "tatami/utils/convert_to_dense.hpp"

#include "tatami_test/tatami_test.hpp"

class ConvertToDenseTest : public ::testing::TestWithParam<std::tuple<int, int> > {
protected:
    size_t NR, NC;

    template<class Param>
    void assemble(const Param& param) {
        NR = std::get<0>(param);
        NC = std::get<1>(param);
    }
};

TEST_P(ConvertToDenseTest, RowToRow) {
    assemble(GetParam());
    auto vec = tatami_test::simulate_dense_vector<double>(NR * NC);
    tatami::DenseMatrix<true, double, int> mat(NR, NC, vec);

    auto converted = tatami::convert_to_dense<true>(&mat);
    EXPECT_TRUE(converted->prefer_rows());
    EXPECT_FALSE(converted->sparse());
    tatami_test::test_simple_row_access(converted.get(), &mat, true, 1);
    tatami_test::test_simple_column_access(converted.get(), &mat, true, 1);

    auto converted2 = tatami::convert_to_dense<true, int, size_t>(&mat); // works for a different type.
    EXPECT_TRUE(converted2->prefer_rows());
    EXPECT_FALSE(converted2->sparse());

    auto wrk2 = converted2->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto start = vec.begin() + i * NC;
        std::vector<int> expected2(start, start + NC);
        EXPECT_EQ(wrk2->fetch(i), expected2);
    }

    // Works in parallel.
    auto convertedP = tatami::convert_to_dense<true>(&mat, 3);
    tatami_test::test_simple_row_access(convertedP.get(), &mat, true, 1);
    tatami_test::test_simple_column_access(convertedP.get(), &mat, true, 1);
}

TEST_P(ConvertToDenseTest, ColumnToColumn) {
    assemble(GetParam());
    auto vec = tatami_test::simulate_dense_vector<double>(NR * NC);
    tatami::DenseMatrix<false, double, int> mat(NR, NC, vec);

    auto converted = tatami::convert_to_dense<false>(&mat);
    EXPECT_FALSE(converted->prefer_rows());
    EXPECT_FALSE(converted->sparse());
    tatami_test::test_simple_row_access(converted.get(), &mat, true, 1);
    tatami_test::test_simple_column_access(converted.get(), &mat, true, 1);

    auto converted2 = tatami::convert_to_dense<false, int, size_t>(&mat); // works for a different type.
    EXPECT_FALSE(converted2->prefer_rows());
    EXPECT_FALSE(converted2->sparse());

    auto wrk2 = converted2->dense_column();
    for (size_t i = 0; i < NC; ++i) {
        auto start = vec.begin() + i * NR;
        std::vector<int> expected2(start, start + NR);
        EXPECT_EQ(wrk2->fetch(i), expected2);
    }
}

TEST_P(ConvertToDenseTest, RowToColumn) {
    assemble(GetParam());
    auto vec = tatami_test::simulate_dense_vector<double>(NR * NC);
    tatami::DenseMatrix<true, double, int> mat(NR, NC, vec);

    auto converted = tatami::convert_to_dense<false>(&mat);
    EXPECT_FALSE(converted->prefer_rows());
    EXPECT_FALSE(converted->sparse());
    tatami_test::test_simple_row_access(converted.get(), &mat, true, 1);
    tatami_test::test_simple_column_access(converted.get(), &mat, true, 1);

    auto converted2 = tatami::convert_to_dense<false, int, size_t>(&mat); // works for a different type.
    EXPECT_FALSE(converted2->prefer_rows());
    EXPECT_FALSE(converted2->sparse());

    auto wrk2 = converted2->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto start = vec.begin() + i * NC;
        std::vector<int> expected2(start, start + NC);
        EXPECT_EQ(wrk2->fetch(i), expected2);
    }

    // Works in parallel.
    auto convertedP = tatami::convert_to_dense<false>(&mat, 3);
    tatami_test::test_simple_row_access(convertedP.get(), &mat, true, 1);
    tatami_test::test_simple_column_access(convertedP.get(), &mat, true, 1);
}

TEST_P(ConvertToDenseTest, ColumnToRow) {
    assemble(GetParam());
    auto vec = tatami_test::simulate_dense_vector<double>(NR * NC);
    tatami::DenseMatrix<false, double, int> mat(NR, NC, vec);

    auto converted = tatami::convert_to_dense<true>(&mat);
    EXPECT_TRUE(converted->prefer_rows());
    EXPECT_FALSE(converted->sparse());
    tatami_test::test_simple_row_access(converted.get(), &mat, true, 1);
    tatami_test::test_simple_column_access(converted.get(), &mat, true, 1);

    auto converted2 = tatami::convert_to_dense<true, int, size_t>(&mat); // works for a different type.
    EXPECT_TRUE(converted2->prefer_rows());
    EXPECT_FALSE(converted2->sparse());

    auto wrk2 = converted2->dense_column();
    for (size_t i = 0; i < NC; ++i) {
        auto start = vec.begin() + i * NR;
        std::vector<int> expected2(start, start + NR);
        EXPECT_EQ(wrk2->fetch(i), expected2);
    }
}

TEST_P(ConvertToDenseTest, Automatic) {
    assemble(GetParam());
    auto vec = tatami_test::simulate_dense_vector<double>(NR * NC);

    {
        tatami::DenseMatrix<false, double, int> mat(NR, NC, vec);
        auto converted = tatami::convert_to_dense(&mat, -1);
        EXPECT_FALSE(converted->prefer_rows());
    }

    {
        tatami::DenseMatrix<true, double, int> mat(NR, NC, vec);
        auto converted = tatami::convert_to_dense(&mat, -1);
        EXPECT_TRUE(converted->prefer_rows());
    }
}

TEST_P(ConvertToDenseTest, FromSparse) {
    assemble(GetParam());

    // From a row-major sparse matrix.
    {
        auto vec = tatami_test::simulate_sparse_compressed<double>(NR, NC, 0.2);
        tatami::CompressedSparseRowMatrix<double, int> smat(NR, NC, std::move(vec.value), std::move(vec.index), std::move(vec.ptr));

        {
            auto converted = tatami::convert_to_dense<true>(&smat);
            EXPECT_TRUE(converted->prefer_rows());
            EXPECT_FALSE(converted->sparse());
            tatami_test::test_simple_row_access(converted.get(), &smat, true, 1);
            tatami_test::test_simple_column_access(converted.get(), &smat, true, 1);
        }

        {
            auto converted = tatami::convert_to_dense<false>(&smat);
            EXPECT_FALSE(converted->prefer_rows());
            EXPECT_FALSE(converted->sparse());
            tatami_test::test_simple_row_access(converted.get(), &smat, true, 1);
            tatami_test::test_simple_column_access(converted.get(), &smat, true, 1);
        }
    }

    // From a row-major sparse matrix.
    {
        auto vec = tatami_test::simulate_sparse_compressed<double>(NC, NR, 0.2);
        tatami::CompressedSparseColumnMatrix<double, int> smat(NR, NC, std::move(vec.value), std::move(vec.index), std::move(vec.ptr));

        {
            auto converted = tatami::convert_to_dense<true>(&smat);
            EXPECT_TRUE(converted->prefer_rows());
            EXPECT_FALSE(converted->sparse());
            tatami_test::test_simple_row_access(converted.get(), &smat, true, 1);
            tatami_test::test_simple_column_access(converted.get(), &smat, true, 1);
        }

        {
            auto converted = tatami::convert_to_dense<false>(&smat);
            EXPECT_FALSE(converted->prefer_rows());
            EXPECT_FALSE(converted->sparse());
            tatami_test::test_simple_row_access(converted.get(), &smat, true, 1);
            tatami_test::test_simple_column_access(converted.get(), &smat, true, 1);
        }
    }
}


INSTANTIATE_TEST_CASE_P(
    ConvertToDense,
    ConvertToDenseTest,
    ::testing::Combine(
        ::testing::Values(10, 50, 100), // number of rows
        ::testing::Values(10, 50, 100)  // number of columns
    )
);
