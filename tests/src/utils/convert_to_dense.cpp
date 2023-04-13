#include <gtest/gtest.h>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/CompressedSparseMatrix.hpp"
#include "tatami/utils/convert_to_dense.hpp"

#include "../_tests/simulate_vector.h"
#include "../_tests/test_row_access.h"
#include "../_tests/test_column_access.h"

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
    auto vec = simulate_dense_vector<double>(NR * NC);
    tatami::DenseMatrix<true, double, int> mat(NR, NC, vec);

    auto converted = tatami::convert_to_dense<true>(&mat);
    EXPECT_TRUE(converted->prefer_rows());
    EXPECT_FALSE(converted->sparse());
    test_simple_row_access(converted.get(), &mat);
    test_simple_column_access(converted.get(), &mat);

    auto converted2 = tatami::convert_to_dense<true, int, size_t>(&mat); // works for a different type.
    EXPECT_TRUE(converted2->prefer_rows());
    EXPECT_FALSE(converted2->sparse());

    auto wrk2 = converted2->dense_row_workspace();
    for (size_t i = 0; i < NR; ++i) {
        auto start = vec.begin() + i * NC;
        std::vector<int> expected2(start, start + NC);
        EXPECT_EQ(converted2->row(i, wrk2.get()), expected2);
    }
}

TEST_P(ConvertToDenseTest, ColumnToColumn) {
    assemble(GetParam());
    auto vec = simulate_dense_vector<double>(NR * NC);
    tatami::DenseMatrix<false, double, int> mat(NR, NC, vec);

    auto converted = tatami::convert_to_dense<false>(&mat);
    EXPECT_FALSE(converted->prefer_rows());
    EXPECT_FALSE(converted->sparse());
    test_simple_row_access(converted.get(), &mat);
    test_simple_column_access(converted.get(), &mat);

    auto converted2 = tatami::convert_to_dense<false, int, size_t>(&mat); // works for a different type.
    EXPECT_FALSE(converted2->prefer_rows());
    EXPECT_FALSE(converted2->sparse());

    auto wrk2 = converted2->dense_column_workspace();
    for (size_t i = 0; i < NC; ++i) {
        auto start = vec.begin() + i * NR;
        std::vector<int> expected2(start, start + NR);
        EXPECT_EQ(converted2->column(i, wrk2.get()), expected2);
    }
}

TEST_P(ConvertToDenseTest, RowToColumn) {
    assemble(GetParam());
    auto vec = simulate_dense_vector<double>(NR * NC);
    tatami::DenseMatrix<true, double, int> mat(NR, NC, vec);

    auto converted = tatami::convert_to_dense<false>(&mat);
    EXPECT_FALSE(converted->prefer_rows());
    EXPECT_FALSE(converted->sparse());
    test_simple_row_access(converted.get(), &mat);
    test_simple_column_access(converted.get(), &mat);

    auto converted2 = tatami::convert_to_dense<false, int, size_t>(&mat); // works for a different type.
    EXPECT_FALSE(converted2->prefer_rows());
    EXPECT_FALSE(converted2->sparse());

    auto wrk2 = converted2->dense_row_workspace();
    for (size_t i = 0; i < NR; ++i) {
        auto start = vec.begin() + i * NC;
        std::vector<int> expected2(start, start + NC);
        EXPECT_EQ(converted2->row(i, wrk2.get()), expected2);
    }
}

TEST_P(ConvertToDenseTest, ColumnToRow) {
    assemble(GetParam());
    auto vec = simulate_dense_vector<double>(NR * NC);
    tatami::DenseMatrix<false, double, int> mat(NR, NC, vec);

    auto converted = tatami::convert_to_dense<true>(&mat);
    EXPECT_TRUE(converted->prefer_rows());
    EXPECT_FALSE(converted->sparse());
    test_simple_row_access(converted.get(), &mat);
    test_simple_column_access(converted.get(), &mat);

    auto converted2 = tatami::convert_to_dense<true, int, size_t>(&mat); // works for a different type.
    EXPECT_TRUE(converted2->prefer_rows());
    EXPECT_FALSE(converted2->sparse());

    auto wrk2 = converted2->dense_column_workspace();
    for (size_t i = 0; i < NC; ++i) {
        auto start = vec.begin() + i * NR;
        std::vector<int> expected2(start, start + NR);
        EXPECT_EQ(converted2->column(i, wrk2.get()), expected2);
    }
}

TEST_P(ConvertToDenseTest, Automatic) {
    assemble(GetParam());
    auto vec = simulate_dense_vector<double>(NR * NC);

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
        auto vec = simulate_sparse_compressed<double>(NR, NC, 0.2);
        tatami::CompressedSparseRowMatrix<double, int> smat(NR, NC, std::move(vec.value), std::move(vec.index), std::move(vec.ptr));

        {
            auto converted = tatami::convert_to_dense<true>(&smat);
            EXPECT_TRUE(converted->prefer_rows());
            EXPECT_FALSE(converted->sparse());
            test_simple_row_access(converted.get(), &smat);
            test_simple_column_access(converted.get(), &smat);
        }

        {
            auto converted = tatami::convert_to_dense<false>(&smat);
            EXPECT_FALSE(converted->prefer_rows());
            EXPECT_FALSE(converted->sparse());
            test_simple_row_access(converted.get(), &smat);
            test_simple_column_access(converted.get(), &smat);
        }
    }

    // From a row-major sparse matrix.
    {
        auto vec = simulate_sparse_compressed<double>(NC, NR, 0.2);
        tatami::CompressedSparseColumnMatrix<double, int> smat(NR, NC, std::move(vec.value), std::move(vec.index), std::move(vec.ptr));

        {
            auto converted = tatami::convert_to_dense<true>(&smat);
            EXPECT_TRUE(converted->prefer_rows());
            EXPECT_FALSE(converted->sparse());
            test_simple_row_access(converted.get(), &smat);
            test_simple_column_access(converted.get(), &smat);
        }

        {
            auto converted = tatami::convert_to_dense<false>(&smat);
            EXPECT_FALSE(converted->prefer_rows());
            EXPECT_FALSE(converted->sparse());
            test_simple_row_access(converted.get(), &smat);
            test_simple_column_access(converted.get(), &smat);
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
