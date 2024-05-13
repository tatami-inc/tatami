#include <gtest/gtest.h>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/CompressedSparseMatrix.hpp"
#include "tatami/dense/convert_to_dense.hpp"

#include "tatami_test/tatami_test.hpp"

class ConvertToDenseTest : public ::testing::TestWithParam<std::tuple<int, int, bool, int> > {
protected:
    size_t NR, NC;
    bool row;
    int threads;

    template<class Param>
    void assemble(const Param& param) {
        NR = std::get<0>(param);
        NC = std::get<1>(param);
        row = std::get<2>(param);
        threads = std::get<3>(param);
    }
};

TEST_P(ConvertToDenseTest, FromDenseRowMajor) {
    assemble(GetParam());
    auto vec = tatami_test::simulate_dense_vector<double>(NR * NC);
    auto mat = std::make_shared<tatami::DenseMatrix<double, int> >(NR, NC, vec, true);

    auto converted = tatami::convert_to_dense(mat.get(), row, threads);
    EXPECT_EQ(converted->prefer_rows(), row);
    EXPECT_FALSE(converted->sparse());

    tatami_test::test_simple_row_access(converted.get(), mat.get());
    tatami_test::test_simple_column_access(converted.get(), mat.get());

    auto converted2 = tatami::convert_to_dense<int, size_t>(mat.get(), row, threads); // works for a different type.
    EXPECT_EQ(converted2->prefer_rows(), row);
    EXPECT_FALSE(converted2->sparse());

    auto wrk2 = converted2->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto start = vec.begin() + i * NC;
        std::vector<int> expected2(start, start + NC);
        EXPECT_EQ(tatami_test::fetch(wrk2.get(), i, NC), expected2);
    }
}

TEST_P(ConvertToDenseTest, FromColumnRowMajor) {
    assemble(GetParam());
    auto vec = tatami_test::simulate_dense_vector<double>(NR * NC);
    auto mat = std::make_shared<tatami::DenseMatrix<double, int> >(NR, NC, vec, false);

    auto converted = tatami::convert_to_dense(mat.get(), row, threads);
    EXPECT_EQ(converted->prefer_rows(), row);
    EXPECT_FALSE(converted->sparse());

    tatami_test::test_simple_row_access(converted.get(), mat.get());
    tatami_test::test_simple_column_access(converted.get(), mat.get());

    auto converted2 = tatami::convert_to_dense<int, size_t>(mat.get(), row, threads); // works for a different type.
    EXPECT_EQ(converted2->prefer_rows(), row);
    EXPECT_FALSE(converted2->sparse());

    auto wrk2 = converted2->dense_column();
    for (size_t i = 0; i < NC; ++i) {
        auto start = vec.begin() + i * NR;
        std::vector<int> expected2(start, start + NR);
        EXPECT_EQ(tatami_test::fetch(wrk2.get(), i, NR), expected2);
    }
}

TEST_P(ConvertToDenseTest, Automatic) {
    assemble(GetParam());
    auto vec = tatami_test::simulate_dense_vector<double>(NR * NC);

    if (row) {
        tatami::DenseMatrix<double, int> mat(NR, NC, vec, true);
        auto converted = tatami::convert_to_dense(&mat, -1, threads);
        EXPECT_TRUE(converted->prefer_rows());
        tatami_test::test_simple_row_access(converted.get(), &mat);
    } else {
        tatami::DenseMatrix<double, int> mat(NR, NC, vec, false);
        auto converted = tatami::convert_to_dense(&mat, -1, threads);
        EXPECT_FALSE(converted->prefer_rows());
        tatami_test::test_simple_column_access(converted.get(), &mat);
    }
}

TEST_P(ConvertToDenseTest, FromSparseRowMajor) {
    assemble(GetParam());
    auto vec = tatami_test::simulate_sparse_compressed<double>(NR, NC, 0.2);
    tatami::CompressedSparseRowMatrix<double, int> smat(NR, NC, std::move(vec.value), std::move(vec.index), std::move(vec.ptr));

    auto converted = tatami::convert_to_dense(&smat, row, threads);
    EXPECT_EQ(converted->prefer_rows(), row);
    EXPECT_FALSE(converted->sparse());
    tatami_test::test_simple_row_access(converted.get(), &smat);
    tatami_test::test_simple_column_access(converted.get(), &smat);
}

TEST_P(ConvertToDenseTest, FromSparseColumnMajor) {
    auto vec = tatami_test::simulate_sparse_compressed<double>(NC, NR, 0.2);
    tatami::CompressedSparseColumnMatrix<double, int> smat(NR, NC, std::move(vec.value), std::move(vec.index), std::move(vec.ptr));

    auto converted = tatami::convert_to_dense(&smat, row, threads);
    EXPECT_EQ(converted->prefer_rows(), row);
    EXPECT_FALSE(converted->sparse());
    tatami_test::test_simple_row_access(converted.get(), &smat);
    tatami_test::test_simple_column_access(converted.get(), &smat);
}

INSTANTIATE_TEST_SUITE_P(
    ConvertToDense,
    ConvertToDenseTest,
    ::testing::Combine(
        ::testing::Values(10, 50, 100), // number of rows
        ::testing::Values(10, 50, 100), // number of columns
        ::testing::Values(true, false), // row major?
        ::testing::Values(1, 3)         // number of threads
    )
);
