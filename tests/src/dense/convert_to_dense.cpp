#include <gtest/gtest.h>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/CompressedSparseMatrix.hpp"
#include "tatami/dense/convert_to_dense.hpp"

#include "tatami_test/tatami_test.hpp"

class ConvertToDenseTest : public ::testing::TestWithParam<std::tuple<int, int, bool, bool, int> > {
protected:
    size_t NR, NC;
    bool from_row, to_row;
    int threads;

    template<class Param>
    void assemble(const Param& param) {
        NR = std::get<0>(param);
        NC = std::get<1>(param);
        from_row = std::get<2>(param);
        to_row = std::get<3>(param);
        threads = std::get<4>(param);
    }
};

TEST_P(ConvertToDenseTest, FromDense) {
    assemble(GetParam());
    auto vec = tatami_test::simulate_dense_vector<double>(NR * NC);
    auto mat = std::make_shared<tatami::DenseMatrix<double, int> >(NR, NC, vec, from_row);

    auto converted = tatami::convert_to_dense(mat.get(), to_row, threads);
    EXPECT_EQ(converted->prefer_rows(), to_row);
    EXPECT_FALSE(converted->sparse());

    tatami_test::test_simple_row_access(converted.get(), mat.get());
    tatami_test::test_simple_column_access(converted.get(), mat.get());

    auto converted2 = tatami::convert_to_dense<int, size_t>(mat.get(), to_row, threads); // works for a different type.
    EXPECT_EQ(converted2->prefer_rows(), to_row);
    EXPECT_FALSE(converted2->sparse());

    auto wrk2 = converted2->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto start = vec.begin() + i * NC;
        std::vector<int> expected2(start, start + NC);
        EXPECT_EQ(tatami_test::fetch(wrk2.get(), i, NC), expected2);
    }
}

TEST_P(ConvertToDenseTest, FromSparse) {
    assemble(GetParam());
    auto vec = tatami_test::simulate_sparse_compressed<double>(NR, NC, 0.2);
    tatami::CompressedSparseMatrix<double, int> smat(NR, NC, std::move(vec.value), std::move(vec.index), std::move(vec.ptr), from_row);

    auto converted = tatami::convert_to_dense(&smat, to_row, threads);
    EXPECT_EQ(converted->prefer_rows(), to_row);
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
        ::testing::Values(true, false), // from row major?
        ::testing::Values(true, false), // to row major?
        ::testing::Values(1, 3)         // number of threads
    )
);

TEST(ConvertToDenseTest, Automatic) {
    size_t NR = 50, NC = 20;
    auto vec = tatami_test::simulate_dense_vector<double>(NR * NC);

    {
        tatami::DenseMatrix<double, int> mat(NR, NC, vec, true);
        auto converted = tatami::convert_to_dense(&mat, -1);
        EXPECT_TRUE(converted->prefer_rows());
        tatami_test::test_simple_row_access(converted.get(), &mat);
    }
    
    {
        tatami::DenseMatrix<double, int> mat(NR, NC, vec, false);
        auto converted = tatami::convert_to_dense(&mat, -1);
        EXPECT_FALSE(converted->prefer_rows());
        tatami_test::test_simple_column_access(converted.get(), &mat);
    }
}
