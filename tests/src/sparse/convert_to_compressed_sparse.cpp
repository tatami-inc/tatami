#include <gtest/gtest.h>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

class ConvertToCompressedSparseTest : public ::testing::TestWithParam<std::tuple<int, int, bool, bool, bool, int> > {
protected:
    size_t NR, NC;
    bool from_row, to_row;
    bool two_pass;
    int nthreads;

    template<class Param>
    void assemble(const Param& param) {
        NR = std::get<0>(param);
        NC = std::get<1>(param);
        from_row = std::get<2>(param);
        to_row = std::get<3>(param);
        two_pass = std::get<4>(param);
        nthreads = std::get<5>(param);
    }
};

TEST_P(ConvertToCompressedSparseTest, FromDense) {
    assemble(GetParam());
    auto vec = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1);
    auto mat = std::make_shared<tatami::DenseMatrix<double, int> >(NR, NC, vec, from_row);

    auto converted = tatami::convert_to_compressed_sparse<double, int>(mat.get(), to_row, two_pass, nthreads);
    EXPECT_EQ(converted->nrow(), NR);
    EXPECT_EQ(converted->ncol(), NC);
    EXPECT_TRUE(converted->sparse());
    EXPECT_EQ(converted->prefer_rows(), to_row);
    tatami_test::test_simple_row_access(converted.get(), mat.get());
    tatami_test::test_simple_column_access(converted.get(), mat.get());

    auto converted2 = tatami::convert_to_compressed_sparse<int, size_t>(mat.get(), to_row, two_pass, nthreads); // works for a different type.
    EXPECT_TRUE(converted2->sparse());
    EXPECT_EQ(converted2->prefer_rows(), to_row);

    auto wrk2 = converted2->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto start = vec.begin() + i * NC;
        std::vector<int> expected2(start, start + NC);
        EXPECT_EQ(tatami_test::fetch(wrk2.get(), i, NC), expected2);
    }
}

TEST_P(ConvertToCompressedSparseTest, FromSparse) {
    assemble(GetParam());
    auto trip = tatami_test::simulate_sparse_compressed<double>(NC, NR, 0.15); 
    auto mat = std::make_shared<tatami::CompressedSparseMatrix<double, int> >(NR, NC, trip.value, trip.index, trip.ptr, from_row);

    auto converted = tatami::convert_to_compressed_sparse<double, int>(mat.get(), to_row, two_pass, nthreads);
    EXPECT_EQ(converted->nrow(), NR);
    EXPECT_EQ(converted->ncol(), NC);
    EXPECT_TRUE(converted->sparse());
    EXPECT_EQ(converted->prefer_rows(), to_row);
    tatami_test::test_simple_row_access(converted.get(), mat.get());
    tatami_test::test_simple_column_access(converted.get(), mat.get());

    auto converted2 = tatami::convert_to_compressed_sparse<int, size_t>(mat.get(), to_row, two_pass, nthreads); // works for a different type.
    EXPECT_TRUE(converted2->sparse());
    EXPECT_EQ(converted2->prefer_rows(), to_row);

    auto wrk = mat->dense_column();
    auto wrk2 = converted2->dense_column();
    for (size_t i = 0; i < NC; ++i) {
        auto expected = tatami_test::fetch(wrk.get(), static_cast<int>(i), NR);
        std::vector<int> expected2(expected.begin(), expected.end());
        EXPECT_EQ(tatami_test::fetch(wrk2.get(), i, NR), expected2);
    }
}

INSTANTIATE_TEST_SUITE_P(
    ConvertToCompressedSparse,
    ConvertToCompressedSparseTest,
    ::testing::Combine(
        ::testing::Values(10, 50, 100), // number of rows
        ::testing::Values(10, 50, 100), // number of columns
        ::testing::Values(true, false), // from row major?
        ::testing::Values(true, false), // to row major?
        ::testing::Values(true, false), // two-pass? 
        ::testing::Values(1, 3)         // number of threads
    )
);

TEST(ConvertToCompressedSparseTest, Automatic) {
    size_t NR = 70, NC = 50;
    auto vec = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.23);

    {
        tatami::DenseMatrix<double, int> mat(NR, NC, vec, false);
        auto converted = tatami::convert_to_compressed_sparse(&mat, -1);
        EXPECT_FALSE(converted->prefer_rows());
    }

    {
        tatami::DenseMatrix<double, int> mat(NR, NC, vec, true);
        auto converted = tatami::convert_to_compressed_sparse(&mat, -1);
        EXPECT_TRUE(converted->prefer_rows());
    }
}
