#include <gtest/gtest.h>
#include "../custom_parallel.h"

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

    auto vec = tatami_test::simulate_vector<double>(NR * NC, [&]{
        tatami_test::SimulateVectorOptions opt;
        opt.seed = NR * 10 + NC + static_cast<size_t>(from_row) * 7 + static_cast<size_t>(to_row) * 13 + two_pass * nthreads;
        opt.density = 0.1;
        opt.seed = 23093469;
        return opt;
    }());

    tatami::DenseMatrix<double, int, decltype(vec)> mat(NR, NC, std::move(vec), from_row);
    auto converted = tatami::convert_to_compressed_sparse<double, int>(&mat, to_row, two_pass, nthreads);

    EXPECT_EQ(converted->nrow(), NR);
    EXPECT_EQ(converted->ncol(), NC);
    EXPECT_TRUE(converted->is_sparse());
    EXPECT_EQ(converted->prefer_rows(), to_row);
    tatami_test::test_simple_row_access(*converted, mat);
    tatami_test::test_simple_column_access(*converted.get(), mat);

    auto converted2 = tatami::convert_to_compressed_sparse<int, size_t>(&mat, to_row, two_pass, nthreads); // works for a different type.
    EXPECT_TRUE(converted2->is_sparse());
    EXPECT_EQ(converted2->prefer_rows(), to_row);

    auto old = mat.dense_row();
    std::vector<double> buffer(NC);
    auto wrk2 = converted2->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto ptr = old->fetch(i, buffer.data());
        std::vector<int> expected(ptr, ptr + NC);
        EXPECT_EQ(tatami_test::fetch(*wrk2, i, NC), expected);
    }
}

TEST_P(ConvertToCompressedSparseTest, FromSparse) {
    assemble(GetParam());

    auto trip = tatami_test::simulate_compressed_sparse<double, int>((from_row ? NR : NC), (from_row ? NC : NR), [&]{
        tatami_test::SimulateCompressedSparseOptions opt;
        opt.seed = NR * 10 + NC + static_cast<size_t>(from_row) * 7 + static_cast<size_t>(to_row) * 13 + two_pass * nthreads;
        opt.density = 0.15;
        opt.seed = 3890793;
        return opt;
    }());

    tatami::CompressedSparseMatrix<
        double,
        int,
        decltype(decltype(trip)::data),
        decltype(decltype(trip)::index),
        decltype(decltype(trip)::indptr)
    > mat(
        NR,
        NC,
        std::move(trip.data),
        std::move(trip.index),
        std::move(trip.indptr),
        from_row
    );
    auto converted = tatami::convert_to_compressed_sparse<double, int>(&mat, to_row, two_pass, nthreads);

    EXPECT_EQ(converted->nrow(), NR);
    EXPECT_EQ(converted->ncol(), NC);
    EXPECT_TRUE(converted->is_sparse());
    EXPECT_EQ(converted->prefer_rows(), to_row);
    tatami_test::test_simple_row_access(*converted, mat);
    tatami_test::test_simple_column_access(*converted, mat);

    auto converted2 = tatami::convert_to_compressed_sparse<int, size_t>(&mat, to_row, two_pass, nthreads); // works for a different type.
    EXPECT_TRUE(converted2->is_sparse());
    EXPECT_EQ(converted2->prefer_rows(), to_row);

    auto wrk = mat.dense_column();
    auto wrk2 = converted2->dense_column();
    for (size_t i = 0; i < NC; ++i) {
        auto expected = tatami_test::fetch(*wrk, static_cast<int>(i), NR);
        std::vector<int> expected2(expected.begin(), expected.end());
        EXPECT_EQ(tatami_test::fetch(*wrk2, i, NR), expected2);
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

class ConvertToCompressedSparseManualTest : public ::testing::Test {
protected:
    inline static size_t NR = 120, NC = 50;
    inline static std::shared_ptr<tatami::Matrix<double, int> > mat;

    static void SetUpTestSuite() {
        auto vec = tatami_test::simulate_vector<double>(NR * NC, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.1;
            opt.seed = 92810823;
            return opt;
        }());
        mat.reset(new tatami::DenseRowMatrix<double, int>(NR, NC, std::move(vec)));
    }
};

TEST_F(ConvertToCompressedSparseManualTest, Consistent) {
    std::vector<size_t> pointers(NR + 1);
    tatami::count_compressed_sparse_non_zeros(mat.get(), true, pointers.data() + 1, 1);

    {
        auto wrk = mat->dense_row();
        for (size_t i = 0; i < NR; ++i) {
            auto row = tatami_test::fetch(*wrk, static_cast<int>(i), NC);
            size_t expected = 0;
            for (auto x : row) {
                expected += (x != 0);
            }
            EXPECT_EQ(expected, pointers[i + 1]);
        }
    }

    for (size_t i = 1; i <= NR; ++i) {
        pointers[i] += pointers[i - 1];
    }
    size_t nonzeros = pointers.back();
    std::vector<double> values(nonzeros);
    std::vector<int> indices(nonzeros);
    tatami::fill_compressed_sparse_contents(mat.get(), true, pointers.data(), values.data(), indices.data(), 1);

    tatami::CompressedSparseMatrix<
        double,
        int,
        decltype(values),
        decltype(indices),
        decltype(pointers)
    > spmat(
        mat->nrow(), 
        mat->ncol(), 
        std::move(values), 
        std::move(indices), 
        std::move(pointers),
        true 
    );

    {
        auto wrk = mat->dense_row();
        auto spwrk = spmat.dense_row();
        for (size_t i = 0; i < NR; ++i) {
            auto expected = tatami_test::fetch(*wrk, static_cast<int>(i), NC);
            auto observed = tatami_test::fetch(*spwrk, static_cast<int>(i), NC);
            EXPECT_EQ(expected, observed);
        }
    }
}

TEST_F(ConvertToCompressedSparseManualTest, Inconsistent) {
    std::vector<size_t> pointers(NC + 1);
    tatami::count_compressed_sparse_non_zeros(mat.get(), false, pointers.data() + 1, 1);

    {
        auto wrk = mat->dense_column();
        for (size_t i = 0; i < NC; ++i) {
            auto column = tatami_test::fetch(*wrk, static_cast<int>(i), NR);
            size_t expected = 0;
            for (auto x : column) {
                expected += (x != 0);
            }
            EXPECT_EQ(expected, pointers[i + 1]);
        }
    }

    for (size_t i = 1; i <= NC; ++i) {
        pointers[i] += pointers[i - 1];
    }
    size_t nonzeros = pointers.back();
    std::vector<double> values(nonzeros);
    std::vector<int> indices(nonzeros);
    tatami::fill_compressed_sparse_contents(mat.get(), false, pointers.data(), values.data(), indices.data(), 1);

    tatami::CompressedSparseMatrix<
        double,
        int,
        decltype(values),
        decltype(indices),
        decltype(pointers)
    > spmat( 
        mat->nrow(), 
        mat->ncol(), 
        std::move(values), 
        std::move(indices), 
        std::move(pointers),
        false 
    );

    {
        auto wrk = mat->dense_column();
        auto spwrk = spmat.dense_column();
        for (size_t i = 0; i < NC; ++i) {
            auto expected= tatami_test::fetch(*wrk, static_cast<int>(i), NR);
            auto observed = tatami_test::fetch(*spwrk, static_cast<int>(i), NR);
            EXPECT_EQ(expected, observed);
        }
    }
}
