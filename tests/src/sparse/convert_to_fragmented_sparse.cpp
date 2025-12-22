#include <gtest/gtest.h>
#include "../custom_parallel.h"

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/convert_to_fragmented_sparse.hpp"
#include "tatami/sparse/CompressedSparseMatrix.hpp"

#include "tatami_test/tatami_test.hpp"

class ConvertToFragmentedSparseTest : public ::testing::TestWithParam<std::tuple<int, int, bool, bool, int> > {
protected:
    int NR, NC;
    bool from_row, to_row;
    int nthreads;

    template<class Param>
    void assemble(const Param& param) {
        NR = std::get<0>(param);
        NC = std::get<1>(param);
        from_row = std::get<2>(param);
        to_row = std::get<3>(param);
        nthreads = std::get<4>(param);
    }
};

TEST_P(ConvertToFragmentedSparseTest, FromDense) {
    assemble(GetParam());

    auto vec = tatami_test::simulate_vector<double>(NR * NC, [&]{
        tatami_test::SimulateVectorOptions opt;
        opt.seed = NR * 10 + NC + static_cast<decltype(opt.seed)>(from_row) * 17 + static_cast<decltype(opt.seed)>(to_row) * 3 + nthreads;
        opt.density = 0.1;
        opt.seed = 786823;
        return opt;
    }());

    tatami::DenseMatrix<double, int, decltype(vec)> mat(NR, NC, std::move(vec), from_row);
    auto converted = tatami::convert_to_fragmented_sparse<double, int>(&mat, to_row, nthreads);

    EXPECT_EQ(converted->nrow(), NR);
    EXPECT_EQ(converted->ncol(), NC);
    EXPECT_TRUE(converted->is_sparse());
    EXPECT_EQ(converted->prefer_rows(), to_row);
    tatami_test::test_simple_row_access(*converted, mat);
    tatami_test::test_simple_column_access(*converted, mat);

    auto converted2 = tatami::convert_to_fragmented_sparse<int, std::size_t>(&mat, to_row, nthreads); // works for a different type.
    EXPECT_TRUE(converted2->is_sparse());
    EXPECT_EQ(converted2->prefer_rows(), to_row);

    auto old = mat.dense_row();
    auto conv = converted2->dense_row();
    std::vector<double> obuffer(NC);
    std::vector<int> obuffer_i, cbuffer_i(NC);
    for (int i = 0; i < NR; ++i) {
        auto optr = old->fetch(i, obuffer.data());
        std::copy_n(optr, NC, obuffer_i.data());
        auto cptr = conv->fetch(i, cbuffer_i.data());
        tatami::copy_n(cptr, NC, cbuffer_i.data());
        EXPECT_EQ(obuffer_i, cbuffer_i);
    }
}

TEST_P(ConvertToFragmentedSparseTest, ColumnToColumn) {
    assemble(GetParam());

    auto trip = tatami_test::simulate_compressed_sparse<double, int>((from_row ? NR : NC), (from_row ? NC : NR), [&]{
        tatami_test::SimulateCompressedSparseOptions opt;
        opt.seed = NR * 10 + NC + static_cast<decltype(opt.seed)>(from_row) * 17 + static_cast<decltype(opt.seed)>(to_row) * 3 + nthreads;
        opt.density = 0.1;
        opt.seed = 123847;
        return opt;
    }());

    tatami::CompressedSparseMatrix<
        double,
        int,
        decltype(decltype(trip)::data),
        decltype(decltype(trip)::index),
        decltype(decltype(trip)::indptr)
    > spmat(
        NR,
        NC,
        std::move(trip.data),
        std::move(trip.index),
        std::move(trip.indptr),
        from_row
    );
    auto converted = tatami::convert_to_fragmented_sparse<double, int>(&spmat, to_row, nthreads);

    EXPECT_EQ(converted->nrow(), NR);
    EXPECT_EQ(converted->ncol(), NC);
    EXPECT_TRUE(converted->is_sparse());
    EXPECT_EQ(converted->prefer_rows(), to_row);
    tatami_test::test_simple_row_access(*converted, spmat);
    tatami_test::test_simple_column_access(*converted, spmat);

    auto converted2 = tatami::convert_to_fragmented_sparse<int, std::size_t>(&spmat, to_row, nthreads); // works for a different type.
    EXPECT_TRUE(converted2->is_sparse());
    EXPECT_EQ(converted2->prefer_rows(), to_row);

    auto wrk = spmat.dense_column();
    auto wrk2 = converted2->dense_column();
    std::vector<double> buffer(NR);
    std::vector<int> buffer_i(NR), buffer2_i(NR);
    for (int i = 0; i < NC; ++i) {
        auto ptr = wrk->fetch(i, buffer.data());
        std::copy_n(ptr, NR, buffer_i.data());
        auto ptr2 = wrk2->fetch(i, buffer2_i.data());
        tatami::copy_n(ptr2, NR, buffer2_i.data());
        EXPECT_EQ(buffer_i, buffer2_i);
    }
}

INSTANTIATE_TEST_SUITE_P(
    ConvertToFragmentedSparse,
    ConvertToFragmentedSparseTest,
    ::testing::Combine(
        ::testing::Values(10, 50, 100), // number of rows
        ::testing::Values(10, 50, 100), // number of columns
        ::testing::Values(true, false), // from row major?
        ::testing::Values(true, false), // to row major?
        ::testing::Values(1, 3)         // number of threads
    )
);
