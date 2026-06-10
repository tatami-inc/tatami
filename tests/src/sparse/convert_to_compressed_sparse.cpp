#include <gtest/gtest.h>

#ifdef CONVERT_COMPRESSED_SPARSE_CUSTOM_PARALLEL_TEST

#include <thread>
#include <vector>

#include "sanisizer/sanisizer.hpp"

// convert_to_compressed_sparse() has an unusual work-sharing arrangement when filling the output arrays.
// We need to check that we are not making assumptions about the first range being assigned to thread 0,
// This code is largely copied from subpar::parallelize_range but the first worker gets the last range.
template<class Function_, typename Task_>
int weird_parallelize(Function_ run_task_range, const Task_ num_tasks, const int num_workers) {
    if (num_tasks <= 0) {
        return 0;
    }

    if (num_workers <= 1 || num_tasks == 1) {
        run_task_range(0, 0, num_tasks);
        return 1;
    }

    // All workers with indices below 'remainder' get an extra task to fill up the remainder.
    Task_ tasks_per_worker = 1;
    int remainder = 0;
    if (sanisizer::is_greater_than_or_equal(num_workers, num_tasks)) {
        num_workers = num_tasks;
    } else {
        tasks_per_worker = num_tasks / num_workers;
        remainder = num_tasks % num_workers;
    }

    std::vector<std::thread> workers;
    sanisizer::reserve(workers, num_workers); // preallocate to ensure we don't get alloc errors during emplace_back().

    // Worker 0 has the last range, etc., and the last worker gets the first range.
    for (int w = 0; w < num_workers; ++w) {
        const auto range_id = num_workers - w - 1;
        const Task_ start = range_id * tasks_per_worker + (range_id < remainder ? range_id : remainder);
        const Task_ length = tasks_per_worker + (range_id < remainder); 
        workers.emplace_back(run_task_range, w, start, length);
    }

    for (auto& wrk : workers) {
        wrk.join();
    }

    return num_workers;
}

// Overriding the default tatami parallelization.
#define TATAMI_CUSTOM_PARALLELIZE(fun, tasks, workers) weird_parallelize(fun, tasks, workers)

#else
#include "../custom_parallel.h"
#endif

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

class ConvertToCompressedSparseTest : public ::testing::TestWithParam<std::tuple<int, int, bool, bool, bool, int> > {
protected:
    int NR, NC;
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
        opt.seed = NR * 10 + NC + static_cast<decltype(opt.seed)>(from_row) * 7 + static_cast<decltype(opt.seed)>(to_row) * 13 + two_pass * nthreads;
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

    auto converted2 = tatami::convert_to_compressed_sparse<int, std::size_t>(&mat, to_row, two_pass, nthreads); // works for a different type.
    EXPECT_TRUE(converted2->is_sparse());
    EXPECT_EQ(converted2->prefer_rows(), to_row);

    auto old = mat.dense_row();
    auto conv = converted2->dense_row();
    std::vector<double> obuffer(NC);
    std::vector<int> obuffer_i(NC), cbuffer_i(NC);
    for (int i = 0; i < NR; ++i) {
        auto optr = old->fetch(i, obuffer.data());
        std::copy_n(optr, NC, obuffer_i.data());
        auto cptr = conv->fetch(i, cbuffer_i.data());
        tatami::copy_n(cptr, NC, cbuffer_i.data());
        EXPECT_EQ(obuffer_i, cbuffer_i);
    }
}

TEST_P(ConvertToCompressedSparseTest, FromSparse) {
    assemble(GetParam());

    auto trip = tatami_test::simulate_compressed_sparse<double, int>((from_row ? NR : NC), (from_row ? NC : NR), [&]{
        tatami_test::SimulateCompressedSparseOptions opt;
        opt.seed = NR * 10 + NC + static_cast<decltype(opt.seed)>(from_row) * 7 + static_cast<decltype(opt.seed)>(to_row) * 13 + two_pass * nthreads;
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

    auto converted2 = tatami::convert_to_compressed_sparse<int, std::size_t>(&mat, to_row, two_pass, nthreads); // works for a different type.
    EXPECT_TRUE(converted2->is_sparse());
    EXPECT_EQ(converted2->prefer_rows(), to_row);

    auto wrk = mat.dense_column();
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
    inline static int NR = 120, NC = 50;
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
    std::vector<std::size_t> pointers(NR + 1);
    tatami::count_compressed_sparse_non_zeros(mat.get(), true, pointers.data() + 1, 1);

    {
        auto wrk = mat->dense_row();
        std::vector<double> buffer(NC);
        for (int i = 0; i < NR; ++i) {
            auto ptr = wrk->fetch(i, buffer.data());
            std::size_t expected = 0;
            for (int c = 0; c < NC; ++c) {
                expected += (ptr[c] != 0);
            }
            EXPECT_EQ(expected, pointers[i + 1]);
        }
    }

    for (int i = 1; i <= NR; ++i) {
        pointers[i] += pointers[i - 1];
    }
    const auto nonzeros = pointers.back();
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
        std::vector<double> buffer(NC), spbuffer(NC);
        for (int i = 0; i < NR; ++i) {
            auto ptr = wrk->fetch(i, buffer.data());
            tatami::copy_n(ptr, NC, buffer.data());
            auto spptr = spwrk->fetch(i, spbuffer.data());
            tatami::copy_n(spptr, NC, spbuffer.data());
            EXPECT_EQ(buffer, spbuffer);
        }
    }
}

TEST_F(ConvertToCompressedSparseManualTest, Inconsistent) {
    std::vector<std::size_t> pointers(NC + 1);
    tatami::count_compressed_sparse_non_zeros(mat.get(), false, pointers.data() + 1, 1);

    {
        auto wrk = mat->dense_column();
        std::vector<double> buffer(NR);
        for (int i = 0; i < NC; ++i) {
            auto ptr = wrk->fetch(i, buffer.data());
            std::size_t expected = 0;
            for (int r = 0; r < NR; ++r) {
                expected += (ptr[r] != 0);
            }
            EXPECT_EQ(expected, pointers[i + 1]);
        }
    }

    for (int i = 1; i <= NC; ++i) {
        pointers[i] += pointers[i - 1];
    }
    const auto nonzeros = pointers.back();
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
        std::vector<double> buffer(NR), spbuffer(NR);
        for (int i = 0; i < NC; ++i) {
            auto ptr = wrk->fetch(i, buffer.data());
            tatami::copy_n(ptr, NR, buffer.data());
            auto spptr = spwrk->fetch(i, spbuffer.data());
            tatami::copy_n(spptr, NR, spbuffer.data());
            EXPECT_EQ(buffer, spbuffer);
        }
    }
}

class ConvertToCompressedSparseEmptyTest : public ::testing::TestWithParam<std::tuple<std::pair<int, int>, bool, bool, bool> > {};

TEST_P(ConvertToCompressedSparseEmptyTest, Empty) {
    auto params = GetParam();
    auto dims = std::get<0>(params);
    bool row_major_source = std::get<1>(params);
    bool row_major_target = std::get<2>(params);
    bool row_major_target2 = std::get<3>(params);

    tatami::DenseMatrix<double, int, std::vector<double> > mat(dims.first, dims.second, std::vector<double>(), row_major_source);
    auto spmat = tatami::convert_to_compressed_sparse<double, int>(mat, row_major_target, {});
    EXPECT_EQ(spmat->nrow(), dims.first);
    EXPECT_EQ(spmat->ncol(), dims.second);

    auto spmat2 = tatami::convert_to_compressed_sparse<double, int>(mat, row_major_target2, {});
    EXPECT_EQ(spmat2->nrow(), dims.first);
    EXPECT_EQ(spmat2->ncol(), dims.second);
}

INSTANTIATE_TEST_SUITE_P(
    ConvertToCompressedSparse,
    ConvertToCompressedSparseEmptyTest,
    ::testing::Combine(
        ::testing::Values(
            std::make_pair(10, 0),
            std::make_pair(0, 10)
        ),
        ::testing::Values(true, false), // whether the input dense matrix is row-major.
        ::testing::Values(true, false), // whether the first converted sparse matrix is row-major.
        ::testing::Values(true, false) // whether the second re-converted sparse matrix is row-major.
    )
);
