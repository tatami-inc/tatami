#include <gtest/gtest.h>

#include "tatami/base/dense/DenseMatrix.hpp"
#include "tatami/base/sparse/CompressedSparseMatrix.hpp"
#include "tatami/ext/layered.hpp"

#include "mock_layered_sparse_data.h"

typedef std::vector<int> IntVec;

class ConvertToLayeredSparseTest : public ::testing::TestWithParam<std::tuple<int, int, IntVec, IntVec, IntVec> > {
protected:
    size_t NR, NC;
    IntVec rows, cols, vals;

    template<class PARAM>
    void dump(PARAM param) {
        NR = std::get<0>(param);
        NC = std::get<1>(param);
        rows = std::get<2>(param);
        cols = std::get<3>(param);
        vals = std::get<4>(param);
    }
};

TEST_P(ConvertToLayeredSparseTest, FromCSC) {
    dump(GetParam());

    // Checking against a reference.
    auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
    typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(vals), decltype(rows), decltype(indptrs)> SparseMat; 
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new SparseMat(NR, NC, std::move(vals), std::move(rows), std::move(indptrs))); 

    auto out = tatami::convert_to_layered_sparse(ref.get());

    auto rwrk = ref->dense_row();
    auto owrk = out.matrix->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto stuff = owrk->fetch(out.permutation[i]);
        EXPECT_EQ(stuff, rwrk->fetch(i));
    }
}

TEST_P(ConvertToLayeredSparseTest, FromCSR) {
    dump(GetParam());

    auto indptrs = tatami::compress_sparse_triplets<true>(NR, NC, vals, rows, cols);
    typedef tatami::CompressedSparseRowMatrix<double, int, decltype(vals), decltype(cols), decltype(indptrs)> SparseMat; 
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new SparseMat(NR, NC, std::move(vals), std::move(cols), std::move(indptrs)));

    auto out = tatami::convert_to_layered_sparse(ref.get());

    auto rwrk = ref->dense_row();
    auto owrk = out.matrix->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto stuff = owrk->fetch(out.permutation[i]);
        EXPECT_EQ(stuff, rwrk->fetch(i));
    }
}

TEST_P(ConvertToLayeredSparseTest, FromDenseColumn) {
    dump(GetParam());

    std::vector<int> full(NR * NC);
    for (size_t i = 0; i < vals.size(); ++i) {
        full[rows[i] + cols[i] * NR] = vals[i];
    }
    typedef tatami::DenseColumnMatrix<double, int, decltype(full)> DenseMat; 
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new DenseMat(NR, NC, std::move(full)));

    auto out = tatami::convert_to_layered_sparse(ref.get());

    auto rwrk = ref->dense_row();
    auto owrk = out.matrix->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto stuff = owrk->fetch(out.permutation[i]);
        EXPECT_EQ(stuff, rwrk->fetch(i));
    }
}

TEST_P(ConvertToLayeredSparseTest, FromDenseRow) {
    dump(GetParam());

    std::vector<int> full(NR * NC);
    for (size_t i = 0; i < vals.size(); ++i) {
        full[rows[i] * NC + cols[i]] = vals[i];
    }
    typedef tatami::DenseRowMatrix<double, int, decltype(full)> DenseMat; 
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new DenseMat(NR, NC, std::move(full)));

    auto out = tatami::convert_to_layered_sparse(ref.get());

    auto rwrk = ref->dense_row();
    auto owrk = out.matrix->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto stuff = owrk->fetch(out.permutation[i]);
        EXPECT_EQ(stuff, rwrk->fetch(i));
    }
}

INSTANTIATE_TEST_CASE_P(
    ConvertToLayeredSparse,
    ConvertToLayeredSparseTest,
    ::testing::Values(
        std::make_tuple(
            // this example guarantees a few rows in each chunk.
            10,
            5,
            IntVec{ 1, 5, 8, 2, 9, 0, 4 }, 
            IntVec{ 2, 3, 1, 0, 2, 2, 4 },
            IntVec{ 0, 1, 10, 100, 1000, 10000, 100000 }
        ),
        std::make_tuple(
            // this example checks that the maximum category is used for each row.
            5,
            8,
            IntVec{ 1, 1, 2, 2, 3, 3, 4, 4 },
            IntVec{ 1, 7, 2, 5, 4, 3, 4, 0 },
            IntVec{ 10, 1, 10, 1000, 10000, 100000, 1, 100000 }
        ),
        std::make_tuple(
            // this example checks that we handle missing categories; in this case, only uint8.
            10,
            9,
            IntVec{ 1, 9, 7, 5, 3, 1, 3, 3, 7, 9 },
            IntVec{ 2, 4, 8, 8, 4, 6, 8, 6, 0, 0 },
            IntVec{ 1, 3, 2, 1, 7, 8, 9, 1, 1, 3 }
        ),
        std::make_tuple(
            // this example checks that we handle missing categories; in this case, only uint16.
            20,
            15,
            IntVec{ 15, 0, 4, 14, 0, 19, 19, 8, 11, 18, 2,  3, 6,  4, 9,  3, 16,  4, 13, 12 },
            IntVec{  3, 3, 4,  2, 8, 12,  3, 6,  2,  3, 2, 11, 1, 11, 5, 12,  7, 12,  5,  0 },
            IntVec{ 1000, 3000, 2000, 1000, 7000, 10000, 9000, 1000, 600, 500, 382, 826, 992, 244, 138, 852, 400, 542, 980, 116 }
        ),
        std::make_tuple(
            // this example checks that we handle missing categories; in this case, only uint32.
            100,
            20,
            IntVec{ 27, 83, 85, 60, 17, 45, 62, 30, 98, 47 },
            IntVec{ 0, 12, 17, 3, 17, 0, 8, 3, 8, 8 },
            IntVec{ 130875, 673886, 405953, 989598, 981526, 794394, 680144, 553105, 277529, 540959 }
        ),
        std::make_tuple(
            // this example checks that we handle empties.
            10,
            9,
            IntVec{},
            IntVec{},
            IntVec{}
        )
    )
);

class ConvertToLayeredSparseHardTest : public ::testing::TestWithParam<int> {};

TEST_P(ConvertToLayeredSparseHardTest, Complex) {
    size_t NR = GetParam();
    size_t NC = 10;

    // Checking against a reference.
    {
        std::vector<size_t> rows, cols;
        std::vector<int> vals;
        mock_layered_sparse_data<false>(NR, NC, rows, cols, vals);

        auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
        typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(vals), decltype(rows), decltype(indptrs)> SparseMat; 
        auto ref = std::shared_ptr<tatami::NumericMatrix>(new SparseMat(NR, NC, std::move(vals), std::move(rows), std::move(indptrs))); 

        auto out = tatami::convert_to_layered_sparse(ref.get());

        auto rwrk = ref->dense_row();
        auto owrk = out.matrix->dense_row();
        for (size_t i = 0; i < NR; ++i) {
            auto stuff = owrk->fetch(out.permutation[i]);
            EXPECT_EQ(stuff, rwrk->fetch(i));
        }
    }

    // Checking in the other orientation. 
    {
        std::vector<size_t> rows, cols;
        std::vector<int> vals;
        mock_layered_sparse_data<true>(NR, NC, rows, cols, vals);

        auto indptrs = tatami::compress_sparse_triplets<true>(NR, NC, vals, rows, cols);
        typedef tatami::CompressedSparseRowMatrix<double, int, decltype(vals), decltype(cols), decltype(indptrs)> SparseMat; 
        auto ref = std::shared_ptr<tatami::NumericMatrix>(new SparseMat(NR, NC, std::move(vals), std::move(cols), std::move(indptrs))); 

        auto out = tatami::convert_to_layered_sparse(ref.get());

        auto rwrk = ref->dense_row();
        auto owrk = out.matrix->dense_row();
        for (size_t i = 0; i < NR; ++i) {
            auto stuff = owrk->fetch(out.permutation[i]);
            EXPECT_EQ(stuff, rwrk->fetch(i));
        }
    }
}

INSTANTIATE_TEST_CASE_P(
    ConvertToLayeredSparseHard,
    ConvertToLayeredSparseHardTest,
    ::testing::Values(1000, 20000)
);

