#include <gtest/gtest.h>
#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/CompressedSparseMatrix.hpp"
#include "tatami/chunked/CustomChunkedMatrix.hpp"
#include "tatami_test/tatami_test.hpp"

#include "mock_chunk.h"

class CustomChunkedMatrixMethods {
protected:
    std::unique_ptr<tatami::Matrix<double, int> > ref, mat;

    typedef tatami::SimpleDenseChunkWrapper<MockDenseChunk<true> > DChunk;
    typedef tatami::SimpleSparseChunkWrapper<MockSparseChunk<false> > SChunk;

    void dense_assemble(std::pair<int, int> matdim, std::pair<int, int> chunkdim, bool rowmajor, int cache_size) {
        auto full = tatami_test::simulate_dense_vector<double>(matdim.first * matdim.second, -10, 10, 
            /* seed = */ matdim.first * matdim.second + chunkdim.first * chunkdim.second + rowmajor);
        ref.reset(new tatami::DenseRowMatrix<double, int>(matdim.first, matdim.second, std::move(full)));

        int chunk_nrow = chunkdim.first;
        int chunk_ncol = chunkdim.second;
        int num_chunks_per_row = (matdim.second + chunkdim.second - 1) / chunkdim.second;
        int num_chunks_per_column = (matdim.first + chunkdim.first - 1) / chunkdim.first;
        std::vector<DChunk> chunks(num_chunks_per_row * num_chunks_per_column);

        for (int r = 0; r < num_chunks_per_column; ++r) {
            for (int c = 0; c < num_chunks_per_row; ++c) {
                auto cstart = c * chunk_ncol;
                auto cend = std::min(cstart + chunk_ncol, static_cast<int>(matdim.second));
                auto clen = cend - cstart;

                auto rstart = r * chunk_nrow;
                auto rend = std::min(rstart + chunk_nrow, static_cast<int>(matdim.first));
                auto rlen = rend - rstart;

                std::vector<double> contents(chunkdim.first * chunkdim.second);
                auto ccptr = contents.data();
                auto ext = ref->dense_row(cstart, clen);
                for (int r2 = 0; r2 < rlen; ++r2) {
                    ext->fetch_copy(r2 + rstart, ccptr);
                    ccptr += chunkdim.second;
                }

                chunks[rowmajor ? r * num_chunks_per_row + c : c * num_chunks_per_column + r] = DChunk(MockDenseChunk<true>(chunkdim.first, chunkdim.second, std::move(contents)));
            }
        }

        tatami::CustomChunkedOptions opt;
        opt.maximum_cache_size = cache_size;
        opt.require_minimum_cache = false;

        mat.reset(new tatami::CustomChunkedDenseMatrix<double, int, DChunk>(
            matdim.first,
            matdim.second,
            chunk_nrow,
            chunk_ncol,
            std::move(chunks),
            rowmajor,
            opt
        ));
    }

    void sparse_assemble(std::pair<int, int> matdim, std::pair<int, int> chunkdim, bool rowmajor, int cache_size) {
        auto full = tatami_test::simulate_sparse_compressed<double>(matdim.second, matdim.first, 0.1, -10, 10,
            /* seed = */ matdim.first * matdim.second + chunkdim.first * chunkdim.second + rowmajor);
        ref.reset(new tatami::CompressedSparseColumnMatrix<double, int>(matdim.first, matdim.second, std::move(full.value), std::move(full.index), std::move(full.ptr)));

        int chunk_nrow = chunkdim.first;
        int chunk_ncol = chunkdim.second;
        int num_chunks_per_row = (matdim.second + chunkdim.second - 1) / chunkdim.second;
        int num_chunks_per_column = (matdim.first + chunkdim.first - 1) / chunkdim.first;
        std::vector<SChunk> chunks(num_chunks_per_row * num_chunks_per_column);

        for (int r = 0; r < num_chunks_per_column; ++r) {
            for (int c = 0; c < num_chunks_per_row; ++c) {
                auto cstart = c * chunk_ncol;
                auto cend = std::min(cstart + chunk_ncol, static_cast<int>(matdim.second));
                auto clen = cend - cstart;

                auto rstart = r * chunk_nrow;
                auto rend = std::min(rstart + chunk_nrow, static_cast<int>(matdim.first));
                auto rlen = rend - rstart;

                std::vector<double> vcontents;
                std::vector<int> icontents;
                std::vector<size_t> pcontents(1);

                auto ext = ref->sparse_column(rstart, rlen);
                std::vector<double> vbuffer(rlen);
                std::vector<int> ibuffer(rlen);

                for (int c2 = 0; c2 < clen; ++c2) {
                    auto range = ext->fetch(c2 + cstart, vbuffer.data(), ibuffer.data());
                    vcontents.insert(vcontents.end(), range.value, range.value + range.number);
                    for (int i = 0; i < range.number; ++i) {
                        icontents.push_back(range.index[i] - rstart);
                    }
                    pcontents.push_back(pcontents.back() + range.number);
                }

                chunks[rowmajor ? r * num_chunks_per_row + c : c * num_chunks_per_column + r] = SChunk(MockSparseChunk<false>(
                    chunkdim.first, 
                    chunkdim.second, 
                    std::move(vcontents),
                    std::move(icontents),
                    std::move(pcontents)
                ));
            }
        }

        tatami::CustomChunkedOptions opt;
        opt.maximum_cache_size = cache_size;
        opt.require_minimum_cache = false;

        mat.reset(new tatami::CustomChunkedSparseMatrix<double, int, SChunk>(
            matdim.first,
            matdim.second,
            chunk_nrow,
            chunk_ncol,
            std::move(chunks),
            rowmajor,
            opt
        ));
    }
};

/*******************************************************/

class CustomChunkedMatrixFullTest :
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool, int, bool, int> >, 
    public CustomChunkedMatrixMethods {};

TEST_P(CustomChunkedMatrixFullTest, Column) {
    auto param = GetParam();
    auto matdim = std::get<0>(param);
    auto chunkdim = std::get<1>(param);
    bool rowmajor = std::get<2>(param);
    bool sparse = std::get<3>(param);
    auto cache_size = std::get<4>(param);
    bool FORWARD = std::get<5>(param);
    size_t JUMP = std::get<6>(param);

    if (sparse) {
        sparse_assemble(matdim, chunkdim, rowmajor, cache_size);
    } else {
        dense_assemble(matdim, chunkdim, rowmajor, cache_size);
    }

    tatami_test::test_simple_column_access(mat.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(mat.get(), ref.get(), FORWARD, JUMP);

    if (cache_size && FORWARD && JUMP == 1) {
        tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ true);
        tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ false);
    }
}

TEST_P(CustomChunkedMatrixFullTest, Row) {
    auto param = GetParam();
    auto matdim = std::get<0>(param);
    auto chunkdim = std::get<1>(param);
    bool rowmajor = std::get<2>(param);
    bool sparse = std::get<3>(param);
    auto cache_size = std::get<4>(param);
    bool FORWARD = std::get<5>(param);
    size_t JUMP = std::get<6>(param);

    if (sparse) {
        sparse_assemble(matdim, chunkdim, rowmajor, cache_size);
    } else {
        dense_assemble(matdim, chunkdim, rowmajor, cache_size);
    }

    tatami_test::test_simple_row_access(mat.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(mat.get(), ref.get(), FORWARD, JUMP);

    if (cache_size && FORWARD && JUMP == 1) {
        tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ true);
        tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ false);
    }
}

INSTANTIATE_TEST_SUITE_P(
    CustomChunkedMatrix,
    CustomChunkedMatrixFullTest,
    ::testing::Combine(
        ::testing::Values( // matrix dimensions
            std::make_pair(200, 50),
            std::make_pair(100, 300)
        ),

        ::testing::Values( // chunk dimensions
            std::make_pair(1, 20),
            std::make_pair(20, 1),
            std::make_pair(11, 13) // odd numbers
        ),

        ::testing::Values(true, false), // row major
        ::testing::Values(false, true), // sparse chunks
        ::testing::Values(0, 1000, 10000), // cache size

        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 4) // jump, to test the workspace's memory.
    )
);

/*******************************************************/

class CustomChunkedMatrixBlockTest :
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool, int, std::pair<double, double> > >, 
    public CustomChunkedMatrixMethods {};

TEST_P(CustomChunkedMatrixBlockTest, Row) {
    auto param = GetParam();
    auto matdim = std::get<0>(param);
    auto chunkdim = std::get<1>(param);
    bool rowmajor = std::get<2>(param);
    bool sparse = std::get<3>(param);
    auto cache_size = std::get<4>(param);
    auto bounds = std::get<5>(param);

    if (sparse) {
        sparse_assemble(matdim, chunkdim, rowmajor, cache_size);
    } else {
        dense_assemble(matdim, chunkdim, rowmajor, cache_size);
    }

    bool FORWARD = true;
    size_t JUMP = 1;
    int FIRST = bounds.first * ref->ncol();
    int LAST = bounds.second * ref->ncol();

    tatami_test::test_sliced_row_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    tatami_test::test_sliced_row_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);

    if (cache_size) {
        tatami_test::test_oracle_row_access(mat.get(), ref.get(), /* random = */ true, FIRST, LAST - FIRST);
        tatami_test::test_oracle_row_access(mat.get(), ref.get(), /* random = */ false, FIRST, LAST - FIRST);
    }
}

TEST_P(CustomChunkedMatrixBlockTest, Column) {
    auto param = GetParam();
    auto matdim = std::get<0>(param);
    auto chunkdim = std::get<1>(param);
    bool rowmajor = std::get<2>(param);
    bool sparse = std::get<3>(param);
    auto cache_size = std::get<4>(param);
    auto bounds = std::get<5>(param);

    if (sparse) {
        sparse_assemble(matdim, chunkdim, rowmajor, cache_size);
    } else {
        dense_assemble(matdim, chunkdim, rowmajor, cache_size);
    }

    bool FORWARD = true;
    size_t JUMP = 1;
    int FIRST = bounds.first * ref->nrow();
    int LAST = bounds.second * ref->nrow();

    tatami_test::test_sliced_column_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    tatami_test::test_sliced_column_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);

    if (cache_size) {
        tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ true, FIRST, LAST - FIRST);
        tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ false, FIRST, LAST - FIRST);
    }
}

INSTANTIATE_TEST_SUITE_P(
    CustomChunkedMatrix,
    CustomChunkedMatrixBlockTest,
    ::testing::Combine(
        ::testing::Values( // matrix dimensions
            std::make_pair(201, 67),
            std::make_pair(123, 372)
        ),

        ::testing::Values( // chunk dimensions
            std::make_pair(1, 20),
            std::make_pair(20, 1),
            std::make_pair(10, 10)
        ),

        ::testing::Values(true, false), // row major
        ::testing::Values(false, true), // sparse chunks
        ::testing::Values(0, 1000, 10000), // cache size

        ::testing::Values( // block boundaries
            std::make_pair(0.0, 0.35),
            std::make_pair(0.15, 0.87),
            std::make_pair(0.38, 1.0)
        )
    )
);

/*******************************************************/

class CustomChunkedMatrixIndexTest :
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool, int, std::pair<double, double> > >, 
    public CustomChunkedMatrixMethods {};

TEST_P(CustomChunkedMatrixIndexTest, Row) {
    auto param = GetParam();
    auto matdim = std::get<0>(param);
    auto chunkdim = std::get<1>(param);
    bool rowmajor = std::get<2>(param);
    bool sparse = std::get<3>(param);
    auto cache_size = std::get<4>(param);
    auto bounds = std::get<5>(param);

    if (sparse) {
        sparse_assemble(matdim, chunkdim, rowmajor, cache_size);
    } else {
        dense_assemble(matdim, chunkdim, rowmajor, cache_size);
    }

    bool FORWARD = true;
    size_t JUMP = 1;
    int FIRST = bounds.first * ref->ncol(), STEP = bounds.second;

    tatami_test::test_indexed_row_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    tatami_test::test_indexed_row_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);

    if (cache_size) {
        int NC = ref->ncol();
        std::vector<int> indices;
        {
            int counter = FIRST;
            while (counter < NC) {
                indices.push_back(counter);
                counter += STEP;
            }
        }

        tatami_test::test_oracle_row_access(mat.get(), ref.get(), /* random = */ true, indices);
        tatami_test::test_oracle_row_access(mat.get(), ref.get(), /* random = */ false, indices);
    }
}

TEST_P(CustomChunkedMatrixIndexTest, Column) {
    auto param = GetParam();
    auto matdim = std::get<0>(param);
    auto chunkdim = std::get<1>(param);
    bool rowmajor = std::get<2>(param);
    bool sparse = std::get<3>(param);
    auto cache_size = std::get<4>(param);
    auto bounds = std::get<5>(param);

    if (sparse) {
        sparse_assemble(matdim, chunkdim, rowmajor, cache_size);
    } else {
        dense_assemble(matdim, chunkdim, rowmajor, cache_size);
    }

    bool FORWARD = true;
    size_t JUMP = 1;
    int FIRST = bounds.first * ref->nrow(), STEP = bounds.second;

    tatami_test::test_indexed_column_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    tatami_test::test_indexed_column_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);

    if (cache_size) {
        int NR = ref->nrow();
        std::vector<int> indices;
        {
            int counter = FIRST;
            while (counter < NR) {
                indices.push_back(counter);
                counter += STEP;
            }
        }

        tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ true, indices);
        tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ false, indices);
    }
}

INSTANTIATE_TEST_SUITE_P(
    CustomChunkedMatrix,
    CustomChunkedMatrixIndexTest,
    ::testing::Combine(
        ::testing::Values( // matrix dimensions
            std::make_pair(198, 67),
            std::make_pair(187, 300)
        ),

        ::testing::Values( // chunk dimensions
            std::make_pair(1, 20),
            std::make_pair(20, 1),
            std::make_pair(7, 13)
        ),

        ::testing::Values(true, false), // row major 
        ::testing::Values(false, true), // sparse chunks
        ::testing::Values(0, 1000, 10000), // cache size

        ::testing::Values( // index information.
            std::make_pair(0.0, 10),
            std::make_pair(0.2, 5),
            std::make_pair(0.7, 3)
        )
    )
);
