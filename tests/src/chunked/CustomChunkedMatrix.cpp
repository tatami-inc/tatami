#include <gtest/gtest.h>
#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/CompressedSparseMatrix.hpp"
#include "tatami/chunked/CustomChunkedMatrix.hpp"
#include "tatami_test/tatami_test.hpp"

#include "mock_chunk.h"

class DenseCustomChunkedMatrixMethods {
protected:
    std::unique_ptr<tatami::Matrix<double, int> > ref, mat;

    typedef tatami::SimpleDenseChunkWrapper<MockDenseChunk<true> > Chunk;

    void assemble(std::pair<int, int> matdim, std::pair<int, int> chunkdim, bool rowmajor, int cache_size) {
        auto full = tatami_test::simulate_dense_vector<double>(matdim.first * matdim.second, -10, 10, 
            /* seed = */ matdim.first * matdim.second + chunkdim.first * chunkdim.second + rowmajor);
        ref.reset(new tatami::DenseRowMatrix<double, int>(matdim.first, matdim.second, std::move(full)));

        int chunk_nrow = chunkdim.first;
        int chunk_ncol = chunkdim.second;
        int num_chunks_per_row = (matdim.second + chunkdim.second - 1) / chunkdim.second;
        int num_chunks_per_column = (matdim.first + chunkdim.first - 1) / chunkdim.first;
        std::vector<Chunk> chunks(num_chunks_per_row * num_chunks_per_column);

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

                chunks[rowmajor ? r * num_chunks_per_row + c : c * num_chunks_per_column + r] = Chunk(MockDenseChunk<true>(chunkdim.first, chunkdim.second, std::move(contents)));
            }
        }

        tatami::CustomChunkedOptions opt;
        opt.maximum_cache_size = cache_size;
        opt.require_minimum_cache = false;

        mat.reset(new tatami::CustomChunkedDenseMatrix<double, int, Chunk>(
            matdim.first,
            matdim.second,
            chunk_nrow,
            chunk_ncol,
            num_chunks_per_column,
            num_chunks_per_row,
            std::move(chunks),
            rowmajor,
            opt
        ));
    }
};

/*******************************************************/

class DenseCustomChunkedMatrixFullTest :
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, int, bool, int> >, 
    public DenseCustomChunkedMatrixMethods {};

TEST_P(DenseCustomChunkedMatrixFullTest, Column) {
    auto param = GetParam();
    auto cache_size = std::get<3>(param);
    assemble(std::get<0>(param), std::get<1>(param), std::get<2>(param), cache_size);

    bool FORWARD = std::get<4>(param);
    size_t JUMP = std::get<5>(param);

    tatami_test::test_simple_column_access(mat.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(mat.get(), ref.get(), FORWARD, JUMP);

    if (cache_size && FORWARD && JUMP == 1) {
        tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ true);
        tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ false);
    }
}

TEST_P(DenseCustomChunkedMatrixFullTest, Row) {
    auto param = GetParam();
    auto cache_size = std::get<3>(param);
    assemble(std::get<0>(param), std::get<1>(param), std::get<2>(param), cache_size);

    bool FORWARD = std::get<4>(param);
    size_t JUMP = std::get<5>(param);

    tatami_test::test_simple_row_access(mat.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(mat.get(), ref.get(), FORWARD, JUMP);

    if (cache_size && FORWARD && JUMP == 1) {
        tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ true);
        tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ false);
    }
}

INSTANTIATE_TEST_SUITE_P(
    CustomChunkedMatrix,
    DenseCustomChunkedMatrixFullTest,
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
        ::testing::Values(0, 1000, 10000), // cache size

        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 4) // jump, to test the workspace's memory.
    )
);

/*******************************************************/

class DenseCustomChunkedMatrixBlockTest :
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, int, std::pair<double, double> > >, 
    public DenseCustomChunkedMatrixMethods {};

TEST_P(DenseCustomChunkedMatrixBlockTest, Row) {
    auto param = GetParam();
    auto cache_size = std::get<3>(param);
    assemble(std::get<0>(param), std::get<1>(param), std::get<2>(param), cache_size);

    bool FORWARD = true;
    size_t JUMP = 1;
    auto bounds = std::get<4>(param);
    int FIRST = bounds.first * ref->ncol();
    int LAST = bounds.second * ref->ncol();

    tatami_test::test_sliced_row_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    tatami_test::test_sliced_row_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);

    if (cache_size) {
        tatami_test::test_oracle_row_access(mat.get(), ref.get(), /* random = */ true, FIRST, LAST - FIRST);
        tatami_test::test_oracle_row_access(mat.get(), ref.get(), /* random = */ false, FIRST, LAST - FIRST);
    }
}

TEST_P(DenseCustomChunkedMatrixBlockTest, Column) {
    auto param = GetParam();
    auto cache_size = std::get<3>(param);
    assemble(std::get<0>(param), std::get<1>(param), std::get<2>(param), cache_size);

    bool FORWARD = true;
    size_t JUMP = 1;
    auto bounds = std::get<4>(param);
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
    DenseCustomChunkedMatrixBlockTest,
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
        ::testing::Values(0, 1000, 10000), // cache size

        ::testing::Values( // block boundaries
            std::make_pair(0.0, 0.35),
            std::make_pair(0.15, 0.87),
            std::make_pair(0.38, 1.0)
        )
    )
);

/*******************************************************/

class DenseCustomChunkedMatrixIndexTest :
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, int, std::pair<double, double> > >, 
    public DenseCustomChunkedMatrixMethods {};

TEST_P(DenseCustomChunkedMatrixIndexTest, Row) {
    auto param = GetParam();
    auto cache_size = std::get<3>(param);
    assemble(std::get<0>(param), std::get<1>(param), std::get<2>(param), cache_size);

    bool FORWARD = true;
    size_t JUMP = 1;
    auto bounds = std::get<4>(param);
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

TEST_P(DenseCustomChunkedMatrixIndexTest, Column) {
    auto param = GetParam();
    auto cache_size = std::get<3>(param);
    assemble(std::get<0>(param), std::get<1>(param), std::get<2>(param), cache_size);

    bool FORWARD = true;
    size_t JUMP = 1;
    auto bounds = std::get<4>(param);
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
    DenseCustomChunkedMatrixIndexTest,
    ::testing::Combine(
        ::testing::Values( // matrix dimensions
            std::make_pair(198, 67),
            std::make_pair(187, 300),
            std::make_pair(152, 211)
        ),

        ::testing::Values( // chunk dimensions
            std::make_pair(1, 20),
            std::make_pair(20, 1),
            std::make_pair(10, 10)
        ),

        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 4), // jump, to test the workspace's memory.

        ::testing::Values( // index information.
            std::make_pair(0.0, 10),
            std::make_pair(0.2, 5),
            std::make_pair(0.7, 3)
        )
    )
);

/*******************************************************/

//class CustomSparseChunkManagerMethods {
//protected:
//    struct Chunk {
//        static constexpr bool sparse = true;
//        typedef int index_type;
//        typedef double value_type;
//
//        bool row_major;
//        std::vector<double> vcontents;
//        std::vector<int> icontents;
//        std::vector<size_t> pcontents;
//
//        void inflate(std::vector<double>& vbuffer, std::vector<int>& ibuffer, std::vector<size_t>& pbuffer) const {
//            vbuffer.resize(vcontents.size());
//            std::copy(vcontents.begin(), vcontents.end(), vbuffer.begin());
//            ibuffer.resize(icontents.size());
//            std::copy(icontents.begin(), icontents.end(), ibuffer.begin());
//            pbuffer.resize(pcontents.size());
//            std::copy(pcontents.begin(), pcontents.end(), pbuffer.begin());
//        }
//    };
//
//    std::unique_ptr<tatami::Matrix<double, int> > ref;
//    tatami::CustomChunkManager<Chunk> manager;
//
//    void assemble(std::pair<int, int> matdim, std::pair<int, int> chunkdim, bool rowmajor, bool chunkrowmajor) {
//        auto full = tatami_test::simulate_sparse_compressed<double>(matdim.first, matdim.second, 0.1, -10, 10,
//            /* seed = */ matdim.first * matdim.second + chunkdim.first * chunkdim.second + rowmajor - chunkrowmajor);
//        ref.reset(new tatami::CompressedSparseRowMatrix<double, int>(matdim.first, matdim.second, std::move(full.value), std::move(full.index), std::move(full.ptr)));
//
//        manager.chunk_nrow = chunkdim.first;
//        manager.chunk_ncol = chunkdim.second;
//        manager.num_chunks_per_row = (matdim.second + chunkdim.second - 1) / chunkdim.second;
//        manager.num_chunks_per_column = (matdim.first + chunkdim.first - 1) / chunkdim.first;
//        manager.row_major = rowmajor;
//        manager.chunks.resize(manager.num_chunks_per_row * manager.num_chunks_per_column);
//
//        for (int r = 0; r < manager.num_chunks_per_column; ++r) {
//            for (int c = 0; c < manager.num_chunks_per_row; ++c) {
//                auto cstart = c * manager.chunk_ncol;
//                auto cend = std::min(cstart + manager.chunk_ncol, static_cast<size_t>(matdim.second));
//                auto clen = cend - cstart;
//
//                auto rstart = r * manager.chunk_nrow;
//                auto rend = std::min(rstart + manager.chunk_nrow, static_cast<size_t>(matdim.first));
//                auto rlen = rend - rstart;
//
//                auto& current = manager.chunks[rowmajor ? r * manager.num_chunks_per_row + c : c * manager.num_chunks_per_column + r];
//                current.row_major = chunkrowmajor;
//                current.pcontents.push_back(0);
//
//                if (chunkrowmajor) {
//                    auto ext = ref->sparse_row(cstart, clen);
//                    std::vector<double> vbuffer(clen);
//                    std::vector<int> ibuffer(clen);
//
//                    for (int r2 = 0; r2 < rlen; ++r2) {
//                        auto range = ext->fetch(r2 + rstart, vbuffer.data(), ibuffer.data());
//                        current.vcontents.insert(current.vcontents.end(), range.value, range.value + range.number);
//                        for (int i = 0; i < range.number; ++i) {
//                            current.icontents.push_back(range.index[i] - cstart);
//                        }
//                        current.pcontents.push_back(current.pcontents.back() + range.number);
//                    }
//
//                } else {
//                    auto ext = ref->sparse_column(rstart, rlen);
//                    std::vector<double> vbuffer(rlen);
//                    std::vector<int> ibuffer(rlen);
//
//                    for (int c2 = 0; c2 < clen; ++c2) {
//                        auto range = ext->fetch(c2 + cstart, vbuffer.data(), ibuffer.data());
//                        current.vcontents.insert(current.vcontents.end(), range.value, range.value + range.number);
//                        for (int i = 0; i < range.number; ++i) {
//                            current.icontents.push_back(range.index[i] - rstart);
//                        }
//                        current.pcontents.push_back(current.pcontents.back() + range.number);
//                    }
//                }
//            }
//        }
//    }
//};
//
///*******************************************************/
//
//class CustomSparseChunkManagerFullTest : 
//    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool> >, 
//    public CustomSparseChunkManagerMethods {};
//
//TEST_P(CustomSparseChunkManagerFullTest, Row) {
//    auto param = GetParam();
//    assemble(std::get<0>(param), std::get<1>(param), std::get<2>(param), std::get<3>(param));
//
//    auto slab = manager.create_slab<true, false>(ref->ncol());
//    auto eslab = manager.create_slab<true, true>(ref->ncol());
//    auto work = manager.create_workspace();
//
//    auto ref_ext = ref->sparse_row();
//    int lastr = -1;
//
//    for (int r = 0; r < ref->nrow(); ++r) {
//        int requiredr = r / manager.chunk_nrow;
//        int offsetr = r % manager.chunk_nrow;
//        if (requiredr != lastr) {
//            manager.extract<true, false>(requiredr, offsetr, ref->nrow(), 0, ref->ncol(), slab, work);
//            lastr = requiredr;
//        }
//
//        auto ref_range = ref_ext->fetch(r);
//        EXPECT_EQ(ref_range.value, slab.values[offsetr]);
//        EXPECT_EQ(ref_range.index, slab.indices[offsetr]);
//
//        // Testing with exact.
//        manager.extract<true, true>(requiredr, offsetr, ref->nrow(), 0, ref->ncol(), eslab, work);
//        EXPECT_EQ(ref_range.value, eslab.values[0]);
//        EXPECT_EQ(ref_range.index, eslab.indices[0]);
//    }
//}
//
//TEST_P(CustomSparseChunkManagerFullTest, Column) {
//    auto param = GetParam();
//    assemble(std::get<0>(param), std::get<1>(param), std::get<2>(param), std::get<3>(param));
//
//    auto slab = manager.create_slab<false, false>(ref->nrow());
//    auto eslab = manager.create_slab<false, true>(ref->nrow());
//    auto work = manager.create_workspace();
//
//    auto ref_ext = ref->sparse_column();
//    int lastc = -1;
//
//    for (int c = 0; c < ref->ncol(); ++c) {
//        int requiredc = c / manager.chunk_ncol;
//        int offsetc = c % manager.chunk_ncol;
//        if (requiredc != lastc) {
//            manager.extract<false, false>(requiredc, offsetc, ref->ncol(), 0, ref->nrow(), slab, work);
//            lastc = requiredc;
//        }
//
//        auto ref_range = ref_ext->fetch(c);
//        EXPECT_EQ(ref_range.value, slab.values[offsetc]);
//        EXPECT_EQ(ref_range.index, slab.indices[offsetc]);
//
//        // Testing with exact.
//        manager.extract<false, true>(requiredc, offsetc, ref->ncol(), 0, ref->nrow(), eslab, work);
//        EXPECT_EQ(ref_range.value, eslab.values[0]);
//        EXPECT_EQ(ref_range.index, eslab.indices[0]);
//    }
//}
//
//INSTANTIATE_TEST_SUITE_P(
//    CustomChunkManager,
//    CustomSparseChunkManagerFullTest,
//    ::testing::Combine(
//        ::testing::Values( // matrix dimensions
//            std::make_pair(200, 50),
//            std::make_pair(100, 300),
//            std::make_pair(152, 211),
//            std::make_pair(512, 32)
//        ),
//        ::testing::Values( // chunk dimensions
//            std::make_pair(1, 20),
//            std::make_pair(20, 1),
//            std::make_pair(10, 10),
//            std::make_pair(11, 13) // odd numbers
//        ),
//        ::testing::Values(true, false), // row major
//        ::testing::Values(true, false) // chunk is row major
//    )
//);
//
///*******************************************************/
//
//class CustomSparseChunkManagerBlockTest : 
//    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, std::pair<double, double>, bool, bool> >, 
//    public CustomSparseChunkManagerMethods {};
//
//TEST_P(CustomSparseChunkManagerBlockTest, Row) {
//    auto param = GetParam();
//    assemble(std::get<0>(param), std::get<1>(param), std::get<3>(param), std::get<4>(param));
//
//    auto bounds = std::get<2>(param);
//    int start = bounds.first * ref->ncol();
//    int len = bounds.second * ref->ncol() - start;
//
//    auto slab = manager.create_slab<true, false>(len);
//    auto eslab = manager.create_slab<true, true>(len);
//    auto work = manager.create_workspace();
//
//    auto ref_ext = ref->sparse_row(start, len);
//    int lastr = -1;
//
//    for (int r = 0; r < ref->nrow(); ++r) {
//        int requiredr = r / manager.chunk_nrow;
//        int offsetr = r % manager.chunk_nrow;
//        if (requiredr != lastr) {
//            manager.extract<true, false>(requiredr, offsetr, ref->nrow(), start, len, slab, work);
//            lastr = requiredr;
//        }
//
//        auto ref_range = ref_ext->fetch(r);
//        EXPECT_EQ(ref_range.value, slab.values[offsetr]);
//        EXPECT_EQ(ref_range.index, slab.indices[offsetr]);
//
//        // Testing with exact.
//        manager.extract<true, true>(requiredr, offsetr, ref->nrow(), start, len, eslab, work);
//        EXPECT_EQ(ref_range.value, eslab.values[0]);
//        EXPECT_EQ(ref_range.index, eslab.indices[0]);
//    }
//}
//
//TEST_P(CustomSparseChunkManagerBlockTest, Column) {
//    auto param = GetParam();
//    assemble(std::get<0>(param), std::get<1>(param), std::get<3>(param), std::get<4>(param));
//
//    auto bounds = std::get<2>(param);
//    int start = bounds.first * ref->nrow();
//    int len = bounds.second * ref->nrow() - start;
//
//    auto slab = manager.create_slab<false, false>(len);
//    auto eslab = manager.create_slab<false, true>(len);
//    auto work = manager.create_workspace();
//
//    auto ref_ext = ref->sparse_column(start, len);
//    int lastc = -1;
//
//    for (int c = 0; c < ref->ncol(); ++c) {
//        int requiredc = c / manager.chunk_ncol;
//        int offsetc = c % manager.chunk_ncol;
//        if (requiredc != lastc) {
//            manager.extract<false, false>(requiredc, offsetc, ref->ncol(), start, len, slab, work);
//            lastc = requiredc;
//        }
//
//        auto ref_range = ref_ext->fetch(c);
//        EXPECT_EQ(ref_range.value, slab.values[offsetc]);
//        EXPECT_EQ(ref_range.index, slab.indices[offsetc]);
//
//        // Testing with exact.
//        manager.extract<false, true>(requiredc, offsetc, ref->ncol(), start, len, eslab, work);
//        EXPECT_EQ(ref_range.value, eslab.values[0]);
//        EXPECT_EQ(ref_range.index, eslab.indices[0]);
//    }
//}
//
//INSTANTIATE_TEST_SUITE_P(
//    CustomChunkManager,
//    CustomSparseChunkManagerBlockTest,
//    ::testing::Combine(
//        ::testing::Values( // matrix dimensions
//            std::make_pair(200, 50),
//            std::make_pair(100, 300),
//            std::make_pair(152, 211)
//        ),
//        ::testing::Values( // chunk dimensions
//            std::make_pair(1, 20),
//            std::make_pair(20, 1),
//            std::make_pair(10, 10)
//        ),
//        ::testing::Values( // block boundaries
//            std::make_pair(0.0, 0.35),
//            std::make_pair(0.15, 0.87),
//            std::make_pair(0.38, 1.0)
//        ),
//        ::testing::Values(true, false), // row major
//        ::testing::Values(true, false) // chunk is row major
//    )
//);
//
///*******************************************************/
//
//class CustomSparseChunkManagerIndexTest : 
//    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, std::pair<double, double>, bool, bool> >, 
//    public CustomSparseChunkManagerMethods 
//{
//protected:
//    static std::vector<int> get_indices(std::pair<double, double> bounds, int range) {
//        int start = bounds.first * range;
//        int jump = bounds.second;
//        std::vector<int> indices;
//        while (start < range) {
//            indices.push_back(start);
//            start += jump;
//        }
//        return indices;
//    }
//};
//
//TEST_P(CustomSparseChunkManagerIndexTest, Row) {
//    auto param = GetParam();
//    assemble(std::get<0>(param), std::get<1>(param), std::get<3>(param), std::get<4>(param));
//    auto indices = get_indices(std::get<2>(param), ref->ncol());
//
//    auto slab = manager.create_slab<true, false>(indices.size());
//    auto eslab = manager.create_slab<true, true>(indices.size());
//    auto work = manager.create_workspace();
//
//    auto ref_ext = ref->sparse_row(indices);
//    int lastr = -1;
//
//    for (int r = 0; r < ref->nrow(); ++r) {
//        int requiredr = r / manager.chunk_nrow;
//        int offsetr = r % manager.chunk_nrow;
//        if (requiredr != lastr) {
//            manager.extract<true, false>(requiredr, offsetr, ref->nrow(), indices, slab, work);
//            lastr = requiredr;
//        }
//
//        auto ref_range = ref_ext->fetch(r);
//        EXPECT_EQ(ref_range.value, slab.values[offsetr]);
//        EXPECT_EQ(ref_range.index, slab.indices[offsetr]);
//
//        // Testing with exact.
//        manager.extract<true, true>(requiredr, offsetr, ref->nrow(), indices, eslab, work);
//        EXPECT_EQ(ref_range.value, eslab.values[0]);
//        EXPECT_EQ(ref_range.index, eslab.indices[0]);
//    }
//}
//
//TEST_P(CustomSparseChunkManagerIndexTest, Column) {
//    auto param = GetParam();
//    assemble(std::get<0>(param), std::get<1>(param), std::get<3>(param), std::get<4>(param));
//    auto indices = get_indices(std::get<2>(param), ref->nrow());
//
//    auto slab = manager.create_slab<false, false>(indices.size());
//    auto eslab = manager.create_slab<false, true>(indices.size());
//    auto work = manager.create_workspace();
//
//    std::vector<double> tmp1(indices.size()), tmp2(indices.size());
//    auto ref_ext = ref->sparse_column(indices);
//    int lastc = -1;
//
//    for (int c = 0; c < ref->ncol(); ++c) {
//        int requiredc = c / manager.chunk_ncol;
//        int offsetc = c % manager.chunk_ncol;
//        if (requiredc != lastc) {
//            manager.extract<false, false>(requiredc, offsetc, ref->ncol(), indices, slab, work);
//            lastc = requiredc;
//        }
//
//        auto ref_range = ref_ext->fetch(c);
//        EXPECT_EQ(ref_range.value, slab.values[offsetc]);
//        EXPECT_EQ(ref_range.index, slab.indices[offsetc]);
//
//        // Testing with exact.
//        manager.extract<false, true>(requiredc, offsetc, ref->ncol(), indices, eslab, work);
//        EXPECT_EQ(ref_range.value, eslab.values[0]);
//        EXPECT_EQ(ref_range.index, eslab.indices[0]);
//    }
//}
//
//INSTANTIATE_TEST_SUITE_P(
//    CustomChunkManager,
//    CustomSparseChunkManagerIndexTest,
//    ::testing::Combine(
//        ::testing::Values( // matrix dimensions
//            std::make_pair(200, 50),
//            std::make_pair(100, 300),
//            std::make_pair(152, 211)
//        ),
//        ::testing::Values( // chunk dimensions
//            std::make_pair(1, 20),
//            std::make_pair(20, 1),
//            std::make_pair(10, 10)
//        ),
//        ::testing::Values( // index information.
//            std::make_pair(0.0, 10),
//            std::make_pair(0.2, 5),
//            std::make_pair(0.7, 3)
//        ),
//        ::testing::Values(true, false), // row major
//        ::testing::Values(true, false) // chunk is row major
//    )
//);
//
