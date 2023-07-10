#include <gtest/gtest.h>
#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/CompressedSparseMatrix.hpp"
#include "tatami/chunked/CustomChunkManager.hpp"
#include "tatami_test/tatami_test.hpp"

class CustomDenseChunkManagerMethods {
protected:
    struct Chunk {
        static constexpr bool sparse = false;
        typedef int index_type;
        typedef double value_type;

        bool row_major;
        std::vector<double> contents;

        void inflate(std::vector<double>& buffer) const {
            buffer.resize(contents.size());
            std::copy(contents.begin(), contents.end(), buffer.begin());
        }
    };

    std::unique_ptr<tatami::Matrix<double, int> > ref;
    tatami::CustomChunkManager<Chunk> manager;

    void assemble(std::pair<int, int> matdim, std::pair<int, int> chunkdim, bool rowmajor, bool chunkrowmajor) {
        auto full = tatami_test::simulate_dense_vector<double>(matdim.first * matdim.second, -10, 10, 
            /* seed = */ matdim.first * matdim.second + chunkdim.first * chunkdim.second + rowmajor - chunkrowmajor);
        ref.reset(new tatami::DenseRowMatrix<double, int>(matdim.first, matdim.second, std::move(full)));

        manager.chunk_nrow = chunkdim.first;
        manager.chunk_ncol = chunkdim.second;
        manager.num_chunks_per_row = (matdim.second + chunkdim.second - 1) / chunkdim.second;
        manager.num_chunks_per_column = (matdim.first + chunkdim.first - 1) / chunkdim.first;
        manager.row_major = rowmajor;
        manager.chunks.resize(manager.num_chunks_per_row * manager.num_chunks_per_column);

        for (int r = 0; r < manager.num_chunks_per_column; ++r) {
            for (int c = 0; c < manager.num_chunks_per_row; ++c) {
                auto cstart = c * manager.chunk_ncol;
                auto cend = std::min(cstart + manager.chunk_ncol, static_cast<size_t>(matdim.second));
                auto clen = cend - cstart;

                auto rstart = r * manager.chunk_nrow;
                auto rend = std::min(rstart + manager.chunk_nrow, static_cast<size_t>(matdim.first));
                auto rlen = rend - rstart;

                auto& current = manager.chunks[rowmajor ? r * manager.num_chunks_per_row + c : c * manager.num_chunks_per_column + r];
                current.row_major = chunkrowmajor;
                current.contents.resize(chunkdim.first * chunkdim.second);

                if (chunkrowmajor) {
                    auto ext = ref->dense_row(cstart, clen);
                    auto ccptr = current.contents.data();
                    for (int r2 = 0; r2 < rlen; ++r2) {
                        ext->fetch_copy(r2 + rstart, ccptr);
                        ccptr += chunkdim.second;
                    }
                } else {
                    auto ext = ref->dense_column(rstart, rlen);
                    auto ccptr = current.contents.data();
                    for (int c2 = 0; c2 < clen; ++c2) {
                        ext->fetch_copy(c2 + cstart, ccptr);
                        ccptr += chunkdim.first;
                    }
                }
            }
        }
    }
};

/*******************************************************/

class CustomDenseChunkManagerFullTest : 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool> >, 
    public CustomDenseChunkManagerMethods {};

TEST_P(CustomDenseChunkManagerFullTest, Row) {
    auto param = GetParam();
    assemble(std::get<0>(param), std::get<1>(param), std::get<2>(param), std::get<3>(param));

    auto cache = manager.create_chunk_cache<true, false>(ref->ncol());
    auto ecache = manager.create_chunk_cache<true, true>(ref->ncol());
    std::vector<double> tmp1(ref->ncol()), tmp2(ref->ncol());
    auto ref_ext = ref->dense_row();
    int lastr = -1;
    const double* ccptr = NULL;

    for (int r = 0; r < ref->nrow(); ++r) {
        int requiredr = r / manager.chunk_nrow;
        if (requiredr != lastr) {
            manager.extract<true, false>(requiredr, r % manager.chunk_nrow, ref->nrow(), 0, ref->ncol(), cache);
            lastr = requiredr;
            ccptr = cache.cache.data();
        }

        ref_ext->fetch_copy(r, tmp1.data());
        std::copy(ccptr, ccptr + tmp2.size(), tmp2.data());
        EXPECT_EQ(tmp1, tmp2);
        ccptr += ref->ncol();

        // Testing with exact.
        manager.extract<true, true>(requiredr, r % manager.chunk_nrow, ref->nrow(), 0, ref->ncol(), ecache);
        EXPECT_EQ(tmp1, ecache.cache);
    }
}

TEST_P(CustomDenseChunkManagerFullTest, Column) {
    auto param = GetParam();
    assemble(std::get<0>(param), std::get<1>(param), std::get<2>(param), std::get<3>(param));

    auto cache = manager.create_chunk_cache<false, false>(ref->nrow());
    auto ecache = manager.create_chunk_cache<false, true>(ref->nrow());
    std::vector<double> tmp1(ref->nrow()), tmp2(ref->nrow());
    auto ref_ext = ref->dense_column();
    int lastc = -1;
    const double* ccptr = NULL;

    for (int c = 0; c < ref->ncol(); ++c) {
        int requiredc = c / manager.chunk_ncol;
        if (requiredc != lastc) {
            manager.extract<false, false>(requiredc, c % manager.chunk_ncol, ref->ncol(), 0, ref->nrow(), cache);
            lastc = requiredc;
            ccptr = cache.cache.data();
        }

        ref_ext->fetch_copy(c, tmp1.data());
        std::copy(ccptr, ccptr + tmp2.size(), tmp2.data());
        EXPECT_EQ(tmp1, tmp2);
        ccptr += ref->nrow();

        // Testing with exact.
        manager.extract<false, true>(requiredc, c % manager.chunk_ncol, ref->ncol(), 0, ref->nrow(), ecache);
        EXPECT_EQ(tmp1, ecache.cache);
    }
}

INSTANTIATE_TEST_SUITE_P(
    CustomChunkManager,
    CustomDenseChunkManagerFullTest,
    ::testing::Combine(
        ::testing::Values( // matrix dimensions
            std::make_pair(200, 50),
            std::make_pair(100, 300),
            std::make_pair(152, 211),
            std::make_pair(512, 32)
        ),
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 20),
            std::make_pair(20, 1),
            std::make_pair(10, 10),
            std::make_pair(11, 13) // odd numbers
        ),
        ::testing::Values(true, false), // row major
        ::testing::Values(true, false) // chunk is row major
    )
);

/*******************************************************/

class CustomDenseChunkManagerBlockTest : 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, std::pair<double, double>, bool, bool> >, 
    public CustomDenseChunkManagerMethods {};

TEST_P(CustomDenseChunkManagerBlockTest, Row) {
    auto param = GetParam();
    assemble(std::get<0>(param), std::get<1>(param), std::get<3>(param), std::get<4>(param));

    auto bounds = std::get<2>(param);
    int start = bounds.first * ref->ncol();
    int len = bounds.second * ref->ncol() - start;

    auto cache = manager.create_chunk_cache<true, false>(len);
    auto ecache = manager.create_chunk_cache<true, true>(len);
    std::vector<double> tmp1(len), tmp2(len);
    auto ref_ext = ref->dense_row(start, len);
    int lastr = -1;
    const double* ccptr = NULL;

    for (int r = 0; r < ref->nrow(); ++r) {
        int requiredr = r / manager.chunk_nrow;
        if (requiredr != lastr) {
            manager.extract<true, false>(requiredr, r % manager.chunk_nrow, ref->nrow(), start, len, cache);
            lastr = requiredr;
            ccptr = cache.cache.data();
        }

        ref_ext->fetch_copy(r, tmp1.data());
        std::copy(ccptr, ccptr + tmp2.size(), tmp2.data());
        EXPECT_EQ(tmp1, tmp2);
        ccptr += len;

        // Testing with exact.
        manager.extract<true, true>(requiredr, r % manager.chunk_nrow, ref->nrow(), start, len, ecache);
        EXPECT_EQ(tmp1, ecache.cache);
    }
}

TEST_P(CustomDenseChunkManagerBlockTest, Column) {
    auto param = GetParam();
    assemble(std::get<0>(param), std::get<1>(param), std::get<3>(param), std::get<4>(param));

    auto bounds = std::get<2>(param);
    int start = bounds.first * ref->nrow();
    int len = bounds.second * ref->nrow() - start;

    auto cache = manager.create_chunk_cache<false, false>(len);
    auto ecache = manager.create_chunk_cache<false, true>(len);
    std::vector<double> tmp1(len), tmp2(len);
    auto ref_ext = ref->dense_column(start, len);
    int lastc = -1;
    const double* ccptr = NULL;

    for (int c = 0; c < ref->ncol(); ++c) {
        int requiredc = c / manager.chunk_ncol;
        if (requiredc != lastc) {
            manager.extract<false, false>(requiredc, c % manager.chunk_ncol, ref->ncol(), start, len, cache);
            lastc = requiredc;
            ccptr = cache.cache.data();
        }

        ref_ext->fetch_copy(c, tmp1.data());
        std::copy(ccptr, ccptr + tmp2.size(), tmp2.data());
        EXPECT_EQ(tmp1, tmp2);
        ccptr += len;

        // Testing with exact.
        manager.extract<false, true>(requiredc, c % manager.chunk_ncol, ref->ncol(), start, len, ecache);
        EXPECT_EQ(tmp1, ecache.cache);
    }
}

INSTANTIATE_TEST_SUITE_P(
    CustomChunkManager,
    CustomDenseChunkManagerBlockTest,
    ::testing::Combine(
        ::testing::Values( // matrix dimensions
            std::make_pair(200, 50),
            std::make_pair(100, 300),
            std::make_pair(152, 211)
        ),
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 20),
            std::make_pair(20, 1),
            std::make_pair(10, 10)
        ),
        ::testing::Values( // block boundaries
            std::make_pair(0.0, 0.35),
            std::make_pair(0.15, 0.87),
            std::make_pair(0.38, 1.0)
        ),
        ::testing::Values(true, false), // row major
        ::testing::Values(true, false) // chunk is row major
    )
);

/*******************************************************/

class CustomDenseChunkManagerIndexTest : 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, std::pair<double, double>, bool, bool> >, 
    public CustomDenseChunkManagerMethods 
{
protected:
    static std::vector<int> get_indices(std::pair<double, double> bounds, int range) {
        int start = bounds.first * range;
        int jump = bounds.second;
        std::vector<int> indices;
        while (start < range) {
            indices.push_back(start);
            start += jump;
        }
        return indices;
    }
};

TEST_P(CustomDenseChunkManagerIndexTest, Row) {
    auto param = GetParam();
    assemble(std::get<0>(param), std::get<1>(param), std::get<3>(param), std::get<4>(param));
    auto indices = get_indices(std::get<2>(param), ref->ncol());

    auto cache = manager.create_chunk_cache<true, false>(indices.size());
    auto ecache = manager.create_chunk_cache<true, true>(indices.size());
    std::vector<double> tmp1(indices.size()), tmp2(indices.size());
    auto ref_ext = ref->dense_row(indices);
    int lastr = -1;
    const double* ccptr = NULL;

    for (int r = 0; r < ref->nrow(); ++r) {
        int requiredr = r / manager.chunk_nrow;
        if (requiredr != lastr) {
            manager.extract<true, false>(requiredr, r % manager.chunk_nrow, ref->nrow(), indices, cache);
            lastr = requiredr;
            ccptr = cache.cache.data();
        }

        ref_ext->fetch_copy(r, tmp1.data());
        std::copy(ccptr, ccptr + tmp2.size(), tmp2.data());
        EXPECT_EQ(tmp1, tmp2);
        ccptr += indices.size();

        // Testing with exact.
        manager.extract<true, true>(requiredr, r % manager.chunk_nrow, ref->nrow(), indices, ecache);
        EXPECT_EQ(tmp1, ecache.cache);
    }
}

TEST_P(CustomDenseChunkManagerIndexTest, Column) {
    auto param = GetParam();
    assemble(std::get<0>(param), std::get<1>(param), std::get<3>(param), std::get<4>(param));
    auto indices = get_indices(std::get<2>(param), ref->nrow());

    auto cache = manager.create_chunk_cache<false, false>(indices.size());
    auto ecache = manager.create_chunk_cache<false, true>(indices.size());
    std::vector<double> tmp1(indices.size()), tmp2(indices.size());
    auto ref_ext = ref->dense_column(indices);
    int lastc = -1;
    const double* ccptr = NULL;

    for (int c = 0; c < ref->ncol(); ++c) {
        int requiredc = c / manager.chunk_ncol;
        if (requiredc != lastc) {
            manager.extract<false, false>(requiredc, c % manager.chunk_ncol, ref->ncol(), indices, cache);
            lastc = requiredc;
            ccptr = cache.cache.data();
        }

        ref_ext->fetch_copy(c, tmp1.data());
        std::copy(ccptr, ccptr + tmp2.size(), tmp2.data());
        EXPECT_EQ(tmp1, tmp2);
        ccptr += indices.size();

        // Testing with exact.
        manager.extract<false, true>(requiredc, c % manager.chunk_ncol, ref->ncol(), indices, ecache);
        EXPECT_EQ(tmp1, ecache.cache);
    }
}

INSTANTIATE_TEST_SUITE_P(
    CustomChunkManager,
    CustomDenseChunkManagerIndexTest,
    ::testing::Combine(
        ::testing::Values( // matrix dimensions
            std::make_pair(200, 50),
            std::make_pair(100, 300),
            std::make_pair(152, 211)
        ),
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 20),
            std::make_pair(20, 1),
            std::make_pair(10, 10)
        ),
        ::testing::Values( // index information.
            std::make_pair(0.0, 10),
            std::make_pair(0.2, 5),
            std::make_pair(0.7, 3)
        ),
        ::testing::Values(true, false), // row major
        ::testing::Values(true, false) // chunk is row major
    )
);

/*******************************************************/

class CustomSparseChunkManagerMethods {
protected:
    struct Chunk {
        static constexpr bool sparse = true;
        typedef int index_type;
        typedef double value_type;

        bool row_major;
        std::vector<double> vcontents;
        std::vector<int> icontents;
        std::vector<size_t> pcontents;

        void inflate(std::vector<double>& vbuffer, std::vector<int>& ibuffer, std::vector<size_t>& pbuffer) const {
            vbuffer.resize(vcontents.size());
            std::copy(vcontents.begin(), vcontents.end(), vbuffer.begin());
            ibuffer.resize(icontents.size());
            std::copy(icontents.begin(), icontents.end(), ibuffer.begin());
            pbuffer.resize(pcontents.size());
            std::copy(pcontents.begin(), pcontents.end(), pbuffer.begin());
        }
    };

    std::unique_ptr<tatami::Matrix<double, int> > ref;
    tatami::CustomChunkManager<Chunk> manager;

    void assemble(std::pair<int, int> matdim, std::pair<int, int> chunkdim, bool rowmajor, bool chunkrowmajor) {
        auto full = tatami_test::simulate_sparse_compressed<double>(matdim.first, matdim.second, 0.1, -10, 10,
            /* seed = */ matdim.first * matdim.second + chunkdim.first * chunkdim.second + rowmajor - chunkrowmajor);
        ref.reset(new tatami::CompressedSparseRowMatrix<double, int>(matdim.first, matdim.second, std::move(full.value), std::move(full.index), std::move(full.ptr)));

        manager.chunk_nrow = chunkdim.first;
        manager.chunk_ncol = chunkdim.second;
        manager.num_chunks_per_row = (matdim.second + chunkdim.second - 1) / chunkdim.second;
        manager.num_chunks_per_column = (matdim.first + chunkdim.first - 1) / chunkdim.first;
        manager.row_major = rowmajor;
        manager.chunks.resize(manager.num_chunks_per_row * manager.num_chunks_per_column);

        for (int r = 0; r < manager.num_chunks_per_column; ++r) {
            for (int c = 0; c < manager.num_chunks_per_row; ++c) {
                auto cstart = c * manager.chunk_ncol;
                auto cend = std::min(cstart + manager.chunk_ncol, static_cast<size_t>(matdim.second));
                auto clen = cend - cstart;

                auto rstart = r * manager.chunk_nrow;
                auto rend = std::min(rstart + manager.chunk_nrow, static_cast<size_t>(matdim.first));
                auto rlen = rend - rstart;

                auto& current = manager.chunks[rowmajor ? r * manager.num_chunks_per_row + c : c * manager.num_chunks_per_column + r];
                current.row_major = chunkrowmajor;
                current.pcontents.push_back(0);

                if (chunkrowmajor) {
                    auto ext = ref->sparse_row(cstart, clen);
                    std::vector<double> vbuffer(clen);
                    std::vector<int> ibuffer(clen);

                    for (int r2 = 0; r2 < rlen; ++r2) {
                        auto range = ext->fetch(r2 + rstart, vbuffer.data(), ibuffer.data());
                        current.vcontents.insert(current.vcontents.end(), range.value, range.value + range.number);
                        for (int i = 0; i < range.number; ++i) {
                            current.icontents.push_back(range.index[i] - cstart);
                        }
                        current.pcontents.push_back(current.pcontents.back() + range.number);
                    }

                } else {
                    auto ext = ref->sparse_column(rstart, rlen);
                    std::vector<double> vbuffer(rlen);
                    std::vector<int> ibuffer(rlen);

                    for (int c2 = 0; c2 < clen; ++c2) {
                        auto range = ext->fetch(c2 + cstart, vbuffer.data(), ibuffer.data());
                        current.vcontents.insert(current.vcontents.end(), range.value, range.value + range.number);
                        for (int i = 0; i < range.number; ++i) {
                            current.icontents.push_back(range.index[i] - rstart);
                        }
                        current.pcontents.push_back(current.pcontents.back() + range.number);
                    }
                }
            }
        }
    }
};

/*******************************************************/

class CustomSparseChunkManagerFullTest : 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool> >, 
    public CustomSparseChunkManagerMethods {};

TEST_P(CustomSparseChunkManagerFullTest, Row) {
    auto param = GetParam();
    assemble(std::get<0>(param), std::get<1>(param), std::get<2>(param), std::get<3>(param));

    auto cache = manager.create_chunk_cache<true, false>(ref->ncol());
    auto ecache = manager.create_chunk_cache<true, true>(ref->ncol());
    auto ref_ext = ref->sparse_row();
    int lastr = -1;

    for (int r = 0; r < ref->nrow(); ++r) {
        int requiredr = r / manager.chunk_nrow;
        int offsetr = r % manager.chunk_nrow;
        if (requiredr != lastr) {
            manager.extract<true, false>(requiredr, offsetr, ref->nrow(), 0, ref->ncol(), cache);
            lastr = requiredr;
        }

        auto ref_range = ref_ext->fetch(r);
        EXPECT_EQ(ref_range.value, cache.cache_values[offsetr]);
        EXPECT_EQ(ref_range.index, cache.cache_indices[offsetr]);

        // Testing with exact.
        manager.extract<true, true>(requiredr, offsetr, ref->nrow(), 0, ref->ncol(), ecache);
        EXPECT_EQ(ref_range.value, ecache.cache_values[0]);
        EXPECT_EQ(ref_range.index, ecache.cache_indices[0]);
    }
}

TEST_P(CustomSparseChunkManagerFullTest, Column) {
    auto param = GetParam();
    assemble(std::get<0>(param), std::get<1>(param), std::get<2>(param), std::get<3>(param));

    auto cache = manager.create_chunk_cache<false, false>(ref->nrow());
    auto ecache = manager.create_chunk_cache<false, true>(ref->nrow());
    auto ref_ext = ref->sparse_column();
    int lastc = -1;

    for (int c = 0; c < ref->ncol(); ++c) {
        int requiredc = c / manager.chunk_ncol;
        int offsetc = c % manager.chunk_ncol;
        if (requiredc != lastc) {
            manager.extract<false, false>(requiredc, offsetc, ref->ncol(), 0, ref->nrow(), cache);
            lastc = requiredc;
        }

        auto ref_range = ref_ext->fetch(c);
        EXPECT_EQ(ref_range.value, cache.cache_values[offsetc]);
        EXPECT_EQ(ref_range.index, cache.cache_indices[offsetc]);

        // Testing with exact.
        manager.extract<false, true>(requiredc, offsetc, ref->ncol(), 0, ref->nrow(), ecache);
        EXPECT_EQ(ref_range.value, ecache.cache_values[0]);
        EXPECT_EQ(ref_range.index, ecache.cache_indices[0]);
    }
}

INSTANTIATE_TEST_SUITE_P(
    CustomChunkManager,
    CustomSparseChunkManagerFullTest,
    ::testing::Combine(
        ::testing::Values( // matrix dimensions
            std::make_pair(200, 50),
            std::make_pair(100, 300),
            std::make_pair(152, 211),
            std::make_pair(512, 32)
        ),
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 20),
            std::make_pair(20, 1),
            std::make_pair(10, 10),
            std::make_pair(11, 13) // odd numbers
        ),
        ::testing::Values(true, false), // row major
        ::testing::Values(true, false) // chunk is row major
    )
);

/*******************************************************/

class CustomSparseChunkManagerBlockTest : 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, std::pair<double, double>, bool, bool> >, 
    public CustomSparseChunkManagerMethods {};

TEST_P(CustomSparseChunkManagerBlockTest, Row) {
    auto param = GetParam();
    assemble(std::get<0>(param), std::get<1>(param), std::get<3>(param), std::get<4>(param));

    auto bounds = std::get<2>(param);
    int start = bounds.first * ref->ncol();
    int len = bounds.second * ref->ncol() - start;

    auto cache = manager.create_chunk_cache<true, false>(len);
    auto ecache = manager.create_chunk_cache<true, true>(len);
    auto ref_ext = ref->sparse_row(start, len);
    int lastr = -1;

    for (int r = 0; r < ref->nrow(); ++r) {
        int requiredr = r / manager.chunk_nrow;
        int offsetr = r % manager.chunk_nrow;
        if (requiredr != lastr) {
            manager.extract<true, false>(requiredr, offsetr, ref->nrow(), start, len, cache);
            lastr = requiredr;
        }

        auto ref_range = ref_ext->fetch(r);
        EXPECT_EQ(ref_range.value, cache.cache_values[offsetr]);
        EXPECT_EQ(ref_range.index, cache.cache_indices[offsetr]);

        // Testing with exact.
        manager.extract<true, true>(requiredr, offsetr, ref->nrow(), start, len, ecache);
        EXPECT_EQ(ref_range.value, ecache.cache_values[0]);
        EXPECT_EQ(ref_range.index, ecache.cache_indices[0]);
    }
}

TEST_P(CustomSparseChunkManagerBlockTest, Column) {
    auto param = GetParam();
    assemble(std::get<0>(param), std::get<1>(param), std::get<3>(param), std::get<4>(param));

    auto bounds = std::get<2>(param);
    int start = bounds.first * ref->nrow();
    int len = bounds.second * ref->nrow() - start;

    auto cache = manager.create_chunk_cache<false, false>(len);
    auto ecache = manager.create_chunk_cache<false, true>(len);
    auto ref_ext = ref->sparse_column(start, len);
    int lastc = -1;

    for (int c = 0; c < ref->ncol(); ++c) {
        int requiredc = c / manager.chunk_ncol;
        int offsetc = c % manager.chunk_ncol;
        if (requiredc != lastc) {
            manager.extract<false, false>(requiredc, offsetc, ref->ncol(), start, len, cache);
            lastc = requiredc;
        }

        auto ref_range = ref_ext->fetch(c);
        EXPECT_EQ(ref_range.value, cache.cache_values[offsetc]);
        EXPECT_EQ(ref_range.index, cache.cache_indices[offsetc]);

        // Testing with exact.
        manager.extract<false, true>(requiredc, offsetc, ref->ncol(), start, len, ecache);
        EXPECT_EQ(ref_range.value, ecache.cache_values[0]);
        EXPECT_EQ(ref_range.index, ecache.cache_indices[0]);
    }
}

INSTANTIATE_TEST_SUITE_P(
    CustomChunkManager,
    CustomSparseChunkManagerBlockTest,
    ::testing::Combine(
        ::testing::Values( // matrix dimensions
            std::make_pair(200, 50),
            std::make_pair(100, 300),
            std::make_pair(152, 211)
        ),
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 20),
            std::make_pair(20, 1),
            std::make_pair(10, 10)
        ),
        ::testing::Values( // block boundaries
            std::make_pair(0.0, 0.35),
            std::make_pair(0.15, 0.87),
            std::make_pair(0.38, 1.0)
        ),
        ::testing::Values(true, false), // row major
        ::testing::Values(true, false) // chunk is row major
    )
);

/*******************************************************/

class CustomSparseChunkManagerIndexTest : 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, std::pair<double, double>, bool, bool> >, 
    public CustomSparseChunkManagerMethods 
{
protected:
    static std::vector<int> get_indices(std::pair<double, double> bounds, int range) {
        int start = bounds.first * range;
        int jump = bounds.second;
        std::vector<int> indices;
        while (start < range) {
            indices.push_back(start);
            start += jump;
        }
        return indices;
    }
};

TEST_P(CustomSparseChunkManagerIndexTest, Row) {
    auto param = GetParam();
    assemble(std::get<0>(param), std::get<1>(param), std::get<3>(param), std::get<4>(param));
    auto indices = get_indices(std::get<2>(param), ref->ncol());

    auto cache = manager.create_chunk_cache<true, false>(indices.size());
    auto ecache = manager.create_chunk_cache<true, true>(indices.size());
    auto ref_ext = ref->sparse_row(indices);
    int lastr = -1;

    for (int r = 0; r < ref->nrow(); ++r) {
        int requiredr = r / manager.chunk_nrow;
        int offsetr = r % manager.chunk_nrow;
        if (requiredr != lastr) {
            manager.extract<true, false>(requiredr, offsetr, ref->nrow(), indices, cache);
            lastr = requiredr;
        }

        auto ref_range = ref_ext->fetch(r);
        EXPECT_EQ(ref_range.value, cache.cache_values[offsetr]);
        EXPECT_EQ(ref_range.index, cache.cache_indices[offsetr]);

        // Testing with exact.
        manager.extract<true, true>(requiredr, offsetr, ref->nrow(), indices, ecache);
        EXPECT_EQ(ref_range.value, ecache.cache_values[0]);
        EXPECT_EQ(ref_range.index, ecache.cache_indices[0]);
    }
}

TEST_P(CustomSparseChunkManagerIndexTest, Column) {
    auto param = GetParam();
    assemble(std::get<0>(param), std::get<1>(param), std::get<3>(param), std::get<4>(param));
    auto indices = get_indices(std::get<2>(param), ref->nrow());

    auto cache = manager.create_chunk_cache<false, false>(indices.size());
    auto ecache = manager.create_chunk_cache<false, true>(indices.size());
    std::vector<double> tmp1(indices.size()), tmp2(indices.size());
    auto ref_ext = ref->sparse_column(indices);
    int lastc = -1;

    for (int c = 0; c < ref->ncol(); ++c) {
        int requiredc = c / manager.chunk_ncol;
        int offsetc = c % manager.chunk_ncol;
        if (requiredc != lastc) {
            manager.extract<false, false>(requiredc, offsetc, ref->ncol(), indices, cache);
            lastc = requiredc;
        }

        auto ref_range = ref_ext->fetch(c);
        EXPECT_EQ(ref_range.value, cache.cache_values[offsetc]);
        EXPECT_EQ(ref_range.index, cache.cache_indices[offsetc]);

        // Testing with exact.
        manager.extract<false, true>(requiredc, offsetc, ref->ncol(), indices, ecache);
        EXPECT_EQ(ref_range.value, ecache.cache_values[0]);
        EXPECT_EQ(ref_range.index, ecache.cache_indices[0]);
    }
}

INSTANTIATE_TEST_SUITE_P(
    CustomChunkManager,
    CustomSparseChunkManagerIndexTest,
    ::testing::Combine(
        ::testing::Values( // matrix dimensions
            std::make_pair(200, 50),
            std::make_pair(100, 300),
            std::make_pair(152, 211)
        ),
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 20),
            std::make_pair(20, 1),
            std::make_pair(10, 10)
        ),
        ::testing::Values( // index information.
            std::make_pair(0.0, 10),
            std::make_pair(0.2, 5),
            std::make_pair(0.7, 3)
        ),
        ::testing::Values(true, false), // row major
        ::testing::Values(true, false) // chunk is row major
    )
);

