#include <gtest/gtest.h>
#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/CompressedSparseMatrix.hpp"
#include "tatami/chunked/utils.hpp"
#include "tatami_test/tatami_test.hpp"
#include "mock_chunk.h"

class SimpleChunkWrapperTestMethods {
protected:
    std::unique_ptr<tatami::Matrix<double, int> > ref;

protected:
    tatami::SimpleDenseChunkWrapper<MockDenseChunk<true> > drm_chunk;
    tatami::SimpleDenseChunkWrapper<MockDenseChunk<false> > dcm_chunk;

    void dense_assemble(const std::pair<int, int>& dim) {
        auto full = tatami_test::simulate_dense_vector<double>(dim.first * dim.second, -10, 10, /* seed = */ dim.first * dim.second);
        drm_chunk = tatami::SimpleDenseChunkWrapper<MockDenseChunk<true> >(MockDenseChunk<true>(dim.first, dim.second, full));
        ref.reset(new tatami::DenseRowMatrix<double, int>(dim.first, dim.second, std::move(full)));

        std::vector<double> replacement(dim.first * dim.second);
        auto ext = ref->dense_column();
        for (int c = 0; c < dim.second; ++c) {
            ext->fetch_copy(c, replacement.data() + c * dim.first);
        }
        dcm_chunk = tatami::SimpleDenseChunkWrapper<MockDenseChunk<false> >(MockDenseChunk<false>(dim.first, dim.second, std::move(replacement)));
    }

protected:
    tatami::SimpleSparseChunkWrapper<MockSparseChunk<true> > srm_chunk;
    tatami::SimpleSparseChunkWrapper<MockSparseChunk<false> > scm_chunk;

    void sparse_assemble(const std::pair<int, int>& dim) {
        auto full = tatami_test::simulate_sparse_compressed<double>(dim.first, dim.second, 0.1, -10, 10, /* seed = */ dim.first * dim.second * 2);
        srm_chunk = tatami::SimpleSparseChunkWrapper<MockSparseChunk<true> >(MockSparseChunk<true>(dim.first, dim.second, full.value, full.index, full.ptr));
        ref.reset(new tatami::CompressedSparseRowMatrix<double, int>(dim.first, dim.second, std::move(full.value), std::move(full.index), std::move(full.ptr)));

        std::vector<double> replacement_values;
        std::vector<int> replacement_indices;
        std::vector<size_t> replacement_indptrs(1);
        auto ext = ref->sparse_column();
        for (int c = 0; c < dim.second; ++c) {
            auto out = ext->fetch(c);
            replacement_values.insert(replacement_values.end(), out.value.begin(), out.value.end());
            replacement_indices.insert(replacement_indices.end(), out.index.begin(), out.index.end());
            replacement_indptrs.push_back(replacement_indptrs.back() + out.number);
        }
        scm_chunk = tatami::SimpleSparseChunkWrapper<MockSparseChunk<false> >(MockSparseChunk<false>(dim.first, dim.second, std::move(replacement_values), std::move(replacement_indices), std::move(replacement_indptrs)));
    }
};

/*************************************/

class SimpleChunkWrapperFullTest : public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, double> > >, public SimpleChunkWrapperTestMethods {
protected:
    template<bool sparse_, typename Chunk_>
    void run_tests(const Chunk_& chunk, const std::pair<int, int>& dim, const std::pair<double, double>& bounds) const {
        auto work = chunk.workspace();
        int r_first = bounds.first * dim.first, r_last = bounds.second * dim.first, r_len = r_last - r_first;
        int c_first = bounds.first * dim.second, c_last = bounds.second * dim.second, c_len = c_last - c_first;

        // Row-major row extraction.
        if (r_len) {
            if constexpr(sparse_) {
                std::vector<std::vector<double> > out_values(r_len);
                std::vector<std::vector<int> > out_indices(r_len);
                chunk.template extract<true>(r_first, r_len, 0, dim.second, work, out_values, out_indices, 10);

                auto ext = this->ref->sparse_row();
                for (int r = r_first; r < r_last; ++r) {
                    auto refrow = ext->fetch(r);
                    EXPECT_EQ(refrow.value, out_values[r - r_first]);
                    for (auto& i : refrow.index) {
                        i += 10;
                    }
                    EXPECT_EQ(refrow.index, out_indices[r - r_first]);
                }

            } else {
                int stride = 100;
                std::vector<double> buffer(r_len * stride);
                chunk.template extract<true>(r_first, r_len, 0, dim.second, work, buffer.data(), stride);

                auto ext = this->ref->dense_row();
                for (int r = r_first; r < r_last; ++r) {
                    auto refrow = ext->fetch(r);
                    auto bptr = buffer.data() + static_cast<size_t>(r - r_first) * stride;
                    EXPECT_EQ(refrow, std::vector<double>(bptr, bptr + dim.second));
                }
            }
        }

        // Row-major column extraction. 
        if (c_len) {
            if constexpr(sparse_) {
                std::vector<std::vector<double> > out_values(c_len);
                std::vector<std::vector<int> > out_indices(c_len);
                chunk.template extract<false>(c_first, c_len, 0, dim.first, work, out_values, out_indices, 5);

                auto ext = this->ref->sparse_column();
                for (int c = c_first; c < c_last; ++c) {
                    auto refcol = ext->fetch(c);
                    EXPECT_EQ(refcol.value, out_values[c - c_first]);
                    for (auto& i : refcol.index) {
                        i += 5;
                    }
                    EXPECT_EQ(refcol.index, out_indices[c - c_first]);
                }

            } else {
                int stride = 98;
                std::vector<double> buffer(c_len * stride);
                chunk.template extract<false>(c_first, c_len, 0, dim.first, work, buffer.data(), stride);

                auto ext = this->ref->dense_column();
                for (int c = c_first; c < c_last; ++c) {
                    auto refcol = ext->fetch(c);
                    auto bptr = buffer.data() + static_cast<size_t>(c - c_first) * stride;
                    EXPECT_EQ(refcol, std::vector<double>(bptr, bptr + dim.first));
                }
            }
        }
    }
};

TEST_P(SimpleChunkWrapperFullTest, Dense) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    dense_assemble(dim);

    auto bounds = std::get<1>(params);
    run_tests<false>(drm_chunk, dim, bounds);
    run_tests<false>(dcm_chunk, dim, bounds);
}

TEST_P(SimpleChunkWrapperFullTest, Sparse) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    sparse_assemble(dim);

    auto bounds = std::get<1>(params);
    run_tests<true>(srm_chunk, dim, bounds);
    run_tests<true>(scm_chunk, dim, bounds);
}

INSTANTIATE_TEST_SUITE_P(
    SimpleChunkWrapper,
    SimpleChunkWrapperFullTest,
    ::testing::Combine(
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 50),
            std::make_pair(50, 1),
            std::make_pair(13, 42),
            std::make_pair(27, 23) 
        ),
        ::testing::Values( // primary range
            std::make_pair(0.0, 1.0),
            std::make_pair(0.2, 0.8),
            std::make_pair(0.4, 0.6)
        )
    )
);

/*************************************/

class SimpleChunkWrapperBlockTest : public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, double>, std::pair<double, double> > >, public SimpleChunkWrapperTestMethods {
protected:
    template<bool sparse_, typename Chunk_>
    void run_tests(const Chunk_& chunk, const std::pair<int, int>& dim, const std::pair<double, double>& bounds, const std::pair<double, double>& block) const {
        int r_first = bounds.first * dim.first,  r_last = bounds.second * dim.first,  r_len = r_last - r_first;
        int c_first = bounds.first * dim.second, c_last = bounds.second * dim.second, c_len = c_last - c_first;

        int r2_first = block.first * dim.first,  r2_last = block.second * dim.first,  r2_len = r2_last - r2_first;
        int c2_first = block.first * dim.second, c2_last = block.second * dim.second, c2_len = c2_last - c2_first;

        auto work = chunk.workspace();

        // Row-major row extraction.
        if (r_len && c2_len) {
            if constexpr(sparse_) {
                std::vector<std::vector<double> > out_values(r_len);
                std::vector<std::vector<int> > out_indices(r_len);
                chunk.template extract<true>(r_first, r_len, c2_first, c2_len, work, out_values, out_indices, 7);

                auto ext = this->ref->sparse_row(c2_first, c2_len);
                for (int r = r_first; r < r_last; ++r) {
                    auto refrow = ext->fetch(r);
                    EXPECT_EQ(refrow.value, out_values[r - r_first]);
                    for (auto& i : refrow.index) {
                        i += 7;
                    }
                    EXPECT_EQ(refrow.index, out_indices[r - r_first]);
                }

            } else {
                int stride = 100;
                std::vector<double> buffer(r_len * stride);
                chunk.template extract<true>(r_first, r_len, c2_first, c2_len, work, buffer.data(), stride);

                auto ext = ref->dense_row(c2_first, c2_len);
                for (int r = r_first; r < r_last; ++r) {
                    auto refrow = ext->fetch(r);
                    auto bptr = buffer.data() + static_cast<size_t>(r - r_first) * stride;
                    EXPECT_EQ(refrow, std::vector<double>(bptr, bptr + c2_len));
                }
            }
        }

        // Row-major column extraction. 
        if (c_len && r2_len) {
            if constexpr(sparse_) {
                std::vector<std::vector<double> > out_values(c_len);
                std::vector<std::vector<int> > out_indices(c_len);
                chunk.template extract<false>(c_first, c_len, r2_first, r2_len, work, out_values, out_indices, 8);

                auto ext = this->ref->sparse_column(r2_first, r2_len);
                for (int c = c_first; c < c_last; ++c) {
                    auto refcol = ext->fetch(c);
                    EXPECT_EQ(refcol.value, out_values[c - c_first]);
                    for (auto& i : refcol.index) {
                        i += 8;
                    }
                    EXPECT_EQ(refcol.index, out_indices[c - c_first]);
                }

            } else {
                int stride = 98;
                std::vector<double> buffer(c_len * stride);
                chunk.template extract<false>(c_first, c_len, r2_first, r2_len, work, buffer.data(), stride);

                auto ext = ref->dense_column(r2_first, r2_len);
                for (int c = c_first; c < c_last; ++c) {
                    auto refcol = ext->fetch(c);
                    auto bptr = buffer.data() + static_cast<size_t>(c - c_first) * stride;
                    EXPECT_EQ(refcol, std::vector<double>(bptr, bptr + r2_len));
                }
            }
        }
    }
};

TEST_P(SimpleChunkWrapperBlockTest, Dense) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    dense_assemble(dim);

    auto bounds = std::get<1>(params);
    auto block = std::get<2>(params);
    run_tests<false>(drm_chunk, dim, bounds, block);
    run_tests<false>(dcm_chunk, dim, bounds, block);
}

TEST_P(SimpleChunkWrapperBlockTest, Sparse) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    sparse_assemble(dim);

    auto bounds = std::get<1>(params);
    auto block = std::get<2>(params);
    run_tests<true>(srm_chunk, dim, bounds, block);
    run_tests<true>(scm_chunk, dim, bounds, block);
}

INSTANTIATE_TEST_SUITE_P(
    SimpleChunkWrapper,
    SimpleChunkWrapperBlockTest,
    ::testing::Combine(
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 50),
            std::make_pair(50, 1),
            std::make_pair(13, 42),
            std::make_pair(27, 23) 
        ),
        ::testing::Values( // primary range
            std::make_pair(0.0, 1.0),
            std::make_pair(0.2, 0.8),
            std::make_pair(0.4, 0.6)
        ),
        ::testing::Values( // secondary block
            std::make_pair(0.0, 1.0),
            std::make_pair(0.25, 0.75),
            std::make_pair(0.33, 0.66)
        )
    )
);

/*************************************/

class SimpleChunkWrapperIndexTest : public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, double>, std::pair<double, double> > >, public SimpleChunkWrapperTestMethods {
protected:
    template<bool sparse_, typename Chunk_>
    void run_tests(const Chunk_& chunk, const std::pair<int, int>& dim, const std::pair<double, double>& bounds, const std::pair<double, double>& block) const {
        int r_first = bounds.first * dim.first,  r_last = bounds.second * dim.first,  r_len = r_last - r_first;
        int c_first = bounds.first * dim.second, c_last = bounds.second * dim.second, c_len = c_last - c_first;

        std::vector<int> r2_indices;
        int r2_start = block.first * dim.first, r2_step = block.second;
        while (r2_start < dim.first) {
            r2_indices.push_back(r2_start);
            r2_start += r2_step;
        } 

        std::vector<int> c2_indices;
        int c2_start = block.first * dim.second, c2_step = block.second;
        while (c2_start < dim.second) {
            c2_indices.push_back(c2_start);
            c2_start += c2_step;
        } 

        auto work = chunk.workspace();

        // Row-major row extraction.
        if (r_len && c2_indices.size()) {
            if constexpr(sparse_) {
                std::vector<std::vector<double> > out_values(r_len);
                std::vector<std::vector<int> > out_indices(r_len);
                chunk.template extract<true>(r_first, r_len, c2_indices, work, out_values, out_indices, 2);

                auto ext = this->ref->sparse_row(c2_indices);
                for (int r = r_first; r < r_last; ++r) {
                    auto refrow = ext->fetch(r);
                    EXPECT_EQ(refrow.value, out_values[r - r_first]);
                    for (auto& i : refrow.index) {
                        i += 2;
                    }
                    EXPECT_EQ(refrow.index, out_indices[r - r_first]);
                }

            } else {
                int stride = 100;
                std::vector<double> buffer(r_len * stride);
                chunk.template extract<true>(r_first, r_len, c2_indices, work, buffer.data(), stride);

                auto ext = ref->dense_row(c2_indices);
                for (int r = r_first; r < r_last; ++r) {
                    auto refrow = ext->fetch(r);
                    auto bptr = buffer.data() + static_cast<size_t>(r - r_first) * stride;
                    EXPECT_EQ(refrow, std::vector<double>(bptr, bptr + c2_indices.size()));
                }
            }
        }

        // Row-major column extraction. 
        if (c_len && r2_indices.size()) {
            if constexpr(sparse_) {
                std::vector<std::vector<double> > out_values(c_len);
                std::vector<std::vector<int> > out_indices(c_len);
                chunk.template extract<false>(c_first, c_len, r2_indices, work, out_values, out_indices, 13);

                auto ext = this->ref->sparse_column(r2_indices);
                for (int c = c_first; c < c_last; ++c) {
                    auto refcol = ext->fetch(c);
                    EXPECT_EQ(refcol.value, out_values[c - c_first]);
                    for (auto& i : refcol.index) {
                        i += 13;
                    }
                    EXPECT_EQ(refcol.index, out_indices[c - c_first]);
                }

            } else {
                int stride = 98;
                std::vector<double> buffer(c_len * stride);
                chunk.template extract<false>(c_first, c_len, r2_indices, work, buffer.data(), stride);

                auto ext = ref->dense_column(r2_indices);
                for (int c = c_first; c < c_last; ++c) {
                    auto refcol = ext->fetch(c);
                    auto bptr = buffer.data() + static_cast<size_t>(c - c_first) * stride;
                    EXPECT_EQ(refcol, std::vector<double>(bptr, bptr + r2_indices.size()));
                }
            }
        }
    }
};

TEST_P(SimpleChunkWrapperIndexTest, Dense) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    dense_assemble(dim);

    auto bounds = std::get<1>(params);
    auto block = std::get<2>(params);
    run_tests<false>(drm_chunk, dim, bounds, block);
    run_tests<false>(dcm_chunk, dim, bounds, block);
}

TEST_P(SimpleChunkWrapperIndexTest, Sparse) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    sparse_assemble(dim);

    auto bounds = std::get<1>(params);
    auto block = std::get<2>(params);
    run_tests<true>(srm_chunk, dim, bounds, block);
    run_tests<true>(scm_chunk, dim, bounds, block);
}

INSTANTIATE_TEST_SUITE_P(
    SimpleChunkWrapper,
    SimpleChunkWrapperIndexTest,
    ::testing::Combine(
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 50),
            std::make_pair(50, 1),
            std::make_pair(13, 42),
            std::make_pair(27, 23) 
        ),
        ::testing::Values( // primary range
            std::make_pair(0.0, 1.0),
            std::make_pair(0.2, 0.8),
            std::make_pair(0.4, 0.6)
        ),
        ::testing::Values( // secondary indices
            std::make_pair(0.0, 10),
            std::make_pair(0.25, 5),
            std::make_pair(0.33, 3)
        )
    )
);

