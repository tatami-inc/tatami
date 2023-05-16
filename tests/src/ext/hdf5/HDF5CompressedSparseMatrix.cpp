#include <gtest/gtest.h>

#ifdef TEST_CUSTOM_PARALLEL // make sure this is included before tatami::apply.
#include "../../stats/custom_parallel.h"
#include "hdf5_custom_lock.h"
#endif

#include "H5Cpp.h"
#include "tatami/base/sparse/CompressedSparseMatrix.hpp"
#include "tatami/base/other/DelayedTranspose.hpp"
#include "tatami/ext/hdf5/HDF5CompressedSparseMatrix.hpp"
#include "tatami/stats/sums.hpp"

#include "../temp_file_path.h"
#include <vector>
#include <random>

#include "../../_tests/test_column_access.h"
#include "../../_tests/test_row_access.h"
#include "../../_tests/test_oracle_access.h"
#include "../../_tests/simulate_vector.h"

class HDF5SparseMatrixTestMethods {
protected:
    std::vector<double> values;
    std::string fpath;
    std::string name;
    CompressedSparseDetails<double> triplets;

    void dump(int chunk, size_t NR, size_t NC) {
        fpath = temp_file_path("tatami-sparse-test.h5");
        H5::H5File fhandle(fpath, H5F_ACC_TRUNC);
        name = "stuff";
        auto ghandle = fhandle.createGroup(name);

        triplets = simulate_sparse_compressed<double>(NR, NC, 0.05, 0, 100);
        for (auto& v : triplets.value) {
            v = std::round(v);
        }

        H5::DSetCreatPropList plist(H5::DSetCreatPropList::DEFAULT.getId());
        if (chunk == 0) {
            plist.setLayout(H5D_CONTIGUOUS);
        } else {
            plist.setLayout(H5D_CHUNKED);
            hsize_t chunkdim = std::min(triplets.value.size(), static_cast<size_t>(chunk));
            plist.setChunk(1, &chunkdim);
        }

        hsize_t dims = triplets.value.size();
        H5::DataSpace dspace(1, &dims);
        {
            H5::DataType dtype(H5::PredType::NATIVE_UINT8);
            auto dhandle = ghandle.createDataSet("data", dtype, dspace, plist);
            dhandle.write(triplets.value.data(), H5::PredType::NATIVE_DOUBLE);
        }

        {
            H5::DataType dtype(H5::PredType::NATIVE_UINT16);
            auto dhandle = ghandle.createDataSet("index", dtype, dspace, plist);
            dhandle.write(triplets.index.data(), H5::PredType::NATIVE_INT);
        }

        {
            hsize_t ncp1 = triplets.ptr.size();
            H5::DataSpace dspace(1, &ncp1);
            H5::DataType dtype(H5::PredType::NATIVE_UINT64);
            auto dhandle = ghandle.createDataSet("indptr", dtype, dspace);
            dhandle.write(triplets.ptr.data(), H5::PredType::NATIVE_LONG);
        }

        return;
    }

    void dump(size_t NR, size_t NC) {
        dump(50, NR, NC); // just using any chunk size here.
    }

    static size_t compute_cache_size(size_t NR, size_t NC, double fraction) {
        return static_cast<double>(NR * NC) * fraction * static_cast<double>(sizeof(double) + sizeof(int));
    }

    template<bool row_>
    auto create_matrix(size_t NR, size_t NC, size_t cache_size) const {
        return tatami::HDF5CompressedSparseMatrix<row_, double, int>(
            row_ ? NR : NC, 
            row_ ? NC : NR, 
            fpath, 
            name + "/data", 
            name + "/index", 
            name + "/indptr", 
            cache_size
        ); 
    }

    template<bool row_>
    auto create_reference(size_t NR, size_t NC) const {
        return tatami::CompressedSparseMatrix<
            row_, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        >(
            row_ ? NR : NC, 
            row_ ? NC : NR, 
            triplets.value, 
            triplets.index, 
            triplets.ptr
        );
    }
};

/*************************************
 *************************************/

class HDF5SparseUtilsTest : public ::testing::Test, public HDF5SparseMatrixTestMethods {};

TEST_F(HDF5SparseUtilsTest, Basic) {
    const size_t NR = 200, NC = 100;
    dump(NR, NC);

    {
        tatami::HDF5CompressedSparseMatrix<true, double, int> mat(NR, NC, fpath, name + "/data", name + "/index", name + "/indptr");
        EXPECT_EQ(mat.nrow(), NR);
        EXPECT_EQ(mat.ncol(), NC);
        EXPECT_TRUE(mat.sparse());
        EXPECT_EQ(mat.sparse_proportion(), 1);
        EXPECT_TRUE(mat.prefer_rows());
        EXPECT_EQ(mat.prefer_rows_proportion(), 1);
    }

    {
        tatami::HDF5CompressedSparseMatrix<false, double, int> mat(NC, NR, fpath, name + "/data", name + "/index", name + "/indptr");
        EXPECT_EQ(mat.nrow(), NC);
        EXPECT_EQ(mat.ncol(), NR);
        EXPECT_TRUE(mat.sparse());
        EXPECT_EQ(mat.sparse_proportion(), 1);
        EXPECT_FALSE(mat.prefer_rows());
        EXPECT_EQ(mat.prefer_rows_proportion(), 0);
    }
}

/*************************************
 *************************************/

class HDF5SparseAccessTest : public ::testing::TestWithParam<std::tuple<bool, int, int, double> >, public HDF5SparseMatrixTestMethods {};

TEST_P(HDF5SparseAccessTest, Primary) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

    const size_t NR = 200, NC = 100;
    dump(std::get<2>(param), NR, NC);
    int cache_size = compute_cache_size(NR, NC, std::get<3>(param)); // We limit the cache size to ensure that the cache management is not trivial.

    {
        auto mat = create_matrix<true>(NR, NC, cache_size);
        auto ref = create_reference<true>(NR, NC);
        test_simple_row_access(&mat, &ref, FORWARD, JUMP);
    }

    {
        auto mat = create_matrix<false>(NR, NC, cache_size);
        auto ref = create_reference<false>(NR, NC);
        test_simple_column_access(&mat, &ref, FORWARD, JUMP);
    }
}

TEST_P(HDF5SparseAccessTest, Secondary) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

    const size_t NR = 50, NC = 10; // much smaller for the secondary dimension.
    dump(std::get<2>(param), NR, NC);
    int cache_size = compute_cache_size(NR, NC, std::get<3>(param));

    {
        auto mat = create_matrix<true>(NR, NC, cache_size);
        auto ref = create_reference<true>(NR, NC);
        test_simple_column_access(&mat, &ref, FORWARD, JUMP);
    }

    {
        auto mat = create_matrix<false>(NR, NC, cache_size);
        auto ref = create_reference<false>(NR, NC);
        test_simple_row_access(&mat, &ref, FORWARD, JUMP);
    }
}

INSTANTIATE_TEST_CASE_P(
    HDF5SparseMatrix,
    HDF5SparseAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false),
        ::testing::Values(1, 3),
        ::testing::Values(0, 100), // chunk size
        ::testing::Values(0, 0.001, 0.01) // cache size
    )
);

/*************************************
 *************************************/

class HDF5SparseSlicedTest : public ::testing::TestWithParam<std::tuple<bool, int, std::vector<double>, int, double> >, public HDF5SparseMatrixTestMethods {};

TEST_P(HDF5SparseSlicedTest, Primary) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    const size_t NR = 128, NC = 256;
    dump(std::get<3>(param), NR, NC);
    int cache_size = compute_cache_size(NR, NC, std::get<4>(param));

    {
        auto mat = create_matrix<true>(NR, NC, cache_size);
        auto ref = create_reference<true>(NR, NC);
        size_t FIRST = interval_info[0] * NC, LAST = interval_info[1] * NC;
        test_sliced_row_access(&mat, &ref, FORWARD, JUMP, FIRST, LAST);
    }

    {
        auto mat = create_matrix<false>(NR, NC, cache_size);
        auto ref = create_reference<false>(NR, NC);
        size_t FIRST = interval_info[0] * NC, LAST = interval_info[1] * NC; // NC is deliberate, due to the transposition.
        test_sliced_column_access(&mat, &ref, FORWARD, JUMP, FIRST, LAST);
    }
}

TEST_P(HDF5SparseSlicedTest, Secondary) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    const size_t NR = 50, NC = 10; // much smaller for the secondary dimension.
    dump(std::get<3>(param), NR, NC);
    int cache_size = compute_cache_size(NR, NC, std::get<4>(param));

    {
        auto mat = create_matrix<true>(NR, NC, cache_size);
        auto ref = create_reference<true>(NR, NC);
        size_t FIRST = interval_info[0] * NR, LAST = interval_info[1] * NR;
        test_sliced_column_access(&mat, &ref, FORWARD, JUMP, FIRST, LAST);
    }

    {
        auto mat = create_matrix<false>(NR, NC, cache_size);
        auto ref = create_reference<false>(NR, NC);
        size_t FIRST = interval_info[0] * NC, LAST = interval_info[1] * NC;
        test_sliced_row_access(&mat, &ref, FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    HDF5SparseMatrix,
    HDF5SparseSlicedTest,
    ::testing::Combine(
        ::testing::Values(true, false),
        ::testing::Values(1, 3),
        ::testing::Values(
            std::vector<double>({ 0, 0.333 }), 
            std::vector<double>({ 0.222, 0.888 }), 
            std::vector<double>({ 0.555, 1 })
        ),
        ::testing::Values(0, 100), // chunk size
        ::testing::Values(0, 0.001, 0.01, 0.02) // cache size
    )
);

/*************************************
 *************************************/

class HDF5SparseIndexedTest : public ::testing::TestWithParam<std::tuple<bool, int, std::vector<double>, int, double> >, public HDF5SparseMatrixTestMethods {};

TEST_P(HDF5SparseIndexedTest, Primary) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    const size_t NR = 200, NC = 100;
    dump(std::get<3>(param), NR, NC);
    int cache_size = compute_cache_size(NR, NC, std::get<4>(param));

    {
        auto mat = create_matrix<true>(NR, NC, cache_size);
        auto ref = create_reference<true>(NR, NC);
        size_t FIRST = interval_info[0] * NC, STEP = interval_info[1];
        test_indexed_row_access(&mat, &ref, FORWARD, JUMP, FIRST, STEP);
    }

    {
        auto mat = create_matrix<false>(NR, NC, cache_size);
        auto ref = create_reference<false>(NR, NC);
        size_t FIRST = interval_info[0] * NC, STEP = interval_info[1]; // NC is deliberate, due to the transposition.
        test_indexed_column_access(&mat, &ref, FORWARD, JUMP, FIRST, STEP);
    }
}

TEST_P(HDF5SparseIndexedTest, Secondary) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    const size_t NR = 50, NC = 10; // much smaller for the secondary dimension.
    dump(std::get<3>(param), NR, NC);
    int cache_size = compute_cache_size(NR, NC, std::get<4>(param));

    {
        auto mat = create_matrix<true>(NR, NC, cache_size);
        auto ref = create_reference<true>(NR, NC);
        size_t FIRST = interval_info[0] * NR, STEP = interval_info[1];
        test_indexed_column_access(&mat, &ref, FORWARD, JUMP, FIRST, STEP);
    }

    {
        auto mat = create_matrix<false>(NR, NC, cache_size);
        auto ref = create_reference<false>(NR, NC);
        size_t FIRST = interval_info[0] * NC, STEP = interval_info[1];
        test_indexed_row_access(&mat, &ref, FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    HDF5SparseMatrix,
    HDF5SparseIndexedTest,
    ::testing::Combine(
        ::testing::Values(true, false),
        ::testing::Values(1, 3),
        ::testing::Values(
            std::vector<double>({ 0.3, 5 }), 
            std::vector<double>({ 0.322, 8 }), 
            std::vector<double>({ 0.455, 9 })
        ),
        ::testing::Values(0, 100), // chunk size
        ::testing::Values(0, 0.001, 0.01, 0.02) // cache size
    )
);

/*************************************
 *************************************/

class HDF5SparseBasicCacheTest : public ::testing::TestWithParam<double>, public HDF5SparseMatrixTestMethods {};

TEST_P(HDF5SparseBasicCacheTest, LruRandomized) {
    // Check that the LRU cache works as expected with totally random access.
    const size_t NR = 100, NC = 150; 
    dump(NR, NC);
    int cache_size = compute_cache_size(NR, NC, GetParam());

    {
        auto mat = create_matrix<true>(NR, NC, cache_size);
        auto ref = create_reference<true>(NR, NC);

        std::mt19937_64 rng(cache_size * 123);
        auto m_ext = mat.dense_row();
        auto r_ext = ref.dense_row();
        for (size_t r0 = 0, end = NR * 10; r0 < end; ++r0) {
            auto r = rng() % NR;
            EXPECT_EQ(m_ext->fetch(r), r_ext->fetch(r));
        }
    }

    {
        auto mat = create_matrix<false>(NR, NC, cache_size);
        auto ref = create_reference<false>(NR, NC);

        std::mt19937_64 rng(cache_size * 456);
        auto m_ext = mat.dense_column();
        auto r_ext = ref.dense_column();
        for (size_t c0 = 0, end = NR * 10; c0 < end; ++c0) {
            auto c = rng() % NR;
            EXPECT_EQ(m_ext->fetch(c), r_ext->fetch(c));
        }
    }
}

TEST_P(HDF5SparseBasicCacheTest, SimpleOracle) {
    // Checking that access with an oracle behaves as expected.
    const size_t NR = 189, NC = 123; 
    dump(NR, NC);
    int cache_size = compute_cache_size(NR, NC, GetParam());

    {
        auto mat = create_matrix<true>(NR, NC, cache_size);
        auto ref = create_reference<true>(NR, NC);

        test_oracle_row_access<tatami::NumericMatrix>(&mat, &ref, false); // consecutive
        test_oracle_row_access<tatami::NumericMatrix>(&mat, &ref, false, 0.3 * NR, 0.5 * NR); // consecutive with bounds

        test_oracle_row_access<tatami::NumericMatrix>(&mat, &ref, true); // randomized
        test_oracle_row_access<tatami::NumericMatrix>(&mat, &ref, true, 0.2 * NR, 0.6 * NR); // randomized with bounds

        // Oracle-based extraction still works if we turn off value extraction.
        tatami::Options opt;
        opt.sparse_extract_value = false;
        auto mwork = mat.sparse_row(opt);
        mwork->set_oracle(std::make_unique<tatami::ConsecutiveOracle<int> >(0, NR));
        auto rwork = ref.sparse_row(opt);
        for (size_t r = 0; r < NR; ++r) {
            auto mout = mwork->fetch(r);
            auto rout = rwork->fetch(r);
            EXPECT_EQ(mout.index, rout.index);
            EXPECT_EQ(mout.value.size(), 0);
        }
    }

    {
        auto mat = create_matrix<false>(NR, NC, cache_size);
        auto ref = create_reference<false>(NR, NC);

        test_oracle_column_access<tatami::NumericMatrix>(&mat, &ref, false); // consecutive
        test_oracle_column_access<tatami::NumericMatrix>(&mat, &ref, false, 0.1 * NR, 0.7 * NR); // consecutive with bounds

        test_oracle_column_access<tatami::NumericMatrix>(&mat, &ref, true); // randomized
        test_oracle_column_access<tatami::NumericMatrix>(&mat, &ref, true, 0.25 * NR, 0.6 * NR); // randomized with bounds
    }
}

TEST_P(HDF5SparseBasicCacheTest, Repeated) {
    size_t NR = 199, NC = 288;
    dump(NR, NC);
    int cache_size = compute_cache_size(NR, NC, GetParam());

    // Check that we re-use the cache effectively when no new elements are
    // requested; no additional extractions from file should occur.
    std::vector<int> predictions;
    int counter = 0;
    for (size_t i = 0; i < NR * 10; ++i) {
        predictions.push_back(i % 2);
    }

    auto mat = create_matrix<true>(NR, NC, cache_size);
    auto ref = create_reference<true>(NR, NC);

    auto rwork = ref.dense_row();
    auto mwork = mat.dense_row();
    auto mwork_o = mat.dense_row();
    mwork_o->set_oracle(std::make_unique<tatami::FixedOracle<int> >(predictions.data(), predictions.size()));

    for (auto i : predictions) {
        auto expected = rwork->fetch(i);
        EXPECT_EQ(mwork->fetch(i), expected);
        EXPECT_EQ(mwork_o->fetch(i), expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    HDF5SparseMatrix,
    HDF5SparseBasicCacheTest,
    ::testing::Values(0, 0.001, 0.01, 0.1) 
);

/*************************************
 *************************************/

class HDF5SparseReuseCacheTest : public ::testing::TestWithParam<std::tuple<double, int, int> >, public HDF5SparseMatrixTestMethods {
protected:
    size_t NR = 150, NC = 200; 
    std::shared_ptr<tatami::NumericMatrix> mat, ref;
    std::vector<int> predictions;

    template<class Params_>
    void assemble(const Params_& params) {
        dump(NR, NC);

        double cache_multiplier = std::get<0>(params);
        int interval_jump = std::get<1>(params);
        int interval_size = std::get<2>(params);
        int cache_size = compute_cache_size(NR, NC, cache_multiplier);

        mat.reset(new tatami::HDF5CompressedSparseMatrix<true, double, int>(NR, NC, fpath, name + "/data", name + "/index", name + "/indptr", cache_size));
        ref.reset(new tatami::CompressedSparseMatrix<
            true, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        >(NR, NC, triplets.value, triplets.index, triplets.ptr));

        // Repeated scans over the same area, to check for correct re-use of
        // cache elements. We scramble the interval to check that the
        // reordering of elements is done correctly in oracle mode.
        std::vector<int> interval(interval_size);
        std::iota(interval.begin(), interval.end(), 0);
        std::mt19937_64 rng(cache_size + interval_size);

        for (size_t r0 = 0; r0 < NR; r0 += interval_jump) {
            std::shuffle(interval.begin(), interval.end(), rng);
            for (auto i : interval) {
                if (i + r0 < NR) {
                    predictions.push_back(i + r0);
                }
            }
        }
    }
};

TEST_P(HDF5SparseReuseCacheTest, FullExtent) {
    assemble(GetParam());

    auto rwork = ref->dense_row();
    auto mwork = mat->dense_row();
    auto mwork_o = mat->dense_row();
    mwork_o->set_oracle(std::make_unique<tatami::FixedOracle<int> >(predictions.data(), predictions.size()));

    for (auto i : predictions) {
        auto expected = rwork->fetch(i);
        EXPECT_EQ(mwork->fetch(i), expected);
        EXPECT_EQ(mwork_o->fetch(i), expected);
    }
}

TEST_P(HDF5SparseReuseCacheTest, SlicedBounds) {
    assemble(GetParam());

    // Testing that the extraction bounds are actually re-used during extraction.
    // This requires that the elements leave the cache and are then requested again;
    // if it's still in the cache, the non-bounded cache element will just be used.
    auto cstart = NC * 0.25, clen = NC * 0.5;
    tatami::Options opt;
    opt.cache_for_reuse = true;

    auto rwork = ref->dense_row(cstart, clen, opt);
    auto mwork = mat->dense_row(cstart, clen, opt);
    auto mwork_o = mat->dense_row(cstart, clen, opt);

    // Doing one final scan across all rows to re-request elements that have
    // likely left the cache, to force them to be reloaded with bounds.
    for (size_t r = 0; r < NR; ++r) {
        predictions.push_back(r);
    }
    mwork_o->set_oracle(std::make_unique<tatami::FixedOracle<int> >(predictions.data(), predictions.size()));

    for (auto i : predictions) {
        auto expected = rwork->fetch(i);
        EXPECT_EQ(mwork->fetch(i), expected);
        EXPECT_EQ(mwork_o->fetch(i), expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    HDF5SparseMatrix,
    HDF5SparseReuseCacheTest,
    ::testing::Combine(
        ::testing::Values(0, 0.001, 0.01, 0.1), // cache size multiplier
        ::testing::Values(1, 3), // jump between intervals
        ::testing::Values(5, 10, 20) // reuse interval size
    )
);

/*************************************
 *************************************/

class HDF5SparseApplyTest : public ::testing::TestWithParam<double>, public HDF5SparseMatrixTestMethods {};

TEST_P(HDF5SparseApplyTest, Basic) {
    // Just putting it through its paces for correct oracular parallelization via apply.
    size_t NR = 500;
    size_t NC = 200;
    dump(NR, NC);
    int cache_size = compute_cache_size(NR, NC, GetParam());

    auto mat = create_matrix<true>(NR, NC, cache_size);
    auto ref = create_reference<true>(NR, NC);

    EXPECT_EQ(tatami::row_sums(&mat), tatami::row_sums(&ref));
    EXPECT_EQ(tatami::column_sums(&mat), tatami::column_sums(&ref));
}

INSTANTIATE_TEST_CASE_P(
    HDF5SparseMatrix,
    HDF5SparseApplyTest,
    ::testing::Values(0, 0.001, 0.01, 0.1) // cache size multiplier
);
