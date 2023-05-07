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

    void dump(const int& caching, size_t NR, size_t NC) {
        fpath = temp_file_path("tatami-sparse-test.h5");
        H5::H5File fhandle(fpath, H5F_ACC_TRUNC);
        name = "stuff";
        auto ghandle = fhandle.createGroup(name);

        triplets = simulate_sparse_compressed<double>(NR, NC, 0.05, 0, 100);
        for (auto& v : triplets.value) {
            v = std::round(v);
        }

        H5::DSetCreatPropList plist(H5::DSetCreatPropList::DEFAULT.getId());
        if (caching == 0) {
            plist.setLayout(H5D_CONTIGUOUS);
        } else {
            plist.setLayout(H5D_CHUNKED);
            hsize_t chunkdim = std::min(triplets.value.size(), static_cast<size_t>(caching));
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
};

/*************************************
 *************************************/

class HDF5SparseUtilsTest : public ::testing::Test, public HDF5SparseMatrixTestMethods {};

TEST_F(HDF5SparseUtilsTest, Basic) {
    const size_t NR = 200, NC = 100;
    dump(50, NR, NC);

    {
        tatami::HDF5CompressedSparseMatrix<true, double, int> mat(NR, NC, fpath, name + "/data", name + "/index", name + "/indptr");
        EXPECT_EQ(mat.nrow(), NR);
        EXPECT_EQ(mat.ncol(), NC);
        EXPECT_TRUE(mat.sparse());
    }

    {
        tatami::HDF5CompressedSparseMatrix<false, double, int> mat(NC, NR, fpath, name + "/data", name + "/index", name + "/indptr");
        EXPECT_EQ(mat.nrow(), NC);
        EXPECT_EQ(mat.ncol(), NR);
        EXPECT_TRUE(mat.sparse());
    }
}

/*************************************
 *************************************/

class HDF5SparseAccessTest : public ::testing::TestWithParam<std::tuple<bool, int, int> >, public HDF5SparseMatrixTestMethods {};

TEST_P(HDF5SparseAccessTest, Primary) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

    auto caching = std::get<2>(param);
    const size_t NR = 200, NC = 100;
    dump(caching, NR, NC);

    {
        // We limit the cache size to ensure that the cache management is not trivial.
        tatami::HDF5CompressedSparseMatrix<true, double, int> mat(NR, NC, fpath, name + "/data", name + "/index", name + "/indptr", NR * 2); 
        tatami::CompressedSparseMatrix<
            true, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NR, NC, triplets.value, triplets.index, triplets.ptr);

        test_simple_row_access(&mat, &ref, FORWARD, JUMP);
    }

    {
        tatami::HDF5CompressedSparseMatrix<false, double, int> mat(NC, NR, fpath, name + "/data", name + "/index", name + "/indptr", NC * 3);
        tatami::CompressedSparseMatrix<
            false, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NC, NR, triplets.value, triplets.index, triplets.ptr);

        test_simple_column_access(&mat, &ref, FORWARD, JUMP);
    }
}

TEST_P(HDF5SparseAccessTest, Secondary) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

    auto caching = std::get<2>(param);
    const size_t NR = 50, NC = 10; // much smaller for the secondary dimension.
    dump(caching, NR, NC);

    {
        tatami::HDF5CompressedSparseMatrix<true, double, int> mat(NR, NC, fpath, name + "/data", name + "/index", name + "/indptr", NC * 2);
        tatami::CompressedSparseMatrix<
            true, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NR, NC, triplets.value, triplets.index, triplets.ptr);

        test_simple_column_access(&mat, &ref, FORWARD, JUMP);
    }

    {
        tatami::HDF5CompressedSparseMatrix<false, double, int> mat(NC, NR, fpath, name + "/data", name + "/index", name + "/indptr", NR * 1.5);
        tatami::CompressedSparseMatrix<
            false, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NC, NR, triplets.value, triplets.index, triplets.ptr);

        test_simple_column_access(&mat, &ref, FORWARD, JUMP);
    }
}

TEST_P(HDF5SparseAccessTest, Apply) {
    // Just putting it through its paces for correct parallelization via apply.
    size_t NR = 500;
    size_t NC = 200;
    dump(10, NR, NC);

    tatami::HDF5CompressedSparseMatrix<true, double, int> mat(NR, NC, fpath, name + "/data", name + "/index", name + "/indptr", NR * 1.5);
    tatami::CompressedSparseMatrix<
        true, 
        double, 
        int, 
        decltype(triplets.value), 
        decltype(triplets.index), 
        decltype(triplets.ptr)
    > ref(NR, NC, triplets.value, triplets.index, triplets.ptr);

    EXPECT_EQ(tatami::row_sums(&mat), tatami::row_sums(&ref));
    EXPECT_EQ(tatami::column_sums(&mat), tatami::column_sums(&ref));
}

INSTANTIATE_TEST_CASE_P(
    HDF5SparseMatrix,
    HDF5SparseAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false),
        ::testing::Values(1, 3),
        ::testing::Values(0, 10, 100)
    )
);

/*************************************
 *************************************/

class HDF5SparseSlicedTest : public ::testing::TestWithParam<std::tuple<bool, int, std::vector<double>, int> >, public HDF5SparseMatrixTestMethods {};

TEST_P(HDF5SparseSlicedTest, Primary) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    auto caching = std::get<3>(param);
    const size_t NR = 200, NC = 100;
    dump(caching, NR, NC);

    {
        tatami::HDF5CompressedSparseMatrix<true, double, int> mat(NR, NC, fpath, name + "/data", name + "/index", name + "/indptr", NR * 1.5);
        tatami::CompressedSparseMatrix<
            true, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NR, NC, triplets.value, triplets.index, triplets.ptr);

        size_t FIRST = interval_info[0] * NC, LAST = interval_info[1] * NC;
        test_sliced_row_access(&mat, &ref, FORWARD, JUMP, FIRST, LAST);
    }

    {
        tatami::HDF5CompressedSparseMatrix<false, double, int> mat(NC, NR, fpath, name + "/data", name + "/index", name + "/indptr", NC * 1);
        tatami::CompressedSparseMatrix<
            false, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NC, NR, triplets.value, triplets.index, triplets.ptr);

        size_t FIRST = interval_info[0] * NC, LAST = interval_info[1] * NC; // NC is deliberate, due to the transposition.
        test_sliced_column_access(&mat, &ref, FORWARD, JUMP, FIRST, LAST);
    }
}

TEST_P(HDF5SparseSlicedTest, Secondary) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    auto caching = std::get<3>(param);
    const size_t NR = 50, NC = 10; // much smaller for the secondary dimension.
    dump(caching, NR, NC);

    {
        tatami::HDF5CompressedSparseMatrix<true, double, int> mat(NR, NC, fpath, name + "/data", name + "/index", name + "/indptr", NC * 1.5);
        tatami::CompressedSparseMatrix<
            true, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NR, NC, triplets.value, triplets.index, triplets.ptr);

        size_t FIRST = interval_info[0] * NR, LAST = interval_info[1] * NR;
        test_sliced_column_access(&mat, &ref, FORWARD, JUMP, FIRST, LAST);
    }

    {
        tatami::HDF5CompressedSparseMatrix<false, double, int> mat(NC, NR, fpath, name + "/data", name + "/index", name + "/indptr", NC * 2);
        tatami::CompressedSparseMatrix<
            false, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NC, NR, triplets.value, triplets.index, triplets.ptr);

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
        ::testing::Values(0, 10, 100) // chunk size
    )
);

/*************************************
 *************************************/

class HDF5SparseIndexedTest : public ::testing::TestWithParam<std::tuple<bool, int, std::vector<double>, int> >, public HDF5SparseMatrixTestMethods {};

TEST_P(HDF5SparseIndexedTest, Primary) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    auto caching = std::get<3>(param);
    const size_t NR = 200, NC = 100;
    dump(caching, NR, NC);

    {
        tatami::HDF5CompressedSparseMatrix<true, double, int> mat(NR, NC, fpath, name + "/data", name + "/index", name + "/indptr", NR * 1.5);
        tatami::CompressedSparseMatrix<
            true, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NR, NC, triplets.value, triplets.index, triplets.ptr);

        size_t FIRST = interval_info[0] * NC, STEP = interval_info[1];
        test_indexed_row_access(&mat, &ref, FORWARD, JUMP, FIRST, STEP);
    }

    {
        tatami::HDF5CompressedSparseMatrix<false, double, int> mat(NC, NR, fpath, name + "/data", name + "/index", name + "/indptr", NC * 1);
        tatami::CompressedSparseMatrix<
            false, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NC, NR, triplets.value, triplets.index, triplets.ptr);

        size_t FIRST = interval_info[0] * NC, STEP = interval_info[1]; // NC is deliberate, due to the transposition.
        test_indexed_column_access(&mat, &ref, FORWARD, JUMP, FIRST, STEP);
    }
}

TEST_P(HDF5SparseIndexedTest, Secondary) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    auto caching = std::get<3>(param);
    const size_t NR = 50, NC = 10; // much smaller for the secondary dimension.
    dump(caching, NR, NC);

    {
        tatami::HDF5CompressedSparseMatrix<true, double, int> mat(NR, NC, fpath, name + "/data", name + "/index", name + "/indptr", NC * 1.5);
        tatami::CompressedSparseMatrix<
            true, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NR, NC, triplets.value, triplets.index, triplets.ptr);

        size_t FIRST = interval_info[0] * NR, STEP = interval_info[1];
        test_indexed_column_access(&mat, &ref, FORWARD, JUMP, FIRST, STEP);
    }

    {
        tatami::HDF5CompressedSparseMatrix<false, double, int> mat(NC, NR, fpath, name + "/data", name + "/index", name + "/indptr", NC * 2);
        tatami::CompressedSparseMatrix<
            false, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NC, NR, triplets.value, triplets.index, triplets.ptr);

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
        ::testing::Values(0, 10, 100) // chunk size
    )
);

/*************************************
 *************************************/

class HDF5SparseBasicCacheTest : public ::testing::TestWithParam<double>, public HDF5SparseMatrixTestMethods {};

TEST_P(HDF5SparseBasicCacheTest, LruRandomized) {
    // Check that the LRU cache works as expected with totally random access.
    const size_t NR = 100, NC = 150; 
    dump(100, NR, NC);

    int cache_multiplier = GetParam();
    int cache_size = static_cast<double>(NR * NC) * cache_multiplier * static_cast<double>(sizeof(double) + sizeof(int));

    {
        tatami::HDF5CompressedSparseMatrix<true, double, int> mat(NR, NC, fpath, name + "/data", name + "/index", name + "/indptr", cache_size);
        tatami::CompressedSparseMatrix<
            true, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NR, NC, triplets.value, triplets.index, triplets.ptr);

        std::mt19937_64 rng(cache_size * 123);
        auto m_ext = mat.dense_row();
        auto r_ext = ref.dense_row();
        for (size_t r0 = 0, end = NR * 10; r0 < end; ++r0) {
            auto r = rng() % NR;
            EXPECT_EQ(m_ext->fetch(r), r_ext->fetch(r));
        }
    }

    {
        tatami::HDF5CompressedSparseMatrix<false, double, int> mat(NC, NR, fpath, name + "/data", name + "/index", name + "/indptr", cache_size);
        tatami::CompressedSparseMatrix<
            false, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NC, NR, triplets.value, triplets.index, triplets.ptr);

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
    dump(100, NR, NC);

    int cache_multiplier = GetParam();
    int cache_size = static_cast<double>(NR * NC) * cache_multiplier * static_cast<double>(sizeof(double) + sizeof(int));

    {
        tatami::HDF5CompressedSparseMatrix<true, double, int> mat(NR, NC, fpath, name + "/data", name + "/index", name + "/indptr", cache_size);
        tatami::CompressedSparseMatrix<
            true, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NR, NC, triplets.value, triplets.index, triplets.ptr);

        test_oracle_row_access<tatami::NumericMatrix>(&mat, &ref, false); // consecutive
        test_oracle_row_access<tatami::NumericMatrix>(&mat, &ref, true); // randomized
    }

    {
        tatami::HDF5CompressedSparseMatrix<false, double, int> mat(NC, NR, fpath, name + "/data", name + "/index", name + "/indptr", cache_size);
        tatami::CompressedSparseMatrix<
            false, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NC, NR, triplets.value, triplets.index, triplets.ptr);

        test_oracle_column_access<tatami::NumericMatrix>(&mat, &ref, false); // consecutive
        test_oracle_column_access<tatami::NumericMatrix>(&mat, &ref, true); // randomized
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
        dump(100, NR, NC);

        double cache_multiplier = std::get<0>(params);
        int interval_jump = std::get<1>(params);
        int interval_size = std::get<2>(params);
        int cache_size = static_cast<double>(NR * NC) * cache_multiplier * static_cast<double>(sizeof(double) + sizeof(int));

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
