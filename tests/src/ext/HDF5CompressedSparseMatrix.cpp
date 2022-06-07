#include <gtest/gtest.h>

#ifdef TEST_CUSTOM_PARALLEL // make sure this is included before tatami::apply.
#include "../stats/custom_parallel.h"
#include "hdf5_custom_lock.h"
#endif

#include "H5Cpp.h"
#include "tatami/base/CompressedSparseMatrix.hpp"
#include "tatami/base/DelayedTranspose.hpp"
#include "tatami/ext/HDF5CompressedSparseMatrix.hpp"
#include "tatami/stats/sums.hpp"

#include "temp_file_path.h"
#include <vector>
#include <random>

#include "../_tests/test_column_access.h"
#include "../_tests/test_row_access.h"
#include "../_tests/simulate_vector.h"

class HDF5SparseMatrixTestMethods {
protected:
    std::vector<double> values;
    std::string fpath;
    std::string name;
    SparseDetails<double> triplets;

    void dump(const int& caching, size_t NR, size_t NC) {
        fpath = temp_file_path("tatami-sparse-test.h5");
        H5::H5File fhandle(fpath, H5F_ACC_TRUNC);
        name = "stuff";
        auto ghandle = fhandle.createGroup(name);

        triplets = simulate_sparse_triplets<double>(NR, NC, 0.05, 0, 100);
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

class HDF5SparseSlicedTest : public ::testing::TestWithParam<std::tuple<bool, int, std::vector<int>, int> >, public HDF5SparseMatrixTestMethods {};

TEST_P(HDF5SparseSlicedTest, Primary) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

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

        test_sliced_row_access(&mat, &ref, FORWARD, JUMP, FIRST, LEN, SHIFT);
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

        test_sliced_column_access(&mat, &ref, FORWARD, JUMP, FIRST, LEN, SHIFT);
    }
}

TEST_P(HDF5SparseSlicedTest, Secondary) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

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

        test_sliced_column_access(&mat, &ref, FORWARD, JUMP, FIRST, LEN, SHIFT);
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

        test_sliced_row_access(&mat, &ref, FORWARD, JUMP, FIRST, LEN, SHIFT);
    }
}

INSTANTIATE_TEST_CASE_P(
    HDF5SparseMatrix,
    HDF5SparseSlicedTest,
    ::testing::Combine(
        ::testing::Values(true, false),
        ::testing::Values(1, 3),
        ::testing::Values(
            std::vector<int>({ 0, 8, 3 }), // overlapping shifts
            std::vector<int>({ 1, 4, 4 }), // non-overlapping shifts
            std::vector<int>({ 3, 10, 0 })
        ),
        ::testing::Values(0, 10, 100) // chunk size
    )
);
