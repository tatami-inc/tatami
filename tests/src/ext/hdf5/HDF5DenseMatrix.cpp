#include <gtest/gtest.h>

#ifdef TEST_CUSTOM_PARALLEL // make sure this is included before tatami::apply.
#include "../../stats/custom_parallel.h"
#include "hdf5_custom_lock.h"
std::mutex hdf5_lock; // declared in hdf5_custom_lock.h, defined here.
#endif

#include "H5Cpp.h"
#include "tatami/base/dense/DenseMatrix.hpp"
#include "tatami/base/other/DelayedTranspose.hpp"
#include "tatami/ext/hdf5/HDF5DenseMatrix.hpp"
#include "tatami/stats/sums.hpp"

#include "../temp_file_path.h"
#include <vector>
#include <random>

#include "../../_tests/test_column_access.h"
#include "../../_tests/test_row_access.h"
#include "../../_tests/simulate_vector.h"

const size_t NR = 200, NC = 100;

class HDF5DenseMatrixTestMethods {
protected:
    std::vector<double> values;
    std::string fpath;
    std::string name;

    void dump(const std::pair<int, int>& chunk_sizes) {
        fpath = temp_file_path("tatami-dense-test.h5");
        name = "stuff";
        H5::H5File fhandle(fpath, H5F_ACC_TRUNC);

        hsize_t dims[2];
        dims[0] = NR;
        dims[1] = NC;
        H5::DataSpace dspace(2, dims);
        H5::DataType dtype(H5::PredType::NATIVE_UINT8);

        H5::DSetCreatPropList plist(H5::DSetCreatPropList::DEFAULT.getId());
        if (chunk_sizes.first == 0) {
            plist.setLayout(H5D_CONTIGUOUS);
        } else {
            plist.setLayout(H5D_CHUNKED);
            hsize_t chunkdims[2];
            chunkdims[0] = chunk_sizes.first;
            chunkdims[1] = chunk_sizes.second;
            plist.setChunk(2, chunkdims);
        }

        auto dhandle = fhandle.createDataSet(name, dtype, dspace, plist);
        values = simulate_dense_vector<double>(NR * NC, 0, 100);
        for (auto& v : values) {
            v = std::round(v);
        }

        dhandle.write(values.data(), H5::PredType::NATIVE_DOUBLE);
        return;
    }
};

/*************************************
 *************************************/

class HDF5DenseUtilsTest : public ::testing::Test, public HDF5DenseMatrixTestMethods {};

TEST_F(HDF5DenseUtilsTest, Basic) {
    dump(std::make_pair<int, int>(10, 10));
    tatami::HDF5DenseMatrix<double, int> mat(fpath, name);
    EXPECT_EQ(mat.nrow(), NR);
    EXPECT_EQ(mat.ncol(), NC);
    EXPECT_FALSE(mat.sparse());
}

TEST_F(HDF5DenseUtilsTest, Preference) {
    {
        dump(std::make_pair<int, int>(10, 10));
        tatami::HDF5DenseMatrix<double, int> mat(fpath, name);
        EXPECT_TRUE(mat.prefer_rows());

        tatami::HDF5DenseMatrix<double, int, true> tmat(fpath, name);
        EXPECT_FALSE(tmat.prefer_rows());
    }

    {
        // First dimension is compromised, switching to the second dimension.
        dump(std::make_pair<int, int>(NR, 1));
        tatami::HDF5DenseMatrix<double, int> mat(fpath, name, NR * sizeof(double));
        EXPECT_FALSE(mat.prefer_rows());

        tatami::HDF5DenseMatrix<double, int, true> tmat(fpath, name, NR * sizeof(double));
        EXPECT_TRUE(tmat.prefer_rows());
    }

    {
        // Second dimension is compromised, but we just use the first anyway.
        dump(std::make_pair<int, int>(1, NC));
        tatami::HDF5DenseMatrix<double, int> mat(fpath, name, NC * sizeof(double));
        EXPECT_TRUE(mat.prefer_rows());

        tatami::HDF5DenseMatrix<double, int, true> tmat(fpath, name, NC * sizeof(double));
        EXPECT_FALSE(tmat.prefer_rows());
    }

    {
        // Both are compromised.
        dump(std::make_pair<int, int>(10, 10));
        tatami::HDF5DenseMatrix<double, int> mat(fpath, name, 0);
        EXPECT_TRUE(mat.prefer_rows());

        tatami::HDF5DenseMatrix<double, int, true> tmat(fpath, name, 0);
        EXPECT_FALSE(tmat.prefer_rows());
    }
}

/*************************************
 *************************************/

class HDF5DenseAccessTest : public ::testing::TestWithParam<std::tuple<bool, int, std::pair<int, int>, bool> >, public HDF5DenseMatrixTestMethods {};

TEST_P(HDF5DenseAccessTest, Basic) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

    auto chunk_sizes = std::get<2>(param);
    dump(chunk_sizes);

    // Make sure the cache size is smaller than the dataset, to get some more interesting things happening.
    auto cache_size = std::get<3>(param) ? NR * 10 : 0;
    tatami::HDF5DenseMatrix<double, int> mat(fpath, name, cache_size); 
    tatami::DenseRowMatrix<double, int> ref(NR, NC, values);

    test_simple_column_access(&mat, &ref, FORWARD, JUMP);
    test_simple_row_access(&mat, &ref, FORWARD, JUMP);
}

TEST_P(HDF5DenseAccessTest, Transposed) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

    auto chunk_sizes = std::get<2>(param);
    dump(chunk_sizes);

    auto cache_size = std::get<3>(param) ? NC * 5 : 0;
    tatami::HDF5DenseMatrix<double, int, true> mat(fpath, name, cache_size);
    std::shared_ptr<tatami::Matrix<double, int> > ptr(new tatami::DenseRowMatrix<double, int>(NR, NC, values));
    tatami::DelayedTranspose<double, int> ref(std::move(ptr));

    test_simple_column_access(&mat, &ref, FORWARD, JUMP);
    test_simple_row_access(&mat, &ref, FORWARD, JUMP);
}

TEST_P(HDF5DenseAccessTest, Apply) {
    // Putting it through its paces for correct parallelization via apply.
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

    auto chunk_sizes = std::get<2>(param);
    dump(chunk_sizes);

    auto cache_size = std::get<3>(param) ? NC * 4 : 0;
    tatami::HDF5DenseMatrix<double, int> mat(fpath, name, cache_size);
    tatami::DenseRowMatrix<double, int> ref(NR, NC, values);

    EXPECT_EQ(tatami::row_sums(&mat), tatami::row_sums(&ref));
    EXPECT_EQ(tatami::column_sums(&mat), tatami::column_sums(&ref));
}

INSTANTIATE_TEST_CASE_P(
    HDF5DenseMatrix,
    HDF5DenseAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false),
        ::testing::Values(1, 3),
        ::testing::Values(
            std::make_pair(NR, 1),
            std::make_pair(1, NC),
            std::make_pair(7, 13), // using chunk sizes that are a little odd to check for off-by-one errors.
            std::make_pair(13, 7),
            std::make_pair(11, 11),
            std::make_pair(0, 0)
        ),
        ::testing::Values(true, false) // Whether to cache or not.
    )
);

/*************************************
 *************************************/

class HDF5DenseSlicedTest : public ::testing::TestWithParam<std::tuple<bool, size_t, std::vector<double>, std::pair<int, int>, bool> >, public HDF5DenseMatrixTestMethods {};

TEST_P(HDF5DenseSlicedTest, Basic) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    auto chunk_sizes = std::get<3>(param);
    dump(chunk_sizes);

    auto cache_size = std::get<4>(param) ? NR * 5 : 0;
    tatami::HDF5DenseMatrix<double, int> mat(fpath, name, cache_size);
    tatami::DenseRowMatrix<double, int> ref(NR, NC, values);

    {
        size_t FIRST = interval_info[0] * NC, LAST = interval_info[1] * NC;
        test_sliced_row_access(&mat, &ref, FORWARD, JUMP, FIRST, LAST); 
    }

    {
        size_t FIRST = interval_info[0] * NR, LAST = interval_info[1] * NR;
        test_sliced_column_access(&mat, &ref, FORWARD, JUMP, FIRST, LAST);
    }
}

TEST_P(HDF5DenseSlicedTest, Transposed) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    auto chunk_sizes = std::get<3>(param);
    dump(chunk_sizes);

    auto cache_size = std::get<4>(param) ? NC * 10 : 0;
    tatami::HDF5DenseMatrix<double, int, true> mat(fpath, name, cache_size);
    std::shared_ptr<tatami::Matrix<double, int> > ptr(new tatami::DenseRowMatrix<double, int>(NR, NC, values));
    tatami::DelayedTranspose<double, int> ref(std::move(ptr));

    {
        size_t FIRST = interval_info[0] * NR, LAST = interval_info[1] * NR; // NR is deliberate here, it's transposed.
        test_sliced_row_access(&mat, &ref, FORWARD, JUMP, FIRST, LAST); 
    }

    {
        size_t FIRST = interval_info[0] * NC, LAST = interval_info[1] * NC; // NC is deliberate here, it's transposed.
        test_sliced_column_access(&mat, &ref, FORWARD, JUMP, FIRST, LAST);
    }
}

INSTANTIATE_TEST_CASE_P(
    HDF5DenseMatrix,
    HDF5DenseSlicedTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.5 }), 
            std::vector<double>({ 0.25, 0.75 }), 
            std::vector<double>({ 0.51, 1 })
        ),
        ::testing::Values(
            std::make_pair(7, 13), // using chunk sizes that are a little odd to check for off-by-one errors.
            std::make_pair(13, 7),
            std::make_pair(11, 11)
        ),
        ::testing::Values(true, false) // Whether to cache or not.
    )
);

/*************************************
 *************************************/

class HDF5DenseIndexedTest : public ::testing::TestWithParam<std::tuple<bool, size_t, std::vector<double>, std::pair<int, int>, bool> >, public HDF5DenseMatrixTestMethods {};

TEST_P(HDF5DenseIndexedTest, Basic) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    auto chunk_sizes = std::get<3>(param);
    dump(chunk_sizes);

    auto cache_size = std::get<4>(param) ? NR * 3 : 0;
    tatami::HDF5DenseMatrix<double, int> mat(fpath, name, cache_size);
    tatami::DenseRowMatrix<double, int> ref(NR, NC, values);

    {
        size_t FIRST = interval_info[0] * NC, STEP = interval_info[1];
        test_indexed_row_access(&mat, &ref, FORWARD, JUMP, FIRST, STEP); 
    }

    {
        size_t FIRST = interval_info[0] * NR, STEP = interval_info[1];
        test_indexed_column_access(&mat, &ref, FORWARD, JUMP, FIRST, STEP);
    }
}

TEST_P(HDF5DenseIndexedTest, Transposed) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);

    auto chunk_sizes = std::get<3>(param);
    dump(chunk_sizes);

    auto cache_size = std::get<4>(param) ? NC * 8 : 0;
    tatami::HDF5DenseMatrix<double, int, true> mat(fpath, name, cache_size);
    std::shared_ptr<tatami::Matrix<double, int> > ptr(new tatami::DenseRowMatrix<double, int>(NR, NC, values));
    tatami::DelayedTranspose<double, int> ref(std::move(ptr));

    {
        size_t FIRST = interval_info[0] * NR, STEP = interval_info[1]; // NR is deliberate here, it's transposed.
        test_indexed_row_access(&mat, &ref, FORWARD, JUMP, FIRST, STEP); 
    }

    {
        size_t FIRST = interval_info[0] * NC, STEP = interval_info[1]; // NC is deliberate here, it's transposed.
        test_indexed_column_access(&mat, &ref, FORWARD, JUMP, FIRST, STEP);
    }
}

INSTANTIATE_TEST_CASE_P(
    HDF5DenseMatrix,
    HDF5DenseIndexedTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0.3, 5 }), 
            std::vector<double>({ 0.11, 9 }), 
            std::vector<double>({ 0.4, 7 })
        ),
        ::testing::Values(
            std::make_pair(7, 13), // using chunk sizes that are a little odd to check for off-by-one errors.
            std::make_pair(13, 7),
            std::make_pair(11, 11)
        ),
        ::testing::Values(true, false) // Whether to cache or not.
    )
);
