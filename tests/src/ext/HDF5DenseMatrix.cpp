#include <gtest/gtest.h>

#include "H5Cpp.h"
#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedTranspose.hpp"
#include "tatami/ext/HDF5DenseMatrix.hpp"

#include "temp_file_path.h"
#include <vector>
#include <random>

const size_t NR = 200, NC = 100;

class HDF5DenseMatrixTestMethods {
protected:
    std::vector<double> values;
    std::string fpath;
    std::string name;

    void dump(const std::pair<int, int>& caching) {
        fpath = temp_file_path("tatami-dense-test.h5");
        name = "stuff";
        H5::H5File fhandle(fpath, H5F_ACC_TRUNC);

        hsize_t dims[2];
        dims[0] = NR;
        dims[1] = NC;
        H5::DataSpace dspace(2, dims);
        H5::DataType dtype(H5::PredType::NATIVE_UINT8);

        H5::DSetCreatPropList plist(H5::DSetCreatPropList::DEFAULT.getId());
        if (caching.first == 0) {
            plist.setLayout(H5D_CONTIGUOUS);
        } else {
            plist.setLayout(H5D_CHUNKED);
            hsize_t chunkdims[2];
            chunkdims[0] = caching.first;
            chunkdims[1] = caching.second;
            plist.setChunk(2, chunkdims);
        }

        auto dhandle = fhandle.createDataSet(name, dtype, dspace, plist);

        // Generating the dataset.
        std::mt19937_64 rng(1234567890);
        std::uniform_real_distribution<> unif(0.0, 1.0);
        values.resize(NR * NC);
        for (auto& v : values) {
            v = std::round(unif(rng) * 100);
        }

        dhandle.write(values.data(), H5::PredType::NATIVE_DOUBLE);
        return;
    }

    static std::pair<size_t, size_t> wrap_intervals(size_t first, size_t last, size_t max) {
        size_t diff = last - first;
        first %= max;
        last = std::min(max, first + diff);
        return std::make_pair(first, last);
    }
};

/*************************************
 *************************************/

class HDF5DenseMatrixTest : public ::testing::TestWithParam<std::pair<int, int> >, public HDF5DenseMatrixTestMethods {};

TEST_P(HDF5DenseMatrixTest, SimpleRowAccess) {
    auto caching = GetParam();
    dump(caching);

    tatami::HDF5DenseMatrix<double, int> mat(fpath, name);
    EXPECT_EQ(mat.nrow(), NR);
    EXPECT_EQ(mat.ncol(), NC);

    if (caching.first) {
        EXPECT_EQ(mat.prefer_rows(), NR/caching.first > NC/caching.second);
    } else {
        EXPECT_TRUE(mat.prefer_rows());
    }

    // No workspace.
    tatami::DenseRowMatrix<double, int> ref(NR, NC, values);
    for (size_t i = 0; i < NR; ++i) {
        EXPECT_EQ(mat.row(i), ref.row(i));
    }

    // With workspace.
    auto wrk = mat.new_workspace(true);
    for (size_t i = 0; i < NR; ++i) {
        EXPECT_EQ(mat.row(i, wrk.get()), ref.row(i));
    }
}

TEST_P(HDF5DenseMatrixTest, SimpleColumnAccess) {
    dump(GetParam());
    tatami::HDF5DenseMatrix<double, int> mat(fpath, name);
    tatami::DenseRowMatrix<double, int> ref(NR, NC, values);

    // No workspace.
    for (size_t i = 0; i < NC; ++i) {
        EXPECT_EQ(mat.column(i), ref.column(i));
    }

    // With workspace.
    auto wrk = mat.new_workspace(false);
    for (size_t i = 0; i < NC; ++i) {
        EXPECT_EQ(mat.column(i, wrk.get()), ref.column(i));
    }
}

TEST_P(HDF5DenseMatrixTest, TransposedRowAccess) {
    auto caching = GetParam();
    dump(caching);

    tatami::HDF5DenseMatrix<double, int, true> mat(fpath, name);
    EXPECT_EQ(mat.nrow(), NC);
    EXPECT_EQ(mat.ncol(), NR);

    if (caching.first) {
        EXPECT_EQ(mat.prefer_rows(), NR/caching.first < NC/caching.second);
    } else {
        EXPECT_FALSE(mat.prefer_rows());
    }

    // No workspace.
    std::shared_ptr<tatami::Matrix<double, int> > ptr(new tatami::DenseRowMatrix<double, int>(NR, NC, values));
    auto ref = tatami::DelayedTranspose<double, int>(std::move(ptr));
    for (size_t i = 0; i < NC; ++i) {
        EXPECT_EQ(mat.row(i), ref.row(i));
    }

    // With workspace.
    auto wrk = mat.new_workspace(true);
    for (size_t i = 0; i < NC; ++i) {
        EXPECT_EQ(mat.row(i, wrk.get()), ref.row(i));
    }
}

TEST_P(HDF5DenseMatrixTest, TransposedColumnAccess) {
    dump(GetParam());
    tatami::HDF5DenseMatrix<double, int, true> mat(fpath, name);

    // No workspace.
    std::shared_ptr<tatami::Matrix<double, int> > ptr(new tatami::DenseRowMatrix<double, int>(NR, NC, values));
    auto ref = tatami::DelayedTranspose<double, int>(std::move(ptr));
    for (size_t i = 0; i < NR; ++i) {
        EXPECT_EQ(mat.column(i), ref.column(i));
    }

    // With workspace.
    auto wrk = mat.new_workspace(false);
    for (size_t i = 0; i < NR; ++i) {
        EXPECT_EQ(mat.column(i, wrk.get()), ref.column(i));
    }
}

INSTANTIATE_TEST_CASE_P(
    HDF5DenseMatrix,
    HDF5DenseMatrixTest,
    ::testing::Values(
        std::make_pair(NR, 1),
        std::make_pair(1, NC),
        std::make_pair(10, 10),
        std::make_pair(0, 0)
    )
);

/*************************************
 *************************************/

class HDF5DenseMatrixSlicedTest : public ::testing::TestWithParam<std::tuple<bool, size_t, std::vector<size_t> > >, public HDF5DenseMatrixTestMethods {};

TEST_P(HDF5DenseMatrixSlicedTest, ColumnAccess) {
    dump(std::make_pair(10, 10)); // basic square chunks.
    tatami::HDF5DenseMatrix<double, int> mat(fpath, name);
    tatami::DenseRowMatrix<double, int> ref(NR, NC, values);

    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    auto work = mat.new_workspace(false);
    for (size_t i = 0; i < NC; i += JUMP, FIRST += SHIFT) {
        size_t c = (FORWARD ? i : NC - i - 1);
        auto interval = wrap_intervals(FIRST, FIRST + LEN, NR);
        size_t start = interval.first, end = interval.second;

        auto expected = ref.column(i, start, end);
        EXPECT_EQ(expected, mat.column(i, start, end));
        EXPECT_EQ(expected, mat.column(i, start, end, work.get()));
    }
}

TEST_P(HDF5DenseMatrixSlicedTest, RowAccess) {
    dump(std::make_pair(10, 10)); // basic square chunks.
    tatami::HDF5DenseMatrix<double, int> mat(fpath, name);
    tatami::DenseRowMatrix<double, int> ref(NR, NC, values);

    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    auto work = mat.new_workspace(true);
    for (size_t i = 0; i < NR; i += JUMP, FIRST += SHIFT) {
        size_t c = (FORWARD ? i : NR - i - 1);
        auto interval = wrap_intervals(FIRST, FIRST + LEN, NC);
        size_t start = interval.first, end = interval.second;

        auto expected = ref.row(i, start, end);
        EXPECT_EQ(expected, mat.row(i, start, end));
        EXPECT_EQ(expected, mat.row(i, start, end, work.get()));
    }
}

INSTANTIATE_TEST_CASE_P(
    HDF5DenseMatrix,
    HDF5DenseMatrixSlicedTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<size_t>({ 0, 8, 3 }), // overlapping shifts
            std::vector<size_t>({ 1, 4, 4 }), // non-overlapping shifts
            std::vector<size_t>({ 3, 10, 0 })
        )
    )
);
