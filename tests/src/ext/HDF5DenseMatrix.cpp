#include <gtest/gtest.h>

#include "H5Cpp.h"
#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedTranspose.hpp"
#include "tatami/ext/HDF5DenseMatrix.hpp"

#include "temp_file_path.h"
#include <vector>
#include <random>

const size_t NR = 200, NC = 100;

class HDF5DenseMatrixTest : public ::testing::TestWithParam<std::pair<int, int> > {
protected:
    std::vector<double> assemble() {
        std::mt19937_64 rng(1234567890);
        std::uniform_real_distribution<> unif(0.0, 1.0);
        std::vector<double> values(NR * NC);
        for (auto& v : values) {
            v = std::round(unif(rng) * 100);
        }
        return values;
    }

    void dump(std::string fpath, std::string dpath, const std::vector<double>& values, const std::pair<int, int>& caching) {
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

        auto dhandle = fhandle.createDataSet(dpath, dtype, dspace, plist);
        dhandle.write(values.data(), H5::PredType::NATIVE_DOUBLE);
        return;
    }
};

/*************************************
 *************************************/

TEST_P(HDF5DenseMatrixTest, SimpleRowAccess) {
    auto caching = GetParam();

    auto truth = assemble();
    std::string fpath = temp_file_path("tatami-dense-test.h5");
    std::string name = "stuff";
    dump(fpath, name, truth, caching);

    tatami::HDF5DenseMatrix<double, int> mat(fpath, name);
    EXPECT_EQ(mat.nrow(), NR);
    EXPECT_EQ(mat.ncol(), NC);

    if (caching.first) {
        EXPECT_EQ(mat.prefer_rows(), NR/caching.first > NC/caching.second);
    } else {
        EXPECT_TRUE(mat.prefer_rows());
    }

    // No workspace.
    tatami::DenseRowMatrix<double, int> ref(NR, NC, truth);
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
    auto caching = GetParam();

    auto truth = assemble();
    std::string fpath = temp_file_path("tatami-dense-test.h5");
    std::string name = "stuff";
    dump(fpath, name, truth, caching);

    tatami::HDF5DenseMatrix<double, int> mat(fpath, name);
    tatami::DenseRowMatrix<double, int> ref(NR, NC, truth);

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

    auto truth = assemble();
    std::string fpath = temp_file_path("tatami-dense-test.h5");
    std::string name = "stuff";
    dump(fpath, name, truth, caching);

    tatami::HDF5DenseMatrix<double, int, true> mat(fpath, name);
    EXPECT_EQ(mat.nrow(), NC);
    EXPECT_EQ(mat.ncol(), NR);

    if (caching.first) {
        EXPECT_EQ(mat.prefer_rows(), NR/caching.first < NC/caching.second);
    } else {
        EXPECT_FALSE(mat.prefer_rows());
    }

    // No workspace.
    auto ref = tatami::DelayedTranspose<double, int>(std::shared_ptr<tatami::Matrix<double, int> >(new tatami::DenseRowMatrix<double, int>(NR, NC, truth)));
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
    auto caching = GetParam();

    auto truth = assemble();
    std::string fpath = temp_file_path("tatami-dense-test.h5");
    std::string name = "stuff";
    dump(fpath, name, truth, caching);

    tatami::HDF5DenseMatrix<double, int, true> mat(fpath, name);

    // No workspace.
    auto ref = tatami::DelayedTranspose<double, int>(std::shared_ptr<tatami::Matrix<double, int> >(new tatami::DenseRowMatrix<double, int>(NR, NC, truth)));
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

//TEST_P(HDF5DenseMatrixTest, SlicedRowAccess) {
//    auto caching = GetParam();
//
//    auto truth = assemble();
//    std::string fpath = "tatami-dense-test.h5";
//    std::string name = "stuff";
//    dump(fpath, name, truth, caching);
//
//    tatami::HDF5DenseMatrix<double, int> mat(fpath, name);
//    tatami::DenseRowMatrix<double, int> ref(NR, NC, truth);
//
//    size_t JUMP = 7;
//    size_t SHIFT = 12;
//    for (size_t i = 0; i < NC; i += JUMP, FIRST += SHIFT) {
//        size_t c = (FORWARD ? i : NC - i - 1);
//        EXPECT_EQ(mat.column(i, )), ref.column(i));
//    }
//}


