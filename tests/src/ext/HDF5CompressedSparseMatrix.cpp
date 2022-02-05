#include <gtest/gtest.h>

#include "H5Cpp.h"
#include "tatami/base/CompressedSparseMatrix.hpp"
#include "tatami/base/DelayedTranspose.hpp"
#include "tatami/ext/HDF5CompressedSparseMatrix.hpp"

#include "temp_file_path.h"
#include <vector>
#include <random>

#include "../_tests/test_column_access.h"
#include "../_tests/test_row_access.h"
#include "../_tests/simulate_vector.h"

const size_t NR = 200, NC = 100;

class HDF5SparseMatrixTestMethods {
protected:
    std::vector<double> values;
    std::string fpath;
    std::string name;

    void dump(const int& caching) {
        fpath = temp_file_path("tatami-sparse-test.h5");
        H5::H5File fhandle(fpath, H5F_ACC_TRUNC);
        name = "stuff";
        auto ghandle = fhandle.createGroup(name);

        auto triplets = simulate_sparse_triplets<double>(NR, NC, 0.05, 0, 100);
        for (auto& v : triplets.value) {
            v = std::round(v);
        }

        H5::DSetCreatPropList plist(H5::DSetCreatPropList::DEFAULT.getId());
        if (caching == 0) {
            plist.setLayout(H5D_CONTIGUOUS);
        } else {
            plist.setLayout(H5D_CHUNKED);
            hsize_t chunkdim = caching;
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
    dump(50);
    tatami::HDF5CompressedSparseMatrix<double, int, true> mat(NR, NC, fpath, name + "/data", name + "/index", name + "/indptr");
    EXPECT_EQ(mat.nrow(), NR);
    EXPECT_EQ(mat.ncol(), NC);
    EXPECT_TRUE(mat.sparse());
}
