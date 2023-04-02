#include <gtest/gtest.h>

#ifdef TEST_CUSTOM_PARALLEL // make sure this is included before tatami::apply.
#include "../../stats/custom_parallel.h"
#include "hdf5_custom_lock.h"
#endif

#include "H5Cpp.h"
#include "tatami/ext/hdf5/load_hdf5_matrix.hpp"

#include "../temp_file_path.h"
#include <vector>
#include <random>

#include "../../_tests/test_column_access.h"
#include "../../_tests/test_row_access.h"
#include "../../_tests/simulate_vector.h"

TEST(LoadHDF5MatrixTest, Sparse) {
    const size_t NR = 200, NC = 100;

    // Dumping a sparse matrix.
    auto fpath = temp_file_path("tatami-sparse-test.h5");
    std::string name = "stuff";
    auto triplets = simulate_sparse_compressed<double>(NR, NC, 0.05, 0, 100);

    {
        H5::H5File fhandle(fpath, H5F_ACC_TRUNC);
        auto ghandle = fhandle.createGroup(name);

        hsize_t dims = triplets.value.size();
        H5::DataSpace dspace(1, &dims);
        {
            H5::DataType dtype(H5::PredType::NATIVE_DOUBLE);
            auto dhandle = ghandle.createDataSet("data", dtype, dspace);
            dhandle.write(triplets.value.data(), H5::PredType::NATIVE_DOUBLE);
        }

        {
            H5::DataType dtype(H5::PredType::NATIVE_UINT16);
            auto dhandle = ghandle.createDataSet("index", dtype, dspace);
            dhandle.write(triplets.index.data(), H5::PredType::NATIVE_INT);
        }

        {
            hsize_t ncp1 = triplets.ptr.size();
            H5::DataSpace dspace(1, &ncp1);
            H5::DataType dtype(H5::PredType::NATIVE_UINT64);
            auto dhandle = ghandle.createDataSet("indptr", dtype, dspace);
            dhandle.write(triplets.ptr.data(), H5::PredType::NATIVE_LONG);
        }
    }

    // Basic load as a CSR matrix (as rows are the primary dimension in this simulation)
    {
        auto mat = tatami::load_hdf5_compressed_sparse_matrix<true, double, int>(NR, NC, fpath, name + "/data", name + "/index", name + "/indptr");
        tatami::CompressedSparseMatrix<
            true, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NR, NC, triplets.value, triplets.index, triplets.ptr);

        test_simple_row_access(&mat, &ref, true, 1);
    }

    // Pretending it's a CSC matrix.
    {
        auto mat = tatami::load_hdf5_compressed_sparse_matrix<false, double, int>(NC, NR, fpath, name + "/data", name + "/index", name + "/indptr");
        tatami::CompressedSparseMatrix<
            false, 
            double, 
            int, 
            decltype(triplets.value), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NC, NR, triplets.value, triplets.index, triplets.ptr);

        test_simple_column_access(&mat, &ref, true, 1);
    }

    // Trying a variety of storage types.
    {
        auto mat = tatami::load_hdf5_compressed_sparse_matrix<
            true, 
            double, 
            int,
            std::vector<uint16_t>,
            std::vector<uint32_t>,
            std::vector<uint64_t>
        >(NR, NC, fpath, name + "/data", name + "/index", name + "/indptr");

        std::vector<double> truncated = triplets.value;
        for (auto& x : truncated) {
            x = std::trunc(x);
        }

        tatami::CompressedSparseMatrix<
            true, 
            double, 
            int, 
            decltype(truncated), 
            decltype(triplets.index), 
            decltype(triplets.ptr)
        > ref(NR, NC, std::move(truncated), triplets.index, triplets.ptr);

        test_simple_column_access(&mat, &ref, true, 1);
    }
}

TEST(LoadHDF5MatrixTest, Dense) {
    size_t NR = 200, NC = 100;
    auto fpath = temp_file_path("tatami-dense-test.h5");
    std::string name = "stuff";

    std::vector<double> values = simulate_dense_vector<double>(NR * NC, 0, 100);

    {
        H5::H5File fhandle(fpath, H5F_ACC_TRUNC);
        hsize_t dims[2];
        dims[0] = NR;
        dims[1] = NC;
        H5::DataSpace dspace(2, dims);
        H5::DataType dtype(H5::PredType::NATIVE_DOUBLE);
        auto dhandle = fhandle.createDataSet(name, dtype, dspace);
        dhandle.write(values.data(), H5::PredType::NATIVE_DOUBLE);
    }

    // Basic load as a row-major matrix. 
    {
        auto mat = tatami::load_hdf5_dense_matrix<double, int>(fpath, name);
        tatami::DenseRowMatrix<double, int> ref(NR, NC, values);
        test_simple_row_access(&mat, &ref, true, 1);
    }

    // Pretending it's a column-major matrix.
    {
        auto mat = tatami::load_hdf5_dense_matrix<double, int, std::vector<double>, true>(fpath, name);
        tatami::DenseColumnMatrix<double, int> ref(NC, NR, values);
        test_simple_column_access(&mat, &ref, true, 1);
    }

    // Trying a different storage type.
    {
        auto mat = tatami::load_hdf5_dense_matrix<double, int, std::vector<int32_t> >(fpath, name);

        std::vector<double> truncated = values;
        for (auto& x : truncated) {
            x = std::trunc(x);
        }
        tatami::DenseRowMatrix<double, int> ref(NR, NC, std::move(truncated));

        test_simple_column_access(&mat, &ref, true, 1);
    }
}
