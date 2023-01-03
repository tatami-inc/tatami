#include <exception>
#include <iostream>
#include <chrono>

/**
 * INTRODUCTION:
 *
 * This is little utility that demonstrates how to call the  
 * HDF5-backed matrices, given the path to a file and a name
 * of a resource inside the file. This can either be a dense
 * dataset containing a matrix, or a group containing a sparse
 * matrix in the same vein as the 10X matrices.
 */

#include "H5Cpp.h"
#include "tatami/stats/sums.hpp"
#include "tatami/ext/hdf5.hpp"

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "need two arguments: [file] [path]" << std::endl;
        return 1;
    }

    std::string file(argv[1]);
    H5::H5File handle(file, H5F_ACC_RDONLY);
    std::string path(argv[2]);

    std::shared_ptr<tatami::Matrix<double, int> > ptr;
    bool is_dense = (handle.childObjType(path) == H5O_TYPE_DATASET);
    if (is_dense) {
        handle.close();
        ptr.reset(new tatami::HDF5DenseMatrix<double, int, true>(file, path));
    } else {
        std::vector<long> dims(2);
        {
            auto shandle = handle.openDataSet(path + "/shape");
            shandle.read(dims.data(), H5::PredType::NATIVE_LONG);
        }
        handle.close();
        ptr.reset(new tatami::HDF5CompressedSparseMatrix<false, double, int>(dims[0], dims[1], file, path + "/data", path + "/indices", path + "/indptr"));
    }

    auto start = std::chrono::system_clock::now();
    auto sums = tatami::row_sums(ptr.get());
    auto total = std::accumulate(sums.begin(), sums.end(), 0);
    auto end = std::chrono::system_clock::now();

    auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Total sum: " << total << "\n";
    std::cout << "Time taken: " << diff.count() << " ms\n";

    // Comparing the time taken to read the entire file in.
    if (is_dense) {
        H5::H5File handle(file, H5F_ACC_RDONLY);
        auto dhandle = handle.openDataSet(path);
        std::vector<double> values(ptr->nrow() * ptr->ncol());

        auto start = std::chrono::system_clock::now();
        dhandle.read(values.data(), H5::PredType::NATIVE_DOUBLE);
        auto val = std::accumulate(values.begin(), values.end(), 0);
        auto end = std::chrono::system_clock::now();

        auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Total sum: " << val << std::endl;
        std::cout << "Time taken: " << diff.count() << " ms\n";
    } else {
        H5::H5File handle(file, H5F_ACC_RDONLY);
        auto dhandle = handle.openDataSet(path + "/data");
        auto ihandle = handle.openDataSet(path + "/indices");

        hsize_t n;
        dhandle.getSpace().getSimpleExtentDims(&n, NULL);
        std::vector<double> values(n);
        std::vector<int > index(n);

        auto start = std::chrono::system_clock::now();
        dhandle.read(values.data(), H5::PredType::NATIVE_DOUBLE);
        ihandle.read(index.data(), H5::PredType::NATIVE_INT);
        auto val = std::accumulate(values.begin(), values.end(), 0);
        auto end = std::chrono::system_clock::now();

        auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Total sum: " << val << std::endl;
        std::cout << "Time taken: " << diff.count() << " ms\n";
    }

    return 0;
}
