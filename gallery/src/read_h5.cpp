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
#include "tatami/ext/HDF5DenseMatrix.hpp"
#include "tatami/ext/HDF5CompressedSparseMatrix.hpp"

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "need two arguments: [file] [path]" << std::endl;
        return 1;
    }

    std::string file(argv[1]);
    H5::H5File handle(file, H5F_ACC_RDONLY);
    std::string path(argv[2]);

    std::shared_ptr<tatami::Matrix<double, int> > ptr;
    if (handle.childObjType(path) == H5O_TYPE_DATASET) {
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
    auto sums = tatami::column_sums(ptr.get());
    std::cout << "First sum: " << sums.front() << "\n";
    std::cout << "Last sum: " << sums.back() << "\n";

    auto end = std::chrono::system_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time taken: " << diff.count() << " ms\n";
    return 0;
}
