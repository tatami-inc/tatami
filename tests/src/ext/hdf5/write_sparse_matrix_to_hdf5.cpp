#include <gtest/gtest.h>

#ifdef TEST_CUSTOM_PARALLEL // make sure this is included before tatami::apply.
#include "../../stats/custom_parallel.h"
#include "hdf5_custom_lock.h"
#endif

#include "H5Cpp.h"
#include "tatami/base/CompressedSparseMatrix.hpp"
#include "tatami/ext/hdf5/load_hdf5_matrix.hpp"
#include "tatami/ext/hdf5/write_sparse_matrix_to_hdf5.hpp"

#include "../temp_file_path.h"
#include <vector>
#include <random>


#include "../../_tests/test_column_access.h"
#include "../../_tests/test_row_access.h"
#include "../../_tests/simulate_vector.h"

TEST(WriteSparseMatrixToHdf5Test, Sparse) {
    const size_t NR = 200, NC = 100;

    // Mocking up a sparse matrix.
    std::string name = "stuff";
    SparseDetails<double> triplets = simulate_sparse_triplets<double>(NC, NR, 0.05, 0, 100);
    tatami::CompressedSparseMatrix<false, double, int> mat(NR, NC, std::move(triplets.value), std::move(triplets.index), std::move(triplets.ptr));

    // Dumping it.
    auto fpath = temp_file_path("tatami-write-test.h5");
    H5::H5File fhandle(fpath, H5F_ACC_TRUNC);
    auto mhandle = fhandle.createGroup("matrix");
    tatami::write_sparse_matrix_to_hdf5(&mat, mhandle);

    // Roundtripping.
    auto reloaded = tatami::load_hdf5_compressed_sparse_matrix<false, double, int>(NR, NC, fpath, "matrix/data", "matrix/indices", "matrix/indptr");
    for (size_t r = 0; r < NR; ++r) {
        auto matrow = mat.row(r);
        auto relrow = reloaded.row(r);
        EXPECT_EQ(matrow, relrow);
    }
}
