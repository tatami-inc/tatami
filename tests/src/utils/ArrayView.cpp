#include <gtest/gtest.h>
#include "tatami/utils/ArrayView.hpp"
#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/CompressedSparseMatrix.hpp"

#include "tatami_test/tatami_test.hpp"

#include <cstdint>
#include <numeric>

TEST(ArrayView, DenseMatrix) {
    int nr = 20, nc = 10;
    std::vector<int> values(nr * nc);
    std::iota(values.begin(), values.end(), -50);
    tatami::DenseColumnMatrix<double, int, decltype(values)> ref(nr, nc, values);

    tatami::ArrayView<int> arr(values.data(), values.size());
    tatami::DenseColumnMatrix<double, int, decltype(arr)> alt(nr, nc, arr);

    tatami_test::test_simple_row_access(&alt, &ref, true, 1);
    tatami_test::test_simple_column_access(&alt, &ref, true, 1);
}

TEST(ArrayView, SparseMatrix) {
    int nr = 10, nc = 6;
    std::vector<double> values  { 2, 5, 3, 4, 5, 5, 7, 3, 2, 4, 5, 1 };
    std::vector<int> indices { 0, 1, 7, 1, 4, 6, 2, 5, 5, 1, 8, 9 };
    std::vector<size_t> indptrs { 0, 3, 6, 8, 9, 11, 12 };
    tatami::CompressedSparseColumnMatrix<double, int, decltype(values), decltype(indices), decltype(indptrs)> ref(nr, nc, values, indices, indptrs);

    tatami::ArrayView<double> varr  (values.data(), values.size());
    tatami::ArrayView<int>    iarr  (indices.data(), indices.size());
    tatami::ArrayView<size_t> indarr(indptrs.data(), indptrs.size());
    tatami::CompressedSparseColumnMatrix<double, int, decltype(varr), decltype(iarr), decltype(indarr)> alt(nr, nc, varr, iarr, indarr);

    tatami_test::test_simple_row_access(&alt, &ref, true, 1);
    tatami_test::test_simple_column_access(&alt, &ref, true, 1);
}
