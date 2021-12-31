#include <gtest/gtest.h>
#include "tatami/utils/NakedArray.hpp"
#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/CompressedSparseMatrix.hpp"

TEST(NakedArray, Dense) {
    std::vector<double> contents(200);
    double counter = -105;
    for (auto& i : contents) { i = counter++; }

    tatami::DenseColumnMatrix<double> mat(10, 20, contents);
    tatami::DenseColumnMatrix<double, int, tatami::NakedArray<double> > mat2(10, 20, tatami::NakedArray(contents.data(), contents.size()));

    EXPECT_EQ(mat2.nrow(), 10);
    EXPECT_EQ(mat2.ncol(), 20);

    for (size_t i = 0; i < mat.ncol(); ++i) {
        auto col1 = mat.column(i);
        auto col2 = mat2.column(i);
        EXPECT_EQ(col1, col2);
    }
}

TEST(NakedArray, Sparse) {
    std::vector<double> values { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    std::vector<int> indices   { 2, 4, 1, 6, 7, 0, 3, 9, 5 };
    std::vector<size_t> indptrs{ 0, 2, 5, 8, 9 };

    tatami::CompressedSparseColumnMatrix<double, int> mat(10, 4, values, indices, indptrs);
    tatami::CompressedSparseColumnMatrix<double, int, 
        tatami::NakedArray<double>,
        tatami::NakedArray<int>,
        tatami::NakedArray<size_t>
    > mat2(10, 4, 
        tatami::NakedArray(values.data(), values.size()),
        tatami::NakedArray(indices.data(), indices.size()),
        tatami::NakedArray(indptrs.data(), indptrs.size()));

    EXPECT_EQ(mat2.nrow(), 10);
    EXPECT_EQ(mat2.ncol(), 4);

    for (size_t i = 0; i < mat.ncol(); ++i) {
        auto col1 = mat.column(i);
        auto col2 = mat2.column(i);
        EXPECT_EQ(col1, col2);
    }
}
