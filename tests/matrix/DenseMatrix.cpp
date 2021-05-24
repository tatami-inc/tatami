#include <gtest/gtest.h>
#include "matrix/DenseMatrix.hpp"
#include <vector>
#include <deque>
#include <numeric>

TEST(DenseMatrix, Construction) {
    bioc::DenseMatrix<double> mat(10, 20);
    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);
    EXPECT_EQ(mat.type(), bioc::_double);

    // Checks run properly.
    std::vector<double> contents;
    EXPECT_ANY_THROW({
        bioc::DenseMatrix<double> mat(10, 20, contents);
    });
    EXPECT_ANY_THROW({
        bioc::DenseMatrix<double> mat(10, 20, std::move(contents));
    });

    std::deque<double> more_contents(200);
    std::iota(more_contents.begin(), more_contents.end(), 1);
    bioc::DenseMatrix<double, std::deque<double> > mat2(10, 20, more_contents);
    EXPECT_EQ(more_contents.size(), 200);

    bioc::DenseMatrix<double, std::deque<double> > mat3(20, 10, std::move(more_contents)); // can't run many tests for state of more_contents here.
}
