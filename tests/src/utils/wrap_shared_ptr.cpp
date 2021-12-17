#include <gtest/gtest.h>
#include "tatami/utils/wrap_shared_ptr.hpp"
#include "tatami/base/DenseMatrix.hpp"

TEST(WrapSharedPtrTest, Simple) {
    std::vector<double> contents(200);
    double counter = -105;
    for (auto& i : contents) { i = counter++; }
    tatami::DenseColumnMatrix<double> mat(10, 20, contents);
    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);

    {
        auto wrap = tatami::wrap_shared_ptr(&mat);
        EXPECT_EQ(wrap->nrow(), 10);
        EXPECT_EQ(wrap->ncol(), 20);

        EXPECT_EQ(wrap->column(0), mat.column(0));
        EXPECT_EQ(wrap->row(0), mat.row(0));
    }

    // Still runs properly as wrapped pointer deletion is a no-op.
    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);
}
