#include <gtest/gtest.h>
#include "tatami/utils/wrap_shared_ptr.hpp"
#include "tatami/base/DenseMatrix.hpp"

#include "../_tests/test_row_access.h"
#include "../_tests/test_column_access.h"

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

        test_simple_row_access(wrap.get(), &mat);
        test_simple_column_access(wrap.get(), &mat);
    }

    // Still runs properly as wrapped pointer deletion is a no-op.
    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);
}
