#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedBind.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../data/data.h"
#include "TestCore.h"

class BindTest: public TestCore {
protected:
    std::shared_ptr<tatami::numeric_matrix> dense;
    std::shared_ptr<tatami::typed_matrix<double, int> > sparse;
protected:
    void assemble(size_t nr, size_t nc, const std::vector<double>& source) {
        dense = std::shared_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(nr, nc, source));
        sparse = tatami::convert_to_sparse(dense.get(), false); // column-major.
    }

    void SetUp() {
        assemble(sparse_nrow, sparse_ncol, sparse_matrix);
    }
};

TEST_F(BindTest, RowBindFullDenseAccess) {
    auto rbind = tatami::make_DelayedBind<0>(std::vector{dense, sparse});
    EXPECT_EQ(rbind->ncol(), dense->ncol());
    const size_t NR = dense->nrow();
    EXPECT_EQ(rbind->nrow(), 2 * NR);

    auto preference = rbind->dimension_preference();
    EXPECT_EQ(preference.first, preference.second);

    // Column access.
    auto wrk = rbind->new_workspace(false);
    set_sizes(0, 2 * NR);

    for (size_t i = 0; i < rbind->ncol(); ++i) {
        wipe_expected();
        fill_expected(dense->column(i, expected.data()));

        wipe_output();
        fill_output(rbind->column(i, output.data()));
        std::vector<double> ref(expected.begin(), expected.begin() + NR);
        EXPECT_EQ(std::vector<double>(output.begin(), output.begin() + NR), ref);
        EXPECT_EQ(std::vector<double>(output.begin() + NR, output.begin() + 2 * NR), ref);

        auto output2 = output;
        wipe_output();
        fill_output(rbind->column(i, output.data(), wrk.get()));
        EXPECT_EQ(output, output2);
    }

    // Row access.
    wrk = rbind->new_workspace(true);
    set_sizes(0, rbind->ncol());
    for (size_t i = 0; i < NR; ++i) {
        wipe_expected();
        fill_expected(dense->row(i % NR, expected.data()));

        wipe_output();
        fill_output(rbind->row(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(rbind->row(i, expected.data(), wrk.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(BindTest, RowBindSubsetDenseAccess) {
    auto rbind = tatami::make_DelayedBind<0>(std::vector{dense, sparse});
    size_t LEN = 6;
    size_t first = 2;
    const size_t NR = dense->nrow();

    // Column access.
    auto wrk = rbind->new_workspace(false);
    for (size_t i = 0; i < rbind->ncol(); ++i) {
        set_sizes(first, std::min(first + LEN, rbind->nrow()));

        wipe_expected();
        if (first < NR) {
            if (last <= NR) {
                fill_expected(dense->column(i, expected.data(), first, last));
            } else {
                auto ptr1 = dense->column(i, expected.data(), first, NR);
                if (expected.data() != ptr1) {
                    std::copy(ptr1, ptr1 + NR - first, expected.data());
                }

                auto exp = expected.data() + NR - first;
                auto ptr2 = dense->column(i, exp, 0, last - NR);
                if (exp != ptr2) {
                    std::copy(ptr2, ptr2 + last - NR, exp);
                }
            }
        }  else {
            fill_expected(dense->column(i, expected.data(), first - NR, last - NR));
        }

        wipe_output();
        fill_output(rbind->column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(rbind->column(i, output.data(), first, last, wrk.get()));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= rbind->nrow();
    }

    // Row access.
    wrk = rbind->new_workspace(true);
    LEN = 7;
    first = 0;
    for (size_t i = 0; i < rbind->nrow(); ++i) {
        set_sizes(first, std::min(first + LEN, rbind->ncol()));

        wipe_expected();
        fill_expected(dense->row(i % NR, expected.data(), first, last));

        wipe_output();
        fill_output(rbind->row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(rbind->row(i, output.data(), first, last, wrk.get()));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= rbind->ncol();
    }
}

TEST_F(BindTest, RowBindFullSparseAccess) {
    auto rbind = tatami::make_DelayedBind<0>(std::vector{dense, sparse});
    const size_t NR = dense->nrow();

    auto wrk = rbind->new_workspace(false);
    set_sizes(0, 2*NR);
    for (size_t i = 0; i < rbind->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->sparse_column(i, outval.data(), outidx.data()));

        wipe_output();
        fill_output(rbind->sparse_column(i, outval.data(), outidx.data()));

        std::vector<double> ref(expected.begin(), expected.begin() + NR);
        EXPECT_EQ(std::vector<double>(output.begin(), output.begin() + NR), ref);
        EXPECT_EQ(std::vector<double>(output.begin() + NR, output.begin() + 2 * NR), ref);

        auto output2 = output;
        wipe_output();
        fill_output(rbind->sparse_column(i, outval.data(), outidx.data(), wrk.get()));
        EXPECT_EQ(output, output2);
    }

    wrk = rbind->new_workspace(true);
    set_sizes(0, rbind->ncol());
    for (size_t i = 0; i < rbind->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->sparse_row(i % NR, outval.data(), outidx.data()));

        wipe_output();
        fill_output(rbind->sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(rbind->sparse_row(i, outval.data(), outidx.data(), wrk.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(BindTest, RowBindSubsetSparseAccess) {
    auto rbind = tatami::make_DelayedBind<0>(std::vector{dense, sparse});
    size_t LEN = 7;
    size_t first = 3;
    const size_t NR = dense->nrow();

    // Column access.
    auto wrk = rbind->new_workspace(false);
    for (size_t i = 0; i < rbind->ncol(); ++i) {
        set_sizes(first, std::min(first + LEN, rbind->nrow()));

        wipe_expected();
        if (first < NR) {
            if (last <= NR) {
                fill_expected(dense->sparse_column(i, outval.data(), outidx.data(), first, last));
            } else {
                auto range1 = dense->sparse_column(i, outval.data(), outidx.data(), first, NR);
                for (size_t i = 0; i < range1.number; ++i) {
                    expected[range1.index[i] - first] = range1.value[i];
                }

                auto range2 = dense->sparse_column(i, outval.data(), outidx.data(), 0, last - NR);
                for (size_t i = 0; i < range2.number; ++i) {
                    expected[NR - first + range2.index[i]] = range2.value[i];
                }
            }
        }  else {
            auto range = dense->sparse_column(i, outval.data(), outidx.data(), first - NR, last - NR);
            for (size_t i = 0; i < range.number; ++i) {
                expected[range.index[i] - (first - NR)] = range.value[i];
            }
        }

        wipe_output();
        fill_output(rbind->sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(rbind->sparse_column(i, outval.data(), outidx.data(), first, last, wrk.get()));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= rbind->nrow();
    }

    // Row access.
    wrk = rbind->new_workspace(true);
    LEN = 11;
    first = 0;
    for (size_t i = 0; i < rbind->nrow(); ++i) {
        set_sizes(first, std::min(first + LEN, rbind->ncol()));

        wipe_expected();
        fill_expected(dense->sparse_row(i % NR, outval.data(), outidx.data(), first, last));

        wipe_output();
        fill_output(rbind->sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(rbind->sparse_row(i, outval.data(), outidx.data(), first, last, wrk.get()));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= rbind->ncol();
    }
}

TEST_F(BindTest, ColumnBindFullDenseAccess) {
    auto cbind = tatami::make_DelayedBind<1>(std::vector{dense, sparse});
    EXPECT_EQ(cbind->nrow(), dense->nrow());
    const size_t NC = dense->ncol();
    EXPECT_EQ(cbind->ncol(), 2 * NC);

    auto preference = cbind->dimension_preference();
    EXPECT_EQ(preference.first, preference.second);

    // Row access.
    auto wrk = cbind->new_workspace(true);
    set_sizes(0, 2 * NC);

    for (size_t i = 0; i < cbind->nrow(); ++i) {
        wipe_expected();
        fill_expected(dense->row(i, expected.data()));

        wipe_output();
        fill_output(cbind->row(i, output.data()));
        std::vector<double> ref(expected.begin(), expected.begin() + NC);
        EXPECT_EQ(std::vector<double>(output.begin(), output.begin() + NC), ref);
        EXPECT_EQ(std::vector<double>(output.begin() + NC, output.begin() + 2 * NC), ref);

        auto output2 = output;
        wipe_output();
        fill_output(cbind->row(i, output.data(), wrk.get()));
        EXPECT_EQ(output, output2);
    }

    // Column access.
    wrk = cbind->new_workspace(false);
    set_sizes(0, cbind->nrow());
    for (size_t i = 0; i < NC; ++i) {
        wipe_expected();
        fill_expected(dense->column(i % NC, expected.data()));

        wipe_output();
        fill_output(cbind->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(cbind->column(i, expected.data(), wrk.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(BindTest, ColumnBindSubsetDenseAccess) {
    auto cbind = tatami::make_DelayedBind<1>(std::vector{dense, sparse});
    size_t LEN = 6;
    size_t first = 2;
    const size_t NC = dense->ncol();

    // Row access.
    auto wrk = cbind->new_workspace(true);
    for (size_t i = 0; i < cbind->nrow(); ++i) {
        set_sizes(first, std::min(first + LEN, cbind->ncol()));

        wipe_expected();
        if (first < NC) {
            if (last <= NC) {
                fill_expected(dense->row(i, expected.data(), first, last));
            } else {
                auto ptr1 = dense->row(i, expected.data(), first, NC);
                if (expected.data() != ptr1) {
                    std::copy(ptr1, ptr1 + NC - first, expected.data());
                }

                auto exp = expected.data() + NC - first;
                auto ptr2 = dense->row(i, exp, 0, last - NC);
                if (exp != ptr2) {
                    std::copy(ptr2, ptr2 + last - NC, exp);
                }
            }
        }  else {
            fill_expected(dense->row(i, expected.data(), first - NC, last - NC));
        }

        wipe_output();
        fill_output(cbind->row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(cbind->row(i, output.data(), first, last, wrk.get()));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= cbind->ncol();
    }

    // Column access.
    wrk = cbind->new_workspace(false);
    LEN = 7;
    first = 0;
    for (size_t i = 0; i < cbind->nrow(); ++i) {
        set_sizes(first, std::min(first + LEN, cbind->nrow()));

        wipe_expected();
        fill_expected(dense->column(i % NC, expected.data(), first, last));

        wipe_output();
        fill_output(cbind->column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(cbind->column(i, output.data(), first, last, wrk.get()));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= cbind->nrow();
    }
}

TEST_F(BindTest, ColumnBindFullSparseAccess) {
    auto cbind = tatami::make_DelayedBind<1>(std::vector{dense, sparse});
    const size_t NC = dense->ncol();

    // Row access
    auto wrk = cbind->new_workspace(true);
    set_sizes(0, 2*NC);
    for (size_t i = 0; i < cbind->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->sparse_row(i, outval.data(), outidx.data()));

        wipe_output();
        fill_output(cbind->sparse_row(i, outval.data(), outidx.data()));

        std::vector<double> ref(expected.begin(), expected.begin() + NC);
        EXPECT_EQ(std::vector<double>(output.begin(), output.begin() + NC), ref);
        EXPECT_EQ(std::vector<double>(output.begin() + NC, output.begin() + 2 * NC), ref);

        auto output2 = output;
        wipe_output();
        fill_output(cbind->sparse_row(i, outval.data(), outidx.data(), wrk.get()));
        EXPECT_EQ(output, output2);
    }

    // Column access
    wrk = cbind->new_workspace(false);
    set_sizes(0, cbind->ncol());
    for (size_t i = 0; i < cbind->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->sparse_column(i % NC, outval.data(), outidx.data()));

        wipe_output();
        fill_output(cbind->sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(cbind->sparse_column(i, outval.data(), outidx.data(), wrk.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(BindTest, ColumnBindSubsetSparseAccess) {
    auto cbind = tatami::make_DelayedBind<1>(std::vector{dense, sparse});
    size_t LEN = 7;
    size_t first = 3;
    const size_t NC = dense->ncol();

    // Row access.
    auto wrk = cbind->new_workspace(true);
    for (size_t i = 0; i < cbind->nrow(); ++i) {
        set_sizes(first, std::min(first + LEN, cbind->ncol()));

        wipe_expected();
        if (first < NC) {
            if (last <= NC) {
                fill_expected(dense->sparse_row(i, outval.data(), outidx.data(), first, last));
            } else {
                auto range1 = dense->sparse_row(i, outval.data(), outidx.data(), first, NC);
                for (size_t i = 0; i < range1.number; ++i) {
                    expected[range1.index[i] - first] = range1.value[i];
                }

                auto range2 = dense->sparse_row(i, outval.data(), outidx.data(), 0, last - NC);
                for (size_t i = 0; i < range2.number; ++i) {
                    expected[NC - first + range2.index[i]] = range2.value[i];
                }
            }
        }  else {
            auto range = dense->sparse_row(i, outval.data(), outidx.data(), first - NC, last - NC);
            for (size_t i = 0; i < range.number; ++i) {
                expected[range.index[i] - (first - NC)] = range.value[i];
            }
        }

        wipe_output();
        fill_output(cbind->sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(cbind->sparse_row(i, outval.data(), outidx.data(), first, last, wrk.get()));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= cbind->ncol();
    }

    // Column access.
    wrk = cbind->new_workspace(false);
    LEN = 11;
    first = 0;
    for (size_t i = 0; i < cbind->ncol(); ++i) {
        set_sizes(first, std::min(first + LEN, cbind->nrow()));

        wipe_expected();
        fill_expected(dense->sparse_column(i % NC, outval.data(), outidx.data(), first, last));

        wipe_output();
        fill_output(cbind->sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(cbind->sparse_column(i, outval.data(), outidx.data(), first, last, wrk.get()));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= cbind->nrow();
    }
}
