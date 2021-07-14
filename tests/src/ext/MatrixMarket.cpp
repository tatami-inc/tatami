#include <gtest/gtest.h>

#include "tatami/ext/MatrixMarket.hpp"

#include <limits>
#include <random>
#include <cstdio>

template<class U, class V, class W>
auto write_MatrixMarket(size_t nr, size_t nc, const U& vals, const V& rows, const W& cols) {
    std::string path("testing-WHEE.txt"); 
    std::ofstream out(path);

    out << "%%MatrixMarket matrix coordinate integer general\n";
    out << nr << " " << nc << " " << vals.size();

    for (size_t i = 0; i < vals.size(); ++i) {
        out << "\n" << rows[i] << " " << cols[i] << " " << vals[i];
    }
    out << std::endl;
    out.close();
    return path;
}

TEST(MatrixMarket, LayeredLoaderSimple) {
    std::vector<int> rows{ 1, 5, 8, 2, 9, 0, 4 };
    std::vector<int> cols{ 2, 3, 1, 0, 2, 2, 4 };
    std::vector<int> vals{ 0, 1, 10, 100, 1000, 10000, 100000 }; 

    size_t NR = 10, NC = 5;
    auto newfile = write_MatrixMarket(NR, NC, vals, rows, cols);
    auto loaded = tatami::MatrixMarket::load_layered_sparse_matrix(newfile.c_str());
    const auto& out = loaded.matrix;

    auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new tatami::CompressedSparseColumnMatrix<double, int,
                                                                                               decltype(vals),
                                                                                               decltype(rows),
                                                                                               decltype(indptrs)
                                                                                               >(10, 5, std::move(vals), std::move(rows), std::move(indptrs))); 

    EXPECT_EQ(out->nrow(), ref->nrow());
    EXPECT_EQ(out->ncol(), ref->ncol());
    EXPECT_TRUE(out->sparse());
    EXPECT_FALSE(out->prefer_rows());

    for (size_t i = 0; i < NR; ++i) {
        int adjusted = loaded.permutation[i];
        EXPECT_EQ(out->row(adjusted), ref->row(i));
    }
}
