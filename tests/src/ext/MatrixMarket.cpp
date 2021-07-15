#include <gtest/gtest.h>

#include "tatami/ext/MatrixMarket.hpp"

#include <limits>
#include <random>
#include <fstream>
#include <string>
#include <vector>

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

template<class PARAM>
class MatrixMarketTest : public ::testing::TestWithParam<PARAM> {
protected:
    size_t NR, NC;
    std::vector<int> rows, cols, vals;

protected:
    auto extra_assemble(PARAM param) {
        NR = std::get<0>(param);
        NC = std::get<1>(param);
        rows = std::get<2>(param);
        cols = std::get<3>(param);
        vals = std::get<4>(param);
        return write_MatrixMarket(NR, NC, vals, rows, cols);
    }
};

typedef std::vector<int> IntVec;
using MatrixMarketSimpleTest = MatrixMarketTest<std::tuple<int, int, IntVec, IntVec, IntVec> >;

TEST_P(MatrixMarketSimpleTest, LayeredLoaderSimple) {
    auto filepath = extra_assemble(GetParam());
    ASSERT_EQ(rows.size(), cols.size());
    ASSERT_EQ(vals.size(), cols.size());

    // Loading it into an assignment object.
    std::ifstream in(filepath);
    tatami::MatrixMarket::LineAssignments ass;
    std::string line;
    while (std::getline(in, line) && ass.preamble(line.c_str())) {}
    while (std::getline(in, line)) {
        ass.add(line.c_str());
    }
    ass.finish();
    in.close();

    EXPECT_EQ(std::accumulate(ass.lines_per_category.begin(), ass.lines_per_category.end(), 0), rows.size());
    EXPECT_EQ(ass.lines_per_category.size(), 3);

    EXPECT_EQ(std::accumulate(ass.rows_per_category.begin(), ass.rows_per_category.end(), 0), NR);
    EXPECT_EQ(ass.rows_per_category.size(), 3);

    auto copy = ass.permutation;
    EXPECT_EQ(copy.size(), NR);
    std::sort(copy.begin(), copy.end()); // 0->(NR-1);
    for (size_t i = 0; i < copy.size(); ++i) {
        EXPECT_EQ(i, copy[i]);
    }

    // Checking against a reference.
    auto loaded = tatami::MatrixMarket::load_layered_sparse_matrix(filepath.c_str());
    const auto& out = loaded.matrix;

    auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new tatami::CompressedSparseColumnMatrix<double, int,
                                                                                               decltype(vals),
                                                                                               decltype(rows),
                                                                                               decltype(indptrs)
                                                                                              >(NR, NC, std::move(vals), std::move(rows), std::move(indptrs))); 
    EXPECT_EQ(out->nrow(), ref->nrow());
    EXPECT_EQ(out->ncol(), ref->ncol());
    EXPECT_TRUE(out->sparse());
    EXPECT_FALSE(out->prefer_rows());

    for (size_t i = 0; i < NR; ++i) {
        int adjusted = loaded.permutation[i];
        auto stuff = out->row(adjusted);
        EXPECT_EQ(stuff, ref->row(i));

        auto maxed = *std::max_element(stuff.begin(), stuff.end());
        if (ass.category[i] == 2) {
            EXPECT_TRUE(maxed > std::numeric_limits<uint16_t>::max());
        } else if (ass.category[i] == 0) {
            EXPECT_TRUE(maxed <= std::numeric_limits<uint8_t>::max());
        } else if (ass.category[i] == 1) {
            EXPECT_TRUE(maxed <= std::numeric_limits<uint16_t>::max());
            EXPECT_TRUE(maxed > std::numeric_limits<uint8_t>::max());
        }
    }
}

INSTANTIATE_TEST_CASE_P(
    MatrixMarket,
    MatrixMarketSimpleTest,
    ::testing::Values(
        std::make_tuple(
            // this example guarantees a few rows in each chunk.
            10,
            5,
            IntVec{ 1, 5, 8, 2, 9, 0, 4 }, 
            IntVec{ 2, 3, 1, 0, 2, 2, 4 },
            IntVec{ 0, 1, 10, 100, 1000, 10000, 100000 }
        ),
        std::make_tuple(
            // this example checks that the maximum category is used for each row.
            5,
            8,
            IntVec{ 1, 1, 2, 2, 3, 3, 4, 4 },
            IntVec{ 1, 7, 2, 5, 4, 3, 4, 0 },
            IntVec{ 10, 1, 10, 1000, 10000, 100000, 1, 100000 }
        ),
        std::make_tuple(
            // this example checks that we handle missing categories; in this case, only shorts.
            10,
            9,
            IntVec{ 1, 9, 7, 5, 3, 1, 3, 3, 7, 9 },
            IntVec{ 2, 4, 8, 8, 4, 6, 8, 6, 0, 0 },
            IntVec{ 1, 3, 2, 1, 7, 8, 9, 1, 1, 3 }
        ),
        std::make_tuple(
            // this example checks that we handle missing categories; in this case, only uint16.
            20,
            15,
            IntVec{ 15, 0, 4, 14, 0, 19, 19, 8, 11, 18, 2,  3, 6,  4, 9,  3, 16,  4, 13, 12 },
            IntVec{  3, 3, 4,  2, 8, 12,  3, 6,  2,  3, 2, 11, 1, 11, 5, 12,  7, 12,  5,  0 },
            IntVec{ 1000, 3000, 2000, 1000, 7000, 10000, 9000, 1000, 600, 500, 382, 826, 992, 244, 138, 852, 400, 542, 980, 116 }
        ),
        std::make_tuple(
            // this example checks that we handle missing categories; in this case, only uint32.
            100,
            20,
            IntVec{ 27, 83, 85, 60, 17, 45, 62, 30, 98, 47 },
            IntVec{ 0, 12, 17, 3, 17, 0, 8, 3, 8, 8 },
            IntVec{ 130875, 673886, 405953, 989598, 981526, 794394, 680144, 553105, 277529, 540959 }
        ),
        std::make_tuple(
            // this example checks that we handle empties.
            10,
            9,
            IntVec{},
            IntVec{},
            IntVec{}
        )
    )
);


