#include "tatami_test/tatami_test.hpp"
#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/CompressedSparseMatrix.hpp"
#include "tatami/utils/make_Wrapper.hpp"

class MakeDenseWrapperTest : public ::testing::TestWithParam<bool> {
protected:
    static inline const int NR = 10, NC = 20;
    static inline std::vector<double> contents;

    static void SetUpTestSuite() {
        contents.resize(NR * NC);
        double counter = -105;
        for (auto& i : contents) {
            i = counter++;
        }
    }
};

TEST_P(MakeDenseWrapperTest, Myopic) {
    auto row = GetParam();
    tatami::DenseMatrix<double, int, std::vector<double> > mat(NR, NC, contents, row);

    const auto target = (row ? NR : NC); 
    const auto nontarget = (row ? NC : NR); 

    {
        auto wrapped = tatami::make_DenseWrapper(mat, row, tatami::Options());
        EXPECT_EQ(wrapped.length(), nontarget);

        const int chosen = target / 2;
        auto res = wrapped.fetch(chosen);
        for (int i = 0; i < nontarget; ++i) {
            EXPECT_EQ(res[i], contents[i + chosen * nontarget]);
        }
    }

    {
        const int bstart = nontarget / 4, blength = (nontarget * 3) / 4;
        auto wrapped = tatami::make_DenseWrapper(mat, row, bstart, blength, tatami::Options());
        EXPECT_EQ(wrapped.length(), blength);

        const int chosen = target / 4;
        auto res = wrapped.fetch(chosen);
        for (int i = 0; i < blength; ++i) {
            EXPECT_EQ(res[i], contents[i + bstart + chosen * nontarget]);
        }
    }

    {
        auto indices = tatami_test::define_non_target_indices(nontarget, 0.1, 0.8, 0.2, 12345);
        auto wrapped = tatami::make_DenseWrapper(mat, row, indices, tatami::Options());
        EXPECT_EQ(wrapped.length(), indices->size());

        const int chosen = (target * 4) / 5 ;
        auto res = wrapped.fetch(chosen);
        for (std::size_t i = 0, end = indices->size(); i < end; ++i) {
            EXPECT_EQ(res[i], contents[(*indices)[i] + chosen * nontarget]);
        }
    }
}

TEST_P(MakeDenseWrapperTest, Oracular) {
    auto row = GetParam();
    tatami::DenseMatrix<double, int, std::vector<double> > mat(NR, NC, contents, row);

    const auto target = (row ? NR : NC); 
    const auto nontarget = (row ? NC : NR); 
    auto oracle = tatami::make_ConsecutiveOracle<int>(target / 3, target / 3);

    {
        auto wrapped = tatami::make_DenseWrapper(mat, row, oracle, tatami::Options());
        EXPECT_EQ(wrapped.length(), nontarget);

        auto res = wrapped.next();
        const auto chosen = oracle->get(0);
        for (int i = 0; i < nontarget; ++i) {
            EXPECT_EQ(res[i], contents[i + chosen * nontarget]);
        }
    }

    {
        const int bstart = nontarget / 4, blength = (nontarget * 3) / 4;
        auto wrapped = tatami::make_DenseWrapper(mat, row, oracle, bstart, blength, tatami::Options());
        EXPECT_EQ(wrapped.length(), blength);

        wrapped.next();
        auto res = wrapped.next();
        const auto chosen = oracle->get(1);
        for (int i = 0; i < blength; ++i) {
            EXPECT_EQ(res[i], contents[i + bstart + chosen * nontarget]);
        }
    }

    {
        auto indices = tatami_test::define_non_target_indices(nontarget, 0.1, 0.8, 0.2, 12345);
        auto wrapped = tatami::make_DenseWrapper(mat, row, oracle, indices, tatami::Options());
        EXPECT_EQ(wrapped.length(), indices->size());

        wrapped.next();
        wrapped.next();
        auto res = wrapped.next();
        const auto chosen = oracle->get(2);
        for (std::size_t i = 0, end = indices->size(); i < end; ++i) {
            EXPECT_EQ(res[i], contents[(*indices)[i] + chosen * nontarget]);
        }
    }
}

INSTANTIATE_TEST_SUITE_P(
    MakeWrapper,
    MakeDenseWrapperTest,
    ::testing::Values(false, true)
);

/**************************************/

class MakeSparseWrapperTest : public ::testing::TestWithParam<bool> {
protected:
    static inline const int target = 100, nontarget = 200;
    static inline tatami_test::SimulateCompressedSparseResult<double, int> simulated;

    static void SetUpTestSuite() {
        simulated = tatami_test::simulate_compressed_sparse<double, int>(target, nontarget, {});
    }

    static std::pair<std::vector<double>, int> create_block_mapping(int chosen, int bstart, int blen) {
        std::vector<double> expected(nontarget, std::numeric_limits<double>::quiet_NaN());
        int found = 0;
        const int bend = bstart + blen;
        for (std::size_t i = simulated.indptr[chosen], end = simulated.indptr[chosen + 1]; i < end; ++i) {
            if (simulated.index[i] >= bstart && simulated.index[i] < bend) {
                ++found;
                expected[simulated.index[i]] = simulated.data[i];
            }
        }
        return std::make_pair(expected, found);
    }

    static std::pair<std::vector<double>, int> create_index_mapping(int chosen, const std::vector<int>& indices) {
        std::vector<char> required(nontarget);
        for (auto i : indices) {
            required[i] = 1;
        }

        std::vector<double> expected(nontarget, std::numeric_limits<double>::quiet_NaN());
        int found = 0;
        for (std::size_t i = simulated.indptr[chosen], end = simulated.indptr[chosen + 1]; i < end; ++i) {
            if (required[simulated.index[i]]) {
                ++found;
                expected[simulated.index[i]] = simulated.data[i];
            }
        }
        return std::make_pair(expected, found);
    }
};


TEST_P(MakeSparseWrapperTest, Myopic) {
    auto row = GetParam();
    const auto NR = (row ? target : nontarget); 
    const auto NC = (row ? nontarget : target); 

    tatami::CompressedSparseMatrix<double, int, std::vector<double>, std::vector<int>, std::vector<std::size_t> > mat(
        NR,
        NC,
        simulated.data,
        simulated.index,
        simulated.indptr,
        row
    );

    {
        auto wrapped = tatami::make_SparseWrapper(mat, row, tatami::Options());
        EXPECT_EQ(wrapped.length(), nontarget);

        const int chosen = target / 2;
        auto res = wrapped.fetch(chosen);
        EXPECT_EQ(res.number, simulated.indptr[chosen + 1] - simulated.indptr[chosen]);
        for (int i = 0; i < res.number; ++i) {
            EXPECT_EQ(res.value[i], simulated.data[i + simulated.indptr[chosen]]);
            EXPECT_EQ(res.index[i], simulated.index[i + simulated.indptr[chosen]]);
        }

        tatami::Options opt;
        opt.sparse_extract_index = false;
        opt.sparse_extract_value = false;
        auto wrapped2 = tatami::make_SparseWrapper(mat, row, opt);
        auto res2 = wrapped2.fetch(chosen);
        EXPECT_EQ(res.number, res2.number);
        EXPECT_TRUE(res2.value == NULL);
        EXPECT_TRUE(res2.index == NULL);
    }

    {
        const int bstart = nontarget / 4, blength = nontarget / 2;
        auto wrapped = tatami::make_SparseWrapper(mat, row, bstart, blength, tatami::Options());
        EXPECT_EQ(wrapped.length(), blength);

        const int chosen = (target * 3) / 4;
        auto mapping = create_block_mapping(chosen, bstart, blength);
        auto res = wrapped.fetch(chosen);
        for (int i = 0; i < res.number; ++i) {
            EXPECT_EQ(res.value[i], mapping.first[res.index[i]]);
        }
        EXPECT_EQ(res.number, mapping.second);
    }

    {
        auto indices = tatami_test::define_non_target_indices(nontarget, 0.1, 0.8, 0.2, 12345);
        auto wrapped = tatami::make_SparseWrapper(mat, row, indices, tatami::Options());
        EXPECT_EQ(wrapped.length(), indices->size());

        const int chosen = target / 4;
        auto mapping = create_index_mapping(chosen, *indices);
        auto res = wrapped.fetch(chosen);
        for (int i = 0; i < res.number; ++i) {
            EXPECT_EQ(res.value[i], mapping.first[res.index[i]]);
        }
        EXPECT_EQ(res.number, mapping.second);
    }
}

TEST_P(MakeSparseWrapperTest, Oracular) {
    auto row = GetParam();
    const auto NR = (row ? target : nontarget); 
    const auto NC = (row ? nontarget : target); 

    tatami::CompressedSparseMatrix<double, int, std::vector<double>, std::vector<int>, std::vector<std::size_t> > mat(
        NR,
        NC,
        simulated.data,
        simulated.index,
        simulated.indptr,
        row
    );

    auto oracle = tatami::make_ConsecutiveOracle<int>(target / 4, target / 2);

    {
        auto wrapped = tatami::make_SparseWrapper(mat, row, oracle, tatami::Options());
        EXPECT_EQ(wrapped.length(), nontarget);

        wrapped.next();
        auto res = wrapped.next();
        const int chosen = oracle->get(1);
        EXPECT_EQ(res.number, simulated.indptr[chosen + 1] - simulated.indptr[chosen]);
        for (int i = 0; i < res.number; ++i) {
            EXPECT_EQ(res.value[i], simulated.data[i + simulated.indptr[chosen]]);
            EXPECT_EQ(res.index[i], simulated.index[i + simulated.indptr[chosen]]);
        }

        tatami::Options opt;
        opt.sparse_extract_index = false;
        opt.sparse_extract_value = false;
        auto wrapped2 = tatami::make_SparseWrapper(mat, row, oracle, opt);
        wrapped2.next();
        auto res2 = wrapped2.next();
        EXPECT_EQ(res.number, res2.number);
        EXPECT_TRUE(res2.value == NULL);
        EXPECT_TRUE(res2.index == NULL);
    }

    {
        const int bstart = nontarget / 5, blength = (nontarget * 3) / 5;
        auto wrapped = tatami::make_SparseWrapper(mat, row, oracle, bstart, blength, tatami::Options());
        EXPECT_EQ(wrapped.length(), blength);

        auto res = wrapped.next();
        const int chosen = oracle->get(0);
        auto mapping = create_block_mapping(chosen, bstart, blength);
        for (int i = 0; i < res.number; ++i) {
            EXPECT_EQ(res.value[i], mapping.first[res.index[i]]);
        }
        EXPECT_EQ(res.number, mapping.second);
    }

    {
        auto indices = tatami_test::define_non_target_indices(nontarget, 0.2, 0.6, 0.3, 98765);
        auto wrapped = tatami::make_SparseWrapper(mat, row, oracle, indices, tatami::Options());
        EXPECT_EQ(wrapped.length(), indices->size());

        wrapped.next(); wrapped.next();
        auto res = wrapped.next();
        const int chosen = oracle->get(2);
        auto mapping = create_index_mapping(chosen, *indices);
        for (int i = 0; i < res.number; ++i) {
            EXPECT_EQ(res.value[i], mapping.first[res.index[i]]);
        }
        EXPECT_EQ(res.number, mapping.second);
    }
}

INSTANTIATE_TEST_SUITE_P(
    MakeWrapper,
    MakeSparseWrapperTest,
    ::testing::Values(false, true)
);
