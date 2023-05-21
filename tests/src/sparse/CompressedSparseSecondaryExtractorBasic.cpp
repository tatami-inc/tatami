#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/sparse/CompressedSparseSecondaryExtractorBasic.hpp"

class CompressedSparseSecondaryExtractorTest : public ::testing::Test {
protected:
    struct SimpleModifier {
        template<typename StoredIndex_>
        static void increment(size_t& ptr, const std::vector<StoredIndex_>&, const size_t&) { ++ptr; }

        template<typename StoredIndex_>
        static void decrement(size_t& ptr, const std::vector<StoredIndex_>&, const size_t&) { --ptr; }

        static size_t get(size_t ptr) { return ptr; }
        static void set(size_t& ptr, size_t val) { ptr = val; }
    };

    size_t n = 3;
    std::vector<size_t> results;

    void SetUp() {
        results.resize(3);
    }

    static int identity(int i) { 
        return i; 
    };

    static std::vector<size_t> expected(size_t a, size_t b, size_t c) { 
        return std::vector<size_t>{ a, b, c }; 
    };
};

TEST_F(CompressedSparseSecondaryExtractorTest, Increment) {
    std::vector<int> indices {
        0, 5, 6, 9, 10,
        1, 8, 12, 18,
        5, 9, 10, 12, 13, 15, 17
    };

    std::vector<size_t> indptrs {
        0,
        5,
        9,
        16
    };

    auto store_fun = [&](int i, size_t p) { results[i] = p; };
    auto skip_fun = [&](int i) { results[i] = -1; };

    // Checking consecutive or semi-consecutive increments.
    {
        tatami::CompressedSparseSecondaryExtractorBasic<int, int, size_t, SimpleModifier> test(19, indices, indptrs);

        EXPECT_TRUE(test.search(0, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(0, -1, -1));

        EXPECT_TRUE(test.search(1, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, 5, -1));

        EXPECT_TRUE(test.search(3, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, -1, -1));

        EXPECT_TRUE(test.search(5, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(1, -1, 9));

        EXPECT_TRUE(test.search(6, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(2, -1, -1));
    }

    // Checking a big jump that triggers a binary search.
    {
        tatami::CompressedSparseSecondaryExtractorBasic<int, int, size_t, SimpleModifier> test(19, indices, indptrs);

        EXPECT_TRUE(test.search(15, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, -1, 14));

        // subsequent consecutive searches still work.
        EXPECT_TRUE(test.search(16, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, -1, -1));

        EXPECT_TRUE(test.search(17, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, -1, 15));
    }

    // Checking a direct big jump to the end.
    {
        tatami::CompressedSparseSecondaryExtractorBasic<int, int, size_t, SimpleModifier> test(19, indices, indptrs);

        EXPECT_TRUE(test.search(18, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, 8, -1));
    }

    // Short-circuits work correctly.
    {
        tatami::CompressedSparseSecondaryExtractorBasic<int, int, size_t, SimpleModifier> test(19, indices, indptrs);

        EXPECT_TRUE(test.search(1, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, 5, -1));

        EXPECT_TRUE(test.search(2, n, identity, indices, indptrs, store_fun, skip_fun)); 
        EXPECT_EQ(results, expected(-1, -1, -1));

        EXPECT_FALSE(test.search(2, n, identity, indices, indptrs, store_fun, skip_fun)); 
        EXPECT_FALSE(test.search(3, n, identity, indices, indptrs, store_fun, skip_fun)); 
        EXPECT_FALSE(test.search(4, n, identity, indices, indptrs, store_fun, skip_fun)); 
    }

    // Repeated requests are honored.
    {
        tatami::CompressedSparseSecondaryExtractorBasic<int, int, size_t, SimpleModifier> test(19, indices, indptrs);

        EXPECT_TRUE(test.search(5, n, identity, indices, indptrs, store_fun, skip_fun));
        auto ex = expected(1, -1, 9);
        EXPECT_EQ(results, ex);

        EXPECT_TRUE(test.search(5, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, ex);

        EXPECT_TRUE(test.search(18, n, identity, indices, indptrs, store_fun, skip_fun));
        ex = expected(-1, 8, -1);
        EXPECT_EQ(results, ex);

        EXPECT_TRUE(test.search(18, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, ex);
    }
}

TEST_F(CompressedSparseSecondaryExtractorTest, Decrement) {
    std::vector<int> indices {
        2, 8, 11, 12, 18,
        3, 4, 7, 10, 12, 15,  
        0, 1, 4, 8, 14
    };

    std::vector<int> indptrs {
        0,
        5,
        11,
        16
    };

    auto store_fun = [&](int i, size_t p) { results[i] = p; };
    auto skip_fun = [&](int i) { results[i] = -1; };

    // Checking consecutive or semi-consecutive decrements.
    {
        tatami::CompressedSparseSecondaryExtractorBasic<int, int, size_t, SimpleModifier> test(19, indices, indptrs);

        EXPECT_TRUE(test.search(18, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(4, -1, -1));

        EXPECT_TRUE(test.search(17, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, -1, -1));

        EXPECT_TRUE(test.search(15, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, 10, -1));

        EXPECT_TRUE(test.search(14, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, -1, 15));

        EXPECT_TRUE(test.search(12, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(3, 9, -1));
    }

    // Checking the jumps.
    {
        tatami::CompressedSparseSecondaryExtractorBasic<int, int, size_t, SimpleModifier> test(19, indices, indptrs);
        EXPECT_TRUE(test.search(18, n, identity, indices, indptrs, store_fun, skip_fun));

        EXPECT_TRUE(test.search(10, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, 8, -1));

        EXPECT_TRUE(test.search(1, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, -1, 12));
    }

    // Big jump to zero.
    {
        tatami::CompressedSparseSecondaryExtractorBasic<int, int, size_t, SimpleModifier> test(19, indices, indptrs);
        EXPECT_TRUE(test.search(18, n, identity, indices, indptrs, store_fun, skip_fun));

        EXPECT_TRUE(test.search(0, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, -1, 11));
    }

    // Short-circuits work correctly.
    {
        tatami::CompressedSparseSecondaryExtractorBasic<int, int, size_t, SimpleModifier> test(19, indices, indptrs);
        EXPECT_TRUE(test.search(18, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_TRUE(test.search(17, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_FALSE(test.search(16, n, identity, indices, indptrs, store_fun, skip_fun)); 
    }

    // Unsigned index type is handled correctly.
    {
        std::vector<unsigned char> indices2(indices.begin(), indices.end());
        tatami::CompressedSparseSecondaryExtractorBasic<int, unsigned char, size_t, SimpleModifier> test(19, indices2, indptrs);
        EXPECT_TRUE(test.search(18, n, identity, indices2, indptrs, store_fun, skip_fun));
        EXPECT_TRUE(test.search(17, n, identity, indices2, indptrs, store_fun, skip_fun));

        // Short circuits correctly.
        EXPECT_FALSE(test.search(16, n, identity, indices2, indptrs, store_fun, skip_fun));

        // Then continues working.
        EXPECT_TRUE(test.search(15, n, identity, indices2, indptrs, store_fun, skip_fun)); 
        EXPECT_EQ(results, expected(-1, 10, -1));

        // Jump and then more decrements.
        EXPECT_TRUE(test.search(7, n, identity, indices2, indptrs, store_fun, skip_fun)); 
        EXPECT_EQ(results, expected(-1, 7, -1));

        EXPECT_FALSE(test.search(6, n, identity, indices2, indptrs, store_fun, skip_fun)); 
        EXPECT_FALSE(test.search(5, n, identity, indices2, indptrs, store_fun, skip_fun));

        EXPECT_TRUE(test.search(4, n, identity, indices2, indptrs, store_fun, skip_fun)); 
        EXPECT_EQ(results, expected(-1, 6, 13));

        EXPECT_TRUE(test.search(2, n, identity, indices2, indptrs, store_fun, skip_fun)); 
        EXPECT_EQ(results, expected(0, -1, -1));

        EXPECT_TRUE(test.search(0, n, identity, indices2, indptrs, store_fun, skip_fun)); 
        EXPECT_EQ(results, expected(-1, -1, 11));
    }
}

TEST_F(CompressedSparseSecondaryExtractorTest,  Alternating) {
    std::vector<int> indices {
        0, 1, 3, 5, 6, 18,
        2, 6, 8, 9, 10, 11, 17, 18,
        4, 6, 7, 8, 10, 11, 14, 15, 16
    };

    std::vector<int> indptrs {
        0,
        6,
        14,
        23
    };

    auto store_fun = [&](int i, size_t p) { results[i] = p; };
    auto skip_fun = [&](int i) { results[i] = -1; };

    {
        tatami::CompressedSparseSecondaryExtractorBasic<int, int, size_t, SimpleModifier> test(19, indices, indptrs);

        // Jump up, followed by decrements.
        EXPECT_TRUE(test.search(15, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, -1, 21));

        EXPECT_TRUE(test.search(14, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, -1, 20));

        EXPECT_TRUE(test.search(11, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, 11, 19));

        // Jumps down, followed by increments
        EXPECT_TRUE(test.search(4, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, -1, 14));

        EXPECT_TRUE(test.search(6, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(4, 7, 15));

        EXPECT_TRUE(test.search(8, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, 8, 17));
    }

    {
        tatami::CompressedSparseSecondaryExtractorBasic<int, int, size_t, SimpleModifier> test(19, indices, indptrs);

        // Jump to end, followed by decrements.
        EXPECT_TRUE(test.search(18, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(5, 13, -1));

        EXPECT_TRUE(test.search(17, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, 12, -1));

        // Jump to start, followed by increments.
        EXPECT_TRUE(test.search(0, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(0, -1, -1));

        EXPECT_TRUE(test.search(1, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(1, -1, -1));

        EXPECT_TRUE(test.search(2, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, 6, -1));

        // Jump to end, followed by decrements, and then a repeated element (which is processed using the increment functionality).
        EXPECT_TRUE(test.search(18, n, identity, indices, indptrs, store_fun, skip_fun));

        EXPECT_TRUE(test.search(17, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, 12, -1));

        EXPECT_TRUE(test.search(17, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, 12, -1));
    }
}

TEST_F(CompressedSparseSecondaryExtractorTest,  Duplicated) {
    struct DuplicatedModifier {
        static void increment(size_t& ptr, const std::vector<int>& indices, size_t limit) {
            auto current = indices[ptr];
            while (ptr < limit && indices[ptr] == current) { ++ptr; }
        }

        static void decrement(size_t& ptr, const std::vector<int>& indices, size_t limit) {
            if (ptr != limit) {
                --ptr;
                auto current = indices[ptr];
                while (ptr > limit && indices[ptr - 1] == current) {
                    --ptr;
                }
            }
        }

        static size_t get(const size_t& ptr) { return ptr; }

        static void set(size_t& ptr, size_t val) { ptr = val; }
    };

    std::vector<int> indices {
        0, 0, 0, 2, 6, 6, 9, 10, 15, 15, 15, 
        1, 1, 3, 3, 3, 5, 8, 8, 8, 11, 14, 18, 18, 18,
        2, 4, 7, 7, 8, 8, 12, 13, 15
    };

    std::vector<size_t> indptrs {
        0,
        11,
        25,
        34
    };

    auto store_fun = [&](int i, size_t p) { results[i] = p; };
    auto skip_fun = [&](int i) { results[i] = -1; };

    // Increments.
    {
        tatami::CompressedSparseSecondaryExtractorBasic<int, int, size_t, DuplicatedModifier> test(19, indices, indptrs);

        EXPECT_TRUE(test.search(0, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(0, -1, -1));

        EXPECT_TRUE(test.search(1, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(-1, 11, -1));

        EXPECT_TRUE(test.search(2, n, identity, indices, indptrs, store_fun, skip_fun));
        EXPECT_EQ(results, expected(3, -1, 25));

        EXPECT_TRUE(test.search(7, n, identity, indices, indptrs, store_fun, skip_fun)); // big jump.
        EXPECT_EQ(results, expected(-1, -1, 27));

        EXPECT_TRUE(test.search(15, n, identity, indices, indptrs, store_fun, skip_fun)); // another big jump.
        EXPECT_EQ(results, expected(8, -1, 33));
    }

    {
        tatami::CompressedSparseSecondaryExtractorBasic<int, int, size_t, DuplicatedModifier> test(19, indices, indptrs);

        EXPECT_TRUE(test.search(10, n, identity, indices, indptrs, store_fun, skip_fun)); // jumps work correctly.
        EXPECT_EQ(results, expected(7, -1, -1));

        EXPECT_TRUE(test.search(18, n, identity, indices, indptrs, store_fun, skip_fun)); // jump to end works correctly.
        EXPECT_EQ(results, expected(-1, 22, -1));
    }

    // Decrement works correctly.
    {
        tatami::CompressedSparseSecondaryExtractorBasic<int, int, size_t, DuplicatedModifier> test(19, indices, indptrs);
        EXPECT_TRUE(test.search(18, n, identity, indices, indptrs, store_fun, skip_fun)); // jump to end works correctly.
        EXPECT_EQ(results, expected(-1, 22, -1));

        EXPECT_TRUE(test.search(15, n, identity, indices, indptrs, store_fun, skip_fun)); // decrement.
        EXPECT_EQ(results, expected(8, -1, 33));
    }
}
