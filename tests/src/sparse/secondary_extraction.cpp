#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/sparse/secondary_extraction.hpp"
#include "tatami/utils/ElementType.hpp"

class SparseSecondaryExtractionCacheTest : public ::testing::Test {
protected:
    template<typename Index_, class IndexStorage_, class PointerStorage_>
    struct ServeIndices {
        ServeIndices(const IndexStorage_& i, const PointerStorage_& p) : indices(i), indptr(p) {}
        const IndexStorage_& indices;
        const PointerStorage_& indptr;

    public:
        typedef tatami::ElementType<PointerStorage_> pointer_type;

        pointer_type start_offset(Index_ primary) const {
            return indptr[primary];
        }

        pointer_type end_offset(Index_ primary) const {
            return indptr[primary + 1];
        }

        auto raw(Index_) const {
            return indices.begin();
        }
    };

    template<typename Index_, class IndexStorage_, class PointerStorage_>
    auto mock_cache(const IndexStorage_& i, const PointerStorage_& p, Index_ sec) {
        ServeIndices<Index_, IndexStorage_, PointerStorage_> server(i, p);
        return tatami::sparse_utils::FullSecondaryExtractionCache<Index_, decltype(server)>(std::move(server), sec, p.size() - 1);
    }

protected:
    template<typename Index_, class IndexVectorStorage_>
    struct ServeFragmentedIndices {
        ServeFragmentedIndices(const IndexVectorStorage_& i) : indices(i) {}
        const IndexVectorStorage_& indices;

    public:
        typedef size_t pointer_type;

        pointer_type start_offset(Index_ primary) const {
            return 0;
        }

        pointer_type end_offset(Index_ primary) const {
            return indices[primary].size();
        }

        auto raw(Index_ primary) const {
            return indices[primary].begin();
        }
    };

    template<typename Index_, class IndexVectorStorage_>
    auto mock_fragmented_cache(const IndexVectorStorage_& i, Index_ sec) {
        ServeFragmentedIndices<Index_, IndexVectorStorage_> server(i);
        return tatami::sparse_utils::FullSecondaryExtractionCache<Index_, decltype(server)>(std::move(server), sec, i.size());
    }

protected:
    std::vector<size_t> results;

    void SetUp() {
        results.resize(3);
    }

    void reset() {
        std::fill(results.begin(), results.end(), -1);
    }

    static std::vector<size_t> expected(size_t a, size_t b, size_t c) { 
        return std::vector<size_t>{ a, b, c }; 
    };
};

TEST_F(SparseSecondaryExtractionCacheTest, Increment) {
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

    auto store_fun = [&](int, int ip, size_t offset) { results[ip] = offset; };

    // Checking consecutive or semi-consecutive increments.
    {
        auto test = mock_cache<int>(indices, indptrs, 19);

        reset();
        EXPECT_TRUE(test.search(0, store_fun));
        EXPECT_EQ(results, expected(0, -1, -1));

        reset();
        EXPECT_TRUE(test.search(1, store_fun));
        EXPECT_EQ(results, expected(-1, 5, -1));

        reset();
        EXPECT_TRUE(test.search(3, store_fun));
        EXPECT_EQ(results, expected(-1, -1, -1));

        reset();
        EXPECT_TRUE(test.search(5, store_fun));
        EXPECT_EQ(results, expected(1, -1, 9));

        reset();
        EXPECT_TRUE(test.search(6, store_fun));
        EXPECT_EQ(results, expected(2, -1, -1));
    }

    // Checking a big jump that triggers a binary search.
    {
        auto test = mock_cache<int>(indices, indptrs, 19);

        reset();
        EXPECT_TRUE(test.search(15, store_fun));
        EXPECT_EQ(results, expected(-1, -1, 14));

        // subsequent consecutive searches still work.
        reset();
        EXPECT_TRUE(test.search(16, store_fun));
        EXPECT_EQ(results, expected(-1, -1, -1));

        reset();
        EXPECT_TRUE(test.search(17, store_fun));
        EXPECT_EQ(results, expected(-1, -1, 15));
    }

    // Checking a direct big jump to the end.
    {
        auto test = mock_cache<int>(indices, indptrs, 19);

        reset();
        EXPECT_TRUE(test.search(18, store_fun));
        EXPECT_EQ(results, expected(-1, 8, -1));
    }

    // Short-circuits work correctly.
    {
        auto test = mock_cache<int>(indices, indptrs, 19);

        reset();
        EXPECT_TRUE(test.search(1, store_fun));
        EXPECT_EQ(results, expected(-1, 5, -1));

        reset();
        EXPECT_TRUE(test.search(2, store_fun)); 
        EXPECT_EQ(results, expected(-1, -1, -1));

        EXPECT_FALSE(test.search(2, store_fun)); 
        EXPECT_FALSE(test.search(3, store_fun)); 
        EXPECT_FALSE(test.search(4, store_fun)); 
    }

    // Repeated requests are honored.
    {
        auto test = mock_cache<int>(indices, indptrs, 19);

        reset();
        EXPECT_TRUE(test.search(5, store_fun));
        auto ex = expected(1, -1, 9);
        EXPECT_EQ(results, ex);

        reset();
        EXPECT_TRUE(test.search(5, store_fun));
        EXPECT_EQ(results, ex);

        reset();
        EXPECT_TRUE(test.search(18, store_fun));
        ex = expected(-1, 8, -1);
        EXPECT_EQ(results, ex);

        reset();
        EXPECT_TRUE(test.search(18, store_fun));
        EXPECT_EQ(results, ex);
    }
}

TEST_F(SparseSecondaryExtractionCacheTest, Decrement) {
    std::vector<int> indices {
        2, 8, 11, 12, 18,
        3, 4, 7, 10, 12, 15,  
        0, 1, 4, 8, 14
    };

    std::vector<size_t> indptrs {
        0,
        5,
        11,
        16
    };

    auto store_fun = [&](int, int i, size_t p) { results[i] = p; };

    // Checking consecutive or semi-consecutive decrements.
    {
        auto test = mock_cache<int>(indices, indptrs, 19);

        reset();
        EXPECT_TRUE(test.search(18, store_fun));
        EXPECT_EQ(results, expected(4, -1, -1));

        reset();
        EXPECT_TRUE(test.search(17, store_fun));
        EXPECT_EQ(results, expected(-1, -1, -1));

        reset();
        EXPECT_TRUE(test.search(15, store_fun));
        EXPECT_EQ(results, expected(-1, 10, -1));

        reset();
        EXPECT_TRUE(test.search(14, store_fun));
        EXPECT_EQ(results, expected(-1, -1, 15));

        reset();
        EXPECT_TRUE(test.search(12, store_fun));
        EXPECT_EQ(results, expected(3, 9, -1));
    }

    // Checking the jumps.
    {
        auto test = mock_cache<int>(indices, indptrs, 19);
        EXPECT_TRUE(test.search(18, store_fun));

        reset();
        EXPECT_TRUE(test.search(10, store_fun));
        EXPECT_EQ(results, expected(-1, 8, -1));

        reset();
        EXPECT_TRUE(test.search(1, store_fun));
        EXPECT_EQ(results, expected(-1, -1, 12));
    }

    // Big jump to zero.
    {
        auto test = mock_cache<int>(indices, indptrs, 19);
        EXPECT_TRUE(test.search(18, store_fun));

        reset();
        EXPECT_TRUE(test.search(0, store_fun));
        EXPECT_EQ(results, expected(-1, -1, 11));
    }

    // Short-circuits work correctly.
    {
        auto test = mock_cache<int>(indices, indptrs, 19);
        EXPECT_TRUE(test.search(18, store_fun));
        EXPECT_TRUE(test.search(17, store_fun));
        EXPECT_FALSE(test.search(16, store_fun)); 
    }

    // Unsigned index type is handled correctly.
    {
        std::vector<unsigned char> indices2(indices.begin(), indices.end());
        auto test = mock_cache<int>(indices2, indptrs, 19);

        EXPECT_TRUE(test.search(18, store_fun));
        EXPECT_TRUE(test.search(17, store_fun));

        // Short circuits correctly.
        EXPECT_FALSE(test.search(16, store_fun));

        // Then continues working.
        reset();
        EXPECT_TRUE(test.search(15, store_fun)); 
        EXPECT_EQ(results, expected(-1, 10, -1));

        // Jump and then more decrements.
        reset();
        EXPECT_TRUE(test.search(7, store_fun)); 
        EXPECT_EQ(results, expected(-1, 7, -1));

        reset();
        EXPECT_TRUE(test.search(6, store_fun)); 
        EXPECT_EQ(results, expected(-1, -1, -1));

        EXPECT_FALSE(test.search(5, store_fun)); // short-circuit.

        reset();
        EXPECT_TRUE(test.search(4, store_fun)); 
        EXPECT_EQ(results, expected(-1, 6, 13));

        reset();
        EXPECT_TRUE(test.search(2, store_fun)); 
        EXPECT_EQ(results, expected(0, -1, -1));

        reset();
        EXPECT_TRUE(test.search(0, store_fun)); 
        EXPECT_EQ(results, expected(-1, -1, 11));
    }
}

TEST_F(SparseSecondaryExtractionCacheTest, Alternating) {
    std::vector<int> indices {
        0, 1, 3, 5, 6, 18,
        2, 6, 8, 9, 10, 11, 17, 18,
        4, 6, 7, 8, 10, 11, 14, 15, 16
    };

    std::vector<int> indptrs { // using int pointers, for some variety.
        0,
        6,
        14,
        23
    };

    auto store_fun = [&](int, int i, size_t p) { results[i] = p; };

    {
        auto test = mock_cache<int>(indices, indptrs, 19);

        // Jump up, followed by decrements.
        reset();
        EXPECT_TRUE(test.search(15, store_fun));
        EXPECT_EQ(results, expected(-1, -1, 21));

        reset();
        EXPECT_TRUE(test.search(14, store_fun));
        EXPECT_EQ(results, expected(-1, -1, 20));

        reset();
        EXPECT_TRUE(test.search(11, store_fun));
        EXPECT_EQ(results, expected(-1, 11, 19));

        // Jumps down, followed by increments.
        reset();
        EXPECT_TRUE(test.search(4, store_fun));
        EXPECT_EQ(results, expected(-1, -1, 14));

        reset();
        EXPECT_TRUE(test.search(6, store_fun));
        EXPECT_EQ(results, expected(4, 7, 15));

        reset();
        EXPECT_TRUE(test.search(8, store_fun));
        EXPECT_EQ(results, expected(-1, 8, 17));
    }

    {
        auto test = mock_cache<int>(indices, indptrs, 19);

        // Jump to end, followed by decrements.
        reset();
        EXPECT_TRUE(test.search(18, store_fun));
        EXPECT_EQ(results, expected(5, 13, -1));

        reset();
        EXPECT_TRUE(test.search(17, store_fun));
        EXPECT_EQ(results, expected(-1, 12, -1));

        // Jump to start, followed by increments.
        reset();
        EXPECT_TRUE(test.search(0, store_fun));
        EXPECT_EQ(results, expected(0, -1, -1));

        reset();
        EXPECT_TRUE(test.search(1, store_fun));
        EXPECT_EQ(results, expected(1, -1, -1));

        reset();
        EXPECT_TRUE(test.search(2, store_fun));
        EXPECT_EQ(results, expected(-1, 6, -1));

        // Jump to end, followed by decrements, and then a repeated element.
        EXPECT_TRUE(test.search(18, store_fun));

        reset();
        EXPECT_TRUE(test.search(17, store_fun));
        EXPECT_EQ(results, expected(-1, 12, -1));

        reset();
        EXPECT_TRUE(test.search(17, store_fun));
        EXPECT_EQ(results, expected(-1, 12, -1));
    }
}

TEST_F(SparseSecondaryExtractionCacheTest, Empty) {
    std::vector<int> indices {};

    std::vector<size_t> indptrs { 0, 0, 0, 0 };

    auto store_fun = [&](int, int i, size_t p) { results[i] = p; };

    auto test = mock_cache<int>(indices, indptrs, 19);
    EXPECT_FALSE(test.search(0, store_fun));

    // Increments returns false.
    EXPECT_FALSE(test.search(1, store_fun));
    EXPECT_FALSE(test.search(2, store_fun));
    EXPECT_FALSE(test.search(10, store_fun));

    // Big jump to the end.
    EXPECT_FALSE(test.search(18, store_fun));

    // First decrement switches to last_increasing = false, so it can't short-circuit.
    reset();
    EXPECT_TRUE(test.search(17, store_fun));
    EXPECT_EQ(results, expected(-1, -1, -1));

    // Next decrements short-circuit properly.
    EXPECT_FALSE(test.search(16, store_fun));
    EXPECT_FALSE(test.search(15, store_fun));
    EXPECT_FALSE(test.search(8, store_fun));

    // Big jump back to the front.
    EXPECT_FALSE(test.search(0, store_fun));

    // First increment switches back to last_increasing = true and can't short-circuit.
    reset();
    EXPECT_TRUE(test.search(1, store_fun));
    EXPECT_EQ(results, expected(-1, -1, -1));
}

TEST_F(SparseSecondaryExtractionCacheTest, Fragmented) {
    std::vector<std::vector<int> > indices {
        { 1, 2, 7, 9, 11, 15 },
        { 0, 5, 7, 14, 18 },
        { 3, 8, 10, 13, 16 } 
    };

    auto store_fun = [&](int, int i, size_t p) { results[i] = p; };

    // Increments.
    {
        auto test = mock_fragmented_cache(indices, 19);

        reset();
        EXPECT_TRUE(test.search(0, store_fun));
        EXPECT_EQ(results, expected(-1, 0, -1));

        reset();
        EXPECT_TRUE(test.search(1, store_fun));
        EXPECT_EQ(results, expected(0, -1, -1));

        reset();
        EXPECT_TRUE(test.search(2, store_fun));
        EXPECT_EQ(results, expected(1, -1, -1));

        reset();
        EXPECT_TRUE(test.search(7, store_fun)); // big jump.
        EXPECT_EQ(results, expected(2, 2, -1));

        reset();
        EXPECT_TRUE(test.search(15, store_fun)); // another big jump.
        EXPECT_EQ(results, expected(5, -1, -1));
    }

    {
        auto test = mock_fragmented_cache(indices, 19);

        reset();
        EXPECT_TRUE(test.search(10, store_fun)); // jumps work correctly.
        EXPECT_EQ(results, expected(-1, -1, 2));

        reset();
        EXPECT_TRUE(test.search(18, store_fun)); // jump to end works correctly.
        EXPECT_EQ(results, expected(-1, 4, -1));
    }

    // Decrement works correctly.
    {
        auto test = mock_fragmented_cache(indices, 19);

        reset();
        EXPECT_TRUE(test.search(18, store_fun)); // jump to end works correctly.
        EXPECT_EQ(results, expected(-1, 4, -1));

        reset();
        EXPECT_TRUE(test.search(15, store_fun)); // decrement.
        EXPECT_EQ(results, expected(5, -1, -1));
    }
}
