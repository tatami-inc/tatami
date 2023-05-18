#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/sparse/utils.hpp"

struct SimpleModifier {
    static void increment(size_t& ptr, const std::vector<int>&, const size_t&) { ++ptr; }
    static void decrement(size_t& ptr, const std::vector<int>&, const size_t&) { --ptr; }
    static size_t get(size_t ptr) { return ptr; }
    static void set(size_t& ptr, size_t val) { ptr = val; }
};

TEST(SecondaryExtractionWorkspaceBase, Increment) {
    tatami::sparse::SecondaryExtractionWorkspaceBase<int, SimpleModifier> test(19);

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

    {
        int curdex = 0;
        size_t curptr = 0;
        test.search_above(0, 0, indices, indptrs, curdex, curptr); // no-op.
        EXPECT_EQ(curdex, 0);
        EXPECT_EQ(curptr, 0);

        test.search_above(1, 0, indices, indptrs, curdex, curptr); // consecutive increase.
        EXPECT_EQ(curdex, 5);
        EXPECT_EQ(curptr, 1);

        test.search_above(3, 0, indices, indptrs, curdex, curptr); // no-op, curdex is still above the requested secondary index.
        EXPECT_EQ(curdex, 5);
        EXPECT_EQ(curptr, 1);

        test.search_above(5, 0, indices, indptrs, curdex, curptr); // hit!
        EXPECT_EQ(curdex, 5);
        EXPECT_EQ(curptr, 1);

        test.search_above(6, 0, indices, indptrs, curdex, curptr); // consecutive increase.
        EXPECT_EQ(curdex, 6);
        EXPECT_EQ(curptr, 2);
    }

    {
        int curdex = 1;
        size_t curptr = 5;
        test.search_above(0, 1, indices, indptrs, curdex, curptr); // no-op, curdex is still above the requested secondary index.
        EXPECT_EQ(curdex, 1);
        EXPECT_EQ(curptr, 5);

        test.search_above(15, 1, indices, indptrs, curdex, curptr); // big jump, triggers binary search.
        EXPECT_EQ(curdex, 18);
        EXPECT_EQ(curptr, 8);
    }

    {
        int curdex = 1;
        size_t curptr = 5;
        test.search_above(18, 1, indices, indptrs, curdex, curptr); // direct big jump to end.
        EXPECT_EQ(curdex, 18);
        EXPECT_EQ(curptr, 8);
    }

    {
        int curdex = 5;
        size_t curptr = 9;
        test.search_above(11, 2, indices, indptrs, curdex, curptr); // big jump (1)
        EXPECT_EQ(curdex, 12);
        EXPECT_EQ(curptr, 12);

        test.search_above(17, 2, indices, indptrs, curdex, curptr); // big jump (2)
        EXPECT_EQ(curdex, 17);
        EXPECT_EQ(curptr, 15);
    }

    {
        int curdex = 5;
        size_t curptr = 9;
        test.search_above(18, 2, indices, indptrs, curdex, curptr); // big jump to end, but without a match.
        EXPECT_EQ(curdex, 19);
        EXPECT_EQ(curptr, 16);
    }
}

TEST(SecondaryExtractionWorkspaceBase, Decrement) {
    tatami::sparse::SecondaryExtractionWorkspaceBase<int, SimpleModifier> test(19);

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

    {
        int curdex = 18;
        size_t curptr = 4;
        test.search_below(18, 0, indices, indptrs, curdex, curptr); // no-op.
        EXPECT_EQ(curdex, 18);
        EXPECT_EQ(curptr, 4);

        test.search_below(15, 0, indices, indptrs, curdex, curptr); // no-op, curdex is still the lower bound.
        EXPECT_EQ(curdex, 18);
        EXPECT_EQ(curptr, 4);

        test.search_below(12, 0, indices, indptrs, curdex, curptr); // decrease.
        EXPECT_EQ(curdex, 12);
        EXPECT_EQ(curptr, 3);

        test.search_below(11, 0, indices, indptrs, curdex, curptr); // another consecutive decrease.
        EXPECT_EQ(curdex, 11);
        EXPECT_EQ(curptr, 2);

        test.search_below(9, 0, indices, indptrs, curdex, curptr); // no-op, curdex is still the lower bound.
        EXPECT_EQ(curdex, 11);
        EXPECT_EQ(curptr, 2);

        test.search_below(1, 0, indices, indptrs, curdex, curptr); // big jump.
        EXPECT_EQ(curdex, 2);
        EXPECT_EQ(curptr, 0);
    }

    {
        int curdex = 19;
        size_t curptr = 11;
        test.search_below(18, 1, indices, indptrs, curdex, curptr); // no-op, still the lower bound..
        EXPECT_EQ(curdex, 19);
        EXPECT_EQ(curptr, 11);

        test.search_below(9, 1, indices, indptrs, curdex, curptr); // big jump (1).
        EXPECT_EQ(curdex, 10);
        EXPECT_EQ(curptr, 8);

        test.search_below(1, 1, indices, indptrs, curdex, curptr); // big jump (2).
        EXPECT_EQ(curdex, 3);
        EXPECT_EQ(curptr, 5);
    }

    {
        int curdex = 19;
        size_t curptr = 11;
        test.search_below(0, 1, indices, indptrs, curdex, curptr); // big jump to zero but no match.
        EXPECT_EQ(curdex, 3);
        EXPECT_EQ(curptr, 5);
    }

    {
        int curdex = 19;
        size_t curptr = 16;
        test.search_below(0, 2, indices, indptrs, curdex, curptr); // big jump to zero, with a match.
        EXPECT_EQ(curdex, 0);
        EXPECT_EQ(curptr, 11);
    }
}

TEST(SecondaryExtractionWorkspaceBase, Duplicated) {
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

    tatami::sparse::SecondaryExtractionWorkspaceBase<int, DuplicatedModifier> test(19);

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

    // Increment works correctly.
    {
        int curdex = 0;
        size_t curptr = 0;
        test.search_above(0, 0, indices, indptrs, curdex, curptr); // no-op.
        EXPECT_EQ(curdex, 0);
        EXPECT_EQ(curptr, 0);

        test.search_above(1, 0, indices, indptrs, curdex, curptr); // increment to next lower bound.
        EXPECT_EQ(curdex, 2);
        EXPECT_EQ(curptr, 3);

        test.search_above(2, 0, indices, indptrs, curdex, curptr); // no-op.
        EXPECT_EQ(curdex, 2);
        EXPECT_EQ(curptr, 3);

        test.search_above(7, 0, indices, indptrs, curdex, curptr); // big jump.
        EXPECT_EQ(curdex, 9);
        EXPECT_EQ(curptr, 6);

        test.search_above(15, 0, indices, indptrs, curdex, curptr); // another big jump.
        EXPECT_EQ(curdex, 15);
        EXPECT_EQ(curptr, 8);
    }

    {
        int curdex = 0;
        size_t curptr = 11;
        test.search_above(10, 1, indices, indptrs, curdex, curptr); // jumps work correctly.
        EXPECT_EQ(curdex, 11);
        EXPECT_EQ(curptr, 20);

        test.search_above(18, 1, indices, indptrs, curdex, curptr); // jump to end works correctly.
        EXPECT_EQ(curdex, 18);
        EXPECT_EQ(curptr, 22);
    }

    // Decrement works correctly.
    {
        int curdex = 19;
        size_t curptr = 34;
        test.search_below(16, 2, indices, indptrs, curdex, curptr); // no-op.
        EXPECT_EQ(curdex, 19);
        EXPECT_EQ(curptr, 34);

        test.search_below(15, 2, indices, indptrs, curdex, curptr); // decrement.
        EXPECT_EQ(curdex, 15);
        EXPECT_EQ(curptr, 33);

        test.search_below(8, 2, indices, indptrs, curdex, curptr); // big jump.
        EXPECT_EQ(curdex, 8);
        EXPECT_EQ(curptr, 29);

        test.search_below(7, 2, indices, indptrs, curdex, curptr); // decrement.
        EXPECT_EQ(curdex, 7);
        EXPECT_EQ(curptr, 27);

        test.search_below(0, 2, indices, indptrs, curdex, curptr); // jump to zero.
        EXPECT_EQ(curdex, 2);
        EXPECT_EQ(curptr, 25);
    }
}
