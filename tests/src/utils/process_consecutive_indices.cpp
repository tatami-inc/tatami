#include <gtest/gtest.h>
#include "tatami/tatami.hpp"

TEST(ProcessConsecutiveIndices, Basic) {
    {
        std::vector<std::pair<int, int> > stored;

        std::vector<int> values{ 1, 2, 3, 6, 7, 10, 12, 13, 14, 15, 20, 21, 30 };
        tatami::process_consecutive_indices<int>(values.data(), values.size(), [&](int x, int l) -> void { stored.emplace_back(x, l); });

        std::vector<std::pair<int, int> > expected { { 1, 3 }, { 6, 2 }, { 10, 1 }, { 12, 4 }, { 20, 2 }, { 30, 1 }};
        EXPECT_EQ(stored, expected);
    }

    {
        std::vector<std::pair<int, int> > stored;

        std::vector<int> values{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        tatami::process_consecutive_indices<int>(values.data(), values.size(), [&](int x, int l) -> void { stored.emplace_back(x, l); });

        std::vector<std::pair<int, int> > expected{ { 0, 11 } };
        EXPECT_EQ(stored, expected);
    }

    {
        std::vector<std::pair<int, int> > stored;

        std::vector<int> values{ 12, 14, 16, 18, 20 };
        tatami::process_consecutive_indices<int>(values.data(), values.size(), [&](int x, int l) -> void { stored.emplace_back(x, l); });

        std::vector<std::pair<int, int> > expected{ { 12, 1 }, { 14, 1 }, { 16, 1 }, { 18, 1 }, { 20, 1 } };
        EXPECT_EQ(stored, expected);
    }
}
