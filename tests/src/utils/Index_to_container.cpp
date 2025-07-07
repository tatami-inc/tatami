#include <gtest/gtest.h>
#include "tatami/utils/Index_to_container.hpp"
#include <cstdint>

TEST(IndexToContainer, SafeNonNegativeEqual) {
    EXPECT_TRUE(tatami::safe_non_negative_equal(0, 0));
    EXPECT_FALSE(tatami::safe_non_negative_equal(-1, -1));
    EXPECT_FALSE(tatami::safe_non_negative_equal(0, -1));
    EXPECT_FALSE(tatami::safe_non_negative_equal(0, 1));
    EXPECT_TRUE(tatami::safe_non_negative_equal(1, 1));

    // Trying with a mix of signed and unsigned types.
    EXPECT_TRUE(tatami::safe_non_negative_equal(1u, 1));
    EXPECT_FALSE(tatami::safe_non_negative_equal(1u, -1));
    EXPECT_TRUE(tatami::safe_non_negative_equal(1, 1u));
    EXPECT_FALSE(tatami::safe_non_negative_equal(-1, 1u));
    EXPECT_TRUE(tatami::safe_non_negative_equal(1u, 1u));
}

template<typename Size_>
class MockVector {
public:
    MockVector(Size_ size) : my_size(size) {}
    Size_ size() const { return my_size; }
    void resize(Size_ size) { my_size = size; }
private:
    Size_ my_size;
};

TEST(IndexToContainer, Cast) {
    EXPECT_EQ(tatami::can_cast_Index_to_container_size<MockVector<std::uint64_t> >(100), 100); 
    EXPECT_EQ(tatami::can_cast_Index_to_container_size<MockVector<std::uint8_t> >(100), 100); 
    EXPECT_EQ(tatami::cast_Index_to_container_size<MockVector<std::uint64_t> >(100), 100); 
    EXPECT_EQ(tatami::cast_Index_to_container_size<MockVector<std::uint8_t> >(100), 100); 

    std::string msg;
    try {
        tatami::cast_Index_to_container_size<MockVector<std::uint8_t> >(10000);
    } catch (std::exception& e) {
        msg = e.what();
    }
    EXPECT_TRUE(msg.find("cast") != std::string::npos);
}

TEST(IndexToContainer, Create) {
    auto x = tatami::create_container_of_Index_size<MockVector<std::uint64_t> >(100);
    EXPECT_EQ(x.size(), 100);
    tatami::resize_container_to_Index_size(x, 20);
    EXPECT_EQ(x.size(), 20);

    auto y = tatami::create_container_of_Index_size<MockVector<std::uint8_t> >(100);
    EXPECT_EQ(y.size(), 100);
    tatami::resize_container_to_Index_size(y, 20);
    EXPECT_EQ(y.size(), 20);

    std::string msg;
    try {
        tatami::create_container_of_Index_size<MockVector<std::uint8_t> >(10000);
    } catch (std::exception& e) {
        msg = e.what();
    }
    EXPECT_TRUE(msg.find("cast") != std::string::npos);

    msg.clear();
    try {
        MockVector<std::uint8_t> z(0);
        tatami::resize_container_to_Index_size(z, 10000);
    } catch (std::exception& e) {
        msg = e.what();
    }
    EXPECT_TRUE(msg.find("cast") != std::string::npos);
}
