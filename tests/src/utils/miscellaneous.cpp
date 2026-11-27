#include <gtest/gtest.h>
#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/utils/ConsecutiveOracle.hpp"
#include "tatami/utils/deprecated/consecutive_extractor.hpp"

class BackCompatibilityTest : public ::testing::Test {
protected:
    std::shared_ptr<tatami::NumericMatrix> ref;

    void SetUp() {
        int nr = 20, nc = 10;
        std::vector<int> values(nr * nc);
        std::iota(values.begin(), values.end(), -50);
        ref.reset(new tatami::DenseColumnMatrix<double, int, decltype(values)>(nr, nc, values));
    }
};

TEST_F(BackCompatibilityTest, SetOracle) {
    {
        auto ext = ref->dense(true, tatami::Options());
        ext->set_oracle(std::make_unique<tatami::ConsecutiveOracle<int> >(20, 10)); // should be a no-op.
    }

    {
        auto ext = ref->dense(false, tatami::Options());
        ext->set_oracle(std::make_unique<tatami::ConsecutiveOracle<int> >(20, 10)); // should be a no-op.
    }

    {
        auto ext = ref->sparse(true, tatami::Options());
        ext->set_oracle(std::make_unique<tatami::ConsecutiveOracle<int> >(20, 10)); // should be a no-op.
    }

    {
        auto ext = ref->sparse(false, tatami::Options());
        ext->set_oracle(std::make_unique<tatami::ConsecutiveOracle<int> >(20, 10)); // should be a no-op.
    }
}

TEST(HasData, Basic) {
    {
        constexpr bool found = tatami::has_data<int, std::vector<int> >::value;
        EXPECT_TRUE(found);
    }

    {
        constexpr bool found = tatami::has_data<int, const std::vector<int> >::value;
        EXPECT_TRUE(found);
    }

    {
        constexpr bool found = tatami::has_data<double, std::vector<int> >::value;
        EXPECT_FALSE(found);
    }

    {
        constexpr bool found = tatami::has_data<double, std::pair<double, double> >::value;
        EXPECT_FALSE(found);
    }
}
