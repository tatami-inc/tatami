#include <gtest/gtest.h>
#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/utils/ConsecutiveOracle.hpp"
#include "tatami/utils/consecutive_extractor.hpp"

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
        auto ext = ref->dense_row();
        ext->set_oracle(std::make_unique<tatami::ConsecutiveOracle<int> >(20, 10)); // should be a no-op.
    }

    {
        auto ext = ref->dense_column();
        ext->set_oracle(std::make_unique<tatami::ConsecutiveOracle<int> >(20, 10)); // should be a no-op.
    }

    {
        auto ext = ref->sparse_row();
        ext->set_oracle(std::make_unique<tatami::ConsecutiveOracle<int> >(20, 10)); // should be a no-op.
    }

    {
        auto ext = ref->sparse_column();
        ext->set_oracle(std::make_unique<tatami::ConsecutiveOracle<int> >(20, 10)); // should be a no-op.
    }
}

TEST_F(BackCompatibilityTest, NewExtractor) {
    std::vector<double> buffer(ref->ncol());

    {
        auto ext = tatami::new_extractor<true, false>(ref.get());
        auto ptr = ext->fetch(0, buffer.data());
        EXPECT_EQ(ptr[0], -50);
        EXPECT_EQ(ptr[9], 130);
    }

    {
        auto ext = tatami::new_extractor<true, false>(ref.get(), tatami::Options());
        auto ptr = ext->fetch(0, buffer.data());
        EXPECT_EQ(ptr[0], -50);
        EXPECT_EQ(ptr[9], 130);
    }

    {
        auto ext = tatami::new_extractor<true, false>(ref.get(), 2, 5, tatami::Options());
        auto ptr = ext->fetch(0, buffer.data());
        EXPECT_EQ(ptr[0], -10);
        EXPECT_EQ(ptr[6], 70);
    }

    {
        auto ext = tatami::new_extractor<true, false>(ref.get(), { 1, 3, 5 }, tatami::Options());
        auto ptr = ext->fetch(0, buffer.data());
        EXPECT_EQ(ptr[0], -30);
        EXPECT_EQ(ptr[2], 50);
    }
}

TEST_F(BackCompatibilityTest, ConsecutiveExtractor) {
    auto ext = tatami::consecutive_extractor<true, false>(ref.get(), 5, 10);
    std::vector<double> buffer(ref->ncol());
    auto ptr = ext->fetch(5, buffer.data());
    EXPECT_EQ(ptr[0], -45);
    EXPECT_EQ(ptr[9], 135);
}
