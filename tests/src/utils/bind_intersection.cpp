#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/utils/bind_intersection.hpp"
#include "tatami/dense/DenseMatrix.hpp"

#include "tatami_test/tatami_test.hpp"

class BindIntersectionTest : public ::testing::Test {
protected:   
    std::vector<std::shared_ptr<tatami::NumericMatrix> > collected;

    std::vector<const int*> get_id_ptrs(const std::vector<std::vector<int> >& ids) { 
        std::vector<const int*> id_ptr;
        for (const auto& x : ids) {
            id_ptr.push_back(x.data());
        }
        return id_ptr;
    }

    size_t find_chosen(int chosen_id, const std::vector<int>& ids) {
        size_t actual_chosen = -1;
        for (size_t j = 0; j < ids.size(); ++j) {
            if (chosen_id == ids[j]) {
                actual_chosen = j;
                break;
            }
        }
        if (actual_chosen == static_cast<size_t>(-1)) {
            throw std::runtime_error("could not find the chosen ID");
        }
        return actual_chosen;
    }
};

TEST_F(BindIntersectionTest, NoOp) {
    // Simplest case: everything is the same.
    size_t NR = 20;
    size_t NC = 0;

    for (int i = 0; i < 3; ++i) {
        size_t cur_NC = 5 * (i + 1);
        auto vec = tatami_test::simulate_dense_vector<double>(NR * cur_NC, i /** seed **/);
        collected.emplace_back(new tatami::DenseMatrix<true, double, int>(NR, cur_NC, std::move(vec)));
        NC += cur_NC;
    }

    std::vector<std::vector<int> > ids(collected.size());
    for (auto& x : ids) {
        x.resize(NR);
        std::iota(x.begin(), x.end(), 0);
    }
    auto id_ptr = get_id_ptrs(ids);

    // Combining by column.
    auto output = tatami::bind_intersection<1>(collected, id_ptr);
    EXPECT_EQ(output.first->nrow(), NR);
    EXPECT_EQ(output.first->ncol(), NC);
    EXPECT_EQ(output.second.size(), NR);

    std::vector<size_t> expected(NR);
    std::iota(expected.begin(), expected.end(), 0);
    EXPECT_EQ(output.second, expected);

    // Checking that the matrix contains the expected values.
    size_t chosen = 5;
    auto wrk = output.first->dense_row();
    auto vals = wrk->fetch(chosen);

    size_t offset = 0;
    for (int i = 0; i < 3; ++i) {
        auto iwrk = collected[i]->dense_row();
        auto expected = iwrk->fetch(chosen);
        auto sofar = offset;
        offset += collected[i]->ncol();
        std::vector<double> observed(vals.begin() + sofar, vals.begin() + offset);
        EXPECT_EQ(observed, expected);
    }
}

TEST_F(BindIntersectionTest, Shuffled) {
    std::vector<std::shared_ptr<tatami::NumericMatrix> > collected;
    size_t NR = 5;
    size_t NC = 0;

    for (int i = 0; i < 3; ++i) {
        size_t cur_NC = 5 * (i + 1);
        auto vec = tatami_test::simulate_dense_vector<double>(NR * cur_NC, i /** seed **/);
        collected.emplace_back(new tatami::DenseMatrix<true, double, int>(NR, cur_NC, std::move(vec)));
        NC += cur_NC;
    }

    // Shuffling the identifiers.
    std::vector<std::vector<int> > ids(collected.size());
    ids[0] = std::vector<int>{ 4, 3, 2, 1, 0 };
    ids[1] = std::vector<int>{ 0, 1, 2, 3, 4 };
    ids[2] = std::vector<int>{ 0, 2, 4, 3, 1 };

    auto id_ptr = get_id_ptrs(ids);

    // Combining by column.
    auto output = tatami::bind_intersection<1>(collected, id_ptr);
    EXPECT_EQ(output.first->nrow(), NR);
    EXPECT_EQ(output.first->ncol(), NC);
    EXPECT_EQ(output.second.size(), NR);

    std::vector<size_t> expected{ 4, 3, 2, 1, 0 }; // flipped.
    EXPECT_EQ(output.second, expected);

    // Checking that the matrix contains the expected values.
    size_t chosen = 1;
    auto wrk = output.first->dense_row();
    auto vals = wrk->fetch(chosen);
    size_t chosen_id = ids[0][output.second[chosen]];

    size_t offset = 0;
    for (int i = 0; i < 3; ++i) {
        auto actual_chosen = find_chosen(chosen_id, ids[i]);
        auto iwrk = collected[i]->dense_row();
        auto expected = iwrk->fetch(actual_chosen);
        auto sofar = offset;
        offset += collected[i]->ncol();
        std::vector<double> observed(vals.begin() + sofar, vals.begin() + offset);
        EXPECT_EQ(observed, expected);
    }
}

TEST_F(BindIntersectionTest, Uncommon) {
    std::vector<std::shared_ptr<tatami::NumericMatrix> > collected;
    size_t NC = 0;

    for (int i = 0; i < 3; ++i) {
        size_t NR = 5 + i;
        size_t cur_NC = 5 * (i + 1);
        auto vec = tatami_test::simulate_dense_vector<double>(NR * cur_NC, i /** seed **/);
        collected.emplace_back(new tatami::DenseMatrix<true, double, int>(NR, cur_NC, std::move(vec)));
        NC += cur_NC;
    }

    // Shuffling the identifiers.
    std::vector<std::vector<int> > ids(collected.size());
    ids[0] = std::vector<int>{ 0, 2, 4, 6, 8 };
    ids[1] = std::vector<int>{ 6, 5, 4, 3, 2, 1 }; // mixing it up a little
    ids[2] = std::vector<int>{ 3, 5, 7, 2, 4, 6, 8 }; // mixing it up a little

    auto id_ptr = get_id_ptrs(ids);

    // Combining by column.
    auto output = tatami::bind_intersection<1>(collected, id_ptr);
    EXPECT_EQ(output.first->nrow(), 3);
    EXPECT_EQ(output.first->ncol(), NC);
    EXPECT_EQ(output.second.size(), 3);

    std::vector<int> true_ids;
    for (auto i : output.second) {
        true_ids.push_back(ids[0][i]);
    }
    std::vector<int> expected { 2, 4, 6 };
    EXPECT_EQ(true_ids, expected);

    // Checking that the matrix contains the expected values.
    size_t chosen = 0;
    auto wrk = output.first->dense_row();
    auto vals = wrk->fetch(chosen);
    size_t chosen_id = ids[0][output.second[chosen]];

    size_t offset = 0;
    for (int i = 0; i < 3; ++i) {
        auto actual_chosen = find_chosen(chosen_id, ids[i]);
        auto iwrk = collected[i]->dense_row();
        auto expected = iwrk->fetch(actual_chosen);
        auto sofar = offset;
        offset += collected[i]->ncol();
        std::vector<double> observed(vals.begin() + sofar, vals.begin() + offset);
        EXPECT_EQ(observed, expected);
    }
}

TEST_F(BindIntersectionTest, NoCommon) {
    std::vector<std::shared_ptr<tatami::NumericMatrix> > collected;
    size_t NC = 0;
    size_t NR = 5;

    for (int i = 0; i < 2; ++i) {
        size_t cur_NC = 5 * (i + 1);
        auto vec = tatami_test::simulate_dense_vector<double>(NR * cur_NC, i /** seed **/);
        collected.emplace_back(new tatami::DenseMatrix<true, double, int>(NR, cur_NC, std::move(vec)));
        NC += cur_NC;
    }

    // Shuffling the identifiers.
    std::vector<std::vector<int> > ids(collected.size());
    ids[0] = std::vector<int>{ 0, 2, 4, 6, 8 };
    ids[1] = std::vector<int>{ 1, 3, 5, 7, 9 };

    auto id_ptr = get_id_ptrs(ids);

    // Combining by column.
    auto output = tatami::bind_intersection<1>(collected, id_ptr);
    EXPECT_EQ(output.first->nrow(), 0);
    EXPECT_EQ(output.first->ncol(), NC);
    EXPECT_EQ(output.second.size(), 0);
}

TEST_F(BindIntersectionTest, ByRows) {
    std::vector<std::shared_ptr<tatami::NumericMatrix> > collected;
    size_t NR = 0;

    for (int i = 0; i < 3; ++i) {
        size_t NC = 10 - i * 2;
        size_t cur_NR = (i + 1) * 10;
        auto vec = tatami_test::simulate_dense_vector<double>(cur_NR * NC, i /** seed **/);
        collected.emplace_back(new tatami::DenseMatrix<true, double, int>(cur_NR, NC, std::move(vec)));
        NR += cur_NR;
    }

    // Shuffling the identifiers.
    std::vector<std::vector<int> > ids(collected.size());
    ids[0] = std::vector<int>{ 0, 2, 4, 6, 8, 1, 3, 5, 7, 9 };
    ids[1] = std::vector<int>{ 0, 1, 2, 3, 4, 5, 6, 7 };
    ids[2] = std::vector<int>{ 6, 5, 4, 3, 2, 1 };

    auto id_ptr = get_id_ptrs(ids);

    // Combining by column.
    auto output = tatami::bind_intersection<0>(collected, id_ptr);
    EXPECT_EQ(output.first->ncol(), 6);
    EXPECT_EQ(output.first->nrow(), NR);
    EXPECT_EQ(output.second.size(), 6);

    std::vector<int> true_ids;
    for (auto i : output.second) {
        true_ids.push_back(ids[0][i]);
    }
    std::vector<int> expected { 1, 2, 3, 4, 5, 6 };
    EXPECT_EQ(true_ids, expected);

    // Checking that the matrix contains the expected values.
    size_t chosen = 4;
    auto wrk = output.first->dense_column();
    auto vals = wrk->fetch(chosen);
    size_t chosen_id = ids[0][output.second[chosen]];

    size_t offset = 0;
    for (int i = 0; i < 3; ++i) {
        auto actual_chosen = find_chosen(chosen_id, ids[i]);
        auto iwrk = collected[i]->dense_column();
        auto expected = iwrk->fetch(actual_chosen);
        auto sofar = offset;
        offset += collected[i]->nrow();
        std::vector<double> observed(vals.begin() + sofar, vals.begin() + offset);
        EXPECT_EQ(observed, expected);
    }
}


