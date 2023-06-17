#include <gtest/gtest.h>
#include "tatami/utils/SomeNumericArray.hpp"
#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/CompressedSparseMatrix.hpp"

#include "tatami_test/tatami_test.hpp"

#include <cstdint>
#include <numeric>

template<typename T, typename O>
void check_getters(typename tatami::SomeNumericArray<O>::Type t) {
    std::vector<T> stuff(100);
    std::iota(stuff.begin(), stuff.end(), 0);

    // Checking the basic methods.
    tatami::SomeNumericArray<O> vec(static_cast<void*>(stuff.data()), stuff.size(), t);
    EXPECT_EQ(vec.size(), stuff.size());
    for (size_t i = 0; i < stuff.size(); ++i) {
        EXPECT_EQ(stuff[i], vec[i]);
    }

    tatami::SomeNumericArray<O> uncast(stuff.data(), stuff.size()); // alternative constructor.
    EXPECT_EQ(uncast.size(), stuff.size());
    for (size_t i = 0; i < stuff.size(); ++i) {
        EXPECT_EQ(stuff[i], uncast[i]);
    }

    // Checking the iterator's methods.
    {
        auto it = vec.begin();
        EXPECT_EQ(*it, 0);
        EXPECT_EQ(it[1], 1);
    }

    {
        auto it = vec.begin();
        ++it;
        EXPECT_EQ(*it, 1);
        EXPECT_EQ(*(it++), 1);
        EXPECT_EQ(*it, 2);

        it += 3;
        EXPECT_EQ(*it, 5);
    }

    {
        auto it = vec.end();
        --it;
        EXPECT_EQ(*it, 99);
        EXPECT_EQ(*(it--), 99);
        EXPECT_EQ(*it, 98);

        it -= 3;
        EXPECT_EQ(*it, 95);
    }

    {
        auto it = vec.begin();
        EXPECT_EQ(*(it + 1), 1);
        EXPECT_EQ(*(2 + it), 2);

        it = vec.end();
        EXPECT_EQ(*(it - 5), 95);
    }

    {
        auto it = vec.begin();
        auto it2 = vec.end();
        EXPECT_TRUE(it == it);
        EXPECT_TRUE(it != it2);

        EXPECT_TRUE(it < it2);
        EXPECT_TRUE(it2 > it);

        EXPECT_TRUE(it2 >= it);
        EXPECT_TRUE(it <= it2);

        EXPECT_EQ(it2 - it, 100);
    }
}

TEST(SomeNumericArray, Integer) {
    check_getters<uint8_t, int>(tatami::SomeNumericArray<int>::U8);
    check_getters<int8_t, int>(tatami::SomeNumericArray<int>::I8);
    check_getters<uint16_t, int>(tatami::SomeNumericArray<int>::U16);
    check_getters<int16_t, int>(tatami::SomeNumericArray<int>::I16);
    check_getters<uint32_t, int>(tatami::SomeNumericArray<int>::U32);
    check_getters<int32_t, int>(tatami::SomeNumericArray<int>::I32);
    check_getters<uint64_t, int>(tatami::SomeNumericArray<int>::U64);
    check_getters<int64_t, int>(tatami::SomeNumericArray<int>::I64);
    check_getters<float, int>(tatami::SomeNumericArray<int>::F32);
    check_getters<double, int>(tatami::SomeNumericArray<int>::F64);
}

TEST(SomeNumericArray, Double) {
    check_getters<uint8_t, double>(tatami::SomeNumericArray<double>::U8);
    check_getters<int8_t, double>(tatami::SomeNumericArray<double>::I8);
    check_getters<uint16_t, double>(tatami::SomeNumericArray<double>::U16);
    check_getters<int16_t, double>(tatami::SomeNumericArray<double>::I16);
    check_getters<uint32_t, double>(tatami::SomeNumericArray<double>::U32);
    check_getters<int32_t, double>(tatami::SomeNumericArray<double>::I32);
    check_getters<uint64_t, double>(tatami::SomeNumericArray<double>::U64);
    check_getters<int64_t, double>(tatami::SomeNumericArray<double>::I64);
    check_getters<float, double>(tatami::SomeNumericArray<double>::F32);
    check_getters<double, double>(tatami::SomeNumericArray<double>::F64);
}

TEST(SomeNumericArray, DenseMatrix) {
    int nr = 20, nc = 10;
    std::vector<int32_t> values(nr * nc);
    std::iota(values.begin(), values.end(), -50);
    tatami::DenseColumnMatrix<double, int, decltype(values)> ref(nr, nc, values);

    tatami::SomeNumericArray<int> arr(values.data(), values.size());
    tatami::DenseColumnMatrix<double, int, decltype(arr)> alt(nr, nc, arr);

    tatami_test::test_simple_row_access(&alt, &ref, true, 1);
    tatami_test::test_simple_column_access(&alt, &ref, true, 1);
}

TEST(SomeNumericArray, SparseMatrix) {
    int nr = 10, nc = 6;
    std::vector<double>   values  { 2, 5, 3, 4, 5, 5, 7, 3, 2, 4, 5, 1 };
    std::vector<uint8_t>  indices { 0, 1, 7, 1, 4, 6, 2, 5, 5, 1, 8, 9 };
    std::vector<uint32_t> indptrs { 0, 3, 6, 8, 9, 11, 12 };
    tatami::CompressedSparseColumnMatrix<double, int, decltype(values), decltype(indices), decltype(indptrs)> ref(nr, nc, values, indices, indptrs);

    tatami::SomeNumericArray<double> varr  (values.data(), values.size());
    tatami::SomeNumericArray<int>    iarr  (indices.data(), indices.size());
    tatami::SomeNumericArray<size_t> indarr(indptrs.data(), indptrs.size());
    tatami::CompressedSparseColumnMatrix<double, int, decltype(varr), decltype(iarr), decltype(indarr)> alt(nr, nc, varr, iarr, indarr);

    tatami_test::test_simple_row_access(&alt, &ref, true, 1);
    tatami_test::test_simple_column_access(&alt, &ref, true, 1);
}
