#ifndef TATAMI_TEST_FETCH_HPP
#define TATAMI_TEST_FETCH_HPP

#include "../tatami/base/Extractor.hpp"
#include "../tatami/base/SparseRange.hpp"
#include "../tatami/utils/copy.hpp"

#include <vector>
#include <algorithm>

namespace tatami_test {

namespace internal {

template<typename Value_, typename Index_>
void trim_sparse(const tatami::SparseRange<Value_, Index_>& raw, std::vector<Value_>& output_v, std::vector<Index_>& output_i) {
    tatami::copy_n(raw.value, raw.number, output_v.data());
    output_v.resize(raw.number);
    tatami::copy_n(raw.index, raw.number, output_i.data());
    output_i.resize(raw.number);
}

}

template<typename Value_, typename Index_>
std::vector<Value_> fetch(tatami::MyopicDenseExtractor<Value_, Index_>* ext, Index_ i, size_t number) {
    std::vector<Value_> output(number);
    auto raw = ext->fetch(i, output.data());
    tatami::copy_n(raw, output.size(), output.data());
    return output;
}

template<typename Value_, typename Index_>
std::vector<Value_> fetch(tatami::OracularDenseExtractor<Value_, Index_>* ext, size_t number) {
    std::vector<Value_> output(number);
    auto raw = ext->fetch(output.data());
    tatami::copy_n(raw, output.size(), output.data());
    return output;
}

template<typename Value_, typename Index_>
struct SparseVector {
    SparseVector(size_t n) : value(n), index(n) {}
    std::vector<Value_> value;
    std::vector<Index_> index;
};

template<typename Value_, typename Index_>
SparseVector<Value_, Index_> fetch(tatami::MyopicSparseExtractor<Value_, Index_>* ext, Index_ i, size_t number) {
    SparseVector<Value_, Index_> output(number);
    auto raw = ext->fetch(i, output.value.data(), output.index.data());
    internal::trim_sparse(raw, output.value, output.index);
    return output;
}

template<typename Value_, typename Index_>
SparseVector<Value_, Index_> fetch(tatami::OracularSparseExtractor<Value_, Index_>* ext, size_t number) {
    SparseVector<Value_, Index_> output(number);
    auto raw = ext->fetch(output.value.data(), output.index.data());
    internal::trim_sparse(raw, output.value, output.index);
    return output;
}

}

#endif
