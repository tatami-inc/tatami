#ifndef TATAMI_DELAYED_SUBSET_UTILS_HPP
#define TATAMI_DELAYED_SUBSET_UTILS_HPP

#include "../base/Matrix.hpp"
#include "../base/new_extractor.hpp"

#include <vector>
#include <algorithm>
#include <type_traits>

namespace tatami {

namespace subset_utils {

template<typename Value_, typename Index_>
const Value_* remap_dense(const Value_* input, Value_* buffer, const std::vector<Index_>& rmapping) {
    auto temp = buffer;
    for (auto i : rmapping) {
        *temp = input[i];
        ++temp;
    } 
    return buffer;
}

template<typename Index_, class IndexStorage_>
struct SubsetOracle : public Oracle<Index_> {
    SubsetOracle(std::shared_ptr<Oracle<Index_> > ora, const IndexStorage_& ix) : source(std::move(ora)), indices(ix) {}

    Index_ get(size_t i) const {
        return indices[source->get(i)];
    }

    size_t total() const {
        return source->total();
    }

private:
    std::shared_ptr<Oracle<Index_> > source;
    const IndexStorage_& indices;
};

template<typename Value_, typename Index_, class IndexStorage_>
struct MyopicPerpendicularDense : public MyopicDenseExtractor<Value_, Index_> {
    template<bool row_, typename ... Args_>
    MyopicPerpendicularDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& in, std::integral_constant<bool, row_>, Args_&& ... args) : 
        indices(in), internal(new_extractor<row_, false>(mat, std::forward<Args_>(args)...)) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        return internal->fetch(indices[i], buffer);
    }

protected:
    const IndexStorage_& indices;
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > internal;
};

template<typename Value_, typename Index_, class IndexStorage_>
struct MyopicPerpendicularSparse : public MyopicSparseExtractor<Value_, Index_> {
    template<bool row_, typename ... Args_>
    MyopicPerpendicularSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& in, std::integral_constant<bool, row_>, Args_&& ... args) : 
        indices(in), internal(new_extractor<row_, true>(mat, std::forward<Args_>(args)...)) {}

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        return internal->fetch(indices[i], vbuffer, ibuffer);
    }

protected:
    const IndexStorage_& indices;
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > internal;
};

template<typename Value_, typename Index_, class IndexStorage_>
struct OracularPerpendicularDense : public OracularDenseExtractor<Value_, Index_> {
    template<bool row_, typename ... Args_>
    OracularPerpendicularDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& in, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > ora, Args_&& ... args) :
        internal(new_extractor<row_, false>(mat, std::make_shared<SubsetOracle<Index_, IndexStorage_> >(std::move(ora), in), std::forward<Args_>(args)...)) {}

    const Value_* fetch(Value_* buffer) {
        return internal->fetch(buffer);
    }

protected:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > internal;
};

template<typename Value_, typename Index_, class IndexStorage_>
struct OracularPerpendicularSparse : public OracularSparseExtractor<Value_, Index_> {
    template<bool row_, typename ... Args_>
    OracularPerpendicularSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& in, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > ora, Args_&& ... args) :
        internal(new_extractor<row_, true>(mat, std::make_shared<SubsetOracle<Index_, IndexStorage_> >(std::move(ora), in), std::forward<Args_>(args)...)) {}

    SparseRange<Value_, Index_> fetch(Value_* vbuffer, Index_* ibuffer) {
        return internal->fetch(vbuffer, ibuffer);
    }

protected:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > internal;
};

}

}

#endif
