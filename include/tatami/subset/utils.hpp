#ifndef TATAMI_DELAYED_SUBSET_UTILS_HPP
#define TATAMI_DELAYED_SUBSET_UTILS_HPP

#include "../base/Matrix.hpp"
#include "../utils/new_extractor.hpp"

#include <vector>
#include <algorithm>
#include <type_traits>

namespace tatami {

namespace subset_utils {

template<typename Index_, class IndexStorage_>
struct SubsetOracle : public Oracle<Index_> {
    SubsetOracle(std::shared_ptr<const Oracle<Index_> > ora, const IndexStorage_& ix) : source(std::move(ora)), indices(ix) {}

    Index_ get(size_t i) const {
        return indices[source->get(i)];
    }

    size_t total() const {
        return source->total();
    }

private:
    std::shared_ptr<const Oracle<Index_> > source;
    const IndexStorage_& indices;
};

template<typename Value_, typename Index_, class IndexStorage_>
struct MyopicPerpendicularDense : public MyopicDenseExtractor<Value_, Index_> {
    template<typename ... Args_>
    MyopicPerpendicularDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& in, bool row, Args_&& ... args) : 
        indices(in), internal(new_extractor<false, false>(mat, row, false, std::forward<Args_>(args)...)) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        return internal->fetch(indices[i], buffer);
    }

protected:
    const IndexStorage_& indices;
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > internal;
};

template<typename Value_, typename Index_, class IndexStorage_>
struct MyopicPerpendicularSparse : public MyopicSparseExtractor<Value_, Index_> {
    template<typename ... Args_>
    MyopicPerpendicularSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& in, bool row, Args_&& ... args) : 
        indices(in), internal(new_extractor<true, false>(mat, row, false, std::forward<Args_>(args)...)) {}

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        return internal->fetch(indices[i], vbuffer, ibuffer);
    }

protected:
    const IndexStorage_& indices;
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > internal;
};

template<typename Value_, typename Index_>
struct OracularPerpendicularDense : public OracularDenseExtractor<Value_, Index_> {
    template<class IndexStorage_, typename ... Args_>
    OracularPerpendicularDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& in, bool row, std::shared_ptr<const Oracle<Index_> > ora, Args_&& ... args) :
        internal(new_extractor<false, true>(mat, row, std::make_shared<SubsetOracle<Index_, IndexStorage_> >(std::move(ora), in), std::forward<Args_>(args)...)) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        return internal->fetch(i, buffer);
    }

protected:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > internal;
};

template<typename Value_, typename Index_>
struct OracularPerpendicularSparse : public OracularSparseExtractor<Value_, Index_> {
    template<class IndexStorage_, typename ... Args_>
    OracularPerpendicularSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& in, bool row, std::shared_ptr<const Oracle<Index_> > ora, Args_&& ... args) :
        internal(new_extractor<true, true>(mat, row, std::make_shared<SubsetOracle<Index_, IndexStorage_> >(std::move(ora), in), std::forward<Args_>(args)...)) {}

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        return internal->fetch(i, vbuffer, ibuffer);
    }

protected:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > internal;
};

}

}

#endif
