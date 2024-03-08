#ifndef TATAMI_PSEUDO_ORACULAR_EXTRACTOR_HPP
#define TATAMI_PSEUDO_ORACULAR_EXTRACTOR_HPP

#include "../base/Matrix.hpp"
#include "../base/Extractor.hpp"

namespace tatami {

template<typename Value_, typename Index_> 
struct PseudoOracularDenseExtractor : public OracularDenseExtractor<Value_, Index_> {
    PseudoOracularDenseExtractor(std::shared_ptr<Oracle<Index_> > ora, std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > r) :
        oracle(std::move(ora)), raw(std::move(r)) {}

    const Value_* fetch(Index_& i, Value_* buffer) {
        i = oracle->get(used);
        ++used;
        return raw->fetch(i, buffer);
    }

    Index_ number() const {
        return raw->number();
    }

private:
    std::shared_ptr<Oracle<Index_> > oracle;
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > raw;
    size_t used = 0;
};

template<typename Value_, typename Index_> 
struct PseudoOracularSparseExtractor : public OracularSparseExtractor<Value_, Index_> {
    PseudoOracularSparseExtractor(std::shared_ptr<Oracle<Index_> > ora, std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > r) :
        oracle(std::move(ora)), raw(std::move(r)) {}

    SparseRange<Value_, Index_> fetch(Index_& i, Value_* vbuffer, Index_* ibuffer) {
        i = oracle->get(used);
        ++used;
        return raw->fetch(i, vbuffer, ibuffer);
    }

    Index_ number() const {
        return raw->number();
    }

private:
    std::shared_ptr<Oracle<Index_> > oracle;
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > raw;
    size_t used = 0;
};

}

#endif
