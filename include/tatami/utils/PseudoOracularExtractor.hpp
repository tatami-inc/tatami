#ifndef TATAMI_PSEUDO_ORACULAR_EXTRACTOR_HPP
#define TATAMI_PSEUDO_ORACULAR_EXTRACTOR_HPP

#include "../base/Matrix.hpp"
#include "../base/Extractor.hpp"

/**
 * @file PseudoOracularExtractor.hpp
 * @brief Mimic the oracle-aware extractor interface.
 */

namespace tatami {

/**
 * @brief Mimic the `OracularDenseExtractor` interface.
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * This is used to quickly implement the dense extraction methods for a `Matrix` subclass that does not benefit from an `Oracle`.
 * Specifically, the oracle is used to generate a prediction that is passed to an oracle-unaware `MyopicDenseExtractor`.
 * This allows `Matrix` subclasses to satisfy the dense oracle-aware extraction interface.
 */
template<typename Value_, typename Index_> 
struct PseudoOracularDenseExtractor : public OracularDenseExtractor<Value_, Index_> {
    /**
     * @param ora The oracle.
     * @param r A dense extractor to extract dimension elements according to `ora`.
     */
    PseudoOracularDenseExtractor(std::shared_ptr<const Oracle<Index_> > ora, std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > r) :
        oracle(std::move(ora)), raw(std::move(r)) {}

    const Value_* fetch(Index_, Value_* buffer) {
        auto i = oracle->get(used);
        ++used;
        return raw->fetch(i, buffer);
    }

private:
    std::shared_ptr<const Oracle<Index_> > oracle;
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > raw;
    size_t used = 0;
};

/**
 * @brief Mimic the `OracularSparseExtractor` interface.
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * This is used to quickly implement the sparse extraction methods for a `Matrix` subclass that does not benefit from an `Oracle`.
 * Specifically, the oracle is used to generate a prediction that is passed to an oracle-unaware `MyopicSparseExtractor`.
 * This allows `Matrix` subclasses to satisfy the sparse oracle-aware extraction interface without much effort.
 */
template<typename Value_, typename Index_> 
struct PseudoOracularSparseExtractor : public OracularSparseExtractor<Value_, Index_> {
    /**
     * @param ora The oracle.
     * @param r A sparse extractor to extract dimension elements according to `ora`.
     */
    PseudoOracularSparseExtractor(std::shared_ptr<const Oracle<Index_> > ora, std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > r) :
        oracle(std::move(ora)), raw(std::move(r)) {}

    SparseRange<Value_, Index_> fetch(Index_, Value_* vbuffer, Index_* ibuffer) {
        auto i = oracle->get(used);
        ++used;
        return raw->fetch(i, vbuffer, ibuffer);
    }

private:
    std::shared_ptr<const Oracle<Index_> > oracle;
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > raw;
    size_t used = 0;
};

}

#endif
