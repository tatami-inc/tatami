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
struct PseudoOracularDenseExtractor final : public OracularDenseExtractor<Value_, Index_> {
    /**
     * @param oracle The oracle.
     * @param ext A dense extractor to extract dimension elements according to `ora`.
     */
    PseudoOracularDenseExtractor(std::shared_ptr<const Oracle<Index_> > oracle, std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > ext) :
        my_oracle(std::move(oracle)), my_ext(std::move(ext)) {}

    const Value_* fetch(const Index_, Value_* const buffer) {
        auto i = my_oracle->get(my_used++);
        return my_ext->fetch(i, buffer);
    }

private:
    std::shared_ptr<const Oracle<Index_> > my_oracle;
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > my_ext;
    PredictionIndex my_used = 0;
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
struct PseudoOracularSparseExtractor final : public OracularSparseExtractor<Value_, Index_> {
    /**
     * @param oracle The oracle.
     * @param ext A sparse extractor to extract dimension elements according to `ora`.
     */
    PseudoOracularSparseExtractor(std::shared_ptr<const Oracle<Index_> > oracle, std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > ext) :
        my_oracle(std::move(oracle)), my_ext(std::move(ext)) {}

    SparseRange<Value_, Index_> fetch(const Index_, Value_* const value_buffer, Index_* const index_buffer) {
        auto i = my_oracle->get(my_used++);
        return my_ext->fetch(i, value_buffer, index_buffer);
    }

private:
    std::shared_ptr<const Oracle<Index_> > my_oracle;
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > my_ext;
    PredictionIndex my_used = 0;
};

}

#endif
