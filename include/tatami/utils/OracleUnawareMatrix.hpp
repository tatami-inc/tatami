#ifndef TATAMI_ORACLE_UNAWARE_MATRIX_HPP
#define TATAMI_ORACLE_UNAWARE_MATRIX_HPP

#include "../base/Matrix.hpp"
#include "../base/Extractor.hpp"

namespace tatami {

/**
 * @cond
 */
template<DimensionSelectionType selection_, bool sparse_, typename Value_, typename Index_> 
struct DummyOracleAwareExtractor : public OracleAwareExtractor<selection_, sparse_, Value_, Index_> {
    DummyOracleAwareExtractor(std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > ex, std::shared_ptr<Oracle<Index_> > o) : 
        internal_extractor(std::move(ex)), internal_oracle(std::move(o)) 
    {
        this->total_predictions = internal_oracle->total();
        if constexpr(selection_ == tatami::DimensionSelectionType::FULL) {
            this->full_length = internal_extractor->full_length;
        } else if constexpr(selection_ == tatami::DimensionSelectionType::BLOCK) {
            this->block_start = internal_extractor->block_start;
            this->block_length = internal_extractor->block_length;
        } else {
            this->index_length = internal_extractor->index_length;
        }
    }

    const Value_* fetch(Index_& i, Value_* buffer) {
        i = internal_oracle->get(this->used_predictions);
        ++(this->used_predictions);
        return internal_extractor->fetch(i, buffer);
    }

    SparseRange<Value_, Index_> fetch(Index_& i, Value_* vbuffer, Index_* ibuffer) {
        i = internal_oracle->get(this->used_predictions);
        ++(this->used_predictions);
        return internal_extractor->fetch(i, vbuffer, ibuffer);
    }

    const Index_* index_start() const {
        return internal_extractor->index_start();
    }

    const Oracle<Index_>* oracle() const {
        return internal_oracle.get();
    }

private:
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > internal_extractor;
    std::shared_ptr<Oracle<Index_> > internal_oracle;
};
/**
 * @endcond
 */

/**
 * @brief Virtual class for an oracle-unaware matrix.
 * 
 * @tparam Value Data value type, should be numeric.
 * @tparam Index Row/column index type, should be integer.
 *
 * This is intended for use as a base class of a concrete `Matrix` subclass that does not use an `Oracle` during extraction.
 * It sets defaults for all of the `Oracle`-aware methods so that we don't need to repeat this in every single subclass.
 */
template <typename Value_, typename Index_ = int>
class OracleUnawareMatrix : public Matrix<Value_, Index_> {
public:
    /**
     * @cond
     */
    OracleUnawareMatrix() = default;

    using Matrix<Value_, Index_>::dense_row;
    using Matrix<Value_, Index_>::dense_column;
    using Matrix<Value_, Index_>::sparse_row;
    using Matrix<Value_, Index_>::sparse_column;
    /**
     * @endcond
     */

public:
    std::unique_ptr<FullDenseOracleAwareExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return std::make_unique<DummyOracleAwareExtractor<DimensionSelectionType::FULL, false, Value_, Index_> >(this->dense_row(opt), std::move(oracle));
    }

    std::unique_ptr<BlockDenseOracleAwareExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<DummyOracleAwareExtractor<DimensionSelectionType::BLOCK, false, Value_, Index_> >(this->dense_row(block_start, block_length, opt), std::move(oracle));
    }

    std::unique_ptr<IndexDenseOracleAwareExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return std::make_unique<DummyOracleAwareExtractor<DimensionSelectionType::INDEX, false, Value_, Index_> >(this->dense_row(std::move(indices), opt), std::move(oracle));
    }

    std::unique_ptr<FullDenseOracleAwareExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return std::make_unique<DummyOracleAwareExtractor<DimensionSelectionType::FULL, false, Value_, Index_> >(this->dense_column(opt), std::move(oracle));
    }

    std::unique_ptr<BlockDenseOracleAwareExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<DummyOracleAwareExtractor<DimensionSelectionType::BLOCK, false, Value_, Index_> >(this->dense_column(block_start, block_length, opt), std::move(oracle));
    }

    std::unique_ptr<IndexDenseOracleAwareExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return std::make_unique<DummyOracleAwareExtractor<DimensionSelectionType::INDEX, false, Value_, Index_> >(this->dense_column(std::move(indices), opt), std::move(oracle));
    }

public:
    std::unique_ptr<FullSparseOracleAwareExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return std::make_unique<DummyOracleAwareExtractor<DimensionSelectionType::FULL, true, Value_, Index_> >(this->sparse_row(opt), std::move(oracle));
    }

    std::unique_ptr<BlockSparseOracleAwareExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<DummyOracleAwareExtractor<DimensionSelectionType::BLOCK, true, Value_, Index_> >(this->sparse_row(block_start, block_length, opt), std::move(oracle));
    }

    std::unique_ptr<IndexSparseOracleAwareExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return std::make_unique<DummyOracleAwareExtractor<DimensionSelectionType::INDEX, true, Value_, Index_> >(this->sparse_row(std::move(indices), opt), std::move(oracle));
    }

    std::unique_ptr<FullSparseOracleAwareExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return std::make_unique<DummyOracleAwareExtractor<DimensionSelectionType::FULL, true, Value_, Index_> >(this->sparse_column(opt), std::move(oracle));
    }

    std::unique_ptr<BlockSparseOracleAwareExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<DummyOracleAwareExtractor<DimensionSelectionType::BLOCK, true, Value_, Index_> >(this->sparse_column(block_start, block_length, opt), std::move(oracle));
    }

    std::unique_ptr<IndexSparseOracleAwareExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return std::make_unique<DummyOracleAwareExtractor<DimensionSelectionType::INDEX, true, Value_, Index_> >(this->sparse_column(std::move(indices), opt), std::move(oracle));
    }

};

}

#endif
