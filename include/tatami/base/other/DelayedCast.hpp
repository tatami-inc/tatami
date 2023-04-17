#ifndef TATAMI_DELAYED_CAST_HPP
#define TATAMI_DELAYED_CAST_HPP

#include "../Matrix.hpp"
#include "../utils.hpp"
#include <memory>
#include <type_traits>

/**
 * @file DelayedCast.hpp
 *
 * @brief Delayed cast to another interface type.
 */

namespace tatami {

/**
 * @brief Recast a `Matrix` to a different interface type.
 *
 * This performs a delayed cast from one interface type to another.
 * It is useful as a compatibility layer between functions that require `Matrix` objects of different types.
 * Casting is achieved by extracting the requested row/column from the input `Matrix` and transforming it to the output types.
 * Note that this is only done per row/column - the entirety of the original matrix is not copied to the new type.
 *
 * @tparam Value_out_ Data type to cast to.
 * @tparam Index_out_ Index type to cast to.
 * @tparam Value_in_ Data type to cast from.
 * @tparam Index_in_ Index type to cast from.
 */
template<typename Value_out_, typename Index_out_, typename Value_in_, typename Index_in_>
class DelayedCast : public Matrix<Value_out_, Index_out_> {
public:
    /**
     * @param p Pointer to the `Matrix` instance to cast from.
     */
    DelayedCast(std::shared_ptr<Matrix<Value_in_, Index_in_> > p) : ptr(std::move(p)) {}

public:
    Index_out_ nrow() const {
        return ptr->nrow();
    }

    Index_out_ ncol() const {
        return ptr->ncol();
    }

    bool sparse() const {
        return ptr->sparse();
    }

    bool prefer_rows() const { 
        return ptr->prefer_rows();
    }

    std::pair<double, double> dimension_preference () const {
        return ptr->dimension_preference();
    }

private:
    std::shared_ptr<Matrix<Value_in_, Index_in_> > ptr;
    static constexpr bool same_Value_type_ = std::is_same<Value_in_, Value_out_>::value;
    static constexpr bool same_Index_type_ = std::is_same<Index_in_, Index_out_>::value;

    template<bool sparse_, typename IndexBuffer_>
    struct CastExtractor : public Extractor<sparse_, Value_out_, Index_out_> {
        CastExtractor(std::unique_ptr<Extractor<sparse_, Value_in_, Index_in_> > i, IndexBuffer_ original_indices) : 
            internal(std::move(i)), 
            indices_buffer(std::move(original_indices)) 
        {
            this->extracted_selection = internal->extracted_selection;
            this->extracted_length = internal->extracted_length;
            this->extracted_block = internal->extracted_block;
        }

        CastExtractor(std::unique_ptr<Extractor<sparse_, Value_in_, Index_in_> > i) : internal(std::move(i)) {
            this->extracted_selection = internal->extracted_selection;
            this->extracted_length = internal->extracted_length;
            this->extracted_block = internal->extracted_block;
        }

        const Index_out_* extracted_index() const {
            if constexpr(same_Index_type_) {
                return internal->extracted_index();
            } else if constexpr(std::is_same<IndexBuffer_, std::vector<Index_out_> >::value) {
                return indices_buffer.data();
            } else {
                return NULL;
            }
        }
    protected:
        std::unique_ptr<Extractor<sparse_, Value_in_, Index_in_> > internal;
        IndexBuffer_ indices_buffer;
    };

    template<typename IndexBuffer_>
    struct DenseCastExtractor : public CastExtractor<false, IndexBuffer_> {
        DenseCastExtractor(std::unique_ptr<DenseExtractor<Value_in_, Index_in_> > i, IndexBuffer_ original_indices) : 
            CastExtractor<false, IndexBuffer_>(std::move(i), std::move(original_indices))
        {
            if constexpr(!same_Value_type_) {
                internal_buffer.resize(this->extracted_length);
            }
        }

        const Value_out_* fetch(Index_out_ i, Value_out_* buffer) {
            if constexpr(same_Value_type_) {
                return this->internal->fetch(i, buffer);
            } else {
                auto out = this->internal->fetch(i, internal_buffer.data());
                std::copy(out, out + this->extracted_length, buffer);
                return buffer;
            }
        }

    protected:
        typename std::conditional<same_Value_type_, bool, std::vector<Value_in_> >::type internal_buffer;
    };

    template<typename IndexBuffer_>
    struct SparseCastExtractor : public CastExtractor<true, IndexBuffer_> {
        SparseCastExtractor(std::unique_ptr<SparseExtractor<Value_in_, Index_in_> > i, IndexBuffer_ original_indices, bool needs_value, bool needs_index) : 
            CastExtractor<true, IndexBuffer_>(std::move(i), std::move(original_indices)), 
            internal_vbuffer(needs_value && !same_Value_type_ ? this->extracted_length : 0),
            internal_ibuffer(needs_index && !same_Index_type_ ? this->extracted_length : 0)
        {}

        SparseRange<Value_out_, Index_out_> fetch(Index_out_ i, Value_out_* vbuffer, Index_out_* ibuffer) {
            if constexpr(same_Value_type_) {
                if constexpr(same_Index_type_) {
                    return this->internal->fetch(i, vbuffer, ibuffer);

                } else {
                    auto out = this->internal->fetch(i, vbuffer, internal_ibuffer.data());
                    if (out.index) {
                        std::copy(out.index, out.index + out.number, ibuffer);
                    } else {
                        ibuffer = NULL;
                    }
                    return SparseRange<Value_out_, Index_out_>(out.number, out.value, ibuffer);
                }

            } else {
                if constexpr(same_Index_type_) {
                    auto out = this->internal->fetch(i, internal_vbuffer.data(), ibuffer);
                    if (out.value) {
                        std::copy(out.value, out.value + out.number, vbuffer);
                    } else {
                        vbuffer = NULL;
                    }
                    return SparseRange<Value_out_, Index_out_>(out.number, vbuffer, out.index);

                } else {
                    auto out = this->internal->fetch(i, internal_vbuffer.data(), internal_ibuffer.data());
                    if (out.value) {
                        std::copy(out.value, out.value + out.number, vbuffer);
                    } else {
                        vbuffer = NULL;
                    }
                    if (out.index) {
                        std::copy(out.index, out.index + out.number, ibuffer);
                    } else {
                        ibuffer = NULL;
                    }
                    return SparseRange<Value_out_, Index_out_>(out.number, vbuffer, ibuffer);
                }
            }
        }

    protected:
        typename std::conditional<same_Value_type_, bool, std::vector<Value_in_> >::type internal_vbuffer;
        typename std::conditional<same_Index_type_, bool, std::vector<Index_in_> >::type internal_ibuffer;
    };

private:
    template<bool accrow_, bool sparse_, typename IndexBuffer_>
    std::unique_ptr<Extractor<sparse_, Value_out_, Index_out_> > populate(IterationOptions<Index_in_> iopt, ExtractionOptions<Index_in_> eopt, IndexBuffer_ original_indices) const {
        if constexpr(sparse_) {
            // Extract values out to avoid invalidation due to the move.
            bool needs_index = eopt.sparse_extract_index;
            bool needs_value = eopt.sparse_extract_value;

            auto inner = new_extractor<accrow_, true>(ptr.get(), std::move(iopt), std::move(eopt));
            return std::unique_ptr<SparseExtractor<Value_out_, Index_out_> >(new SparseCastExtractor<IndexBuffer_>(std::move(inner), std::move(original_indices), needs_value, needs_index));
        } else {
            auto inner = new_extractor<accrow_, false>(ptr.get(), std::move(iopt), std::move(eopt));
            return std::unique_ptr<DenseExtractor<Value_out_, Index_out_> >(new DenseCastExtractor<IndexBuffer_>(std::move(inner), std::move(original_indices)));
        }
    }

    static void cast_dimension_selection(const DimensionSelection<Index_out_>& src, DimensionSelection<Index_in_>& dest) {
        dest.type = src.type;

        if (dest.type == DimensionSelectionType::BLOCK) {
            dest.block_start = src.block_start;
            dest.block_length = src.block_length;

        } else if (dest.type == DimensionSelectionType::INDEX) {
            if (src.index_start) {
                dest.indices.insert(dest.indices.end(), src.index_start, src.index_start + src.index_length);
            } else {
                dest.indices.insert(dest.indices.end(), src.indices.begin(), src.indices.end());
            }
        }
    }

    template<bool accrow_, bool sparse_>
    std::unique_ptr<Extractor<sparse_, Value_out_, Index_out_> > populate(IterationOptions<Index_out_> iopt, ExtractionOptions<Index_out_> eopt) const {
        if constexpr(same_Index_type_) {
            return populate<accrow_, sparse_>(std::move(iopt), std::move(eopt), false);
        } else {
            // Need to cast the options to the correct integer type.
            IterationOptions<Index_in_> icopy;
            cast_dimension_selection(iopt.selection, icopy.selection);
            icopy.cache_for_reuse = iopt.cache_for_reuse;
            icopy.access_pattern = iopt.access_pattern;
            if (icopy.access_pattern == AccessPattern::SEQUENCE) {
                if (iopt.sequence_start) {
                    icopy.sequence.insert(icopy.sequence.end(), iopt.sequence_start, iopt.sequence_start + iopt.sequence_length);
                } else {
                    icopy.sequence.insert(icopy.sequence.end(), iopt.sequence.begin(), iopt.sequence.end());
                }
            }

            ExtractionOptions<Index_in_> ecopy;
            cast_dimension_selection(eopt.selection, ecopy.selection);
            ecopy.sparse_extract_index = eopt.sparse_extract_index;
            ecopy.sparse_extract_value = eopt.sparse_extract_value;
            ecopy.sparse_ordered_index = eopt.sparse_ordered_index;

            if (eopt.selection.type == DimensionSelectionType::INDEX) {
                // Transferring the correctly-typed vector of selection indices
                // for use in CastExtractor::extracted_index().
                if (eopt.selection.index_start) {
                    eopt.selection.indices.insert(eopt.selection.indices.end(), eopt.selection.index_start, eopt.selection.index_start + eopt.selection.index_length);
                }
                return populate<accrow_, sparse_, std::vector<Index_out_> >(std::move(icopy), std::move(ecopy), std::move(eopt.selection.indices));
            } else {
                return populate<accrow_, sparse_>(std::move(icopy), std::move(ecopy), false);
            }
        }
    }

public:
    std::unique_ptr<DenseExtractor<Value_out_, Index_out_> > dense_row(IterationOptions<Index_out_> iopt, ExtractionOptions<Index_out_> eopt) const {
        return populate<true, false>(std::move(iopt), std::move(eopt));
    }

    std::unique_ptr<DenseExtractor<Value_out_, Index_out_> > dense_column(IterationOptions<Index_out_> iopt, ExtractionOptions<Index_out_> eopt) const {
        return populate<false, false>(std::move(iopt), std::move(eopt));
    }

    std::unique_ptr<SparseExtractor<Value_out_, Index_out_> > sparse_row(IterationOptions<Index_out_> iopt, ExtractionOptions<Index_out_> eopt) const {
        return populate<true, true>(std::move(iopt), std::move(eopt));
    }

    std::unique_ptr<SparseExtractor<Value_out_, Index_out_> > sparse_column(IterationOptions<Index_out_> iopt, ExtractionOptions<Index_out_> eopt) const {
        return populate<false, true>(std::move(iopt), std::move(eopt));
    }
};

/**
 * Recast a `Matrix` to a different interface type.
 *
 * @tparam Value_out_ Data type to cast to.
 * @tparam Index_out_ Index type to cast to.
 * @tparam Matrix_ A realized `Matrix` class, possibly one that is `const`.
 *
 * @param p Pointer to the `Matrix` instance to cast from.
 * @return Pointer to a `Matrix` instance of the desired interface type.
 */
template<typename Value_out_, typename Index_out_, class Matrix_>
std::shared_ptr<Matrix<Value_out_, Index_out_> > make_DelayedCast(std::shared_ptr<Matrix_> p) {
    return std::shared_ptr<Matrix<Value_out_, Index_out_> >(new DelayedCast<Value_out_, Index_out_, typename Matrix_::value_type, typename Matrix_::index_type>(std::move(p)));
}

}

#endif
