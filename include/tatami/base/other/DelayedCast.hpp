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

    template<DimensionSelectionType selection_, bool sparse_>
    struct CastExtractor : public Extractor<selection_, sparse_, Value_out_, Index_out_> {
        CastExtractor(std::unique_ptr<Extractor<selection_, sparse_, Value_in_, Index_in_> > i) : internal(std::move(i)) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = internal->full_length;
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = internal->block_start;
                this->block_length = internal->block_length;
            } else {
                this->index_length = internal->index_length;
                if constexpr(!same_Index_type_) {
                    auto ptr = internal->index_start();
                    indices_buffer = std::vector<Index_out_>(ptr, ptr + this->index_length);
                }
            }
        }

        const Index_out_* index_start() const {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                if constexpr(same_Index_type_) {
                    return internal->index_start();
                } else {
                    return indices_buffer.data();
                }
            } else {
                return NULL;
            }
        }

    protected:
        std::unique_ptr<Extractor<selection_, sparse_, Value_in_, Index_in_> > internal;
        typename std::conditional<same_Index_type_ || selection_ != DimensionSelectionType::INDEX, bool, std::vector<Index_out_> >::type indices_buffer;
    };

    template<DimensionSelectionType selection_>
    struct DenseCastExtractor : public CastExtractor<selection_, false> {
        DenseCastExtractor(std::unique_ptr<Extractor<selection_, false, Value_in_, Index_in_> > i) : CastExtractor<selection_, false>(std::move(i)) {
            if constexpr(!same_Value_type_) {
                internal_buffer.resize(extracted_length<selection_, Index_out_>(*this));
            }
        }

        const Value_out_* fetch(Index_out_ i, Value_out_* buffer) {
            if constexpr(same_Value_type_) {
                return this->internal->fetch(i, buffer);
            } else {
                auto out = this->internal->fetch(i, internal_buffer.data());
                std::copy(out, out + extracted_length<selection_, Index_out_>(*this), buffer);
                return buffer;
            }
        }

    protected:
        typename std::conditional<same_Value_type_, bool, std::vector<Value_in_> >::type internal_buffer;
    };

    template<DimensionSelectionType selection_>
    struct SparseCastExtractor : public CastExtractor<selection_, true> {
        SparseCastExtractor(std::unique_ptr<Extractor<selection_, true, Value_in_, Index_in_> > i) : CastExtractor<selection_, true>(std::move(i)) {
            if constexpr(!same_Value_type_) {
                internal_vbuffer.resize(extracted_length<selection_, Index_out_>(*this));
            }
            if constexpr(!same_Index_type_) {
                internal_ibuffer.resize(extracted_length<selection_, Index_out_>(*this));
            }
        }

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
    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_out_, Index_out_> > populate_core(const Options<Index_in_>& opt, Args_... args) const {
        if constexpr(sparse_) {
            auto inner = new_extractor<accrow_, true>(ptr.get(), args..., opt);
            return std::unique_ptr<Extractor<selection_, true, Value_out_, Index_out_> >(new SparseCastExtractor<selection_>(std::move(inner)));
        } else {
            auto inner = new_extractor<accrow_, false>(ptr.get(), args..., opt);
            return std::unique_ptr<Extractor<selection_, false, Value_out_, Index_out_> >(new DenseCastExtractor<selection_>(std::move(inner)));
        }
    }

    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_out_, Index_out_> > populate_cast(const Options<Index_in_>& opt) const {
        return populate_core<accrow_, selection_, sparse_>(opt);
    }

    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_out_, Index_out_> > populate_cast(const Options<Index_in_>& opt, Index_out_ block_start, Index_out_ block_length) const {
        return populate_core<accrow_, selection_, sparse_>(opt, static_cast<Index_in_>(block_start), static_cast<Index_in_>(block_length));
    }

    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_out_, Index_out_> > populate_cast(const Options<Index_in_>& opt, const Index_out_* index_start, size_t index_length) const {
        std::vector<Index_in_> temp(index_start, index_start + index_length);
        return populate_core<accrow_, selection_, sparse_>(opt, temp.data(), index_length);
    }

    struct CastOracle : SequenceOracle<Index_in_> {
        CastOracle(std::shared_ptr<SequenceOracle<Index_out_> > s) : source(std::move(s)) {}

        std::shared_ptr<SequenceOracle<Index_out_> > source;
        std::vector<Index_in_> buffer;

        std::pair<const Index_in_*, size_t> predict(size_t n) {
            auto raw = source->predict(n);
            buffer.clear();
            buffer.insert(buffer.end(), raw.first, raw.first + raw.second);
            return std::make_pair(buffer.data(), raw.second);
        }
    };

    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_out_, Index_out_> > populate(const Options<Index_out_>& opt, Args_... args) const {
        if constexpr(same_Index_type_) {
            return populate_core<accrow_, selection_, sparse_>(opt, args...);

        } else {
            // Need to copy everything over to a new options type.
            Options<Index_in_> optcopy;
            optcopy.sparse = opt.sparse;

            optcopy.access.cache_for_reuse = opt.access.cache_for_reuse;
            if (opt.access.pattern) {
                optcopy.access.pattern.reset(new CastOracle(opt.access.pattern));
            }

            return populate_cast<accrow_, selection_, sparse_>(optcopy, args...);
        }
    }

public:
    std::unique_ptr<FullDenseExtractor<Value_out_, Index_out_> > dense_row(const Options<Index_out_>& opt) const {
        return populate<true, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_out_, Index_out_> > dense_row(Index_out_ block_start, Index_out_ block_length, const Options<Index_out_>& opt) const {
        return populate<true, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_out_, Index_out_> > dense_row(const Index_out_* index_start, size_t index_length, const Options<Index_out_>& opt) const {
        return populate<true, DimensionSelectionType::INDEX, false>(opt, index_start, index_length);
    }

    std::unique_ptr<FullDenseExtractor<Value_out_, Index_out_> > dense_column(const Options<Index_out_>& opt) const {
        return populate<false, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_out_, Index_out_> > dense_column(Index_out_ block_start, Index_out_ block_length, const Options<Index_out_>& opt) const {
        return populate<false, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_out_, Index_out_> > dense_column(const Index_out_* index_start, size_t index_length, const Options<Index_out_>& opt) const {
        return populate<false, DimensionSelectionType::INDEX, false>(opt, index_start, index_length);
    }

public:
    std::unique_ptr<FullSparseExtractor<Value_out_, Index_out_> > sparse_row(const Options<Index_out_>& opt) const {
        return populate<true, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_out_, Index_out_> > sparse_row(Index_out_ block_start, Index_out_ block_length, const Options<Index_out_>& opt) const {
        return populate<true, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexSparseExtractor<Value_out_, Index_out_> > sparse_row(const Index_out_* index_start, size_t index_length, const Options<Index_out_>& opt) const {
        return populate<true, DimensionSelectionType::INDEX, true>(opt, index_start, index_length);
    }

    std::unique_ptr<FullSparseExtractor<Value_out_, Index_out_> > sparse_column(const Options<Index_out_>& opt) const {
        return populate<false, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_out_, Index_out_> > sparse_column(Index_out_ block_start, Index_out_ block_length, const Options<Index_out_>& opt) const {
        return populate<false, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexSparseExtractor<Value_out_, Index_out_> > sparse_column(const Index_out_* index_start, size_t index_length, const Options<Index_out_>& opt) const {
        return populate<false, DimensionSelectionType::INDEX, true>(opt, index_start, index_length);
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
