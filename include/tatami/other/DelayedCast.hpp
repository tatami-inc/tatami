#ifndef TATAMI_DELAYED_CAST_HPP
#define TATAMI_DELAYED_CAST_HPP

#include "../base/Matrix.hpp"
#include "../base/utils.hpp"
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
    DelayedCast(std::shared_ptr<const Matrix<Value_in_, Index_in_> > p) : ptr(std::move(p)) {}

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

    double sparse_proportion() const {
        return ptr->sparse_proportion();
    }

    bool prefer_rows() const { 
        return ptr->prefer_rows();
    }

    double prefer_rows_proportion() const {
        return ptr->prefer_rows_proportion();
    }

    bool uses_oracle(bool row) const {
        return ptr->uses_oracle(row);
    }

private:
    std::shared_ptr<const Matrix<Value_in_, Index_in_> > ptr;
    static constexpr bool same_Value_type_ = std::is_same<Value_in_, Value_out_>::value;
    static constexpr bool same_Index_type_ = std::is_same<Index_in_, Index_out_>::value;

private:
    template<DimensionSelectionType selection_, bool sparse_>
    struct CastExtractor : public Extractor<selection_, sparse_, Value_out_, Index_out_> {
        CastExtractor(std::unique_ptr<Extractor<selection_, sparse_, Value_in_, Index_in_> > inner) : internal(std::move(inner)) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = internal->full_length;
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = internal->block_start;
                this->block_length = internal->block_length;
            } else if constexpr(selection_ == DimensionSelectionType::INDEX) {
                this->index_length = internal->index_length;
            }
        }

        CastExtractor(std::unique_ptr<Extractor<selection_, sparse_, Value_in_, Index_in_> > inner, std::vector<Index_out_> idx) : CastExtractor(std::move(inner)) {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                this->indices = std::move(idx);
            }
        }

    protected:
        std::unique_ptr<Extractor<selection_, sparse_, Value_in_, Index_in_> > internal;
        typename std::conditional<same_Index_type_ || selection_ != DimensionSelectionType::INDEX, bool, std::vector<Index_out_> >::type indices;

    public:
        const Index_out_* index_start() const {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                if constexpr(same_Index_type_) {
                    return internal->index_start();
                } else {
                    return indices.data();
                }
            } else {
                return NULL;
            }
        }

    private:
        struct CastOracle : public Oracle<Index_in_> {
            CastOracle(std::unique_ptr<Oracle<Index_out_> > s) : source(std::move(s)) {}

            size_t predict(Index_in_* target, size_t length) {
                buffer.resize(length);
                size_t filled = source->predict(buffer.data(), length);
                std::copy(buffer.begin(), buffer.begin() + filled, target);
                return filled;
            }
        private:
            std::unique_ptr<Oracle<Index_out_> > source;
            std::vector<Index_out_> buffer;
        };

    public:
        void set_oracle(std::unique_ptr<Oracle<Index_out_> > o) {
            internal->set_oracle(std::make_unique<CastOracle>(std::move(o)));
        }
    };

    template<DimensionSelectionType selection_>
    struct DenseCastExtractor : public CastExtractor<selection_, false> {
        template<typename ...Args_>
        DenseCastExtractor(std::unique_ptr<Extractor<selection_, false, Value_in_, Index_in_> > inner, Args_&& ... args) : CastExtractor<selection_, false>(std::move(inner), std::forward<Args_>(args)...) {
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
        template<typename ...Args_>
        SparseCastExtractor(std::unique_ptr<Extractor<selection_, true, Value_in_, Index_in_> > inner, Args_&& ... args) : 
            CastExtractor<selection_, true>(std::move(inner), std::forward<Args_>(args)...) 
        {
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
    template<bool accrow_, DimensionSelectionType selection_, bool sparse_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_out_, Index_out_> > populate(const Options& opt) const {
        std::unique_ptr<Extractor<selection_, sparse_, Value_out_, Index_out_> > output;
        auto inner = new_extractor<accrow_, sparse_>(ptr.get(), opt);
        if constexpr(sparse_) {
            output.reset(new SparseCastExtractor<selection_>(std::move(inner)));
        } else {
            output.reset(new DenseCastExtractor<selection_>(std::move(inner)));
        }
        return output;
    }

    template<bool accrow_, DimensionSelectionType selection_, bool sparse_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_out_, Index_out_> > populate(const Options& opt, Index_out_ block_start, Index_out_ block_length) const {
        std::unique_ptr<Extractor<selection_, sparse_, Value_out_, Index_out_> > output;
        auto inner = new_extractor<accrow_, sparse_>(ptr.get(), static_cast<Index_in_>(block_start), static_cast<Index_in_>(block_length), opt);
        if constexpr(sparse_) {
            output.reset(new SparseCastExtractor<selection_>(std::move(inner)));
        } else {
            output.reset(new DenseCastExtractor<selection_>(std::move(inner)));
        }
        return output;
    }

    template<bool accrow_, DimensionSelectionType selection_, bool sparse_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_out_, Index_out_> > populate(const Options& opt, std::vector<Index_out_> indices) const {
        std::unique_ptr<Extractor<selection_, sparse_, Value_out_, Index_out_> > output;

        if constexpr(!same_Index_type_) {
            std::vector<Index_in_> temp(indices.begin(), indices.end());
            auto inner = new_extractor<accrow_, sparse_>(ptr.get(), std::move(temp), opt);
            if constexpr(sparse_) {
                output.reset(new SparseCastExtractor<selection_>(std::move(inner), std::move(indices)));
            } else {
                output.reset(new DenseCastExtractor<selection_>(std::move(inner), std::move(indices)));
            }
        } else {
            auto inner = new_extractor<accrow_, sparse_>(ptr.get(), std::move(indices), opt);
            if constexpr(sparse_) {
                output.reset(new SparseCastExtractor<selection_>(std::move(inner)));
            } else {
                output.reset(new DenseCastExtractor<selection_>(std::move(inner)));
            }
        }

        return output;
    }

public:
    std::unique_ptr<FullDenseExtractor<Value_out_, Index_out_> > dense_row(const Options& opt) const {
        return populate<true, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_out_, Index_out_> > dense_row(Index_out_ block_start, Index_out_ block_length, const Options& opt) const {
        return populate<true, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_out_, Index_out_> > dense_row(std::vector<Index_out_> indices, const Options& opt) const {
        return populate<true, DimensionSelectionType::INDEX, false>(opt, std::move(indices));
    }

    std::unique_ptr<FullDenseExtractor<Value_out_, Index_out_> > dense_column(const Options& opt) const {
        return populate<false, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_out_, Index_out_> > dense_column(Index_out_ block_start, Index_out_ block_length, const Options& opt) const {
        return populate<false, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_out_, Index_out_> > dense_column(std::vector<Index_out_> indices, const Options& opt) const {
        return populate<false, DimensionSelectionType::INDEX, false>(opt, std::move(indices));
    }

public:
    std::unique_ptr<FullSparseExtractor<Value_out_, Index_out_> > sparse_row(const Options& opt) const {
        return populate<true, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_out_, Index_out_> > sparse_row(Index_out_ block_start, Index_out_ block_length, const Options& opt) const {
        return populate<true, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexSparseExtractor<Value_out_, Index_out_> > sparse_row(std::vector<Index_out_> indices, const Options& opt) const {
        return populate<true, DimensionSelectionType::INDEX, true>(opt, std::move(indices));
    }

    std::unique_ptr<FullSparseExtractor<Value_out_, Index_out_> > sparse_column(const Options& opt) const {
        return populate<false, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_out_, Index_out_> > sparse_column(Index_out_ block_start, Index_out_ block_length, const Options& opt) const {
        return populate<false, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexSparseExtractor<Value_out_, Index_out_> > sparse_column(std::vector<Index_out_> indices, const Options& opt) const {
        return populate<false, DimensionSelectionType::INDEX, true>(opt, std::move(indices));
    }
};

/**
 * Recast a `Matrix` to a different interface type.
 *
 * @tparam Value_out_ Data type to cast to.
 * @tparam Index_out_ Index type to cast to.
 * @tparam Value_in_ Data type to cast from.
 * @tparam Index_in_ Index type to cast from.
 *
 * @param p Pointer to the (possbly `const`) `Matrix` instance to cast from.
 * @return Pointer to a `Matrix` instance of the desired interface type.
 */
template<typename Value_out_, typename Index_out_, typename Value_in_, typename Index_in_>
std::shared_ptr<Matrix<Value_out_, Index_out_> > make_DelayedCast(std::shared_ptr<const Matrix<Value_in_, Index_in_> > p) {
    return std::shared_ptr<Matrix<Value_out_, Index_out_> >(new DelayedCast<Value_out_, Index_out_, Value_in_, Index_in_>(std::move(p)));
}

/**
 * @cond
 */
template<typename Value_out_, typename Index_out_, typename Value_in_, typename Index_in_>
std::shared_ptr<Matrix<Value_out_, Index_out_> > make_DelayedCast(std::shared_ptr<Matrix<Value_in_, Index_in_> > p) {
    return std::shared_ptr<Matrix<Value_out_, Index_out_> >(new DelayedCast<Value_out_, Index_out_, Value_in_, Index_in_>(std::move(p)));
}
/**
 * @endcond
 */


}

#endif
