#ifndef TATAMI_DELAYED_UNARY_ISOMETRIC_OP_H
#define TATAMI_DELAYED_UNARY_ISOMETRIC_OP_H

#include <memory>
#include "../../base/Matrix.hpp"
#include "../../base/utils.hpp"

/**
 * @file DelayedUnaryIsometricOp.hpp
 *
 * @brief Delayed unary isometric operations.
 *
 * This is equivalent to the class of the same name in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @brief Delayed isometric operations on a single matrix.
 *
 * Implements any operation that preserves the shape of the matrix and operates on each matrix value independently.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `DenseExtractor::fetch()` or friends.
 * We only consider "unary" operations that involve a single `Matrix` - see `DelayedBinaryIsometricOp` for operations between two `Matrix` instances.
 *
 * The `Operation_` class is expected to provide the following static `constexpr` member variables:
 *
 * - `needs_row`: whether the operation needs the row index, e.g., if the operation varies by row.
 * - `needs_column`: whether the operation needs the column index, e.g., if the operation varies by column.
 * - `always_dense`: whether the operation always yields a dense result when applied on sparse data.
 * - `always_sparse`: whether the operation always yields a sparse result when applied on sparse data.
 *   This should not be `true` if `always_dense` is `true`.
 * 
 * The class should implement the following method:
 *
 * - `void dense<row_>(Index_ i, Index_ start, Index_ length, Value_* buffer) const`: 
 *   This method should apply the operation to all values in `buffer`, which contains a contiguous block of elements from row `i` (when `row_ = true`).
 *   This block spans columns from `[start, start + length)`.
 *   If `row_ = false`, `i` is instead a column and the block contains rows.
 * - `void dense<row_>(Index_ i, const Index_* indices, Index_ length, Value_* buffer) const`: 
 *   This method should apply the operation to all values in `buffer`, which contains to a subset of elements from row `i` (when `row_ = true`).
 *   Column indices of the elements in the subset are defined by `indices`, which is a sorted and unique array of length `length`w.
 *   If `row_ = false`, `i` is instead a column and `indices` contains rows.
 * - `void sparse<row_>(Index_ i, Index_ number, Value_* buffer, const Index_* indices) const`:
 *   This method should apply the operation to all values in `buffer`, which contains `number` elements from row `i` (when `row_ = true`).
 *   The column indices of all elements is contained in `indices` - unless `needs_column = false`, in which case `indices` may be `NULL`.
 *   If `row_ = false`, `i` is instead a column and `indices` contains row indices (or is `NULL`, if `needs_row = false`).
 * 
 * If neither of `always_dense` or `always_sparse` is `true`, the class should also implement:
 *
 * - `bool actual_sparse() const`: whether this particular instance of the operation yields a sparse result when applied on sparse data of type `Value_`.
 *   For example, an addition operation remains sparse if the added value is zero.
 *
 * If `always_sparse` is not `true`, the class should also implement:
 *
 * - `Value_ zero<row_>(Index_ i) const`:
 *   This method should return the result of applying the operation on a zero input for row `i` (when `row_ = true`) or column `i` (otherwise).
 * 
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Class implementing the operation.
 */
template<typename Value_, typename Index_, class Operation_>
class DelayedUnaryIsometricOp : public Matrix<Value_, Index_> {
public:
    /**
     * @param p Pointer to the underlying matrix.
     * @param op Instance of the functor class.
     */
    DelayedUnaryIsometricOp(std::shared_ptr<const Matrix<Value_, Index_> > p, Operation_ op) : mat(std::move(p)), operation(std::move(op)) {}

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    Operation_ operation;
    static_assert(!Operation_::always_dense || !Operation_::always_sparse);

public:
    Index_ nrow() const {
        return mat->nrow();
    }
    
    Index_ ncol() const {
        return mat->ncol();
    }

    /**
     * @return `true` if both the underlying (pre-operation) matrix is sparse and the operation preserves sparsity.
     * Otherwise returns `false`.
     */
    bool sparse() const {
        if constexpr(Operation_::always_sparse) {
            return mat->sparse();
        } else if constexpr(!Operation_::always_dense) {
            if (operation.actual_sparse()) {
                return mat->sparse();
            } 
        }
        return false;
    }

    double sparse_proportion() const {
        if constexpr(Operation_::always_sparse) {
            return mat->sparse_proportion();
        } else if constexpr(!Operation_::always_dense) {
            if (operation.actual_sparse()) {
                return mat->sparse_proportion();
            }
        }
        return 0;
    }

    bool prefer_rows() const { 
        return mat->prefer_rows();
    }

    double prefer_rows_proportion() const { 
        return mat->prefer_rows_proportion();
    }

    bool uses_oracle(bool row) const {
        return mat->uses_oracle(row);
    }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

private:
    template<DimensionSelectionType selection_, bool sparse_, bool inner_sparse_>
    struct IsometricExtractorBase : public Extractor<selection_, sparse_, Value_, Index_> {
        IsometricExtractorBase(const DelayedUnaryIsometricOp* p, std::unique_ptr<Extractor<selection_, inner_sparse_, Value_, Index_> > i) : parent(p), internal(std::move(i)) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = internal->full_length;
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = internal->block_start;
                this->block_length = internal->block_length;
            } else {
                this->index_length = internal->index_length;
            }
        }

        const Index_* index_start() const {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                return internal->index_start();
            } else {
                return NULL;
            }
        }

        void set_oracle(std::unique_ptr<Oracle<Index_> > o) {
            internal->set_oracle(std::move(o));
        }

    protected:
        const DelayedUnaryIsometricOp* parent;

        std::unique_ptr<Extractor<selection_, inner_sparse_, Value_, Index_> > internal;
    };

    /**************************************
     ********** Dense extraction **********
     **************************************/
private:
    /**
     * Used if:
     *
     * - the underlying matrix is dense.
     *
     * OR
     *
     * - the underlying matrix is sparse
     * - the operation discards sparsity in a variable manner.
     */
    template<bool accrow_, DimensionSelectionType selection_> 
    struct DenseIsometricExtractor_Basic : public IsometricExtractorBase<selection_, false, false> {
        template<typename ... Args_>
        DenseIsometricExtractor_Basic(const DelayedUnaryIsometricOp* p, const Options& opt, Args_&& ... args) :
            IsometricExtractorBase<selection_, false, false>(p, new_extractor<accrow_, false>(p->mat.get(), std::forward<Args_>(args)..., opt)) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            this->internal->fetch_copy(i, buffer);

            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->parent->operation.template dense<accrow_>(i, 0, this->full_length, buffer);
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->parent->operation.template dense<accrow_>(i, this->block_start, this->block_length, buffer);
            } else {
                this->parent->operation.template dense<accrow_>(i, this->internal->index_start(), this->index_length, buffer);
            }

            return buffer;
        }
    };

    template<bool accrow_,class Extractor_>
    void populate_index_mapping(Extractor_& ext, std::vector<Index_>& index_mapping) const {
        auto isize = ext->index_length;
        if (isize) {
            auto istart = ext->index_start();
            index_mapping.resize(accrow_ ? mat->ncol() : mat->nrow());
            for (Index_ i = 0; i < isize; ++i) {
                index_mapping[istart[i]] = i;
            }
        }
    }

    /**
     * Used if:
     *
     * - the underlying matrix is sparse
     * - the operation preserves sparsity
     * 
     * OR
     *
     * - the underlying matrix is sparse
     * - the operation discards sparsity in a constant manner.
     */
    template<bool accrow_, DimensionSelectionType selection_> 
    struct DenseIsometricExtractor_FromSparse : public IsometricExtractorBase<selection_, false, true> {
        template<typename ... Args_>
        DenseIsometricExtractor_FromSparse(const DelayedUnaryIsometricOp* p, Options opt, Args_&& ... args) :
            IsometricExtractorBase<selection_, false, true>(p, [&]{
                auto copy = opt;
                copy.sparse_extract_value = true;
                copy.sparse_extract_index = true;
                return new_extractor<accrow_, true>(p->mat.get(), std::forward<Args_>(args)..., copy);
            }()) 
        {
            auto n = extracted_length<selection_, Index_>(*(this->internal));
            internal_vbuffer.resize(n);
            internal_ibuffer.resize(n);
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                p->template populate_index_mapping<accrow_>(this->internal, index_mapping);
            }
        }

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto vbuffer = internal_vbuffer.data();
            auto range = this->internal->fetch(i, vbuffer, internal_ibuffer.data());
            if (range.value != vbuffer) {
                // Don't use fetch_copy as we don't want to copy the indices.
                std::copy(range.value, range.value + range.number, vbuffer);
            }
            this->parent->operation.template sparse<accrow_>(i, range.number, vbuffer, range.index);

            auto full_length = extracted_length<selection_, Index_>(*(this->internal));
            if (range.number < full_length) { // avoid calling zero() if possible, as this might throw zero-related errors in non-IEEE platforms.
                std::fill(buffer, buffer + full_length, [&]{
                    if constexpr(Operation_::always_sparse) {
                        return static_cast<Value_>(0);
                    } else if constexpr(Operation_::always_dense) {
                        return this->parent->operation.template zero<accrow_>(i);
                    } else if (this->parent->operation.actual_sparse()) {
                        return static_cast<Value_>(0);
                    } else {
                        return this->parent->operation.template zero<accrow_>(i);
                    }
                }());
            }

            if constexpr(selection_ == DimensionSelectionType::FULL) {
                for (Index_ i = 0; i < range.number; ++i) {
                    buffer[range.index[i]] = vbuffer[i];
                }
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                auto shift = this->internal->block_start;
                for (Index_ i = 0; i < range.number; ++i) {
                    buffer[range.index[i] - shift] = vbuffer[i];
                }
            } else {
                for (Index_ i = 0; i < range.number; ++i) {
                    buffer[index_mapping[range.index[i]]] = vbuffer[i];
                }
            }

            return buffer;
        }

    private:
        std::vector<Value_> internal_vbuffer;
        std::vector<Index_> internal_ibuffer;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type index_mapping;
    };

    /***************************************
     ********** Sparse extraction **********
     ***************************************/
private:
    template<DimensionSelectionType selection_, class Extractor_>
    static void fill_dense_indices(const Extractor_& ext, Index_* ibuffer) {
        if constexpr(selection_ == DimensionSelectionType::FULL) {
            std::iota(ibuffer, ibuffer + ext->full_length, 0);
        } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
            std::iota(ibuffer, ibuffer + ext->block_length, ext->block_start);
        } else {
            auto xptr = ext->index_start();
            std::copy(xptr, xptr + ext->index_length, ibuffer);
        }
    }

    /**
     * Used if:
     *
     * - the underlying matrix is dense
     *
     * OR 
     *
     * - the underlying matrix is sparse
     * - the operation discards sparsity in a variable manner
     */
    template<bool accrow_, DimensionSelectionType selection_> 
    struct SparseIsometricExtractor_FromDense : public IsometricExtractorBase<selection_, true, false> {
        template<typename ... Args_>
        SparseIsometricExtractor_FromDense(const DelayedUnaryIsometricOp* p, const Options& opt, Args_&& ... args) :
            IsometricExtractorBase<selection_, true, false>(p, new_extractor<accrow_, false>(p->mat.get(), std::forward<Args_>(args)..., opt)),
            needs_value(opt.sparse_extract_value),
            needs_index(opt.sparse_extract_index)
        {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            SparseRange<Value_, Index_> output(extracted_length<selection_, Index_>(*(this->internal)), NULL, NULL);

            if (needs_value) {
                this->internal->fetch_copy(i, vbuffer);
                if constexpr(selection_ == DimensionSelectionType::FULL) {
                    this->parent->operation.template dense<accrow_>(i, 0, this->full_length, vbuffer);
                } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                    this->parent->operation.template dense<accrow_>(i, this->block_start, this->block_length, vbuffer);
                } else {
                    this->parent->operation.template dense<accrow_>(i, this->internal->index_start(), this->index_length, vbuffer);
                }
                output.value = vbuffer;
            }

            if (needs_index) {
                fill_dense_indices<selection_>(this->internal, ibuffer);
                output.index = ibuffer;
            }

            return output;
        }

    private:
        bool needs_value;
        bool needs_index;
    };

    /**
     * Used if:
     *
     * - the underlying matrix is sparse
     * - the operation preserves sparsity
     * - indices are not necessary to perform the operation 
     */
    template<bool accrow_, DimensionSelectionType selection_> 
    struct SparseIsometricExtractor_Simple : public IsometricExtractorBase<selection_, true, true> {
        template<typename ... Args_>
        SparseIsometricExtractor_Simple(const DelayedUnaryIsometricOp* p, const Options& opt, Args_&& ... args) :
            IsometricExtractorBase<selection_, true, true>(p, new_extractor<accrow_, true>(p->mat.get(), std::forward<Args_>(args)..., opt)) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            auto raw = this->internal->fetch(i, vbuffer, ibuffer);

            if (raw.value) {
                if (raw.value != vbuffer) {
                    // Don't use fetch_copy as we don't want to copy the indices.
                    std::copy(raw.value, raw.value + raw.number, vbuffer);
                }

                this->parent->operation.template sparse<accrow_>(i, raw.number, vbuffer, raw.index);
                raw.value = vbuffer;
            }

            return raw;
        }
    };

    /**
     * Used if:
     *
     * - the underlying matrix is sparse
     * - the operation preserves sparsity
     * - indices are necessary to perform the operation 
     */
    template<bool accrow_, DimensionSelectionType selection_> 
    struct SparseIsometricExtractor_NeedsIndices : public IsometricExtractorBase<selection_, true, true> {
        template<typename ... Args_>
        SparseIsometricExtractor_NeedsIndices(const DelayedUnaryIsometricOp* p, const Options& opt, Args_&& ... args) :
            IsometricExtractorBase<selection_, true, true>(p, [&]{
                // The index is only necessary to (i) compute the operation on the values 
                // and then (ii) insert those values into the dense buffer. So, there's 
                // no need to extract the index if we don't even want the values.
                auto copy = opt;
                if (opt.sparse_extract_value) {
                    copy.sparse_extract_index = true;
                }
                return new_extractor<accrow_, true>(p->mat.get(), std::forward<Args_>(args)..., copy);
            }()),
            report_index(opt.sparse_extract_index)
        {
            if (!opt.sparse_extract_index && opt.sparse_extract_value) {
                // We only need an internal ibuffer if the user wants the
                // values but didn't provide enough space to store the indices
                // (which we need to pass to the operation's functor).
                internal_ibuffer.resize(extracted_length<selection_, Index_>(*this));
            }
        }

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            // If the internal ibuffer is empty, we're either extracting the indices
            // directly into the user's ibuffer, or we don't need the indices at all.
            // Either way, it doesn't hurt to use the user's ibuffer.
            Index_* iin = (internal_ibuffer.empty() ? ibuffer : internal_ibuffer.data());

            auto raw = this->internal->fetch(i, vbuffer, iin);

            if (raw.value) {
                if (raw.value != vbuffer) {
                    // Don't use fetch_copy as we don't want to copy the indices.
                    std::copy(raw.value, raw.value + raw.number, vbuffer);
                }

                this->parent->operation.template sparse<accrow_>(i, raw.number, vbuffer, raw.index);
                raw.value = vbuffer;
            }

            if (!report_index) {
                raw.index = NULL;
            } 

            return raw;
        }

    protected:
        std::vector<Index_> internal_ibuffer;
        bool report_index = false;
    };

    /**
     * Used if:
     *
     * - the underlying matrix is sparse
     * - the operation discards sparsity in a constant manner
     */
    template<bool accrow_, DimensionSelectionType selection_> 
    struct SparseIsometricExtractor_ForcedDense : public IsometricExtractorBase<selection_, true, true> {
        template<typename ... Args_>
        SparseIsometricExtractor_ForcedDense(const DelayedUnaryIsometricOp* p, Options opt, Args_&& ... args) :
            IsometricExtractorBase<selection_, true, true>(p, [&]{
                // Same logic as in SparseIsometricExtractor_NeedsIndices.
                auto copy = opt;
                if (opt.sparse_extract_value) {
                    copy.sparse_extract_index = true;
                }
                return new_extractor<accrow_, true>(p->mat.get(), std::forward<Args_>(args)..., copy);
            }()),
            report_index(opt.sparse_extract_index)
        {
            if (opt.sparse_extract_value) {
                auto n = extracted_length<selection_, Index_>(*(this->internal));
                internal_vbuffer.resize(n);
                if (!opt.sparse_extract_index) {
                    // Same logic as in SparseIsometricExtractor_NeedsIndices.
                    internal_ibuffer.resize(n);
                }

                if constexpr(selection_ == DimensionSelectionType::INDEX) {
                    p->template populate_index_mapping<accrow_>(this->internal, index_mapping);
                }
            }
        }

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            // Same logic as in SparseIsometricExtractor_NeedsIndices.
            Index_* iin = (internal_ibuffer.empty() ? ibuffer : internal_ibuffer.data());

            auto tmp_vbuffer = internal_vbuffer.data();
            auto range = this->internal->fetch(i, tmp_vbuffer, iin);
            SparseRange<Value_, Index_> output(extracted_length<selection_, Index_>(*this), NULL, NULL);

            if (range.value) { 
                if (range.value != tmp_vbuffer) {
                    // Don't use fetch_copy as we don't want to copy the indices.
                    std::copy(range.value, range.value + range.number, tmp_vbuffer);
                }
                this->parent->operation.template sparse<accrow_>(i, range.number, tmp_vbuffer, range.index);

                auto N = extracted_length<selection_, Index_>(*(this->internal));
                if (range.number < N) { // only invoking 'zero()' if we really need to, as this could throw various zero-related errors.
                    std::fill(vbuffer, vbuffer + N, [&]{
                        if constexpr(Operation_::always_sparse) {
                            return 0; // this never actually gets called, we just want to protect the zero() from a compile-time requirement.
                        } else {
                            return this->parent->operation.template zero<accrow_>(i);
                        }
                    }());
                }

                if constexpr(selection_ == DimensionSelectionType::FULL) {
                    for (Index_ i = 0; i < range.number; ++i) {
                        vbuffer[range.index[i]] = tmp_vbuffer[i];
                    }
                } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                    auto shift = this->internal->block_start;
                    for (Index_ i = 0; i < range.number; ++i) {
                        vbuffer[range.index[i] - shift] = tmp_vbuffer[i];
                    }
                } else {
                    for (Index_ i = 0; i < range.number; ++i) {
                        vbuffer[index_mapping[range.index[i]]] = tmp_vbuffer[i];
                    }
                }

                output.value = vbuffer;
            }

            if (report_index) {
                fill_dense_indices<selection_>(this->internal, ibuffer);
                output.index = ibuffer;
            }

            return output;
        }

    private:
        std::vector<Value_> internal_vbuffer;
        std::vector<Index_> internal_ibuffer;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type index_mapping;
        bool report_index = false;
    };

    /**********************************************
     ********** Public extractor methods **********
     **********************************************/
private:
    bool preserves_sparsity() const {
        if constexpr(Operation_::always_sparse) {
            return true;
        } else if constexpr(!Operation_::always_dense) {
            return operation.actual_sparse();
        }
        return false;
    }

    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > propagate(const Options& opt, Args_&& ... args) const {
        std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > output;

        if constexpr(!sparse_) {
            if (!(mat->sparse())) {
                output.reset(new DenseIsometricExtractor_Basic<accrow_, selection_>(this, opt, std::forward<Args_>(args)...));

            } else {
                if (preserves_sparsity()) {
                    output.reset(new DenseIsometricExtractor_FromSparse<accrow_, selection_>(this, opt, std::forward<Args_>(args)...));

                } else {
                    if constexpr(accrow_ && !Operation_::needs_column) {      // constant sparsity breaking
                        output.reset(new DenseIsometricExtractor_FromSparse<accrow_, selection_>(this, opt, std::forward<Args_>(args)...));
                    } else if constexpr(!accrow_ && !Operation_::needs_row) { // constant sparsity breaking
                        output.reset(new DenseIsometricExtractor_FromSparse<accrow_, selection_>(this, opt, std::forward<Args_>(args)...));
                    } else {                                                  // variable sparsity breaking
                        output.reset(new DenseIsometricExtractor_Basic<accrow_, selection_>(this, opt, std::forward<Args_>(args)...));
                    }
                }
            }

        } else {
            if (!(mat->sparse())) {
                output.reset(new SparseIsometricExtractor_FromDense<accrow_, selection_>(this, opt, std::forward<Args_>(args)...));

            } else {
                if (preserves_sparsity()) {
                    if constexpr(accrow_ && !Operation_::needs_column) {      // no need for column indices
                        output.reset(new SparseIsometricExtractor_Simple<accrow_, selection_>(this, opt, std::forward<Args_>(args)...));
                    } else if constexpr(!accrow_ && !Operation_::needs_row) { // no need for row indices
                        output.reset(new SparseIsometricExtractor_Simple<accrow_, selection_>(this, opt, std::forward<Args_>(args)...));
                    } else {                                                  // guess we need indices
                        output.reset(new SparseIsometricExtractor_NeedsIndices<accrow_, selection_>(this, opt, std::forward<Args_>(args)...));
                    }

                } else {
                    if constexpr(accrow_ && !Operation_::needs_column) {      // constant sparsity breaking
                        output.reset(new SparseIsometricExtractor_ForcedDense<accrow_, selection_>(this, opt, std::forward<Args_>(args)...));
                    } else if constexpr(!accrow_ && !Operation_::needs_row) { // constant sparsity breaking
                        output.reset(new SparseIsometricExtractor_ForcedDense<accrow_, selection_>(this, opt, std::forward<Args_>(args)...));
                    } else {                                                  // variable sparsity breaking
                        output.reset(new SparseIsometricExtractor_FromDense<accrow_, selection_>(this, opt, std::forward<Args_>(args)...));
                    }
                }
            }
        }

        return output;
    }

public:
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        return propagate<true, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return propagate<true, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }
    
    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        return propagate<true, DimensionSelectionType::INDEX, false>(opt, std::move(indices));
    }

    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        return propagate<false, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return propagate<false, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }
    
    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        return propagate<false, DimensionSelectionType::INDEX, false>(opt, std::move(indices));
    }

public:
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_row(const Options& opt) const {
        return propagate<true, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return propagate<true, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }
    
    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& opt) const {
        return propagate<true, DimensionSelectionType::INDEX, true>(opt, std::move(indices));
    }

    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_column(const Options& opt) const {
        return propagate<false, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return propagate<false, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }
    
    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options& opt) const {
        return propagate<false, DimensionSelectionType::INDEX, true>(opt, std::move(indices));
    }
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Helper class defining the operation.
 *
 * @param p Pointer to a (possibly `const`) `Matrix`.
 * @param op Instance of the operation helper class.
 *
 * @return Instance of a `DelayedUnaryIsometricOp` clas.
 */
template<typename Value_, typename Index_, class Operation_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedUnaryIsometricOp(std::shared_ptr<const Matrix<Value_, Index_> > p, Operation_ op) {
    typedef typename std::remove_reference<Operation_>::type Op_;
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedUnaryIsometricOp<Value_, Index_, Op_>(std::move(p), std::move(op)));
}

/**
 * @cond
 */
// For automatic template deduction with non-const pointers.
template<typename Value_, typename Index_, class Operation_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedUnaryIsometricOp(std::shared_ptr<Matrix<Value_, Index_> > p, Operation_ op) {
    typedef typename std::remove_reference<Operation_>::type Op_;
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedUnaryIsometricOp<Value_, Index_, Op_>(std::move(p), std::move(op)));
}

// For back-compatibility.
template<typename ... Args_>
auto make_DelayedUnaryIsometricOp(Args_&&... args) {
    return make_DelayedUnaryIsometricOp(std::forward<Args_>(args)...);
}

template<typename Value_, typename Index_, class Operation_>
using DelayedIsometricOp = DelayedUnaryIsometricOp<Value_, Index_, Operation_>;
/**
 * @endcond
 */

}

#include "arith_helpers.hpp"

#include "math_helpers.hpp"

#include "compare_helpers.hpp"

#include "boolean_helpers.hpp"

#endif
