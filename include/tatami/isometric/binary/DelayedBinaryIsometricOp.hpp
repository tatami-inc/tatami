#ifndef TATAMI_DELAYED_BINARY_ISOMETRIC_OP_H
#define TATAMI_DELAYED_BINARY_ISOMETRIC_OP_H

#include <memory>
#include <deque>
#include "../../base/Matrix.hpp"
#include "../../base/utils.hpp"

/**
 * @file DelayedBinaryIsometricOp.hpp
 *
 * @brief Delayed binary isometric operations.
 *
 * This is equivalent to the class of the same name in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @brief Delayed isometric operations on two matrices
 *
 * Implements any operation that takes two matrices of the same shape and returns another matrix of that shape.
 * Each entry of the output matrix is a function of the corresponding values in the two input matrices.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `DenseExtractor::fetch()` or friends.
 *
 * The `Operation_` class is expected to provide the following static `constexpr` member variables:
 *
 * - `always_sparse`: whether the operation can be optimized to return a sparse result if both input matrices are sparse.
 * 
 * The class should implement the following method:
 *
 * - `void dense<row_>(Index_ i, Index_ start, Index_ length, Value_* left_buffer, const Value_* right_buffer) const`: 
 *   This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`,
 *   each of which contain a contiguous block of elements from row `i` of the left and right matrices, respectively (when `row_ = true`).
 *   The result of the operation should be stored in `left_buffer`.
 *   The block starts at column `start` and is of length `length`.
 *   If `row_ = false`, `i` is instead a column and the block starts at row `start`.
 * - `void dense<row_>(Index_ i, const Index_* indices, Index_ length, Value_* buffer1, const Value_* buffer2) const`: 
 *   This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`,
 *   each of which contain a subset of elements from row `i` of the left and right matrices, respectively (when `row_ = true`).
 *   The result of the operation should be stored in `left_buffer`.
 *   The subset is defined by column indices in the `indices` array of length `length`.
 *   If `row_ = false`, `i` is instead a column and `indices` contains rows.
 * 
 * If `always_sparse = true`, the class should implement:
 *
 * - `Index_ sparse<row_, needs_value, needs_index>(Index_ i, const SparseRange<Value_, Index_>& left, const SparseRange<Value_, Index_>& right, Value_* value_buffer, Index_* index_buffer) const`:
 *   This method should apply the operation to the sparse values in `left` and `right`, 
 *   consisting of the contents of row `i` from the left and right matrices, respectively (when `row_ = true`).
 *   All non-zero values resulting from the operation should be stored in `value_buffer` if `needs_value = true`, otherwise `value_buffer = NULL` and should be ignored.
 *   The corresponding indices of those values should be stored in `index_buffer` if `needs_index = true`, otherwise `index_buffer = NULL` and should be ignored.
 *   The return value should be the number of structural non-zero elements in the output buffers.
 *   If `row_ = false`, the contents of `left` and `right` are taken from column `i` instead.
 *   Note that all values in `left` and `right` are already sorted by increasing index.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Class implementing the operation.
 */
template<typename Value_, typename Index_, class Operation_>
class DelayedBinaryIsometricOp : public Matrix<Value_, Index_> {
public:
    /**
     * @param l Pointer to the left matrix.
     * @param r Pointer to the right matrix.
     * @param op Instance of the functor class.
     */
    DelayedBinaryIsometricOp(std::shared_ptr<const Matrix<Value_, Index_> > l, std::shared_ptr<const Matrix<Value_, Index_> > r, Operation_ op) : 
        left(std::move(l)), right(std::move(r)), operation(std::move(op)) 
    {
        if (left->nrow() != right->nrow() || left->ncol() != right->ncol()) {
            throw std::runtime_error("shape of the left and right matrices should be the same");
        }

        prefer_rows_proportion_internal = (left->prefer_rows_proportion() + right->prefer_rows_proportion()) / 2;
    }

private:
    std::shared_ptr<const Matrix<Value_, Index_> > left, right;
    Operation_ operation;
    double prefer_rows_proportion_internal;

public:
    Index_ nrow() const {
        return left->nrow();
    }

    Index_ ncol() const {
        return left->ncol();
    }

    /**
     * @return `true` if both underlying (pre-operation) matrices are sparse and the operation preserves sparsity.
     * Otherwise returns `false`.
     */
    bool sparse() const {
        if constexpr(Operation_::always_sparse) {
            return left->sparse() && right->sparse();
        }
        return false;
    }

    double sparse_proportion() const {
        if constexpr(Operation_::always_sparse) {
            // Well, better than nothing.
            return (left->sparse_proportion() + right->sparse_proportion())/2;
        }
        return 0;
    }

    bool prefer_rows() const { 
        return prefer_rows_proportion_internal > 0.5;
    }

    double prefer_rows_proportion() const { 
        return prefer_rows_proportion_internal;
    }

    bool uses_oracle(bool row) const {
        return left->uses_oracle(row) || right->uses_oracle(row);
    }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

private:
    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, bool inner_sparse_ = sparse_>
    struct IsometricExtractorBase : public Extractor<selection_, sparse_, Value_, Index_> {
        IsometricExtractorBase(
            const DelayedBinaryIsometricOp* p, 
            std::unique_ptr<Extractor<selection_, inner_sparse_, Value_, Index_> > l,
            std::unique_ptr<Extractor<selection_, inner_sparse_, Value_, Index_> > r
        ) : 
            parent(p), 
            left_internal(std::move(l)),
            right_internal(std::move(r))
        {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = left_internal->full_length;
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = left_internal->block_start;
                this->block_length = left_internal->block_length;
            } else {
                this->index_length = left_internal->index_length;
            }
        }

        const Index_* index_start() const {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                return left_internal->index_start();
            } else {
                return NULL;
            }
        }

    protected:
        const DelayedBinaryIsometricOp* parent;
        std::unique_ptr<Extractor<selection_, inner_sparse_, Value_, Index_> > left_internal, right_internal;

    private:
        // Need to basically clone the oracle stream.
        struct ParentOracle {
            ParentOracle(std::unique_ptr<Oracle<Index_> > o) : source(std::move(o)) {}

            size_t fill(bool left, Index_* buffer, size_t number) {
                auto& current = (left ? left_counter : right_counter);
                size_t end = current + number;
                size_t available = stream.size();

                if (available >= end) {
                    std::copy(stream.begin() + current, stream.begin() + end, buffer);
                    current = end;
                    return number;
                }

                size_t handled = 0;
                if (current < available) {
                    std::copy(stream.begin() + current, stream.end(), buffer);
                    handled = available - current;
                    buffer += handled;
                    number -= handled;
                }

                size_t filled = source->predict(buffer, number);
                current = available + filled;

                // Try to slim down if the accumulated stream has gotten too big.
                if (stream.size() >= 10000) { 
                    size_t minimum = std::min(left_counter, right_counter);
                    if (minimum) {
                        stream.erase(stream.begin(), stream.begin() + minimum);
                        left_counter -= minimum;
                        right_counter -= minimum;
                    }
                }

                stream.insert(stream.end(), buffer, buffer + filled);
                return filled + handled;
            }
        private:
            std::unique_ptr<Oracle<Index_> > source;
            std::deque<Index_> stream;
            size_t left_counter = 0, right_counter = 0;
        };

        struct ChildOracle : public Oracle<Index_> {
            ChildOracle(ParentOracle* o, bool l) : parent(o), left(l) {}
            size_t predict(Index_* buffer, size_t number) {
                return parent->fill(left, buffer, number);
            }
        private:
            ParentOracle* parent;
            bool left;
        };

        std::unique_ptr<ParentOracle> parent_oracle;

    public:
        void set_oracle(std::unique_ptr<Oracle<Index_> > o) {
            auto left_use = parent->left->uses_oracle(accrow_);
            auto right_use = parent->right->uses_oracle(accrow_);

            if (left_use && right_use) {
                parent_oracle.reset(new ParentOracle(std::move(o)));
                left_internal->set_oracle(std::make_unique<ChildOracle>(parent_oracle.get(), true));
                right_internal->set_oracle(std::make_unique<ChildOracle>(parent_oracle.get(), false));
            } else if (left_use) {
                left_internal->set_oracle(std::move(o));
            } else if (right_use) {
                right_internal->set_oracle(std::move(o));
            }
        }
    };

    /**************************************
     ********** Dense extraction **********
     **************************************/
private:
    template<bool accrow_, DimensionSelectionType selection_> 
    struct DenseIsometricExtractor : public IsometricExtractorBase<accrow_, selection_, false> {
        DenseIsometricExtractor(
            const DelayedBinaryIsometricOp* p, 
            std::unique_ptr<Extractor<selection_, false, Value_, Index_> > l, 
            std::unique_ptr<Extractor<selection_, false, Value_, Index_> > r 
        ) : 
            IsometricExtractorBase<accrow_, selection_, false, false>(p, std::move(l), std::move(r))
        {
            holding_buffer.resize(extracted_length<selection_, Index_>(*this));
        }

        const Value_* fetch(Index_ i, Value_* buffer) {
            this->left_internal->fetch_copy(i, buffer);
            auto rptr = this->right_internal->fetch(i, holding_buffer.data());

            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->parent->operation.template dense<accrow_>(i, 0, this->full_length, buffer, rptr);
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->parent->operation.template dense<accrow_>(i, this->block_start, this->block_length, buffer, rptr);
            } else {
                this->parent->operation.template dense<accrow_>(i, this->left_internal->index_start(), this->index_length, buffer, rptr);
            }

            return buffer;
        }

    private:
        std::vector<Value_> holding_buffer;
    };

    /***************************************
     ********** Sparse extraction **********
     ***************************************/
private:
    template<bool accrow_, DimensionSelectionType selection_> 
    struct RegularSparseIsometricExtractor : public IsometricExtractorBase<accrow_, selection_, true> {
        RegularSparseIsometricExtractor(
            const DelayedBinaryIsometricOp* p, 
            std::unique_ptr<Extractor<selection_, true, Value_, Index_> > l, 
            std::unique_ptr<Extractor<selection_, true, Value_, Index_> > r, 
            bool rv,
            bool ri
        ) : 
            IsometricExtractorBase<accrow_, selection_, true, true>(p, std::move(l), std::move(r)), 
            report_value(rv),
            report_index(ri)
        {
            auto n = extracted_length<selection_, Index_>(*this);
            left_internal_ibuffer.resize(n);
            right_internal_ibuffer.resize(n);

            if (report_value) {
                left_internal_vbuffer.resize(n);
                right_internal_vbuffer.resize(n);
            }
        }

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            auto left_ranges = this->left_internal->fetch(i, left_internal_vbuffer.data(), left_internal_ibuffer.data());
            auto right_ranges = this->right_internal->fetch(i, right_internal_vbuffer.data(), right_internal_ibuffer.data());

            SparseRange<Value_, Index_> output(0, NULL, NULL);

            if (report_value && report_index) {
                output.number = this->parent->operation.template sparse<accrow_, true, true, Value_, Index_>(i, left_ranges, right_ranges, vbuffer, ibuffer);
                output.value = vbuffer;
                output.index = ibuffer;
            } else if (report_value) {
                output.number = this->parent->operation.template sparse<accrow_, true, false, Value_, Index_>(i, left_ranges, right_ranges, vbuffer, NULL);
                output.value = vbuffer;
            } else if (report_index) {
                output.number = this->parent->operation.template sparse<accrow_, false, true, Value_, Index_>(i, left_ranges, right_ranges, NULL, ibuffer);
                output.index = ibuffer;
            } else {
                output.number = this->parent->operation.template sparse<accrow_, false, false, Value_, Index_>(i, left_ranges, right_ranges, NULL, NULL);
            }

            return output;
        }

    protected:
        std::vector<Value_> left_internal_vbuffer, right_internal_vbuffer;
        std::vector<Index_> left_internal_ibuffer, right_internal_ibuffer;
        bool report_value = false;
        bool report_index = false;
    };

    /*******************************************
     ********** Un-sparsed extraction **********
     *******************************************/
private:
    // Technically, we could avoid constructing the internal extractor if
    // we don't want the values, but that's a pretty niche optimization,
    // so we won't bother doing that.
    template<bool accrow_, DimensionSelectionType selection_>
    struct DensifiedSparseIsometricExtractor : public IsometricExtractorBase<accrow_, selection_, true, false> {
        DensifiedSparseIsometricExtractor(
            const DelayedBinaryIsometricOp* p, 
            std::unique_ptr<Extractor<selection_, false, Value_, Index_> > l, 
            std::unique_ptr<Extractor<selection_, false, Value_, Index_> > r,
            bool rv,
            bool ri
        ) :
            IsometricExtractorBase<accrow_, selection_, true, false>(p, std::move(l), std::move(r)), 
            report_value(rv),
            report_index(ri) 
        {
            holding_buffer.resize(extracted_length<selection_, Index_>(*this));
        }

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            SparseRange<Value_, Index_> output(extracted_length<selection_, Index_>(*this), NULL, NULL);

            if (report_value) {
                this->left_internal->fetch_copy(i, vbuffer);
                auto rptr = this->right_internal->fetch(i, holding_buffer.data());

                if constexpr(!Operation_::always_sparse) {
                    if constexpr(selection_ == DimensionSelectionType::FULL) {
                        this->parent->operation.template dense<accrow_>(i, 0, this->full_length, vbuffer, rptr);
                    } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                        this->parent->operation.template dense<accrow_>(i, this->block_start, this->block_length, vbuffer, rptr);
                    } else {
                        this->parent->operation.template dense<accrow_>(i, this->left_internal->index_start(), this->index_length, vbuffer, rptr);
                    }
                }

                output.value = vbuffer;
            }

            if (report_index) {
                if constexpr(selection_ == DimensionSelectionType::FULL) {
                    std::iota(ibuffer, ibuffer + this->full_length, 0);
                } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                    std::iota(ibuffer, ibuffer + this->block_length, this->block_start);
                } else {
                    auto xptr = this->left_internal->index_start();
                    std::copy(xptr, xptr + this->index_length, ibuffer);
                }

                output.index = ibuffer;
            }

            return output;
        }

    protected:
        std::vector<Value_> holding_buffer;
        bool report_value = false;
        bool report_index = false;
    };

    /**********************************************
     ********** Public extractor methods **********
     **********************************************/
private:
    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > propagate(const Options& opt, Args_ ... args) const {
        std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > output;

        if constexpr(!sparse_) {
            auto left_inner = new_extractor<accrow_, false>(left.get(), args..., opt); // Explicit copy of the variadic args here.
            auto right_inner = new_extractor<accrow_, false>(right.get(), std::move(args)..., opt); // Do a move once we don't need them anymore.
            output.reset(new DenseIsometricExtractor<accrow_, selection_>(this, std::move(left_inner), std::move(right_inner)));

        } else if constexpr(Operation_::always_sparse) {
            bool report_value = opt.sparse_extract_value;
            bool report_index = opt.sparse_extract_index;

            auto optcopy = opt;
            optcopy.sparse_extract_index = true; // We need the indices to combine things properly.
            optcopy.sparse_ordered_index = true; // Make life easier for operation implementers.

            auto left_inner = new_extractor<accrow_, true>(left.get(), args..., optcopy);
            auto right_inner = new_extractor<accrow_, true>(right.get(), std::move(args)..., optcopy);
            output.reset(new RegularSparseIsometricExtractor<accrow_, selection_>(this, std::move(left_inner), std::move(right_inner), report_value, report_index));

        } else {
            bool report_value = opt.sparse_extract_value;
            bool report_index = opt.sparse_extract_index;
            auto left_inner = new_extractor<accrow_, false>(left.get(), args..., opt);
            auto right_inner = new_extractor<accrow_, false>(right.get(), std::move(args)..., opt);
            output.reset(new DensifiedSparseIsometricExtractor<accrow_, selection_>(this, std::move(left_inner), std::move(right_inner), report_value, report_index));
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
 * @param left Pointer to a (possibly `const`) `Matrix`.
 * @param right Pointer to a (possibly `const`) `Matrix`.
 * @param op Instance of the operation helper class.
 *
 * @return Instance of a `DelayedBinaryIsometricOp` clas.
 */
template<typename Value_, typename Index_, class Operation_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedBinaryIsometricOp(std::shared_ptr<const Matrix<Value_, Index_> > left, std::shared_ptr<const Matrix<Value_, Index_> > right, Operation_ op) {
    typedef typename std::remove_reference<Operation_>::type Op_;
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedBinaryIsometricOp<Value_, Index_, Op_>(std::move(left), std::move(right), std::move(op)));
}

/**
 * @cond
 */
// For automatic template deduction with non-const pointers.
template<typename Value_, typename Index_, class Operation_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedBinaryIsometricOp(std::shared_ptr<Matrix<Value_, Index_> > left, std::shared_ptr<Matrix<Value_, Index_> > right, Operation_ op) {
    typedef typename std::remove_reference<Operation_>::type Op_;
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedBinaryIsometricOp<Value_, Index_, Op_>(std::move(left), std::move(right), std::move(op)));
}
/**
 * @endcond
 */

}

#include "arith_helpers.hpp"

#include "compare_helpers.hpp"

#include "boolean_helpers.hpp"

#endif
