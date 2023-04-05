#ifndef TATAMI_DELAYED_SUBSET_BLOCK
#define TATAMI_DELAYED_SUBSET_BLOCK

#include "Matrix.hpp"
#include <algorithm>
#include <memory>

/**
 * @file DelayedSubsetBlock.hpp
 *
 * @brief Delayed subsetting to a single contiguous block.
 *
 * This is a specialized implementation that is more efficient than the `tatami::DelayedSubset` class.
 */

namespace tatami {

/**
 * @brief Delayed subsetting to a contiguous block.
 *
 * Implements delayed subsetting (i.e., slicing) of a matrix to a single contiguous block of rows or columns.
 * This is a specialized implementation that is more efficient than the `tatami::DelayedSubset` class.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam MARGIN Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam T Type of matrix value.
 * @tparam IDX Type of index value.
 */
template<int MARGIN, typename T, typename IDX>
class DelayedSubsetBlock : public Matrix<T, IDX> {
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param f Index of the start of the block. This should be a row index if `MARGIN = 0` and a column index otherwise.
     * @param l Index of the one-past-the-end of the block.
     */
    DelayedSubsetBlock(std::shared_ptr<const Matrix<T, IDX> > p, size_t f, size_t l) : mat(p), block_start(f), block_length(l - f) {}

    /**
     * @copydoc DelayedSubsetBlock
     */
    DelayedSubsetBlock(std::shared_ptr<Matrix<T, IDX> > p, size_t f, size_t l) : mat(p), block_start(f), block_length(l - f) {}

private:
    std::shared_ptr<const Matrix<T, IDX> > mat;
    size_t block_start, block_length;

public:
    size_t nrow() const {
        if constexpr(MARGIN==0) {
            return block_length;
        } else {
            return mat->nrow();
        }
    }

    size_t ncol() const {
        if constexpr(MARGIN==0) {
            return mat->ncol();
        } else {
            return block_length;
        }
    }

    bool sparse() const {
        return mat->sparse();
    }

    bool prefer_rows() const {
        return mat->prefer_rows();
    }

    std::pair<double, double> dimension_preference() const {
        return mat->dimension_preference();
    }

    using Matrix<T, IDX>::column;

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::sparse_column;

    using Matrix<T, IDX>::sparse_row;

public:
    /**
     * @cond
     */
    template<bool WORKROW>
    struct AlongWorkspace : public Workspace<WORKROW> {
        AlongWorkspace(std::shared_ptr<BlockWorkspace<WORKROW> > w) : internal(std::move(w)) {}
        std::shared_ptr<BlockWorkspace<WORKROW> > internal;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowWorkspace> new_row_workspace(bool cache = false) const {
        if constexpr(MARGIN == 0) {
            return mat->new_row_workspace(cache);
        } else {
            return std::shared_ptr<RowWorkspace>(new AlongWorkspace<true>(mat->new_row_workspace(block_start, block_length, cache)));
        }
    }

    std::shared_ptr<ColumnWorkspace> new_column_workspace(bool cache = false) const {
        if constexpr(MARGIN == 0) {
            return std::shared_ptr<ColumnWorkspace>(new AlongWorkspace<false>(mat->new_column_workspace(block_start, block_length, cache)));
        } else {
            return mat->new_column_workspace(cache);
        }
    }

    const T* row(size_t r, T* buffer, RowWorkspace* work) const {
        if constexpr(MARGIN == 0) {
            return mat->row(block_start + r, buffer, work);
        } else {
            return mat->row(r, buffer, static_cast<AlongWorkspace<true>*>(work)->internal.get());
        }
    }

    const T* column(size_t c, T* buffer, ColumnWorkspace* work) const {
        if constexpr(MARGIN == 0) {
            return mat->column(c, buffer, static_cast<AlongWorkspace<false>*>(work)->internal.get());
        } else {
            return mat->column(block_start + c, buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(block_start + r, out_values, out_indices, work, sorted);
        } else {
            return subset_sparse<true>(r, out_values, out_indices, static_cast<AlongWorkspace<true>*>(work)->internal.get(), sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return subset_sparse<false>(c, out_values, out_indices, static_cast<AlongWorkspace<false>*>(work)->internal.get(), sorted);
        } else {
            return mat->sparse_column(block_start + c, out_values, out_indices, work, sorted);
        }
    }

private:
    template<bool WORKROW, class InputWorkspace>
    SparseRange<T, IDX> subset_sparse(size_t i, T* out_values, IDX* out_indices, InputWorkspace* work, bool sorted) const {
        SparseRange<T, IDX> output;

        if constexpr(WORKROW) {
            output = mat->sparse_row(i, out_values, out_indices, work, sorted);
        } else {
            output = mat->sparse_column(i, out_values, out_indices, work, sorted);
        }

        if (block_start) {
            if (out_indices != output.index) {
                std::copy(output.index, output.index + output.number, out_indices);
                output.index = out_indices;
            }
            for (size_t i = 0; i < output.number; ++i) {
                out_indices[i] -= block_start;
            }
        }

        return output;
    }

public:
    /**
     * @cond
     */
    template<bool WORKROW>
    struct AlongBlockWorkspace : public BlockWorkspace<WORKROW> {
        AlongBlockWorkspace(size_t s, size_t l, std::shared_ptr<BlockWorkspace<WORKROW> > w) : details(s, l), internal(std::move(w)) {}
        std::pair<size_t, size_t> details;
        const std::pair<size_t, size_t>& block() const { return details; }
        std::shared_ptr<BlockWorkspace<WORKROW> > internal;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowBlockWorkspace> new_row_workspace(size_t start_, size_t length_, bool cache) const {
        if constexpr(MARGIN == 0) {
            return mat->new_row_workspace(start_, length_, cache);
        } else {
            return std::shared_ptr<RowBlockWorkspace>(new AlongBlockWorkspace<true>(start_, length_, mat->new_row_workspace(block_start + start_, length_, cache)));
        }
    }

    std::shared_ptr<ColumnBlockWorkspace> new_column_workspace(size_t start_, size_t length_, bool cache) const {
        if constexpr(MARGIN == 0) {
            return std::shared_ptr<ColumnBlockWorkspace>(new AlongBlockWorkspace<false>(start_, length_, mat->new_column_workspace(block_start + start_, length_, cache)));
        } else {
            return mat->new_column_workspace(start_, length_, cache);
        }
    }

    const T* row(size_t r, T* buffer, RowBlockWorkspace* work) const {
        if constexpr(MARGIN == 0) {
            return mat->row(block_start + r, buffer, work);
        } else {
            return mat->row(r, buffer, static_cast<AlongBlockWorkspace<true>*>(work)->internal.get());
        }
    }

    const T* column(size_t c, T* buffer, ColumnBlockWorkspace* work) const {
        if constexpr(MARGIN == 0) {
            return mat->column(c, buffer, static_cast<AlongBlockWorkspace<false>*>(work)->internal.get());
        } else {
            return mat->column(block_start + c, buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowBlockWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(block_start + r, out_values, out_indices, work, sorted);
        } else {
            return subset_sparse<true>(r, out_values, out_indices, static_cast<AlongBlockWorkspace<true>*>(work)->internal.get(), sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnBlockWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return subset_sparse<false>(c, out_values, out_indices, static_cast<AlongBlockWorkspace<false>*>(work)->internal.get(), sorted);
        } else {
            return mat->sparse_column(block_start + c, out_values, out_indices, work, sorted);
        }
    }

public:
    /**
     * @cond
     */
    template<bool WORKROW>
    struct AlongIndexWorkspace : public IndexWorkspace<IDX, WORKROW> {
        AlongIndexWorkspace(std::vector<IDX> i) : indices_(std::move(i)) {}

        std::vector<IDX> indices_;
        const std::vector<IDX>& indices() const {
            if (indices_.empty()) {
                return this->internal->indices();
            } else {
                return indices_;
            }
        }

        std::shared_ptr<IndexWorkspace<IDX, WORKROW> > internal;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(std::vector<IDX> indices_, bool cache) const {
        return new_workspace<true>(std::move(indices_), cache);
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX> > new_column_workspace(std::vector<IDX> indices_, bool cache) const {
        return new_workspace<false>(std::move(indices_), cache);
    }

    const T* row(size_t r, T* buffer, RowIndexWorkspace<IDX>* work) const {
        if constexpr(MARGIN == 0) {
            return mat->row(block_start + r, buffer, work);
        } else {
            return mat->row(r, buffer, static_cast<AlongIndexWorkspace<true>*>(work)->internal.get());
        }
    }

    const T* column(size_t c, T* buffer, ColumnIndexWorkspace<IDX>* work) const {
        if constexpr(MARGIN == 0) {
            return mat->column(c, buffer, static_cast<AlongIndexWorkspace<false>*>(work)->internal.get());
        } else {
            return mat->column(block_start + c, buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(block_start + r, out_values, out_indices, work, sorted);
        } else {
            return subset_sparse<true>(r, out_values, out_indices, static_cast<AlongIndexWorkspace<true>*>(work)->internal.get(), sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return subset_sparse<false>(c, out_values, out_indices, static_cast<AlongIndexWorkspace<false>*>(work)->internal.get(), sorted);
        } else {
            return mat->sparse_column(block_start + c, out_values, out_indices, work, sorted);
        }
    }

private:
    template<bool WORKROW>
    std::shared_ptr<IndexWorkspace<IDX, WORKROW> > new_workspace(std::vector<IDX> subset, bool cache) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            if constexpr(WORKROW) {
                return mat->new_row_workspace(std::move(subset), cache);
            } else {
                return mat->new_column_workspace(std::move(subset), cache);
            }
        } else {
            auto ptr = new AlongIndexWorkspace<WORKROW>(std::move(subset));
            std::shared_ptr<IndexWorkspace<IDX, WORKROW> > output(ptr);

            if (block_start) {
                auto copy = ptr->indices_;
                for (auto& x : copy) {
                    x += block_start;
                }
                if constexpr(WORKROW) {
                    ptr->internal = mat->new_row_workspace(std::move(copy), cache);
                } else {
                    ptr->internal = mat->new_column_workspace(std::move(copy), cache);
                }
            } else {
                if constexpr(WORKROW) {
                    ptr->internal = mat->new_row_workspace(std::move(ptr->indices_), cache);
                } else {
                    ptr->internal = mat->new_column_workspace(std::move(ptr->indices_), cache);
                }
                ptr->indices_.clear();
            }

            return output;
        }
    }
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam MARGIN Dimension along which the addition is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam MAT A specialized `Matrix`, to be automatically deducted.
 * @tparam V Vector containing the subset indices, to be automatically deducted.
 *
 * @param p Pointer to the underlying (pre-subset) `Matrix`.
 * @param f Index of the start of the block. This should be a row index if `MARGIN = 0` and a column index otherwise.
 * @param l Index of the one-past-the-end of the block.
 *
 * @return A pointer to a `DelayedSubsetBlock` instance.
 */
template<int MARGIN, class MAT>
std::shared_ptr<MAT> make_DelayedSubsetBlock(std::shared_ptr<MAT> p, size_t f, size_t l) {
    return std::shared_ptr<MAT>(
        new DelayedSubsetBlock<MARGIN, typename MAT::data_type, typename MAT::index_type>(
            p,
            f,
            l
        )
    );
}

}

#endif
