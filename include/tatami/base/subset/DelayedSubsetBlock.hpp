#ifndef TATAMI_DELAYED_SUBSET_BLOCK
#define TATAMI_DELAYED_SUBSET_BLOCK

#include "utils.hpp"
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

    using Matrix<T, IDX>::dense_column_workspace;

    using Matrix<T, IDX>::dense_row_workspace;

    using Matrix<T, IDX>::sparse_column_workspace;

    using Matrix<T, IDX>::sparse_row_workspace;

private:
    template<bool WORKROW, bool SPARSE>
    struct AlongWorkspace : public Workspace<WORKROW, SPARSE> {
        AlongWorkspace(std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > w) : internal(std::move(w)) {}
        std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > internal;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<Workspace<WORKROW, SPARSE> > create_new_workspace(const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return new_workspace<WORKROW, SPARSE>(mat.get(), opt);
        } else {
            return std::shared_ptr<Workspace<WORKROW, SPARSE> >(new AlongWorkspace<WORKROW, SPARSE>(new_workspace<WORKROW, SPARSE>(mat.get(), block_start, block_length, opt)));
        }
    } 

    template<bool WORKROW>
    const T* get_dense(size_t i, T* buffer, DenseWorkspace<WORKROW>* work) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return extract_dense<WORKROW>(mat.get(), block_start + i, buffer, work);
        } else {
            return extract_dense<WORKROW>(mat.get(), i, buffer, static_cast<AlongWorkspace<WORKROW, false>*>(work)->internal.get());
        }
    }

    template<class InputWorkspace>
    SparseRange<T, IDX> subset_sparse(size_t i, T* vbuffer, IDX* ibuffer, InputWorkspace* work) const {
        auto output = extract_sparse<InputWorkspace::row>(mat.get(), i, vbuffer, ibuffer, work);

        if (block_start && output.index) {
            for (size_t i = 0; i < output.number; ++i) {
                ibuffer[i] = output.index[i] - block_start;
            }
            output.index = ibuffer;
        }

        return output;
    }

    template<bool WORKROW>
    SparseRange<T, IDX> get_sparse(size_t i, T* vbuffer, IDX* ibuffer, SparseWorkspace<WORKROW>* work) const {
        if constexpr((MARGIN==0) == WORKROW) {
            return extract_sparse<WORKROW>(mat.get(), block_start + i, vbuffer, ibuffer, work);
        } else {
            return subset_sparse(i, vbuffer, ibuffer, static_cast<AlongWorkspace<WORKROW, true>*>(work)->internal.get());
        }
    }

public:
    std::shared_ptr<DenseRowWorkspace> dense_row_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<true, false>(opt);
    }

    std::shared_ptr<DenseColumnWorkspace> dense_column_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<false, false>(opt);
    }

    const T* row(size_t r, T* buffer, DenseRowWorkspace* work) const {
        return get_dense<true>(r, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnWorkspace* work) const {
        return get_dense<false>(c, buffer, work);
    }

    std::shared_ptr<SparseRowWorkspace> sparse_row_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<true, true>(opt);
    }

    std::shared_ptr<SparseColumnWorkspace> sparse_column_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<false, true>(opt);
    }

    SparseRange<T, IDX> row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowWorkspace* work) const {
        return get_sparse<true>(r, vbuffer, ibuffer, work);
    }

    SparseRange<T, IDX> column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnWorkspace* work) const {
        return get_sparse<false>(c, vbuffer, ibuffer, work);
    }

private:
    template<bool WORKROW, bool SPARSE>
    struct AlongBlockWorkspace : public BlockWorkspace<WORKROW, SPARSE> {
        AlongBlockWorkspace(size_t s, std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > w) : 
            BlockWorkspace<WORKROW, SPARSE>(s, w->length), internal(std::move(w)) {}
        std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > internal;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > create_new_workspace(size_t start_, size_t length_, const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return new_workspace<WORKROW, SPARSE>(mat.get(), start_, length_, opt);
        } else {
            auto ptr = new AlongBlockWorkspace<WORKROW, SPARSE>(start_, new_workspace<WORKROW, SPARSE>(mat.get(), block_start + start_, length_, opt));
            return std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> >(ptr);
        }
    }

    template<bool WORKROW>
    const T* get_dense(size_t i, T* buffer, DenseBlockWorkspace<WORKROW>* work) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return extract_dense<WORKROW>(mat.get(), block_start + i, buffer, work);
        } else {
            return extract_dense<WORKROW>(mat.get(), i, buffer, static_cast<AlongBlockWorkspace<WORKROW, false>*>(work)->internal.get());
        }
    }

    template<bool WORKROW>
    SparseRange<T, IDX> get_sparse(size_t i, T* vbuffer, IDX* ibuffer, SparseBlockWorkspace<WORKROW>* work) const {
        if constexpr((MARGIN==0) == WORKROW) {
            return extract_sparse<WORKROW>(mat.get(), block_start + i, vbuffer, ibuffer, work);
        } else {
            return subset_sparse(i, vbuffer, ibuffer, static_cast<AlongBlockWorkspace<WORKROW, true>*>(work)->internal.get());
        }
    }

public:
    std::shared_ptr<DenseRowBlockWorkspace> dense_row_workspace(size_t start_, size_t length_, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, false>(start_, length_, opt);
    }

    std::shared_ptr<DenseColumnBlockWorkspace> dense_column_workspace(size_t start_, size_t length_, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, false>(start_, length_, opt);
    }

    const T* row(size_t r, T* buffer, DenseRowBlockWorkspace* work) const {
        return get_dense<true>(r, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnBlockWorkspace* work) const {
        return get_dense<false>(c, buffer, work);
    }

    std::shared_ptr<SparseRowBlockWorkspace> sparse_row_workspace(size_t start_, size_t length_, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, true>(start_, length_, opt);
    }

    std::shared_ptr<SparseColumnBlockWorkspace> sparse_column_workspace(size_t start_, size_t length_, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, true>(start_, length_, opt);
    }

    SparseRange<T, IDX> row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowBlockWorkspace* work) const {
        return get_sparse<true>(r, vbuffer, ibuffer, work);
    }

    SparseRange<T, IDX> column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnBlockWorkspace* work) const {
        return get_sparse<false>(c, vbuffer, ibuffer, work);
    }

private:
    template<bool WORKROW, bool SPARSE>
    struct AlongIndexWorkspace : public IndexWorkspace<IDX, WORKROW, SPARSE> {
        AlongIndexWorkspace(std::vector<IDX> i) : IndexWorkspace<IDX, WORKROW, SPARSE>(i.size()), indices_(std::move(i)) {}

        std::vector<IDX> indices_;
        const std::vector<IDX>& indices() const {
            if (indices_.empty() && this->internal) {
                return this->internal->indices();
            } else {
                return indices_;
            }
        }

        std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > internal;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > create_new_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return new_workspace<WORKROW, SPARSE>(mat.get(), std::move(subset), opt);

        } else {
            auto ptr = new AlongIndexWorkspace<WORKROW, SPARSE>(std::move(subset));
            std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > output(ptr);

            if (block_start) {
                auto copy = ptr->indices_;
                for (auto& x : copy) {
                    x += block_start;
                }
                ptr->internal = new_workspace<WORKROW, SPARSE>(mat.get(), std::move(copy), opt);
            } else {
                ptr->internal = new_workspace<WORKROW, SPARSE>(mat.get(), std::move(ptr->indices_), opt);
                ptr->indices_.clear();
            }

            return output;
        }
    }

    template<bool WORKROW>
    const T* get_dense(size_t i, T* buffer, DenseIndexWorkspace<IDX, WORKROW>* work) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return extract_dense<WORKROW>(mat.get(), block_start + i, buffer, work);
        } else {
            return extract_dense<WORKROW>(mat.get(), i, buffer, static_cast<AlongIndexWorkspace<WORKROW, false>*>(work)->internal.get());
        }
    }

    template<bool WORKROW>
    SparseRange<T, IDX> get_sparse(size_t i, T* vbuffer, IDX* ibuffer, SparseIndexWorkspace<IDX, WORKROW>* work) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return extract_sparse<WORKROW>(mat.get(), block_start + i, vbuffer, ibuffer, work);
        } else {
            return subset_sparse(i, vbuffer, ibuffer, static_cast<AlongIndexWorkspace<WORKROW, true>*>(work)->internal.get());
        }
    }

public:
    std::shared_ptr<DenseRowIndexWorkspace<IDX> > dense_row_workspace(std::vector<IDX> indices_, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, false>(std::move(indices_), opt);
    }

    std::shared_ptr<DenseColumnIndexWorkspace<IDX> > dense_column_workspace(std::vector<IDX> indices_, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, false>(std::move(indices_), opt);
    }

    const T* row(size_t r, T* buffer, DenseRowIndexWorkspace<IDX>* work) const {
        return get_dense<true>(r, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnIndexWorkspace<IDX>* work) const {
        return get_dense<false>(c, buffer, work);
    }

    std::shared_ptr<SparseRowIndexWorkspace<IDX> > sparse_row_workspace(std::vector<IDX> indices_, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, true>(std::move(indices_), opt);
    }

    std::shared_ptr<SparseColumnIndexWorkspace<IDX> > sparse_column_workspace(std::vector<IDX> indices_, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, true>(std::move(indices_), opt);
    }

    SparseRange<T, IDX> row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowIndexWorkspace<IDX>* work) const {
        return get_sparse<true>(r, vbuffer, ibuffer, work);
    }

    SparseRange<T, IDX> column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnIndexWorkspace<IDX>* work) const {
        return get_sparse<false>(c, vbuffer, ibuffer, work);
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
