#ifndef TATAMI_DELAYED_ISOMETRIC_OP_H
#define TATAMI_DELAYED_ISOMETRIC_OP_H

#include <memory>
#include "Matrix.hpp"
#include "utils.hpp"
#include "Workspace.hpp"

/**
 * @file DelayedIsometricOp.hpp
 *
 * @brief Delayed isometric operations.
 *
 * This is equivalent to the class of the same name in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @brief Delayed isometric operations on a matrix.
 *
 * Implements any operation that preserves the shape of the matrix and operates on each matrix value independently.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam T Type of matrix value.
 * @tparam OP Functor class implementing the operation.
 * This should accept the row index, column index and value, and return the modified value after applying the operation. 
 * @tparam IDX Type of index value.
 */
template<typename T, typename IDX, class OP>
class DelayedIsometricOp : public Matrix<T, IDX> {
public:
    /**
     * @param p Pointer to the underlying matrix.
     * @param op Instance of the functor class.
     */
    DelayedIsometricOp(std::shared_ptr<const Matrix<T, IDX> > p, OP op) : mat(p), operation(std::move(op)) {}

private:
    std::shared_ptr<const Matrix<T, IDX> > mat;
    OP operation;
    static_assert(std::is_same<T, decltype(operation(0, 0, 0))>::value);

public:
    size_t nrow() const {
        return mat->nrow();
    }
    
    size_t ncol() const {
        return mat->ncol();
    }

    /**
     * @return `true` if both the underlying (pre-operation) matrix is sparse and the operation preserves sparsity.
     * Otherwise returns `false`.
     */
    bool sparse() const {
        if constexpr(OP::sparse) {
            return mat->sparse();
        } else {
            return false;
        }
    }

    /**
     * @return `true` if row-wise extraction is preferred by the underlying (pre-operation) matrix, otherwise returns `false`.
     */
    bool prefer_rows() const { 
        return mat->prefer_rows();
    }

    using Matrix<T, IDX>::column;

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::dense_row_workspace;

    using Matrix<T, IDX>::dense_column_workspace;

    using Matrix<T, IDX>::sparse_row_workspace;

    using Matrix<T, IDX>::sparse_column_workspace;

    /********************************************
     ********** Dense full extraction ***********
     ********************************************/
public:
    std::shared_ptr<DenseRowWorkspace> dense_row_workspace(const WorkspaceOptions& opt) const {
        return mat->dense_row_workspace(opt);
    }

    std::shared_ptr<DenseColumnWorkspace> dense_column_workspace(const WorkspaceOptions& opt) const {
        return mat->dense_column_workspace(opt);
    }

    const T* row(size_t r, T* buffer, DenseRowWorkspace* work) const {
        return operate_on_dimension(r, 0, mat->ncol(), buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnWorkspace* work) const {
        return operate_on_dimension(c, 0, mat->nrow(), buffer, work);
    }

private:
    template<bool WORKROW>
    static void mutate(size_t x, size_t start, size_t end, const T* in, T* out, const OP& operation) { 
        for (size_t i = start; i < end; ++i, ++in) {
            if constexpr(WORKROW) {
                out[i - start] = operation(x, i, *in);
            } else {
                out[i - start] = operation(i, x, *in);
            }
        }
    }

    template<class SomeWorkspace>
    const T* operate_on_dimension(size_t x, size_t start, size_t end, T* buffer, SomeWorkspace* work) const {
        const T* raw = extract_dense<SomeWorkspace::row>(mat.get(), x, buffer, work);
        mutate<SomeWorkspace::row>(x, start, end, raw, buffer, operation);
        return buffer;
    }

    /********************************************
     ********** Sparse full extraction **********
     ********************************************/
private:
    template<bool WORKROW>
    static constexpr bool needs_index() {
        if constexpr(WORKROW) {
            return OP::needs_row;
        } else {
            return OP::needs_column;
        }
    }

    template<bool WORKROW>
    static constexpr bool use_simple_workspace() {
        return (!needs_index<!WORKROW>() && OP::sparse);
    }

    template<class WorkspaceFactory>
    std::shared_ptr<typename WorkspaceFactory::Parent> create_sparse_workspace(const WorkspaceOptions& opt, WorkspaceFactory factory) const {
        if constexpr(use_simple_workspace<WorkspaceFactory::Parent::row>()) {
            return factory.inner_sparse_workspace(mat.get(), opt);

        } else {
            std::shared_ptr<typename WorkspaceFactory::Parent> output;

            if constexpr(OP::sparse) {
                auto ptr = factory.intermediate_sparse_workspace(opt);
                output.reset(ptr);

                if (opt.mode == SparseExtractMode::VALUE) {
                    // Only need to allocate a buffer if the indices are needed AND we didn't extract them in the first place.
                    auto copy = opt;
                    copy.mode = SparseExtractMode::BOTH;
                    ptr->inner = factory.inner_sparse_workspace(mat.get(), copy);
                    ptr->ibuffer.resize(factory.length());
                } else {
                    ptr->inner = factory.inner_sparse_workspace(mat.get(), opt);
                }

            } else {
                auto ptr = factory.intermediate_dense_workspace(opt);
                output.reset(ptr);

                // Going for dense extraction, in which case we don't even need indices.
                if (sparse_extract_value(opt.mode)) {
                    ptr->acquire_inner_workspace(mat.get(), opt);
                }
                ptr->report_index = sparse_extract_index(opt.mode);
            }

            return output;
        }
    }

    template<class SomeWorkspace, class OperatorFactory>
    SparseRange<T, IDX> operate_on_dimension(size_t x, T* vbuffer, IDX* ibuffer, SomeWorkspace* work, OperatorFactory factory) const {
        constexpr bool WORKROW = SomeWorkspace::row;

        if constexpr(use_simple_workspace<WORKROW>()) {
            auto raw = extract_sparse<WORKROW>(mat.get(), x, vbuffer, ibuffer, work);
            if (raw.value) {
                for (size_t i = 0; i < raw.number; ++i) {
                    if constexpr(SomeWorkspace::row) {
                        vbuffer[i] = operation(x, 0, raw.value[i]); // no-op value, we don't need the column indices.
                    } else {
                        vbuffer[i] = operation(0, x, raw.value[i]); // no-op value, we don't need the row indices.
                    }
                }
            } 
            return SparseRange<T, IDX>(raw.number, vbuffer, raw.index);

        } else if constexpr(!OP::sparse) {
            auto wptr = static_cast<typename OperatorFactory::dense_workspace_type*>(work);
            SparseRange<T, IDX> raw(factory.length(), NULL, NULL);

            if (wptr->inner) {
                auto found = extract_dense<WORKROW>(mat.get(), x, vbuffer, wptr->inner.get());
                factory.mutate_values(x, found, vbuffer, operation);
                raw.value = vbuffer;
            }

            if (wptr->report_index) {
                factory.copy_indices(ibuffer);
                raw.index = ibuffer;
            }
            return raw;

        } else {
            auto wptr = static_cast<typename OperatorFactory::sparse_workspace_type*>(work);

            // If the workspace's ibuffer is empty, we're either extracting the indices
            // directly into the user's ibuffer, or we don't need the indices at all.
            // Either way, it doesn't hurt to use the user's ibuffer.
            IDX* iin = (wptr->ibuffer.empty() ? ibuffer : wptr->ibuffer.data());

            auto raw = extract_sparse<WORKROW>(mat.get(), x, vbuffer, iin, wptr->inner.get());

            if (raw.value) {
                for (size_t i = 0; i < raw.number; ++i) {
                    if constexpr(WORKROW) {
                        vbuffer[i] = operation(x, raw.index[i], raw.value[i]);
                    } else {
                        vbuffer[i] = operation(raw.index[i], x, raw.value[i]);
                    }
                }
                raw.value = vbuffer;
            }
            return raw;
        }
    }

private:
    template<bool WORKROW>
    struct SimpleWorkspaceFactory {
        typedef SparseWorkspace<WORKROW> Parent;

        struct S : public Parent {
            std::shared_ptr<SparseWorkspace<WORKROW> > inner;
            std::vector<IDX> ibuffer;
        };

        S* intermediate_sparse_workspace(const WorkspaceOptions&) const {
            return new S;
        }

        auto inner_sparse_workspace(const Matrix<T, IDX>* mat, const WorkspaceOptions& opt) const {
            return new_workspace<WORKROW, true>(mat, opt);
        }

        struct D : public Parent {
            std::shared_ptr<DenseWorkspace<WORKROW> > inner;
            bool report_index;

            void acquire_inner_workspace(const Matrix<T, IDX>* mat, const WorkspaceOptions& opt) {
                inner = new_workspace<WORKROW, false>(mat, opt);
            }
        };

        D* intermediate_dense_workspace(const WorkspaceOptions&) const {
            return new D;
        }
    };

    template<bool WORKROW>
    struct SimpleOperatorFactory {
        SimpleOperatorFactory(size_t l) : len(l) {}

        typedef typename SimpleWorkspaceFactory<WORKROW>::S sparse_workspace_type;
        typedef typename SimpleWorkspaceFactory<WORKROW>::D dense_workspace_type;

        size_t len;

        size_t length() const { return len; }

        void copy_indices(IDX* ibuffer) const {
            std::iota(ibuffer, ibuffer + len, static_cast<IDX>(0));
        }

        void mutate_values(size_t x, const T* in, T* out, const OP& operation) const {
            DelayedIsometricOp<T, IDX, OP>::mutate<WORKROW>(x, 0, len, in, out, operation);
        }
    };

    template<bool WORKROW> friend class SimpleOperatorFactory;

public:
    std::shared_ptr<SparseRowWorkspace> sparse_row_workspace(const WorkspaceOptions& opt) const {
        return create_sparse_workspace(opt, SimpleWorkspaceFactory<true>());
    }

    std::shared_ptr<SparseColumnWorkspace> sparse_column_workspace(const WorkspaceOptions& opt) const {
        return create_sparse_workspace(opt, SimpleWorkspaceFactory<false>());
    }

    SparseRange<T, IDX> row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowWorkspace* work) const {
        return operate_on_dimension(r, vbuffer, ibuffer, work, SimpleOperatorFactory<true>(this->ncol()));
    }

    SparseRange<T, IDX> column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnWorkspace* work) const {
        return operate_on_dimension(c, vbuffer, ibuffer, work, SimpleOperatorFactory<false>(this->nrow()));
    }

    /********************************************
     ********** Dense block extraction **********
     ********************************************/
public:
    std::shared_ptr<DenseRowBlockWorkspace> dense_row_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return mat->dense_row_workspace(start, length, opt);
    }

    std::shared_ptr<DenseColumnBlockWorkspace> dense_column_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return mat->dense_column_workspace(start, length, opt);
    }

    const T* row(size_t r, T* buffer, DenseRowBlockWorkspace* work) const {
        return operate_on_dimension(r, work->start, work->start + work->length, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnBlockWorkspace* work) const {
        return operate_on_dimension(c, work->start, work->start + work->length, buffer, work);
    }

    /********************************************
     ********** Sparse block extraction *********
     ********************************************/
private:
    template<bool WORKROW>
    struct BlockWorkspaceFactory {
        BlockWorkspaceFactory(size_t s, size_t l) : start(s), length(l) {}
   
        size_t start, length;
        typedef SparseBlockWorkspace<WORKROW> Parent;

        struct S : public SparseBlockWorkspace<WORKROW> {
            S(size_t s, size_t l) : SparseBlockWorkspace<WORKROW>(s, l) {}
            std::shared_ptr<SparseBlockWorkspace<WORKROW> > inner;
            std::vector<IDX> ibuffer;
        };

        S* intermediate_sparse_workspace(const WorkspaceOptions&) const {
            return new S(start, length);
        }

        struct D : public SparseBlockWorkspace<WORKROW> {
            D(size_t s, size_t l) : SparseBlockWorkspace<WORKROW>(s, l) {}
            std::shared_ptr<DenseBlockWorkspace<WORKROW> > inner;
            bool report_index;

            void acquire_inner_workspace(const Matrix<T, IDX>* mat, const WorkspaceOptions& opt) {
                inner = new_workspace<WORKROW, false>(mat, this->start, this->length, opt);
            }
        };

        D* intermediate_dense_workspace(const WorkspaceOptions&) const {
            return new D(start, length);
        }

        auto inner_sparse_workspace(const Matrix<T, IDX>* mat, const WorkspaceOptions& opt) const {
            return new_workspace<WORKROW, true>(mat, start, length, opt);
        }
    };

    template<bool WORKROW>
    struct BlockOperatorFactory {
        BlockOperatorFactory(size_t s, size_t l) : start(s), len(l) {}

        typedef typename BlockWorkspaceFactory<WORKROW>::S sparse_workspace_type;
        typedef typename BlockWorkspaceFactory<WORKROW>::D dense_workspace_type;

        size_t start, len;

        size_t length() const { return len; }

        void copy_indices(IDX* ibuffer) const {
            std::iota(ibuffer, ibuffer + len, static_cast<IDX>(start));
        }

        void mutate_values(size_t x, const T* in, T* out, const OP& operation) const {
            DelayedIsometricOp<T, IDX, OP>::mutate<WORKROW>(x, start, start + len, in, out, operation);
        }
    };

    template<bool WORKROW> friend class BlockOperatorFactory;

public:
    std::shared_ptr<SparseRowBlockWorkspace> sparse_row_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return create_sparse_workspace(opt, BlockWorkspaceFactory<true>(start, length));
    }

    std::shared_ptr<SparseColumnBlockWorkspace> sparse_column_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return create_sparse_workspace(opt, BlockWorkspaceFactory<false>(start, length));
    }

    SparseRange<T, IDX> row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowBlockWorkspace* work) const {
        return operate_on_dimension(r, vbuffer, ibuffer, work, BlockOperatorFactory<true>(work->start, work->length));
    }

    SparseRange<T, IDX> column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnBlockWorkspace* work) const {
        return operate_on_dimension(c, vbuffer, ibuffer, work, BlockOperatorFactory<false>(work->start, work->length));
    }

    /********************************************
     ********* Dense indexed extraction *********
     ********************************************/
public:
    std::shared_ptr<DenseRowIndexWorkspace<IDX> > dense_row_workspace(std::vector<IDX> indices, const WorkspaceOptions& opt) const {
        return mat->dense_row_workspace(std::move(indices), opt);
    }

    std::shared_ptr<DenseColumnIndexWorkspace<IDX> > dense_column_workspace(std::vector<IDX> indices, const WorkspaceOptions& opt) const {
        return mat->dense_column_workspace(std::move(indices), opt);
    }

    const T* row(size_t r, T* buffer, DenseRowIndexWorkspace<IDX>* work) const {
        return operate_on_dimension_indexed(r, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnIndexWorkspace<IDX>* work) const {
        return operate_on_dimension_indexed(c, buffer, work);
    }

private:
    template<bool WORKROW>
    static void mutate(size_t x, const std::vector<IDX>& indices, const T* in, T* out, const OP& operation) {
        for (size_t i = 0, end = indices.size(); i < end; ++i, ++in) {
            if constexpr(WORKROW) {
                out[i] = operation(x, indices[i], *in);
            } else {
                out[i] = operation(indices[i], x, *in);
            } 
        }
    }

    template<class SomeWorkspace>
    const T* operate_on_dimension_indexed(size_t x, T* buffer, SomeWorkspace* work) const {
        const T* raw = extract_dense<SomeWorkspace::row>(mat.get(), x, buffer, work);
        mutate<SomeWorkspace::row>(x, work->indices(), raw, buffer, operation);
        return buffer;
    }

    /********************************************
     ********* Sparse indexed extraction ********
     ********************************************/
private:
    template<bool WORKROW>
    struct IndexWorkspaceFactory {
        IndexWorkspaceFactory(std::vector<IDX>& i) : iptr(&i) {}
   
        std::vector<IDX>* iptr;
        typedef SparseIndexWorkspace<IDX, WORKROW> Parent;

        struct S : public Parent {
            S(size_t l) : Parent(l) {}
            std::shared_ptr<SparseIndexWorkspace<IDX, WORKROW> > inner;
            std::vector<IDX> ibuffer;

            const std::vector<IDX>& indices() const { return inner->indices(); }
        };

        S* intermediate_sparse_workspace(const WorkspaceOptions&) const {
            return new S(iptr->size());
        }

        auto inner_sparse_workspace(const Matrix<T, IDX>* mat, const WorkspaceOptions& opt) {
            return new_workspace<WORKROW, true>(mat, std::move(*iptr), opt); // called no more than once!
        }

        struct D : public Parent {
            D(std::vector<IDX> i) : Parent(i.size()), indices_(std::move(i)) {}
            std::shared_ptr<DenseIndexWorkspace<IDX, WORKROW> > inner;
            bool report_index;

            std::vector<IDX> indices_;
            const std::vector<IDX>& indices() const { 
                if (inner) {
                    return inner->indices(); 
                } else {
                    return indices_;
                }
            }

            void acquire_inner_workspace(const Matrix<T, IDX>* mat, const WorkspaceOptions& opt) {
                inner = new_workspace<WORKROW, false>(mat, std::move(indices_), opt); // called no more than once!
            }
        };

        D* intermediate_dense_workspace(const WorkspaceOptions&) {
            return new D(std::move(*iptr)); // called no more than once!
        }
    };

    template<bool WORKROW>
    struct IndexOperatorFactory {
        IndexOperatorFactory(const std::vector<IDX>& i) : iptr(&i) {}

        typedef typename IndexWorkspaceFactory<WORKROW>::S sparse_workspace_type;
        typedef typename IndexWorkspaceFactory<WORKROW>::D dense_workspace_type;

        const std::vector<IDX>* iptr;

        size_t length() const { return iptr->size(); }

        void copy_indices(IDX* ibuffer) const {
            // Copying to avoid lifetime issues with IndexWorkspace's indices.
            std::copy(iptr->begin(), iptr->end(), ibuffer);
        }

        void mutate_values(size_t x, const T* in, T* out, const OP& operation) const {
            DelayedIsometricOp<T, IDX, OP>::mutate<WORKROW>(x, *iptr, in, out, operation);
        }
    };

    template<bool WORKROW> friend struct IndexOperatorFactory;

public:
    std::shared_ptr<SparseRowIndexWorkspace<IDX> > sparse_row_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return create_sparse_workspace(opt, IndexWorkspaceFactory<true>(subset));
    }

    std::shared_ptr<SparseColumnIndexWorkspace<IDX> > sparse_column_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return create_sparse_workspace(opt, IndexWorkspaceFactory<false>(subset));
    }

    SparseRange<T, IDX> row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowIndexWorkspace<IDX>* work) const {
        return operate_on_dimension(r, vbuffer, ibuffer, work, IndexOperatorFactory<true>(work->indices()));
    }

    SparseRange<T, IDX> column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnIndexWorkspace<IDX>* work) const {
        return operate_on_dimension(c, vbuffer, ibuffer, work, IndexOperatorFactory<false>(work->indices()));
    }
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam MAT A specialized `Matrix`, to be automatically deducted.
 * @tparam OP Helper class defining the operation.
 *
 * @param p Pointer to a `Matrix`.
 * @param op Instance of the operation helper class.
 */
template<class MAT, class OP>
std::shared_ptr<MAT> make_DelayedIsometricOp(std::shared_ptr<MAT> p, OP op) {
    return std::shared_ptr<MAT>(
        new DelayedIsometricOp<typename MAT::data_type, typename MAT::index_type, typename std::remove_reference<OP>::type>(
            p,
            std::move(op)
        )
    );
}

}

#include "arith_scalar_helpers.hpp"

#include "arith_vector_helpers.hpp"

#include "math_helpers.hpp"

#endif
