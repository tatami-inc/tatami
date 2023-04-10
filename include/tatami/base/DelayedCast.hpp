#ifndef TATAMI_DELAYED_CAST_HPP
#define TATAMI_DELAYED_CAST_HPP

#include "Matrix.hpp"
#include "utils.hpp"
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
 * @tparam T_out Data type to cast to.
 * @tparam IDX_out Index type to cast to.
 * @tparam T_in Data type to cast from.
 * @tparam IDX_in Index type to cast from.
 */
template<typename T_out, typename IDX_out, typename T_in, typename IDX_in>
class DelayedCast : public Matrix<T_out, IDX_out> {
public:
    /**
     * @param p Pointer to the `Matrix` instance to cast from.
     */
    DelayedCast(std::shared_ptr<Matrix<T_in, IDX_in> > p) : ptr(std::move(p)) {}

public:
    size_t nrow() const {
        return ptr->nrow();
    }

    size_t ncol() const {
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
    std::shared_ptr<Matrix<T_in, IDX_in> > ptr;
    static constexpr bool same_T_type = std::is_same<T_in, T_out>::value;
    static constexpr bool same_IDX_type = std::is_same<IDX_in, IDX_out>::value;

    struct DenseBase {
        DenseBase(size_t n, const WorkspaceOptions& opt) : buffer(!same_T_type ? n : 0) {}

        std::vector<T_in> buffer;
    };

    struct SparseBase {
        SparseBase(size_t n, const WorkspaceOptions& opt) : 
            vbuffer(sparse_extract_value(opt.mode) && !same_T_type ? n : 0),
            ibuffer(sparse_extract_index(opt.mode) && !same_IDX_type ? n : 0)
            {}

        std::vector<T_in> vbuffer;
        std::vector<IDX_in> ibuffer;
    };

    template<bool SPARSE>
    using ConditionalBase = typename std::conditional<SPARSE, SparseBase, DenseBase>::type;

private:
    template<bool WORKROW, bool SPARSE>
    struct CastWorkspace : public Workspace<WORKROW, SPARSE>, public ConditionalBase<SPARSE> {
        CastWorkspace(size_t n, const WorkspaceOptions& opt, std::shared_ptr<Workspace<WORKROW, SPARSE> > p) : 
            ConditionalBase<SPARSE>(n, opt), internal(std::move(p)) {}

        std::shared_ptr<Workspace<WORKROW, SPARSE> > internal;
    };

    typedef CastWorkspace<true, false> DenseCastRowWorkspace;
    typedef CastWorkspace<false, false> DenseCastColumnWorkspace;
    typedef CastWorkspace<true, true> SparseCastRowWorkspace;
    typedef CastWorkspace<false, true> SparseCastColumnWorkspace;

public:
    std::shared_ptr<DenseRowWorkspace> dense_row_workspace(const WorkspaceOptions& opt) const { 
        return std::shared_ptr<DenseRowWorkspace>(new DenseCastRowWorkspace(this->ncol(), opt, ptr->dense_row_workspace(opt)));
    }

    std::shared_ptr<DenseColumnWorkspace> dense_column_workspace(const WorkspaceOptions& opt) const { 
        return std::shared_ptr<DenseColumnWorkspace>(new DenseCastColumnWorkspace(this->nrow(), opt, ptr->dense_column_workspace(opt)));
    }

    const T_out* row(size_t r, T_out* buffer, DenseRowWorkspace* work) const {
        return cast_dense<true>(r, buffer, ptr->ncol(), static_cast<DenseCastRowWorkspace*>(work));
    }

    const T_out* column(size_t c, T_out* buffer, DenseColumnWorkspace* work) const {
        return cast_dense<false>(c, buffer, ptr->nrow(), static_cast<DenseCastColumnWorkspace*>(work));
    }

    std::shared_ptr<SparseRowWorkspace> sparse_row_workspace(const WorkspaceOptions& opt) const { 
        return std::shared_ptr<SparseRowWorkspace>(new SparseCastRowWorkspace(this->ncol(), opt, ptr->sparse_row_workspace(opt)));
    }

    std::shared_ptr<SparseColumnWorkspace> sparse_column_workspace(const WorkspaceOptions& opt) const { 
        return std::shared_ptr<SparseColumnWorkspace>(new SparseCastColumnWorkspace(this->nrow(), opt, ptr->sparse_column_workspace(opt)));
    }

    SparseRange<T_out, IDX_out> row(size_t r, T_out* vbuffer, IDX_out* ibuffer, SparseRowWorkspace* work) const {
        return cast_sparse<true>(r, vbuffer, ibuffer, static_cast<SparseCastRowWorkspace*>(work));
    }

    SparseRange<T_out, IDX_out> column(size_t c, T_out* vbuffer, IDX_out* ibuffer, SparseColumnWorkspace* work) const {
        return cast_sparse<false>(c, vbuffer, ibuffer, static_cast<SparseCastColumnWorkspace*>(work));
    }

private:
    template<bool WORKROW, class InputWorkspace>
    const T_out* cast_dense(size_t i, T_out* buffer, size_t length, InputWorkspace* work) const {
        if constexpr(same_T_type) {
            if constexpr(WORKROW) {
                return ptr->row(i, buffer, work->internal.get()); 
            } else {
                return ptr->column(i, buffer, work->internal.get()); 
            }
        } else {
            const T_in* out;
            if constexpr(WORKROW) {
                out = ptr->row(i, work->buffer.data(), work->internal.get());
            } else {
                out = ptr->column(i, work->buffer.data(), work->internal.get());
            }
            std::copy(out, out + length, buffer);
            return buffer;
        }
    }

    template<bool WORKROW, class InputWorkspace>
    SparseRange<T_out, IDX_out> cast_sparse(size_t i, T_out* vbuffer, IDX_out* ibuffer, InputWorkspace* work) const {
        auto inwork = work->internal.get();

        if constexpr(same_T_type) {
            if constexpr(same_IDX_type) {
                return extract_sparse<WORKROW>(ptr.get(), i, vbuffer, ibuffer, inwork);

            } else {
                auto out = extract_sparse<WORKROW>(ptr.get(), i, vbuffer, work->ibuffer.data(), inwork);
                if (out.index) {
                    std::copy(out.index, out.index + out.number, ibuffer);
                } else {
                    ibuffer = NULL;
                }
                return SparseRange<T_out, IDX_out>(out.number, out.value, ibuffer);
            }

        } else {
            if constexpr(same_IDX_type) {
                auto out = extract_sparse<WORKROW>(ptr.get(), i, work->vbuffer.data(), ibuffer, inwork);
                if (out.value) {
                    std::copy(out.value, out.value + out.number, vbuffer);
                } else {
                    vbuffer = NULL;
                }
                return SparseRange<T_out, IDX_out>(out.number, vbuffer, out.index);

            } else {
                auto out = extract_sparse<WORKROW>(ptr.get(), i, work->vbuffer.data(), work->ibuffer.data(), inwork);
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
                return SparseRange<T_out, IDX_out>(out.number, vbuffer, ibuffer);
            }
        }
    }

private:
    template<bool WORKROW, bool SPARSE>
    struct CastBlockWorkspace : public BlockWorkspace<WORKROW, SPARSE>, public ConditionalBase<SPARSE> {
        CastBlockWorkspace(size_t s, size_t l, const WorkspaceOptions& opt, std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > p) : 
            BlockWorkspace<WORKROW, SPARSE>(s, l), ConditionalBase<SPARSE>(l, opt), internal(std::move(p)) {}

        std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > internal;
    };

    typedef CastBlockWorkspace<true, false> DenseCastRowBlockWorkspace;
    typedef CastBlockWorkspace<false, false> DenseCastColumnBlockWorkspace;
    typedef CastBlockWorkspace<true, true> SparseCastRowBlockWorkspace;
    typedef CastBlockWorkspace<false, true> SparseCastColumnBlockWorkspace;

public:
    std::shared_ptr<DenseRowBlockWorkspace> dense_row_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const { 
        return std::shared_ptr<DenseRowBlockWorkspace>(new DenseCastRowBlockWorkspace(start, length, opt, ptr->dense_row_workspace(start, length, opt)));
    }

    std::shared_ptr<DenseColumnBlockWorkspace> dense_column_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const { 
        return std::shared_ptr<DenseColumnBlockWorkspace>(new DenseCastColumnBlockWorkspace(start, length, opt, ptr->dense_column_workspace(start, length, opt)));
    }

    const T_out* row(size_t r, T_out* buffer, DenseRowBlockWorkspace* work) const {
        auto wptr = static_cast<DenseCastRowBlockWorkspace*>(work);
        return cast_dense<true>(r, buffer, wptr->length, wptr);
    }

    const T_out* column(size_t c, T_out* buffer, DenseColumnBlockWorkspace* work) const {
        auto wptr = static_cast<DenseCastColumnBlockWorkspace*>(work);
        return cast_dense<false>(c, buffer, wptr->length, wptr);
    }

    std::shared_ptr<SparseRowBlockWorkspace> sparse_row_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const { 
        return std::shared_ptr<SparseRowBlockWorkspace>(new SparseCastRowBlockWorkspace(start, length, opt, ptr->sparse_row_workspace(start, length, opt)));
    }

    std::shared_ptr<SparseColumnBlockWorkspace> sparse_column_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const { 
        return std::shared_ptr<SparseColumnBlockWorkspace>(new SparseCastColumnBlockWorkspace(start, length, opt, ptr->sparse_column_workspace(start, length, opt)));
    }

    SparseRange<T_out, IDX_out> row(size_t r, T_out* vbuffer, IDX_out* ibuffer, SparseRowBlockWorkspace* work) const {
        auto wptr = static_cast<SparseCastRowBlockWorkspace*>(work);
        return cast_sparse<true>(r, vbuffer, ibuffer, wptr);
    }

    SparseRange<T_out, IDX_out> column(size_t c, T_out* vbuffer, IDX_out* ibuffer, SparseColumnBlockWorkspace* work) const {
        auto wptr = static_cast<SparseCastColumnBlockWorkspace*>(work);
        return cast_sparse<false>(c, vbuffer, ibuffer, wptr);
    }

private:
    template<bool WORKROW, bool SPARSE>
    struct CastIndexWorkspace : public IndexWorkspace<IDX_out, WORKROW, SPARSE>, public ConditionalBase<SPARSE> {
        CastIndexWorkspace(size_t n, const WorkspaceOptions& opt) : 
            IndexWorkspace<IDX_out, WORKROW, SPARSE>(n),
            ConditionalBase<SPARSE>(n, opt)
            {}

        std::vector<IDX_out> indices_;
        const std::vector<IDX_out>& indices() const { 
            if constexpr(same_IDX_type) {
                return internal->indices();
            } else {
                return indices_; 
            }
        }
        
        std::shared_ptr<IndexWorkspace<IDX_in, WORKROW, SPARSE> > internal;
    };

    typedef CastIndexWorkspace<true, false> DenseCastRowIndexWorkspace;
    typedef CastIndexWorkspace<false, false> DenseCastColumnIndexWorkspace;
    typedef CastIndexWorkspace<true, true> SparseCastRowIndexWorkspace;
    typedef CastIndexWorkspace<false, true> SparseCastColumnIndexWorkspace;

public:
    std::shared_ptr<DenseRowIndexWorkspace<IDX_out> > dense_row_workspace(std::vector<IDX_out> indices, const WorkspaceOptions& opt) const { 
        return create_new_workspace<true, false>(std::move(indices), opt);
    }

    std::shared_ptr<DenseColumnIndexWorkspace<IDX_out> > dense_column_workspace(std::vector<IDX_out> indices, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, false>(std::move(indices), opt);
    }

    const T_out* row(size_t r, T_out* buffer, DenseRowIndexWorkspace<IDX_out>* work) const {
        auto wptr = static_cast<DenseCastRowIndexWorkspace*>(work);
        return cast_dense<true>(r, buffer, wptr->length, wptr);
    }

    const T_out* column(size_t c, T_out* buffer, DenseColumnIndexWorkspace<IDX_out>* work) const {
        auto wptr = static_cast<DenseCastColumnIndexWorkspace*>(work);
        return cast_dense<false>(c, buffer, wptr->length, wptr);
    }

    std::shared_ptr<SparseRowIndexWorkspace<IDX_out> > sparse_row_workspace(std::vector<IDX_out> indices, const WorkspaceOptions& opt) const { 
        return create_new_workspace<true, true>(std::move(indices), opt);
    }

    std::shared_ptr<SparseColumnIndexWorkspace<IDX_out> > sparse_column_workspace(std::vector<IDX_out> indices, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, true>(std::move(indices), opt);
    }

    SparseRange<T_out, IDX_out> row(size_t r, T_out* vbuffer, IDX_out* ibuffer, SparseRowIndexWorkspace<IDX_out>* work) const {
        auto wptr = static_cast<SparseCastRowIndexWorkspace*>(work);
        return cast_sparse<true>(r, vbuffer, ibuffer, wptr);
    }

    SparseRange<T_out, IDX_out> column(size_t c, T_out* vbuffer, IDX_out* ibuffer, SparseColumnIndexWorkspace<IDX_out>* work) const {
        auto wptr = static_cast<SparseCastColumnIndexWorkspace*>(work);
        return cast_sparse<false>(c, vbuffer, ibuffer, wptr);
    }

private:
    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<IndexWorkspace<IDX_out, WORKROW, SPARSE> > create_new_workspace(std::vector<IDX_out> indices, const WorkspaceOptions& opt) const { 
        auto wptr = new CastIndexWorkspace<WORKROW, SPARSE>(indices.size(), opt);
        auto output = std::shared_ptr<IndexWorkspace<IDX_out, WORKROW, SPARSE> >(wptr);

        if constexpr(!same_IDX_type) {
            wptr->indices_ = std::move(indices);
            wptr->internal = new_workspace<WORKROW, SPARSE>(ptr.get(), std::vector<IDX_in>(wptr->indices_.begin(), wptr->indices_.end()), opt);
        } else {
            wptr->internal = new_workspace<WORKROW, SPARSE>(ptr.get(), std::move(indices), opt);
        }

        return output;
    }
};

/**
 * Recast a `Matrix` to a different interface type.
 *
 * @tparam T_out Data type to cast to.
 * @tparam IDX_out Index type to cast to.
 * @tparam T_in Data type to cast from.
 * @tparam IDX_in Index type to cast from.
 *
 * @param p Pointer to the `Matrix` instance to cast from.
 * @return Pointer to a `Matrix` instance of the desired interface type.
 */
template<typename T_out, typename IDX_out, typename T_in, typename IDX_in>
std::shared_ptr<Matrix<T_out, IDX_out> > make_DelayedCast(const std::shared_ptr<Matrix<T_in, IDX_in> > p) {
    return std::shared_ptr<Matrix<T_out, IDX_out> >(new DelayedCast<T_out, IDX_out, T_in, IDX_in>(std::move(p)));
}

}

#endif
