#ifndef TATAMI_DELAYED_CAST_HPP
#define TATAMI_DELAYED_CAST_HPP

#include "Matrix.hpp"
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

public:
    /**
     * @cond
     */
    template<bool WORKROW>
    struct CastWorkspace : Workspace<WORKROW> {
        CastWorkspace(size_t nv, std::shared_ptr<Workspace<WORKROW> > p) : vbuffer(nv), internal(std::move(p)) {}
        std::vector<T_in> vbuffer;
        std::vector<IDX_in> ibuffer;
        std::shared_ptr<Workspace<WORKROW> > internal;
    };

    typedef CastWorkspace<true> CastRowWorkspace;
    typedef CastWorkspace<false> CastColumnWorkspace;
    /**
     * @endcond
     */

    std::shared_ptr<RowWorkspace> new_row_workspace() const { 
        size_t nv = 0;
        if constexpr(!same_T_type) {
            nv = this->ncol();
        }
        return std::shared_ptr<RowWorkspace>(new CastRowWorkspace(nv, ptr->new_row_workspace()));
    }

    std::shared_ptr<ColumnWorkspace> new_column_workspace() const { 
        size_t nv = 0;
        if constexpr(!same_T_type) {
            nv = this->nrow();
        }
        return std::shared_ptr<ColumnWorkspace>(new CastColumnWorkspace(nv, ptr->new_column_workspace()));
    }

    const T_out* row(size_t r, T_out* buffer, RowWorkspace* work) const {
        return cast_dense<true>(r, buffer, ptr->ncol(), static_cast<CastRowWorkspace*>(work));
    }

    const T_out* column(size_t c, T_out* buffer, ColumnWorkspace* work) const {
        return cast_dense<false>(c, buffer, ptr->nrow(), static_cast<CastColumnWorkspace*>(work));
    }

    SparseRange<T_out, IDX_out> sparse_row(size_t r, T_out* vbuffer, IDX_out* ibuffer, RowWorkspace* work, bool sorted=true) const {
        return cast_sparse<true>(r, vbuffer, ibuffer, ptr->ncol(), static_cast<CastRowWorkspace*>(work));
    }

    SparseRange<T_out, IDX_out> sparse_column(size_t c, T_out* vbuffer, IDX_out* ibuffer, ColumnWorkspace* work, bool sorted=true) const {
        return cast_sparse<false>(c, vbuffer, ibuffer, ptr->nrow(), static_cast<CastColumnWorkspace*>(work));
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
                out = ptr->row(i, work->vbuffer.data(), work->internal.get());
            } else {
                out = ptr->column(i, work->vbuffer.data(), work->internal.get());
            }
            std::copy(out, out + length, buffer);
            return buffer;
        }
    }

    template<bool WORKROW, class InputWorkspace>
    SparseRange<T_out, IDX_out> cast_sparse(size_t i, T_out* vbuffer, IDX_out* ibuffer, size_t length, InputWorkspace* work) const {
        if constexpr(same_T_type) {
            if constexpr(same_IDX_type) {
                if constexpr(WORKROW) {
                    return ptr->sparse_row(i, vbuffer, ibuffer, work->internal.get());
                } else {
                    return ptr->sparse_column(i, vbuffer, ibuffer, work->internal.get());
                }
            } else {
                if (work->ibuffer.empty()) {
                    work->ibuffer.resize(length);
                }
                SparseRange<T_in, IDX_in> out;
                if constexpr(WORKROW) {
                    out = ptr->sparse_row(i, vbuffer, work->ibuffer.data(), work->internal.get());
                } else {
                    out = ptr->sparse_column(i, vbuffer, work->ibuffer.data(), work->internal.get());
                }
                std::copy(out.index, out.index + out.number, ibuffer);
                return SparseRange<T_out, IDX_out>(out.number, out.value, ibuffer);
            }
        } else {
            if constexpr(same_IDX_type) {
                SparseRange<T_in, IDX_in> out;
                if constexpr(WORKROW) {
                    out = ptr->sparse_row(i, work->vbuffer.data(), ibuffer, work->internal.get());
                } else {
                    out = ptr->sparse_column(i, work->vbuffer.data(), ibuffer, work->internal.get());
                }
                std::copy(out.value, out.value + out.number, vbuffer);
                return SparseRange<T_out, IDX_out>(out.number, vbuffer, out.index);
            } else {
                if (work->ibuffer.empty()) {
                    work->ibuffer.resize(length);
                }
                SparseRange<T_in, IDX_in> out;
                if constexpr(WORKROW) {
                    out = ptr->sparse_row(i, work->vbuffer.data(), work->ibuffer.data(), work->internal.get());
                } else {
                    out = ptr->sparse_column(i, work->vbuffer.data(), work->ibuffer.data(), work->internal.get());
                }
                std::copy(out.value, out.value + out.number, vbuffer);
                std::copy(out.index, out.index + out.number, ibuffer);
                return SparseRange<T_out, IDX_out>(out.number, vbuffer, ibuffer);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<bool WORKROW>
    struct CastBlockWorkspace : BlockWorkspace<WORKROW> {
        CastBlockWorkspace(size_t nv, std::shared_ptr<BlockWorkspace<WORKROW> > p) : vbuffer(nv), internal(std::move(p)) {}

        std::vector<T_in> vbuffer;
        std::vector<IDX_in> ibuffer;
        std::shared_ptr<BlockWorkspace<WORKROW> > internal;

        const std::pair<size_t, size_t>& block() const { return internal->block(); }
    };

    typedef CastBlockWorkspace<true> CastRowBlockWorkspace;
    typedef CastBlockWorkspace<false> CastColumnBlockWorkspace;
    /**
     * @endcond
     */

    std::shared_ptr<RowBlockWorkspace> new_row_workspace(size_t start, size_t length) const { 
        size_t nv = 0;
        if constexpr(!same_T_type) {
            nv = length;
        }
        return std::shared_ptr<RowBlockWorkspace>(new CastRowBlockWorkspace(nv, ptr->new_row_workspace(start, length)));
    }

    std::shared_ptr<ColumnBlockWorkspace> new_column_workspace(size_t start, size_t length) const { 
        size_t nv = 0;
        if constexpr(!same_T_type) {
            nv = length;
        }
        return std::shared_ptr<ColumnBlockWorkspace>(new CastColumnBlockWorkspace(nv, ptr->new_column_workspace(start, length)));
    }

    const T_out* row(size_t r, T_out* buffer, RowBlockWorkspace* work) const {
        auto wptr = static_cast<CastRowBlockWorkspace*>(work);
        return cast_dense<true>(r, buffer, wptr->internal->length(), wptr);
    }

    const T_out* column(size_t c, T_out* buffer, ColumnBlockWorkspace* work) const {
        auto wptr = static_cast<CastColumnBlockWorkspace*>(work);
        return cast_dense<false>(c, buffer, wptr->internal->length(), wptr);
    }

    SparseRange<T_out, IDX_out> sparse_row(size_t r, T_out* vbuffer, IDX_out* ibuffer, RowBlockWorkspace* work, bool sorted=true) const {
        auto wptr = static_cast<CastRowBlockWorkspace*>(work);
        return cast_sparse<true>(r, vbuffer, ibuffer, wptr->internal->length(), wptr);
    }

    SparseRange<T_out, IDX_out> sparse_column(size_t c, T_out* vbuffer, IDX_out* ibuffer, ColumnBlockWorkspace* work, bool sorted=true) const {
        auto wptr = static_cast<CastColumnBlockWorkspace*>(work);
        return cast_sparse<false>(c, vbuffer, ibuffer, wptr->internal->length(), wptr);
    }

public:
    /**
     * @cond
     */
    template<bool WORKROW>
    struct CastIndexWorkspace : IndexWorkspace<IDX_out, WORKROW> {
        CastIndexWorkspace(size_t nv) : vbuffer(nv) {}

        std::vector<IDX_out> indices_;
        const std::vector<IDX_out>& indices() const { 
            if constexpr(same_IDX_type) {
                return internal->indices();
            } else {
                return indices_; 
            }
        }
        
        std::vector<T_in> vbuffer;
        std::vector<IDX_in> ibuffer;

        std::shared_ptr<IndexWorkspace<IDX_in, WORKROW> > internal;
    };

    typedef CastIndexWorkspace<true> CastRowIndexWorkspace;
    typedef CastIndexWorkspace<false> CastColumnIndexWorkspace;
    /**
     * @endcond
     */

    std::shared_ptr<RowIndexWorkspace<IDX_out> > new_row_workspace(std::vector<IDX_out> indices) const { 
        return new_workspace<true>(std::move(indices));
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX_out> > new_column_workspace(std::vector<IDX_out> indices) const {
        return new_workspace<false>(std::move(indices));
    }

    const T_out* row(size_t r, T_out* buffer, RowIndexWorkspace<IDX_out>* work) const {
        auto wptr = static_cast<CastRowIndexWorkspace*>(work);
        return cast_dense<true>(r, buffer, wptr->internal->length(), wptr);
    }

    const T_out* column(size_t c, T_out* buffer, ColumnIndexWorkspace<IDX_out>* work) const {
        auto wptr = static_cast<CastColumnIndexWorkspace*>(work);
        return cast_dense<false>(c, buffer, wptr->internal->length(), wptr);
    }

    SparseRange<T_out, IDX_out> sparse_row(size_t r, T_out* vbuffer, IDX_out* ibuffer, RowIndexWorkspace<IDX_out>* work, bool sorted=true) const {
        auto wptr = static_cast<CastRowIndexWorkspace*>(work);
        return cast_sparse<true>(r, vbuffer, ibuffer, wptr->internal->length(), wptr);
    }

    SparseRange<T_out, IDX_out> sparse_column(size_t c, T_out* vbuffer, IDX_out* ibuffer, ColumnIndexWorkspace<IDX_out>* work, bool sorted=true) const {
        auto wptr = static_cast<CastColumnIndexWorkspace*>(work);
        return cast_sparse<false>(c, vbuffer, ibuffer, wptr->internal->length(), wptr);
    }

private:
    template<bool WORKROW>
    std::shared_ptr<IndexWorkspace<IDX_out, WORKROW> > new_workspace(std::vector<IDX_out> indices) const { 
        size_t nv = 0;
        if constexpr(!same_T_type) {
            nv = indices.size();
        }

        auto wptr = new CastIndexWorkspace<WORKROW>(nv);
        auto output = std::shared_ptr<IndexWorkspace<IDX_out, WORKROW> >(wptr);

        if constexpr(!same_IDX_type) {
            wptr->indices_ = std::move(indices);
            if constexpr(WORKROW) {
                wptr->internal = ptr->new_row_workspace(std::vector<IDX_in>(wptr->indices_.begin(), wptr->indices_.end()));
            } else {
                wptr->internal = ptr->new_column_workspace(std::vector<IDX_in>(wptr->indices_.begin(), wptr->indices_.end()));
            }
        } else {
            if constexpr(WORKROW) {
                wptr->internal = ptr->new_row_workspace(std::move(indices));
            } else {
                wptr->internal = ptr->new_column_workspace(std::move(indices));
            }
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
