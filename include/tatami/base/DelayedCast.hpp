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

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct CastWorkspace : Workspace<ROW> {
        CastWorkspace(size_t nv, std::shared_ptr<Workspace<ROW> > p) : vbuffer(nv), internal(std::move(p)) {}
        std::vector<T_in> vbuffer;
        std::vector<IDX_in> ibuffer;
        std::shared_ptr<Workspace<ROW> > internal;
    };

    typedef CastWorkspace<true> CastRowWorkspace;
    typedef CastWorkspace<false> CastColumnWorkspace;
    /**
     * @endcond
     */

    std::shared_ptr<RowWorkspace> new_row_workspace() const { 
        size_t nv = 0;
        if constexpr(!std::is_same<T_in, T_out>::value) {
            nv = this->ncol();
        }
        return std::shared_ptr<RowWorkspace>(new CastRowWorkspace(nv, ptr->new_row_workspace()));
    }

    std::shared_ptr<ColumnWorkspace> new_column_workspace() const { 
        size_t nv = 0;
        if constexpr(!std::is_same<T_in, T_out>::value) {
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
    template<bool ROW, class InputWorkspace>
    const T_out* cast_dense(size_t i, T_out* buffer, size_t length, InputWorkspace* work) const {
        if constexpr(std::is_same<T_in, T_out>::value) {
            if constexpr(ROW) {
                return ptr->row(i, buffer, work->internal.get()); 
            } else {
                return ptr->column(i, buffer, work->internal.get()); 
            }
        } else {
            const T_in* out;
            if constexpr(ROW) {
                out = ptr->row(i, work->vbuffer.data(), work->internal.get());
            } else {
                out = ptr->column(i, work->vbuffer.data(), work->internal.get());
            }
            std::copy(out, out + length, buffer);
            return buffer;
        }
    }

    template<bool ROW, class InputWorkspace>
    SparseRange<T_out, IDX_out> cast_sparse(size_t i, T_out* vbuffer, IDX_out* ibuffer, size_t length, InputWorkspace* work) const {
        if constexpr(std::is_same<T_in, T_out>::value) {
            if constexpr(std::is_same<IDX_in, IDX_out>::value) {
                if constexpr(ROW) {
                    return ptr->sparse_row(i, vbuffer, ibuffer, work->internal.get());
                } else {
                    return ptr->sparse_column(i, vbuffer, ibuffer, work->internal.get());
                }
            } else {
                if (work->ibuffer.empty()) {
                    work->ibuffer.resize(length);
                }
                SparseRange<T_in, IDX_in> out;
                if constexpr(ROW) {
                    out = ptr->sparse_row(i, vbuffer, work->ibuffer.data(), work->internal.get());
                } else {
                    out = ptr->sparse_column(i, vbuffer, work->ibuffer.data(), work->internal.get());
                }
                std::copy(out.index, out.index + out.number, ibuffer);
                return SparseRange<T_out, IDX_out>(out.number, out.value, ibuffer);
            }
        } else {
            if constexpr(std::is_same<IDX_in, IDX_out>::value) {
                SparseRange<T_in, IDX_in> out;
                if constexpr(ROW) {
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
                if constexpr(ROW) {
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
    template<bool ROW>
    struct CastBlockWorkspace : BlockWorkspace<ROW> {
        CastBlockWorkspace(size_t start, size_t length, size_t nv, std::shared_ptr<BlockWorkspace<ROW> > p) : 
            BlockWorkspace<ROW>(start, length), vbuffer(nv), internal(std::move(p)) {}

        std::vector<T_in> vbuffer;
        std::vector<IDX_in> ibuffer;
        std::shared_ptr<BlockWorkspace<ROW> > internal;
    };

    typedef CastBlockWorkspace<true> CastRowBlockWorkspace;
    typedef CastBlockWorkspace<false> CastColumnBlockWorkspace;
    /**
     * @endcond
     */

    std::shared_ptr<RowBlockWorkspace> new_row_workspace(size_t start, size_t length) const { 
        size_t nv = 0;
        if constexpr(!std::is_same<T_in, T_out>::value) {
            nv = length;
        }
        return std::shared_ptr<RowBlockWorkspace>(new CastRowBlockWorkspace(start, length, nv, ptr->new_row_workspace(start, length)));
    }

    std::shared_ptr<ColumnBlockWorkspace> new_column_workspace(size_t start, size_t length) const { 
        size_t nv = 0;
        if constexpr(!std::is_same<T_in, T_out>::value) {
            nv = length;
        }
        return std::shared_ptr<ColumnBlockWorkspace>(new CastColumnBlockWorkspace(start, length, nv, ptr->new_column_workspace(start, length)));
    }

    const T_out* row(size_t r, T_out* buffer, RowBlockWorkspace* work) const {
        return cast_dense<true>(r, buffer, work->length, static_cast<CastRowBlockWorkspace*>(work));
    }

    const T_out* column(size_t c, T_out* buffer, ColumnBlockWorkspace* work) const {
        return cast_dense<false>(c, buffer, work->length, static_cast<CastColumnBlockWorkspace*>(work));
    }

    SparseRange<T_out, IDX_out> sparse_row(size_t r, T_out* vbuffer, IDX_out* ibuffer, RowBlockWorkspace* work, bool sorted=true) const {
        return cast_sparse<true>(r, vbuffer, ibuffer, work->length, static_cast<CastRowBlockWorkspace*>(work));
    }

    SparseRange<T_out, IDX_out> sparse_column(size_t c, T_out* vbuffer, IDX_out* ibuffer, ColumnBlockWorkspace* work, bool sorted=true) const {
        return cast_sparse<false>(c, vbuffer, ibuffer, work->length, static_cast<CastColumnBlockWorkspace*>(work));
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct CastIndexWorkspace : IndexWorkspace<IDX_out, ROW> {
        CastIndexWorkspace(size_t length, const IDX_out* indices, size_t nv) : IndexWorkspace<IDX_out, ROW>(length, indices), vbuffer(nv) {
            if constexpr(!std::is_same<IDX_in, IDX_out>::value) {
                more_indices.resize(length);
                std::copy(indices, indices + length, more_indices.begin());
            }
        }

        const IDX_in* host_indices() const {
            if constexpr(!std::is_same<IDX_in, IDX_out>::value) {
                return more_indices.data();
            } else {
                return this->indices;
            }
        }

        std::vector<IDX_in> more_indices;
        std::vector<T_in> vbuffer;
        std::vector<IDX_in> ibuffer;
        std::shared_ptr<IndexWorkspace<IDX_in, ROW> > internal;
    };

    typedef CastIndexWorkspace<true> CastRowIndexWorkspace;
    typedef CastIndexWorkspace<false> CastColumnIndexWorkspace;
    /**
     * @endcond
     */

    std::shared_ptr<RowIndexWorkspace<IDX_out> > new_row_workspace(size_t length, const IDX_out* indices) const { 
        size_t nv = 0;
        if constexpr(!std::is_same<T_in, T_out>::value) {
            nv = length;
        }

        auto wptr = new CastRowIndexWorkspace(length, indices, nv);
        auto output = std::shared_ptr<RowIndexWorkspace<IDX_out> >(wptr);
        wptr->internal = ptr->new_row_workspace(length, wptr->host_indices());
        return output;
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX_out> > new_column_workspace(size_t length, const IDX_out* indices) const { 
        size_t nv = 0;
        if constexpr(!std::is_same<T_in, T_out>::value) {
            nv = length;
        }

        auto wptr = new CastColumnIndexWorkspace(length, indices, nv);
        auto output = std::shared_ptr<ColumnIndexWorkspace<IDX_out> >(wptr);
        wptr->internal = ptr->new_column_workspace(length, wptr->host_indices());
        return output;
    }

    const T_out* row(size_t r, T_out* buffer, RowIndexWorkspace<IDX_out>* work) const {
        return cast_dense<true>(r, buffer, work->length, static_cast<CastRowIndexWorkspace*>(work));
    }

    const T_out* column(size_t c, T_out* buffer, ColumnIndexWorkspace<IDX_out>* work) const {
        return cast_dense<false>(c, buffer, work->length, static_cast<CastColumnIndexWorkspace*>(work));
    }

    SparseRange<T_out, IDX_out> sparse_row(size_t r, T_out* vbuffer, IDX_out* ibuffer, RowIndexWorkspace<IDX_out>* work, bool sorted=true) const {
        return cast_sparse<true>(r, vbuffer, ibuffer, work->length, static_cast<CastRowIndexWorkspace*>(work));
    }

    SparseRange<T_out, IDX_out> sparse_column(size_t c, T_out* vbuffer, IDX_out* ibuffer, ColumnIndexWorkspace<IDX_out>* work, bool sorted=true) const {
        return cast_sparse<false>(c, vbuffer, ibuffer, work->length, static_cast<CastColumnIndexWorkspace*>(work));
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
