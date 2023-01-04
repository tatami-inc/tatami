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
    struct CastWorkspace : public Workspace {
        CastWorkspace(size_t nv, size_t ni, std::shared_ptr<Workspace> p) : vbuffer(nv), ibuffer(ni), internal(std::move(p)) {}
        std::vector<T_in> vbuffer;
        std::vector<IDX_in> ibuffer;
        std::shared_ptr<Workspace> internal;
    };

public:
    std::shared_ptr<Workspace> new_workspace(bool row) const { 
        size_t nv = 0, ni = 0;
        if (row) {
            if constexpr(!std::is_same<T_in, T_out>::value) {
                nv = this->ncol();
            }
            if constexpr(!std::is_same<IDX_in, IDX_out>::value) {
                ni = this->ncol();
            }
        } else {
            if constexpr(!std::is_same<T_in, T_out>::value) {
                nv = this->nrow();
            }
            if constexpr(!std::is_same<IDX_in, IDX_out>::value) {
                ni = this->nrow();
            }
        }
        return std::shared_ptr<Workspace>(new CastWorkspace(nv, ni, ptr->new_workspace(row)));
    }

private:
    template<class Function>
    const T_out* cast_dense(T_out* buffer, size_t first, size_t last, Workspace* work, Function fun) const {
        if constexpr(std::is_same<T_in, T_out>::value) {
            if (work == nullptr) {
                return fun(buffer, nullptr);
            } else {
                auto wptr = static_cast<CastWorkspace*>(work);
                return fun(buffer, wptr->internal.get()); 
            }
        } else {
            size_t n = last - first;
            if (work == nullptr) {
                std::vector<T_in> vbuffer0(n);
                auto out = fun(vbuffer0.data(), nullptr);
                std::copy(out, out + n, buffer);
            } else {
                auto wptr = static_cast<CastWorkspace*>(work);
                auto out = fun(wptr->vbuffer.data(), wptr->internal.get());
                std::copy(out, out + n, buffer);
            }
            return buffer;
        }
    }

public:
    const T_out* row(size_t r, T_out* buffer, size_t first, size_t last, Workspace* work=nullptr) const {
        return cast_dense(buffer, first, last, work, [&](T_in* buffer0, Workspace* work0) -> const T_in* {
            return ptr->row(r, buffer0, first, last, work0);
        });
    }

    const T_out* column(size_t c, T_out* buffer, size_t first, size_t last, Workspace* work=nullptr) const {
        return cast_dense(buffer, first, last, work, [&](T_in* buffer0, Workspace* work0) -> const T_in* {
            return ptr->column(c, buffer0, first, last, work0);
        });
    }

private:
    template<class Function>
    SparseRange<T_out, IDX_out> cast_sparse(T_out* vbuffer, IDX_out* ibuffer, size_t first, size_t last, Workspace* work, Function fun) const {
        if constexpr(std::is_same<T_in, T_out>::value) {
            if constexpr(std::is_same<IDX_in, IDX_out>::value) {
                if (work == nullptr) {
                    return fun(vbuffer, ibuffer, nullptr);
                } else {
                    auto wptr = static_cast<CastWorkspace*>(work);
                    return fun(vbuffer, ibuffer, wptr->internal.get());
                }
            } else {
                if (work == nullptr) {
                    std::vector<IDX_in> ibuffer0(last - first);
                    auto out = fun(vbuffer, ibuffer0.data(), nullptr);
                    std::copy(out.index, out.index + out.number, ibuffer);
                    return SparseRange<T_out, IDX_out>(out.number, out.value, ibuffer);
                } else {
                    auto wptr = static_cast<CastWorkspace*>(work);
                    auto out = fun(vbuffer, wptr->ibuffer.data(), wptr->internal.get());
                    std::copy(out.index, out.index + out.number, ibuffer);
                    return SparseRange<T_out, IDX_out>(out.number, out.value, ibuffer);
                }
            }
        } else {
            if constexpr(std::is_same<IDX_in, IDX_out>::value) {
                if (work == nullptr) {
                    std::vector<T_in> vbuffer0(last - first);
                    auto out = fun(vbuffer0.data(), ibuffer, nullptr);
                    std::copy(out.value, out.value + out.number, vbuffer);
                    return SparseRange<T_out, IDX_out>(out.number, vbuffer, out.index);
                } else {
                    auto wptr = static_cast<CastWorkspace*>(work);
                    auto out = fun(wptr->vbuffer.data(), ibuffer, wptr->internal.get());
                    std::copy(out.value, out.value + out.number, vbuffer);
                    return SparseRange<T_out, IDX_out>(out.number, vbuffer, out.index);
                }
            } else {
                if (work == nullptr) {
                    std::vector<IDX_in> ibuffer0(last - first);
                    std::vector<T_in> vbuffer0(last - first);
                    auto out = fun(vbuffer0.data(), ibuffer0.data(), nullptr);
                    std::copy(out.value, out.value + out.number, vbuffer);
                    std::copy(out.index, out.index + out.number, ibuffer);
                    return SparseRange<T_out, IDX_out>(out.number, vbuffer, ibuffer);
                } else {
                    auto wptr = static_cast<CastWorkspace*>(work);
                    auto out = fun(wptr->vbuffer.data(), wptr->ibuffer.data(), wptr->internal.get());
                    std::copy(out.value, out.value + out.number, vbuffer);
                    std::copy(out.index, out.index + out.number, ibuffer);
                    return SparseRange<T_out, IDX_out>(out.number, vbuffer, ibuffer);
                }
            }
        }
    }

public:
    SparseRange<T_out, IDX_out> sparse_row(size_t r, T_out* vbuffer, IDX_out* ibuffer, size_t first, size_t last, Workspace* work=nullptr, bool sorted=true) const {
        return cast_sparse(vbuffer, ibuffer, first, last, work, 
            [&](T_in* vbuffer0, IDX_in* ibuffer0, Workspace* work0) -> SparseRange<T_in, IDX_in> {
                return ptr->sparse_row(r, vbuffer0, ibuffer0, first, last, work0, sorted);
            }
        );
    }

    SparseRange<T_out, IDX_out> sparse_column(size_t c, T_out* vbuffer, IDX_out* ibuffer, size_t first, size_t last, Workspace* work=nullptr, bool sorted=true) const {
        return cast_sparse(vbuffer, ibuffer, first, last, work, 
            [&](T_in* vbuffer0, IDX_in* ibuffer0, Workspace* work0) -> SparseRange<T_in, IDX_in> {
                return ptr->sparse_column(c, vbuffer0, ibuffer0, first, last, work0, sorted);
            }
        );
    }

private:
    std::shared_ptr<Matrix<T_in, IDX_in> > ptr;
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
