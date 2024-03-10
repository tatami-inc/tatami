#ifndef TATAMI_DELAYED_CAST_HPP
#define TATAMI_DELAYED_CAST_HPP

#include "../base/Matrix.hpp"

#include <memory>
#include <type_traits>
#include <algorithm>

/**
 * @file DelayedCast.hpp
 *
 * @brief Delayed cast to another interface type.
 */

namespace tatami {

/**
 * @cond
 */
namespace DelayedCast_internal {

template<typename ValueOut_, typename ValueIn_>
struct Dense {
    Dense(size_t n) : in_buffer(std::is_same<ValueOut_, ValueIn_>::value ? 0 : n) {}

protected:
    template<class Fetch_>
    const ValueOut_* fetch_base(ValueOut_* buffer, Fetch_ f) {
        if constexpr(std::is_same<ValueOut_, ValueIn_>::value) {
            return f(buffer);
        } else {
            auto ptr = f(in_buffer.data());
            std::copy_n(ptr, in_buffer.size(), buffer);
            return buffer;
        }
    }

    std::vector<ValueIn_> in_buffer;
};

template<typename ValueOut_, typename IndexOut_, typename ValueIn_, typename IndexIn_>
struct MyopicDense : public MyopicDenseExtractor<ValueOut_, IndexOut_>, public Dense<ValueOut_, ValueIn_> {
    MyopicDense(std::unique_ptr<MyopicDenseExtractor<ValueIn_, IndexIn_> > i) :
        Dense<ValueOut_, ValueIn_>(i->number()), internal(std::move(i)) {}

    const ValueOut_* fetch(IndexOut_ i, ValueOut_* buffer) {
        return this->fetch_base(buffer, [&](ValueIn_* b) -> auto { return internal->fetch(i, b); });
    }

    IndexOut_ number() const {
        return internal->number();
    }

private:
    std::unique_ptr<MyopicDenseExtractor<ValueIn_, IndexIn_> > internal;
};

template<typename ValueOut_, typename IndexOut_, typename ValueIn_, typename IndexIn_>
struct OracularDense : public OracularDenseExtractor<ValueOut_, IndexOut_>, public Dense<ValueOut_, ValueIn_> {
    OracularDense(std::unique_ptr<OracularDenseExtractor<ValueIn_, IndexIn_> > i) : 
        Dense<ValueOut_, ValueIn_>(i->number()), internal(std::move(i)) {}

    const ValueOut_* fetch(IndexOut_& i, ValueOut_* buffer) {
        IndexIn_ j;
        auto out = this->fetch_base(buffer, [&](ValueIn_* b) -> auto { return internal->fetch(j, b); });
        i = j;
        return out;
    }

    IndexOut_ number() const {
        return internal->number();
    }

private:
    std::unique_ptr<OracularDenseExtractor<ValueIn_, IndexIn_> > internal;
};

template<typename ValueOut_, typename IndexOut_, typename ValueIn_, typename IndexIn_>
struct Sparse {
    Sparse(size_t n, const Options& opt) : 
        in_vbuffer(std::is_same<ValueOut_, ValueIn_>::value || !opt.sparse_extract_value ? 0 : n),
        in_ibuffer(std::is_same<IndexOut_, IndexIn_>::value || !opt.sparse_extract_index ? 0 : n) 
    {}

    template<class Fetch_>
    SparseRange<ValueOut_, IndexOut_> fetch_base(ValueOut_* vbuffer, IndexOut_* ibuffer, Fetch_ f) {
        IndexIn_* iptr;
        if constexpr(std::is_same<IndexOut_, IndexIn_>::value) {
            iptr = ibuffer;
        } else {
            iptr = in_ibuffer.data();
        }

        ValueIn_* vptr;
        if constexpr(std::is_same<ValueOut_, ValueIn_>::value) {
            vptr = vbuffer;
        } else {
            vptr = in_vbuffer.data();
        }

        auto range = f(vptr, iptr);
        SparseRange<ValueOut_, IndexOut_> output(range.number);

        if constexpr(std::is_same<IndexOut_, IndexIn_>::value) {
            output.index = range.index;
        } else if (range.index != NULL) {
            std::copy_n(range.index, range.number, ibuffer);
            output.index = ibuffer;
        }

        if constexpr(std::is_same<ValueOut_, ValueIn_>::value) {
            output.value = range.value;
        } else if (range.value != NULL) {
            std::copy_n(range.value, range.number, vbuffer);
            output.value = vbuffer;
        }

        return output;
    }

    std::vector<ValueIn_> in_vbuffer;
    std::vector<IndexIn_> in_ibuffer;
};

template<typename ValueOut_, typename IndexOut_, typename ValueIn_, typename IndexIn_>
struct MyopicSparse : public MyopicSparseExtractor<ValueOut_, IndexOut_>, public Sparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_> {
    MyopicSparse(std::unique_ptr<MyopicSparseExtractor<ValueIn_, IndexIn_> > i, const Options& opt) : 
        Sparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_>(i->number(), opt), internal(std::move(i)) {}

    SparseRange<ValueOut_, IndexOut_> fetch(IndexOut_ i, ValueOut_* vbuffer, IndexOut_* ibuffer) {
        return this->fetch_base(vbuffer, ibuffer, [&](ValueIn_* vb, IndexIn_* ib) -> auto { return internal->fetch(i, vb, ib); });
    }

    IndexOut_ number() const {
        return internal->number();
    }

private:
    std::unique_ptr<MyopicSparseExtractor<ValueIn_, IndexIn_> > internal;
};

template<typename ValueOut_, typename IndexOut_, typename ValueIn_, typename IndexIn_>
struct OracularSparse : public OracularSparseExtractor<ValueOut_, IndexOut_>, public Sparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_> {
    OracularSparse(std::unique_ptr<OracularSparseExtractor<ValueIn_, IndexIn_> > i, const Options& opt) : 
        Sparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_>(i->number(), opt), internal(std::move(i)) {}

    SparseRange<ValueOut_, IndexOut_> fetch(IndexOut_& i, ValueOut_* vbuffer, IndexOut_* ibuffer) {
        IndexIn_ j;
        auto output = this->fetch_base(vbuffer, ibuffer, [&](ValueIn_* vb, IndexIn_* ib) -> auto { return internal->fetch(j, vb, ib); });
        i = j;
        return output;
    }

    IndexOut_ number() const {
        return internal->number();
    }

private:
    std::unique_ptr<OracularSparseExtractor<ValueIn_, IndexIn_> > internal;
};

template<typename IndexIn_, typename IndexOut_>
struct CastOracle : public Oracle<IndexIn_> {
    CastOracle(std::shared_ptr<Oracle<IndexOut_> > o) : oracle(std::move(o)) {}

    IndexIn_ get(size_t i) const {
        return oracle->get(i);
    }

    size_t total() const {
        return oracle->total();
    }

    std::shared_ptr<Oracle<IndexOut_> > oracle;
};

template<typename IndexIn_, typename IndexOut_>
std::shared_ptr<Oracle<IndexIn_> > convert(std::shared_ptr<Oracle<IndexOut_> > o) {
    if constexpr(std::is_same<IndexIn_, IndexOut_>::value) {
        return o;
    } else {
        return std::make_shared<CastOracle<IndexIn_, IndexOut_> >(std::move(o));
    }
}

template<typename IndexIn_, typename IndexOut_>
std::vector<IndexIn_> convert(std::vector<IndexOut_> i) {
    if constexpr(std::is_same<IndexIn_, IndexOut_>::value) {
        return i;
    } else {
        return std::vector<IndexIn_>(i.begin(), i.end());
    }
}

}
/**
 * @endcond
 */

/**
 * @brief Recast a `Matrix` to a different interface type.
 *
 * This performs a delayed cast from one interface type to another.
 * It is useful as a compatibility layer between functions that require `Matrix` objects of different types.
 * Casting is achieved by extracting the requested row/column from the input `Matrix` and transforming it to the output types.
 * Note that this is only done per row/column - the entirety of the original matrix is not copied to the new type.
 *
 * @tparam ValueOut_ Data type to cast to.
 * @tparam IndexOut_ Index type to cast to.
 * @tparam ValueIn_ Data type to cast from.
 * @tparam IndexIn_ Index type to cast from.
 */
template<typename ValueOut_, typename IndexOut_, typename ValueIn_, typename IndexIn_>
class DelayedCast : public Matrix<ValueOut_, IndexOut_> {
public:
    /**
     * @param p Pointer to the `Matrix` instance to cast from.
     */
    DelayedCast(std::shared_ptr<const Matrix<ValueIn_, IndexIn_> > p) : ptr(std::move(p)) {}

public:
    IndexOut_ nrow() const {
        return ptr->nrow();
    }

    IndexOut_ ncol() const {
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
    std::shared_ptr<const Matrix<ValueIn_, IndexIn_> > ptr;

    /********************
     *** Myopic dense ***
     ********************/
public:
    std::unique_ptr<MyopicDenseExtractor<ValueOut_, IndexOut_> > dense_row(const Options& opt) const {
        return std::make_unique<DelayedCast_internal::MyopicDense<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(ptr->dense_row(opt));
    }

    std::unique_ptr<MyopicDenseExtractor<ValueOut_, IndexOut_> > dense_row(IndexOut_ block_start, IndexOut_ block_length, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::MyopicDense<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->dense_row(block_start, block_length, opt)
        );
    }

    std::unique_ptr<MyopicDenseExtractor<ValueOut_, IndexOut_> > dense_row(std::vector<IndexOut_> indices, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::MyopicDense<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->dense_row(DelayedCast_internal::convert<IndexIn_>(std::move(indices)), opt)
        );
    }

    std::unique_ptr<MyopicDenseExtractor<ValueOut_, IndexOut_> > dense_column(const Options& opt) const {
        return std::make_unique<DelayedCast_internal::MyopicDense<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(ptr->dense_column(opt));
    }

    std::unique_ptr<MyopicDenseExtractor<ValueOut_, IndexOut_> > dense_column(IndexOut_ block_start, IndexOut_ block_length, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::MyopicDense<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->dense_column(block_start, block_length, opt)
        );
    }

    std::unique_ptr<MyopicDenseExtractor<ValueOut_, IndexOut_> > dense_column(std::vector<IndexOut_> indices, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::MyopicDense<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->dense_column(DelayedCast_internal::convert<IndexIn_>(std::move(indices)), opt)
        );
    }

    /*********************
     *** Myopic sparse ***
     *********************/
public:
    std::unique_ptr<MyopicSparseExtractor<ValueOut_, IndexOut_> > sparse_row(const Options& opt) const {
        return std::make_unique<DelayedCast_internal::MyopicSparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->sparse_row(opt),
            opt
        );
    }

    std::unique_ptr<MyopicSparseExtractor<ValueOut_, IndexOut_> > sparse_row(IndexOut_ block_start, IndexOut_ block_length, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::MyopicSparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->sparse_row(block_start, block_length, opt),
            opt
        );
    }

    std::unique_ptr<MyopicSparseExtractor<ValueOut_, IndexOut_> > sparse_row(std::vector<IndexOut_> indices, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::MyopicSparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->sparse_row(DelayedCast_internal::convert<IndexIn_>(std::move(indices)), opt), 
            opt
        );
    }

    std::unique_ptr<MyopicSparseExtractor<ValueOut_, IndexOut_> > sparse_column(const Options& opt) const {
        return std::make_unique<DelayedCast_internal::MyopicSparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->sparse_column(opt),
            opt
        );
    }

    std::unique_ptr<MyopicSparseExtractor<ValueOut_, IndexOut_> > sparse_column(IndexOut_ block_start, IndexOut_ block_length, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::MyopicSparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->sparse_column(block_start, block_length, opt),
            opt
        );
    }

    std::unique_ptr<MyopicSparseExtractor<ValueOut_, IndexOut_> > sparse_column(std::vector<IndexOut_> indices, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::MyopicSparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->sparse_column(DelayedCast_internal::convert<IndexIn_>(std::move(indices)), opt), 
            opt
        );
    }

    /********************
     *** Oracular dense ***
     ********************/
public:
    std::unique_ptr<OracularDenseExtractor<ValueOut_, IndexOut_> > dense_row(std::shared_ptr<Oracle<IndexOut_> > oracle, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::OracularDense<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->dense_row(DelayedCast_internal::convert<IndexIn_>(std::move(oracle)), opt)
        );
    }

    std::unique_ptr<OracularDenseExtractor<ValueOut_, IndexOut_> > dense_row(std::shared_ptr<Oracle<IndexOut_> > oracle, IndexOut_ block_start, IndexOut_ block_length, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::OracularDense<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->dense_row(DelayedCast_internal::convert<IndexIn_>(std::move(oracle)), block_start, block_length, opt)
        );
    }

    std::unique_ptr<OracularDenseExtractor<ValueOut_, IndexOut_> > dense_row(std::shared_ptr<Oracle<IndexOut_> > oracle, std::vector<IndexOut_> indices, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::OracularDense<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->dense_row(
                DelayedCast_internal::convert<IndexIn_>(std::move(oracle)), 
                DelayedCast_internal::convert<IndexIn_>(std::move(indices)), 
                opt
            )
        );
    }

    std::unique_ptr<OracularDenseExtractor<ValueOut_, IndexOut_> > dense_column(std::shared_ptr<Oracle<IndexOut_> > oracle, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::OracularDense<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->dense_column(DelayedCast_internal::convert<IndexIn_>(std::move(oracle)), opt)
        );
    }

    std::unique_ptr<OracularDenseExtractor<ValueOut_, IndexOut_> > dense_column(std::shared_ptr<Oracle<IndexOut_> > oracle, IndexOut_ block_start, IndexOut_ block_length, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::OracularDense<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->dense_column(DelayedCast_internal::convert<IndexIn_>(std::move(oracle)), block_start, block_length, opt)
        );
    }

    std::unique_ptr<OracularDenseExtractor<ValueOut_, IndexOut_> > dense_column(std::shared_ptr<Oracle<IndexOut_> > oracle, std::vector<IndexOut_> indices, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::OracularDense<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->dense_column(
                DelayedCast_internal::convert<IndexIn_>(std::move(oracle)), 
                DelayedCast_internal::convert<IndexIn_>(std::move(indices)), 
                opt
            )
        );
    }

    /*********************
     *** Oracular sparse ***
     *********************/
public:
    std::unique_ptr<OracularSparseExtractor<ValueOut_, IndexOut_> > sparse_row(std::shared_ptr<Oracle<IndexOut_> > oracle, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::OracularSparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->sparse_row(DelayedCast_internal::convert<IndexIn_>(std::move(oracle)), opt), opt
        );
    }

    std::unique_ptr<OracularSparseExtractor<ValueOut_, IndexOut_> > sparse_row(std::shared_ptr<Oracle<IndexOut_> > oracle, IndexOut_ block_start, IndexOut_ block_length, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::OracularSparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->sparse_row(DelayedCast_internal::convert<IndexIn_>(std::move(oracle)), block_start, block_length, opt),
            opt
        );
    }

    std::unique_ptr<OracularSparseExtractor<ValueOut_, IndexOut_> > sparse_row(std::shared_ptr<Oracle<IndexOut_> > oracle, std::vector<IndexOut_> indices, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::OracularSparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->sparse_row(
                DelayedCast_internal::convert<IndexIn_>(std::move(oracle)), 
                DelayedCast_internal::convert<IndexIn_>(std::move(indices)), 
                opt
            ),
            opt
        );
    }

    std::unique_ptr<OracularSparseExtractor<ValueOut_, IndexOut_> > sparse_column(std::shared_ptr<Oracle<IndexOut_> > oracle, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::OracularSparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->sparse_column(DelayedCast_internal::convert<IndexIn_>(std::move(oracle)), opt),
            opt
        );
    }

    std::unique_ptr<OracularSparseExtractor<ValueOut_, IndexOut_> > sparse_column(std::shared_ptr<Oracle<IndexOut_> > oracle, IndexOut_ block_start, IndexOut_ block_length, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::OracularSparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->sparse_column(DelayedCast_internal::convert<IndexIn_>(std::move(oracle)), block_start, block_length, opt),
            opt
        );
    }

    std::unique_ptr<OracularSparseExtractor<ValueOut_, IndexOut_> > sparse_column(std::shared_ptr<Oracle<IndexOut_> > oracle, std::vector<IndexOut_> indices, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::OracularSparse<ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(
            ptr->sparse_column(
                DelayedCast_internal::convert<IndexIn_>(std::move(oracle)), 
                DelayedCast_internal::convert<IndexIn_>(std::move(indices)), 
                opt
            ), 
            opt
        );
    }
};

/**
 * Recast a `Matrix` to a different interface type.
 *
 * @tparam ValueOut_ Data type to cast to.
 * @tparam IndexOut_ Index type to cast to.
 * @tparam ValueIn_ Data type to cast from.
 * @tparam IndexIn_ Index type to cast from.
 *
 * @param p Pointer to the (possbly `const`) `Matrix` instance to cast from.
 * @return Pointer to a `Matrix` instance of the desired interface type.
 */
template<typename ValueOut_, typename IndexOut_, typename ValueIn_, typename IndexIn_>
std::shared_ptr<Matrix<ValueOut_, IndexOut_> > make_DelayedCast(std::shared_ptr<const Matrix<ValueIn_, IndexIn_> > p) {
    return std::shared_ptr<Matrix<ValueOut_, IndexOut_> >(new DelayedCast<ValueOut_, IndexOut_, ValueIn_, IndexIn_>(std::move(p)));
}

/**
 * @cond
 */
template<typename ValueOut_, typename IndexOut_, typename ValueIn_, typename IndexIn_>
std::shared_ptr<Matrix<ValueOut_, IndexOut_> > make_DelayedCast(std::shared_ptr<Matrix<ValueIn_, IndexIn_> > p) {
    return std::shared_ptr<Matrix<ValueOut_, IndexOut_> >(new DelayedCast<ValueOut_, IndexOut_, ValueIn_, IndexIn_>(std::move(p)));
}
/**
 * @endcond
 */

}

#endif
