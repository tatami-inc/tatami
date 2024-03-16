#ifndef TATAMI_DELAYED_CAST_HPP
#define TATAMI_DELAYED_CAST_HPP

#include "../base/Matrix.hpp"
#include "../utils/new_extractor.hpp"

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

template<typename IndexIn_, typename IndexOut_>
struct CastOracle : public Oracle<IndexIn_> {
    CastOracle(std::shared_ptr<const Oracle<IndexOut_> > o) : oracle(std::move(o)) {}

    IndexIn_ get(size_t i) const {
        return oracle->get(i);
    }

    size_t total() const {
        return oracle->total();
    }

    std::shared_ptr<const Oracle<IndexOut_> > oracle;
};

template<bool oracle_, typename IndexIn_, typename IndexOut_>
MaybeOracle<oracle_, IndexIn_> convert(MaybeOracle<oracle_, IndexOut_> o) {
    if constexpr(!oracle_) {
        return false;
    } else if constexpr(std::is_same<IndexIn_, IndexOut_>::value) {
        return o;
    } else {
        return std::make_shared<CastOracle<IndexIn_, IndexOut_> >(std::move(o));
    }
}

template<typename IndexIn_, typename IndexOut_>
VectorPtr<IndexIn_> convert(VectorPtr<IndexOut_> i) {
    if constexpr(std::is_same<IndexIn_, IndexOut_>::value) {
        return i;
    } else {
        const auto& ix = *i;
        return std::make_shared<std::vector<IndexIn_> >(ix.begin(), ix.end());
    }
}

template<bool oracle_, typename ValueOut_, typename IndexOut_, typename ValueIn_, typename IndexIn_>
struct Dense : public DenseExtractor<oracle_, ValueOut_, IndexOut_> {
    Dense(const Matrix<ValueIn_, IndexIn_>* mat, bool row, MaybeOracle<oracle_, IndexOut_> oracle, const Options& opt) {
        allocate(row ? mat->ncol() : mat->nrow());
        internal_ext = new_extractor<false, oracle_>(mat, row, convert<oracle_, IndexIn_, IndexOut_>(std::move(oracle)), opt);
    }

    Dense(const Matrix<ValueIn_, IndexIn_>* mat, bool row, MaybeOracle<oracle_, IndexOut_> oracle, IndexOut_ block_start, IndexOut_ block_length, const Options& opt) {
        allocate(block_length);
        internal_ext = new_extractor<false, oracle_>(mat, row, convert<oracle_, IndexIn_, IndexOut_>(std::move(oracle)), block_start, block_length, opt);
    }

    Dense(const Matrix<ValueIn_, IndexIn_>* mat, bool row, MaybeOracle<oracle_, IndexOut_> oracle, VectorPtr<IndexOut_> indices_ptr, const Options& opt) {
        allocate(indices_ptr->size());
        internal_ext = new_extractor<false, oracle_>(mat, row, convert<oracle_, IndexIn_, IndexOut_>(std::move(oracle)), convert<IndexIn_>(std::move(indices_ptr)), opt);
    }

private:
    void allocate(size_t n) {
        if constexpr(!no_op) {
            internal_buffer.resize(n);
        }
    }

public:
    const ValueOut_* fetch(IndexOut_ i, ValueOut_* buffer) {
        if constexpr(no_op) {
            return internal_ext->fetch(i, buffer);
        } else {
            auto ptr = internal_ext->fetch(i, internal_buffer.data());
            std::copy_n(ptr, internal_buffer.size(), buffer);
            return buffer;
        }
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, ValueIn_, IndexIn_> > internal_ext;
    static constexpr bool no_op = std::is_same<ValueOut_, ValueIn_>::value;
    typename std::conditional<no_op, bool, std::vector<ValueIn_> >::type internal_buffer;
};

template<bool oracle_, typename ValueOut_, typename IndexOut_, typename ValueIn_, typename IndexIn_>
struct Sparse : public SparseExtractor<oracle_, ValueOut_, IndexOut_> {
    Sparse(const Matrix<ValueIn_, IndexIn_>* mat, bool row, MaybeOracle<oracle_, IndexOut_> oracle, const Options& opt) {
        allocate(row ? mat->ncol() : mat->nrow(), opt);
        internal_ext = new_extractor<true, oracle_>(mat, row, convert<oracle_, IndexIn_, IndexOut_>(std::move(oracle)), opt);
    }

    Sparse(const Matrix<ValueIn_, IndexIn_>* mat, bool row, MaybeOracle<oracle_, IndexOut_> oracle, IndexOut_ block_start, IndexOut_ block_length, const Options& opt) {
        allocate(block_length, opt);
        internal_ext = new_extractor<true, oracle_>(mat, row, convert<oracle_, IndexIn_, IndexOut_>(std::move(oracle)), block_start, block_length, opt);
    }

    Sparse(const Matrix<ValueIn_, IndexIn_>* mat, bool row, MaybeOracle<oracle_, IndexOut_> oracle, VectorPtr<IndexOut_> indices_ptr, const Options& opt) {
        allocate(indices_ptr->size(), opt);
        internal_ext = new_extractor<true, oracle_>(mat, row, convert<oracle_, IndexIn_, IndexOut_>(std::move(oracle)), convert<IndexIn_>(std::move(indices_ptr)), opt);
    }

private:
    void allocate(size_t n, const Options& opt) {
        if constexpr(!no_op_value) {
            if (opt.sparse_extract_value) {
                internal_vbuffer.resize(n);
            }
        }
        if constexpr(!no_op_index) {
            if (opt.sparse_extract_index) {
                internal_ibuffer.resize(n);
            }
        }
    }

public:
    SparseRange<ValueOut_, IndexOut_> fetch(IndexOut_ i, ValueOut_* vbuffer, IndexOut_* ibuffer) {
        IndexIn_* iptr = [&]() {
            if constexpr(no_op_index) {
                return ibuffer;
            } else {
                return internal_ibuffer.data();
            }
        }();

        ValueIn_* vptr = [&]() {
            if constexpr(no_op_value) {
                return vbuffer;
            } else {
                return internal_vbuffer.data();
            }
        }();

        auto range = internal_ext->fetch(i, vptr, iptr);
        SparseRange<ValueOut_, IndexOut_> output(range.number);

        if constexpr(no_op_index) {
            output.index = range.index;
        } else if (range.index != NULL) {
            std::copy_n(range.index, range.number, ibuffer);
            output.index = ibuffer;
        }

        if constexpr(no_op_value) {
            output.value = range.value;
        } else if (range.value != NULL) {
            std::copy_n(range.value, range.number, vbuffer);
            output.value = vbuffer;
        }

        return output;
    }

private:
    std::unique_ptr<SparseExtractor<oracle_, ValueIn_, IndexIn_> > internal_ext;
    static constexpr bool no_op_value = std::is_same<ValueOut_, ValueIn_>::value;
    typename std::conditional<no_op_value, bool, std::vector<ValueIn_> >::type internal_vbuffer;
    static constexpr bool no_op_index = std::is_same<IndexOut_, IndexIn_>::value;
    typename std::conditional<no_op_index, bool, std::vector<IndexIn_> >::type internal_ibuffer;
};

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
    std::unique_ptr<MyopicDenseExtractor<ValueOut_, IndexOut_> > dense(bool row, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::Dense<false, ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(ptr.get(), row, false, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<ValueOut_, IndexOut_> > dense(bool row, IndexOut_ block_start, IndexOut_ block_length, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::Dense<false, ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(ptr.get(), row, false, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<ValueOut_, IndexOut_> > dense(bool row, VectorPtr<IndexOut_> indices_ptr, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::Dense<false, ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(ptr.get(), row, false, std::move(indices_ptr), opt);
    }

    /*********************
     *** Myopic sparse ***
     *********************/
public:
    std::unique_ptr<MyopicSparseExtractor<ValueOut_, IndexOut_> > sparse(bool row, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::Sparse<false, ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(ptr.get(), row, false, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<ValueOut_, IndexOut_> > sparse(bool row, IndexOut_ block_start, IndexOut_ block_length, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::Sparse<false, ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(ptr.get(), row, false, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<ValueOut_, IndexOut_> > sparse(bool row, VectorPtr<IndexOut_> indices_ptr, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::Sparse<false, ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(ptr.get(), row, false, std::move(indices_ptr), opt);
    }

    /**********************
     *** Oracular dense ***
     **********************/
public:
    std::unique_ptr<OracularDenseExtractor<ValueOut_, IndexOut_> > dense(bool row, std::shared_ptr<const Oracle<IndexOut_> > oracle, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::Dense<true, ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(ptr.get(), row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularDenseExtractor<ValueOut_, IndexOut_> > dense(bool row, std::shared_ptr<const Oracle<IndexOut_> > oracle, IndexOut_ block_start, IndexOut_ block_length, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::Dense<true, ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(ptr.get(), row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularDenseExtractor<ValueOut_, IndexOut_> > dense(bool row, std::shared_ptr<const Oracle<IndexOut_> > oracle, VectorPtr<IndexOut_> indices_ptr, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::Dense<true, ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(ptr.get(), row, std::move(oracle), std::move(indices_ptr), opt);
    }

    /***********************
     *** Oracular sparse ***
     ***********************/
public:
    std::unique_ptr<OracularSparseExtractor<ValueOut_, IndexOut_> > sparse(bool row, std::shared_ptr<const Oracle<IndexOut_> > oracle, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::Sparse<true, ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(ptr.get(), row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularSparseExtractor<ValueOut_, IndexOut_> > sparse(bool row, std::shared_ptr<const Oracle<IndexOut_> > oracle, IndexOut_ block_start, IndexOut_ block_length, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::Sparse<true, ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(ptr.get(), row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularSparseExtractor<ValueOut_, IndexOut_> > sparse(bool row, std::shared_ptr<const Oracle<IndexOut_> > oracle, VectorPtr<IndexOut_> indices_ptr, const Options& opt) const {
        return std::make_unique<DelayedCast_internal::Sparse<true, ValueOut_, IndexOut_, ValueIn_, IndexIn_> >(ptr.get(), row, std::move(oracle), std::move(indices_ptr), opt);
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
