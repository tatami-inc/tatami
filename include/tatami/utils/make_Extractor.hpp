#ifndef TATAMI_NEW_EXTRACTOR_HPP
#define TATAMI_NEW_EXTRACTOR_HPP

#include <memory>
#include <type_traits>

#include "../base/Matrix.hpp"

/**
 * @file make_Extractor.hpp
 * @brief Templated construction of a new extractor.
 */

namespace tatami {

/**
 * @tparam oracle_ Whether an oracle should be supplied.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * Maybe an `Oracle`, maybe a placeholder boolean, depending on `oracle_`.
 */
template<bool oracle_, typename Index_>
using MaybeOracle = typename std::conditional<oracle_, std::shared_ptr<const Oracle<Index_> >, bool>::type;

/**
 * @cond
 */
template<bool sparse_, bool oracle_, typename Value_, typename Index_, typename ... Args_>
std::unique_ptr<Extractor<sparse_, oracle_, Value_, Index_> > make_Extractor_base(
    const Matrix<Value_, Index_>& matrix,
    const bool row, MaybeOracle<oracle_, Index_> oracle,
    Args_&&... args
) {
    if constexpr(sparse_) {
        if constexpr(oracle_) {
            return matrix.sparse(row, std::move(oracle), std::forward<Args_>(args)...);
        } else {
            return matrix.sparse(row, std::forward<Args_>(args)...);
        }
    } else {
        if constexpr(oracle_) {
            return matrix.dense(row, std::move(oracle), std::forward<Args_>(args)...);
        } else {
            return matrix.dense(row, std::forward<Args_>(args)...);
        }
    }
}
/**
 * @endcond
 */

/**
 * Create a new `Extractor` that fetches the entire extent of each element of the target dimension. 
 * This is intended to be used in templated code that switches between dense/sparse and myopic/oracular extraction.
 *
 * @tparam sparse_ Whether to perform sparse retrieval.
 * @tparam oracle_ Whether an oracle should be supplied.
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column index. 
 *
 * @param[in] matrix A `tatami::Matrix` object to iterate over.
 * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
 * @param oracle Pointer to an oracle if `oracle_ = true`, otherwise a placeholder boolean that is ignored.
 * @param opt Additional arguments to pass to the relevant `Matrix::dense()` or `Matrix::sparse()` overload.
 *
 * @return Pointer to an `Extractor` to access elements of the target dimension of `matrix`.
 */
template<bool sparse_, bool oracle_, typename Value_, typename Index_>
std::unique_ptr<Extractor<sparse_, oracle_, Value_, Index_> > make_Extractor(
    const Matrix<Value_, Index_>& matrix,
    const bool row,
    MaybeOracle<oracle_, Index_> oracle,
    const tatami::Options& opt
) {
    return make_Extractor_base<sparse_, oracle_>(matrix, row, std::move(oracle), opt);
}

/**
 * Create a new `Extractor` that fetches a contiguous block of each element of the target dimension. 
 * This is intended to be used in templated code that switches between dense/sparse and myopic/oracular extraction.
 *
 * @tparam sparse_ Whether to perform sparse retrieval.
 * @tparam oracle_ Whether an oracle should be supplied.
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column index. 
 *
 * @param[in] matrix A `tatami::Matrix` object to iterate over.
 * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
 * @param oracle Pointer to an oracle if `oracle_ = true`, otherwise a placeholder boolean that is ignored.
 * @param block_start,block_length,opt Additional arguments to pass to the relevant `Matrix::dense()` or `Matrix::sparse()` overload.
 *
 * @return Pointer to an `Extractor` to access a block of each element of the target dimension of `matrix`.
 */
template<bool sparse_, bool oracle_, typename Value_, typename Index_>
std::unique_ptr<Extractor<sparse_, oracle_, Value_, Index_> > make_Extractor(
    const Matrix<Value_, Index_>& matrix,
    const bool row,
    MaybeOracle<oracle_, Index_> oracle,
    Index_ block_start,
    Index_ block_length,
    const tatami::Options& opt
) {
    return make_Extractor_base<sparse_, oracle_>(matrix, row, std::move(oracle), block_start, block_length, opt);
}

/**
 * Create a new `Extractor` that fetches an indexed subset of each element of the target dimension. 
 * This is intended to be used in templated code that switches between dense/sparse and myopic/oracular extraction.
 *
 * @tparam sparse_ Whether to perform sparse retrieval.
 * @tparam oracle_ Whether an oracle should be supplied.
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column index. 
 *
 * @param[in] matrix A `tatami::Matrix` object to iterate over.
 * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
 * @param oracle Pointer to an oracle if `oracle_ = true`, otherwise a placeholder boolean that is ignored.
 * @param indices_ptr,opt Additional arguments to pass to the relevant `Matrix::dense()` or `Matrix::sparse()` overload.
 *
 * @return Pointer to an `Extractor` to access an indexed subset of each element of the target dimension of `matrix`.
 */
template<bool sparse_, bool oracle_, typename Value_, typename Index_>
std::unique_ptr<Extractor<sparse_, oracle_, Value_, Index_> > make_Extractor(
    const Matrix<Value_, Index_>& matrix,
    const bool row,
    MaybeOracle<oracle_, Index_> oracle,
    VectorPtr<Index_> indices_ptr,
    const tatami::Options& opt
) {
    return make_Extractor_base<sparse_, oracle_>(matrix, row, std::move(oracle), std::move(indices_ptr), opt);
}

/**
 * @cond
 */
template<bool sparse_, bool oracle_, typename Value_, typename Index_, typename ... Args_>
[[deprecated("use 'make_Extractor()' instead")]]
auto new_extractor(const Matrix<Value_, Index_>& matrix, const bool row, MaybeOracle<oracle_, Index_> oracle, Args_&&... args) {
    return make_Extractor(matrix, row, std::move(oracle), std::forward<Args_>(args)...);
}

template<bool sparse_, bool oracle_, typename Value_, typename Index_, typename ... Args_>
[[deprecated("use 'make_Extractor()' instead")]]
auto new_extractor(const Matrix<Value_, Index_>* ptr, bool row, tatami::MaybeOracle<oracle_, Index_> oracle, Args_&&... args) {
    return new_extractor<sparse_, oracle_, Value_, Index_>(*ptr, row, std::move(oracle), std::forward<Args_>(args)...);
}
/**
 * @endcond
 */

}

#endif
