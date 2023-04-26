#ifndef TATAMI_BASE_UTILS_HPP
#define TATAMI_BASE_UTILS_HPP

#include "Matrix.hpp"
#include "Options.hpp"
#include <memory>

namespace tatami {

/**
 * @cond
 */
// Default to false.
template<typename T, class V, typename = int>
struct has_data {
    static const bool value = false;
};

// Specialization is only run if it _has_ a data method.
template<typename T, class V>
struct has_data<T, V, decltype((void) V().data(), 0)> { 
    static const bool value = std::is_same<T*, decltype(V().data())>::value;
};
/**
 * @endcond
 */

/**
 * @tparam row_ Whether to iterate over rows.
 * @tparam sparse_ Whether to perform sparse retrieval.
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 * @tparam Args_ Further arguments.
 *
 * @param[in] ptr Pointer to a `Matrix` object to iterate over.
 * @param opt Options for iteration and extraction.
 * @param args Zero or more additional arguments to pass to methods like `Matrix::dense_row()`.
 *
 * @return A `Extractor` object to access the requested dimension of `ptr`.
 */
template<bool row_, bool sparse_, typename Value_, typename Index_, typename ... Args_>
auto new_extractor(const Matrix<Value_, Index_>* ptr, const Options& opt, Args_... args) {
    if constexpr(sparse_) {
        if constexpr(row_) {
            return ptr->sparse_row(std::move(args)..., opt);
        } else {
            return ptr->sparse_column(std::move(args)..., opt);
        }
    } else {
        if constexpr(row_) {
            return ptr->dense_row(std::move(args)..., opt);
        } else {
            return ptr->dense_column(std::move(args)..., opt);
        }
    }
}

/**
 * @tparam row_ Whether to iterate over rows.
 * @tparam sparse_ Whether to perform sparse retrieval.
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 * @tparam Args_ Further arguments.
 *
 * @param[in] ptr Pointer to a `Matrix` object to iterate over.
 *
 * @return A `Extractor` object to access the requested dimension of `ptr`.
 */
template<bool row_, bool sparse_, typename Value_, typename Index_, typename ... Args_>
auto new_extractor(const Matrix<Value_, Index_>* ptr) {
    if constexpr(sparse_) {
        if constexpr(row_) {
            return ptr->sparse_row();
        } else {
            return ptr->sparse_column();
        }
    } else {
        if constexpr(row_) {
            return ptr->dense_row();
        } else {
            return ptr->dense_column();
        }
    }
}

}

#endif
