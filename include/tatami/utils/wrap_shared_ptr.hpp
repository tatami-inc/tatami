#ifndef TATAMI_WRAP_SHARED_PTR_HPP
#define TATAMI_WRAP_SHARED_PTR_HPP

#include <memory>
#include "../base/Matrix.hpp"

/**
 * @file wrap_shared_ptr.hpp
 *
 * @brief Wrap a raw `tatami::Matrix` pointer inside a mock shared pointer.
 */

namespace tatami {

/**
 * Wrap a raw pointer inside a `shared_ptr`, typically to enable use of a raw `tatami::Matrix` pointer with delayed operation wrappers.
 * This enables use of delayed operations inside functions that accept a raw pointer to an externally owned `tatami::Matrix`.
 *
 * @tparam Value_ Numeric type for the matrix value.
 * @tparam Index_ Integer type for the row/column indices.
 * 
 * @param ptr A pointer to a `tatami::Matrix` instance.
 * 
 * @return A shared pointer to the same object addressed by `ptr`.
 * The assumption is that `ptr` will always outlive the returned pointer.
 */
template<typename Value_, typename Index_>
std::shared_ptr<const Matrix<Value_, Index_> > wrap_shared_ptr(const Matrix<Value_, Index_>* const ptr) {
    // Using an aliasing constructor, see https://stackoverflow.com/a/36691828
    return std::shared_ptr<const Matrix<Value_, Index_> >(std::shared_ptr<const Matrix<Value_, Index_> >{}, ptr);
}

}

#endif
