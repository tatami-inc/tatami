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
 * @param ptr A pointer to a `tatami::Matrix` instance.
 * 
 * @return A shared pointer to the same object addressed by `ptr`.
 * The assumption is that `ptr` will always outlive the returned pointer.
 */
template<typename T, typename IDX>
std::shared_ptr<const Matrix<T, IDX> > wrap_shared_ptr(const Matrix<T, IDX>* ptr) {
    // Using a no-op deleter in a shared pointer to get it to work with, e.g.,
    // delayed operations without actually taking ownership of 'ptr'. This
    // shared pointer should only be used as long as 'ptr' is alive.
    return std::shared_ptr<const Matrix<T, IDX> >(ptr, [](const Matrix<T, IDX>*){});
}

}

#endif
