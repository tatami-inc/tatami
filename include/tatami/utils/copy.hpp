#ifndef TATAMI_COPY_HPP
#define TATAMI_COPY_HPP

#include <algorithm>

/**
 * @brief copy.hpp
 * @file Copy data from one buffer to another.
 */

namespace tatami {

/**
 * @tparam Value_ Type of value being copied.
 * @param[in] input Pointer to a source array of size `n`.
 * @param n Length of the array.
 * @param[out] output Pointer to a destination array of size `n`.
 *
 * @return Values are copied from `input` to `output`, and `output` is returned.
 * This is a no-op if `input == output`.
 */
template<typename Value_>
Value_* copy_n(const Value_* input, size_t n, Value_* output) {
    if (input != output) {
        std::copy_n(input, n, output);
    }
    return output;
}

}

#endif
