#ifndef TATAMI_COPY_HPP
#define TATAMI_COPY_HPP

#include <algorithm>

/**
 * @file copy.hpp
 * @brief Copy data from one buffer to another.
 */

namespace tatami {

/**
 * @tparam Value_ Type of value being copied.
 * @tparam Size_ Integer type of the array length.
 *
 * @param[in] input Pointer to a source array of size `n`.
 * @param n Length of the array.
 * This should fit into a `std::size_t`, even if `Size_` is of a larger type.
 * @param[out] output Pointer to a destination array of size `n`.
 *
 * @return Values are copied from `input` to `output`, and `output` is returned.
 * This is a no-op if `input == output`.
 */
template<typename Value_, typename Size_>
Value_* copy_n(const Value_* input, Size_ n, Value_* output) {
    if (input != output) {
        std::copy_n(input, n, output);
    }
    return output;
}

}

#endif
