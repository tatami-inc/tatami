#ifndef TATAMI_CONTENT_TYPE_H
#define TATAMI_CONTENT_TYPE_H

/**
 * @file content_type.hpp
 *
 * Helper functions focusing on the types of the matrix values.
 */

namespace tatami {

/**
 * Available content types.
 */
enum content_type { _unknown, _double, _int32_t, _std_string };

/**
 * @tparam T The type.
 *
 * @return A `content_type` corresponding to `T`.
 */
template<typename T> 
content_type determine_content_type() {
    if constexpr(std::is_same<T, double>::value) {
        return _double;
    } else if constexpr(std::is_same<T, int32_t>::value) {
        return _int32_t;
    } else if constexpr(std::is_same<T, std::string>::value) {
        return _std_string;
    }
    return _unknown;
}

}

#endif
