#ifndef TATAMI_TYPES_H
#define TATAMI_TYPES_H

namespace tatami {

enum content_type { _unknown, _double, _int32_t, _std_string };

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
