#ifndef TATAMI_SPARSE_UTILS_HPP
#define TATAMI_SPARSE_UTILS_HPP

#include <type_traits>

namespace tatami {

namespace sparse_utils {

/*
 * If PointerStorage_ is just a bool, then we assume that the indices for each
 * primary element are fragmented into separate vectors. In such cases, we can
 * just call size() to get the upper limit and we set the lower limit to zero.
 *
 * If PointerStorage_ is not a bool, it is assumed to be a vector of pointers
 * in the typical compressed sparse style. In such cases, we actually need to
 * do some work to get the lower/upper bounds of each primary element.
 */

template<class PointerStorage_>
constexpr bool is_fragmented() {
    return std::is_same<PointerStorage_, bool>::value;
}

template<class PointerStorage_, class IndexStorage_, typename Index_>
const auto& get_indices(const IndexStorage_& all_indices, Index_ primary) {
    if constexpr(is_fragmented<PointerStorage_>()) {
        return all_indices[primary];
    } else {
        return all_indices;
    }
}

template<class IndexStorage_, class PointerStorage_, typename Index_>
auto get_upper_limit(const IndexStorage_& indices, const PointerStorage_& indptrs, Index_ primary) {
    if constexpr(is_fragmented<PointerStorage_>()) {
        // just making sure we got passed in a vector of numbers, rather than a fragmented vector of vectors.
        static_assert(std::is_arithmetic<typename std::remove_reference<decltype(indices[0])>::type>::value); 
        return indices.size();
    } else {
        return indptrs[primary + 1];
    }
}

template<class PointerStorage_, typename Index_>
auto get_lower_limit(const PointerStorage_& indptrs, Index_ primary) {
    if constexpr(is_fragmented<PointerStorage_>()) {
        return 0;
    } else {
        return indptrs[primary];
    }
}

}

}

#endif
