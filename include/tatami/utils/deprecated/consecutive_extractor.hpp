#ifndef TATAMI_CONSEUCTIVE_EXTRACTOR_HPP
#define TATAMI_CONSEUCTIVE_EXTRACTOR_HPP

#include <memory>
#include "../make_Extractor.hpp"
#include "../ConsecutiveOracle.hpp"

namespace tatami {

template<bool sparse_, typename Value_, typename Index_, typename ... Args_>
[[deprecated("use 'make_Extractor()' and 'make_ConsecutiveOracle()' instead")]]
auto consecutive_extractor(const Matrix<Value_, Index_>& matrix, const bool row, const Index_ iter_start, const Index_ iter_length, Args_&&... args) {
    return make_Extractor<sparse_, true>(
        matrix,
        row,
        make_ConsecutiveOracle(iter_start, iter_length),
        std::forward<Args_>(args)...
    );
}

}

#endif
