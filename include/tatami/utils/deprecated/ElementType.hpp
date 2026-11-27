#ifndef TATAMI_ELEMENTTYPE_HPP
#define TATAMI_ELEMENTTYPE_HPP

#include "../copy.hpp"

namespace tatami {

template<class Array_>
using ElementType = I<decltype(std::declval<Array_>()[0])>;

}

#endif
