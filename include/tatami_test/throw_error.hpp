#ifndef TATAMI_THROW_ERROR_HPP
#define TATAMI_THROW_ERROR_HPP

#include <string>
#include <gtest/gtest.h>

namespace tatami_test {

template<class Function_>
void throws_error(Function_ fun, const std::string& msg) {
    try {
        fun();
        FAIL() << "expected error message '" << msg << "', got no error";
    } catch (std::exception& e) {
        std::string observed(e.what());
        if (observed.find(msg) == std::string::npos) {
            FAIL() << "expected error message '" << msg << "', got '" << observed << "'";
        }
    }
}

}

#endif
