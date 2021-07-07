#ifndef PRINT_VECTOR_HPP
#define PRINT_VECTOR_HPP

#include <iostream>
#include <numeric>
#include <iomanip>

template<class IT>
void print_vector(IT start, IT end, int width = 6, int precision = 2) {
    bool first = true;
    std::cout << "[ "; 
    for (IT it = start; it != end; ++it) {
        if (!first) {
            std::cout << ", ";
        }
        std::cout << std::setw(width) << std::fixed << std::setprecision(precision) << *it;
        first = false;
    }
    std::cout << " ]" << std::endl;
}

#endif
