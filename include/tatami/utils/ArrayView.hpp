#ifndef TATAMI_ARRAY_VIEW_HPP
#define TATAMI_ARRAY_VIEW_HPP

#include <cstddef>

/**
 * @file ArrayView.hpp
 *
 * @brief Defines a **tatami**-compatible array view.
 */
namespace tatami {

/**
 * @brief View into a pre-allocated array.
 *
 * This allows us to use a pre-existing array in the **tatami** classes, assuming that the array lives longer than the view created on it.
 * We provide some methods to mimic a `std::vector` for use within the **tatami** constructors.
 * This is implemented in lieu of having access to C++20's `std::span` class.
 *
 * @tparam T Array type, usually numeric.
 */
template<typename T>
class ArrayView {
public:
    /**
     * @param[in] ptr Pointer to the start of the array.
     * The lifetime of the array is assumed to exceed that of the constructed `ArrayView` instance.
     * @param number Number of array elements.
     */
    ArrayView(const T* ptr, std::size_t number) : my_ptr(ptr), my_number(number) {}

    /**
     * Default constructor to create a zero-length view.
     */
    ArrayView() : ArrayView(NULL, 0) {}

    /**
     * @return Number of array elements.
     */
    std::size_t size() const { return my_number; }

    /**
     * @return Pointer to the start of the array.
     */
    const T* data() const { return my_ptr; }

    /**
     * @return Pointer to the start of the array, for use in templated functions expecting a `std::vector`.
     */
    const T* begin() const { return my_ptr; }

    /**
     * @return Pointer to the end of the array, for use in templated functions expecting a `std::vector`.
     */
    const T* end() const { return my_ptr + my_number; }

    /**
     * @param i Index of the array.
     * @return Value of the array at element `i`.
     */
    T operator[](std::size_t i) const {
        return my_ptr[i];
    }

private:
    const T* my_ptr;
    std::size_t my_number;
};

}

#endif
