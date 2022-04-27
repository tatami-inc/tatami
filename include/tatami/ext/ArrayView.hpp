#ifndef TATAMI_ARRAY_VIEW_HPP
#define TATAMI_ARRAY_VIEW_HPP

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
     * @param[in] p Pointer to the start of the array.
     * The lifetime of the array is assumed to exceed that of the constructed `ArrayView` instance.
     * @param n Number of array elements.
     */
    ArrayView(const T* p, size_t n) : ptr(p), num(n) {}

    /**
     * @return Number of array elements.
     */
    size_t size() const { return num; }

    /**
     * @return Pointer to the start of the array.
     */
    const T* data() const { return ptr; }

    /**
     * @return Pointer to the start of the array, for use in templated functions expecting a `std::vector`.
     */
    const T* begin() const { return ptr; }

    /**
     * @return Pointer to the end of the array, for use in templated functions expecting a `std::vector`.
     */
    const T* end() const { return ptr + num; }

    /**
     * @param i Index of the array.
     * @return Value of the array at element `i`.
     */
    T operator[](size_t i) const {
        return ptr[i];
    }
private:
    const T* ptr;
    size_t num;
};

}

#endif
