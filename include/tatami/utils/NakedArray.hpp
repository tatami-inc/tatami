#ifndef TATAMI_NAKED_ARRAY_HPP
#define TATAMI_NAKED_ARRAY_HPP

/**
 * @file NakedArray.hpp
 *
 * @brief Vector-like wrapper around a raw pointer to an array.
 */

namespace tatami {

/**
 * @brief Wrapper around a raw pointer to an array.
 *
 * This serves as a drop-in replacement for some templated uses of `std::vector` in **tatami** classes.
 * The idea is to enable use of **tatami** functions while avoiding unnecessary copies,
 * provided that the raw pointer outlives the instance of the **tatami** class created from it.
 *
 * @tparam T Type of the pointer.
 */
template<typename T>
class NakedArray {
public:
    /**
     * Default constructor.
     */
    NakedArray() : ptr(NULL), len(0) {}

    /**
     * @param p Pointer to an array of values.
     * @param n Length of the array pointed to by `p`.
     */
    NakedArray(const T* p, size_t n) : ptr(p), len(n) {}

    /**
     * @return Pointer to the start of the array.
     */
    const T* begin() const { return ptr; }

    /**
     * @return Pointer to the end of the array.
     */
    const T* end() const { return ptr + len; }

    /**
     * @return Size of the array.
     */
    size_t size () const { return len; }

    /**
     * @param i Position on the array.
     * @return The value at the specified position.
     */
    T operator[](size_t i) const {
        return ptr[i];
    }

    /**
     * @return Pointer to the start of the array.
     */
    const T* data() const { return ptr; }
private:
    const T* ptr;
    size_t len;
};

}

#endif
