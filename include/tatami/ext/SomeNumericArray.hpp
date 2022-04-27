#ifndef TATAMI_SOME_NUMERIC_ARRAY_HPP
#define TATAMI_SOME_NUMERIC_ARRAY_HPP

#include <cstdint>
#include <cstddef>
#include <iterator>

/**
 * @file SomeNumericArray.hpp
 *
 * @brief Defines an array class with run-time numeric type.
 */

namespace tatami {

/**
 * @brief Array of some numeric type, determined at runtime.
 *
 * This holds a pointer to an existing array of some numeric type,
 * mimicking the behavior of `std::vector<T>` for **tatami** use cases. 
 * The aim is to support inputs of variable types without multiple template specializations,
 * especially in cases where there are combinations of such arrays (e.g., `CompressedSparseMatrix`).
 * Of course, this comes with a mild performance penalty as the type must be checked upon extracting any value.
 *
 * @tparam T Type to return when values are extracted.
 * This is allowed to differ from the internal storage type. 
 */
template<typename T = double>
struct SomeNumericArray {
    /**
     * A list of supported types.
     * The letters indicate whether it is an integer (I), unsigned integer (U) or a float,
     * while the numbers specify the number of bits for the type.
     * So, for example, `U16` is an unsigned 16-bit integer, while `F64` is a `double`.
     */
    enum Type { I8, U8, I16, U16, I32, U32, I64, U64, F32, F64 };
   
    /**
     * @param[in] x Pointer to the array of interest, of run-time type specified by `t`.
     * The lifetime of the array should exceed that of the constructed `SomeNumericArray` and any of its copies.
     * @param n Length of the array pointed to by `x`.
     * @param t Type of the array, see `Type`.
     */
    SomeNumericArray(void * x, size_t n, Type t) : len(n), type(t) {
        switch (t) {
            case I8:
                i8 = static_cast<const int8_t*>(x);
                break;
            case U8:
                u8 = static_cast<const uint8_t*>(x);
                break;
            case I16:
                i16 = static_cast<const int16_t*>(x);
                break;
            case U16:
                u16 = static_cast<const uint16_t*>(x);
                break;
            case I32:
                i32 = static_cast<const int32_t*>(x);
                break;
            case U32:
                u32 = static_cast<const uint32_t*>(x);
                break;
            case I64:
                i64 = static_cast<const int64_t*>(x);
                break;
            case U64:
                u64 = static_cast<const uint64_t*>(x);
                break;
            case F32:
                f32 = static_cast<const float*>(x);
                break;
            case F64:
                f64 = static_cast<const double*>(x);
                break;
        }
        return;
    }

    /**
     * @param x Pointer to an existing array of `int8_t`s.
     * @param n Length of the array pointed to by `x`.
     */
    SomeNumericArray(const int8_t* x, size_t n) : i8(x), len(n), type(I8) {}

    /**
     * @param x Pointer to an existing array of `uint8_t`s.
     * @param n Length of the array pointed to by `x`.
     */
    SomeNumericArray(const uint8_t* x, size_t n) : u8(x), len(n), type(U8) {}

    /**
     * @param x Pointer to an existing array of `int16_t`s.
     * @param n Length of the array pointed to by `x`.
     */
    SomeNumericArray(const int16_t* x, size_t n) : i16(x), len(n), type(I16) {}

    /**
     * @param x Pointer to an existing array of `uint16_t`s.
     * @param n Length of the array pointed to by `x`.
     */
    SomeNumericArray(const uint16_t* x, size_t n) : u16(x), len(n), type(U16) {}
    
    /**
     * @param x Pointer to an existing array of `int32_t`s.
     * @param n Length of the array pointed to by `x`.
     */
    SomeNumericArray(const int32_t* x, size_t n) : i32(x), len(n), type(I32) {}
    
    /**
     * @param x Pointer to an existing array of `uint32_t`s.
     * @param n Length of the array pointed to by `x`.
     */
    SomeNumericArray(const uint32_t* x, size_t n) : u32(x), len(n), type(U32) {}
    
    /**
     * @param x Pointer to an existing array of `int64_t`s.
     * @param n Length of the array pointed to by `x`.
     */
    SomeNumericArray(const int64_t* x, size_t n) : i64(x), len(n), type(I64) {}
    
    /**
     * @param x Pointer to an existing array of `uint64_t`s.
     * @param n Length of the array pointed to by `x`.
     */
    SomeNumericArray(const uint64_t* x, size_t n) : u64(x), len(n), type(U64) {}
    
    /**
     * @param x Pointer to an existing array of `float`s.
     * @param n Length of the array pointed to by `x`.
     */
    SomeNumericArray(const float* x, size_t n) : f32(x), len(n), type(F32) {}
    
    /**
     * @param x Pointer to an existing array of `double`s.
     * @param n Length of the array pointed to by `x`.
     */
    SomeNumericArray(const double* x, size_t n) : f64(x), len(n), type(F64) {}

private:
    Type type;
    size_t len;

    const int8_t* i8 = NULL;
    const uint8_t* u8 = NULL;
    const int16_t* i16 = NULL;
    const uint16_t* u16 = NULL;
    const int32_t* i32 = NULL;
    const uint32_t* u32 = NULL;
    const int64_t* i64 = NULL;
    const uint64_t* u64 = NULL;
    const float* f32 = NULL;
    const double* f64 = NULL;
public:
    /**
     * @param i Positional index on the array.
     * @return Value of the `i`-th element as a `T`.
     */
    T operator[](size_t i) const {
        switch (type) {
            case I8:
                return i8[i];
            case U8:
                return u8[i];
            case I16:
                return i16[i];
            case U16:
                return u16[i];
            case I32:
                return i32[i];
            case U32:
                return u32[i];
            case I64:
                return i64[i];
            case U64:
                return u64[i];
            case F32:
                return f32[i];
            case F64:
                return f64[i];
        }
    }

    /**
     * @return Length of the array, in terms of the number of elements of the specified type.
     */
    size_t size() const {
        return len;
    }
public:
    /**
     * @brief Random-access iterator class.
     *
     * This mimics the const iterators for `std::vector` types.
     */
    struct Iterator {
        /**
         * Default constructor.
         */
        Iterator() : parent(NULL), index(0) {}

        /**
         * @param p Pointer to the parental `SomeNumericArray` object.
         * @param i Index along the parental array, representing the current position of the iterator.
         *
         * Needless to say, we assume that the parental array outlives the iterator.
         */
        Iterator(const SomeNumericArray* p, size_t i) : parent(p), index(i) {}

    public:
        // Tags to pretend it's an iterator, at least enough for tatami to work;
        // see https://internalpointers.com/post/writing-custom-iterators-modern-cpp for details.

        /**
         * Random access iterator tag.
         */
        using iterator_category = std::random_access_iterator_tag;

        /**
         * Difference type.
         */
        using difference_type = std::ptrdiff_t;

        /**
         * Value type.
         */
        using value_type = T;

        /**
         * Pointer type, note the `const`.
         */
        using pointer = const T*;  

        /**
         * Reference type, note the `const`.
         */
        using reference = const T&; 

    private:
        const SomeNumericArray* parent;
        size_t index;

    public:
        /**
         * @return The value at the current position of the iterator on the parental `SomeNumericArray` object.
         */
        value_type operator*() const {
            return (*parent)[index];
        }

        /**
         * @param i The number of elements to add to the current position of the iterator to obtain a new position.
         * @return The value at the new position on the parental `SomeNumericArray` object.
         */
        value_type operator[](size_t i) const {
            return (*parent)[index + i];
        }

    public:
        /**
         * @param right Another `Iterator` object referencing the same parental array.
         * @return Whether the current iterator and `right` are pointing to the same position.
         */
        bool operator==(const Iterator& right) const {
            return index == right.index;
        }

        /**
         * @param right Another `Iterator` object referencing the same parental array.
         * @return Whether the current iterator and `right` are pointing to different positions.
         */
        bool operator!=(const Iterator& right) const {
            return !(*this == right);
        }

        /**
         * @param right Another `Iterator` object referencing the same parental array.
         * @return Whether the current iterator is pointing to an earlier position than `right`.
         */
        bool operator<(const Iterator& right) const {
            return index < right.index;
        }

        /**
         * @param right Another `Iterator` object referencing the same parental array.
         * @return Whether the current iterator is pointing to an equal or later position than `right`.
         */
        bool operator>=(const Iterator& right) const {
            return !(*this < right);
        }

        /**
         * @param right Another `Iterator` object referencing the same parental array.
         * @return Whether the current iterator is pointing to a later position than `right`.
         */
        bool operator>(const Iterator& right) const {
            return index > right.index;
        }

        /**
         * @param right Another `Iterator` object referencing the same parental array.
         * @return Whether the current iterator is pointing to an equal or earlier position than `right`.
         */
        bool operator<=(const Iterator& right) const {
            return !(*this > right);
        }

    public:
        /**
         * @param n Number of elements to advance the current iterator.
         * @return The iterator's position is moved forward by `n`, and a reference to the iterator is returned.
         */
        Iterator& operator+=(size_t n) {
            index += n;
            return *this;
        }

        /**
         * @return The iterator's position is moved forward by 1, and a reference to the iterator is returned.
         */
        Iterator& operator++() {
            *this += 1;
            return *this;
        }

        /**
         * @return The position of the iterator is moved forward by 1,
         * while a copy of the iterator (pre-increment) is returned.
         */
        Iterator operator++(int) {
            auto copy = *this;
            ++(*this);
            return copy;
        }

        /**
         * @param n Number of elements to move back the current iterator.
         * @return The iterator's position is moved backward by `n`, and a reference to the iterator is returned.
         */
        Iterator& operator-=(size_t n) {
            index -= n;
            return *this;
        }

        /**
         * @return The iterator's position is moved backwards by 1, and a reference to the iterator is returned.
         */
        Iterator& operator--() {
            *this -= 1;
            return *this;
        }

        /**
         * @return The position of the iterator is moved back by 1,
         * while a copy of the iterator (pre-decrement) is returned.
         */
        Iterator operator--(int) {
            auto copy = *this;
            --(*this);
            return copy;
        }

    public:
        /**
         * @param n Number of elements to advance the iterator.
         * @return A new iterator is returned at the position of the current iterator plus `n`.
         */
        Iterator operator+(size_t n) const {
            return Iterator(parent, index + n);
        }

        /**
         * @param n Number of elements to move back the iterator.
         * @return A new iterator is returned at the position of the current iterator minus `n`.
         */
        Iterator operator-(size_t n) const {
            return Iterator(parent, index - n);
        }

        /**
         * @param n Number of elements to advance the iterator.
         * @param it An existing `Iterator`.
         * @return A new iterator is returned at the position of `it` plus `n`.
         */
        friend Iterator operator+(size_t n, const Iterator& it) {
            return Iterator(it.parent, it.index + n);
        }

        /**
         * @param right Another `Iterator` object referencing the same parental array.
         * @return The difference in positions of the two iterators.
         */
        std::ptrdiff_t operator-(const Iterator& right) const {
            std::ptrdiff_t out;
            if (right.index > index) {
                out = right.index - index;
                out *= -1;
            } else {
                out = index - right.index;
            }
            return out;
        }
    };

    /**
     * @return An iterator to the start of the `SomeNumericArray`.
     */
    Iterator begin() const {
        return Iterator(this, 0);
    }

    /**
     * @return An iterator to the end of the `SomeNumericArray`.
     */
    Iterator end() const {
        return Iterator(this, this->size());
    }
};

}

#endif
