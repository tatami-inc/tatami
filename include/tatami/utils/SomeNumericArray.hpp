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
 * Types supported in `SomeNumericArray`.
 * The letters indicate whether it is an integer (I), unsigned integer (U) or a float,
 * while the numbers specify the number of bits for the type.
 * So, for example, `U16` is an unsigned 16-bit integer, while `F64` is a `double`.
 */
enum SomeNumericType { I8, U8, I16, U16, I32, U32, I64, U64, F32, F64 };

/**
 * @brief Array of some numeric type, determined at runtime.
 *
 * This holds a pointer to an existing array of some numeric type,
 * mimicking the behavior of `std::vector<Value_>` for **tatami** use cases. 
 * The aim is to support inputs of variable types without multiple template specializations,
 * especially in cases where there are combinations of such arrays (e.g., `CompressedSparseMatrix`).
 * Of course, this comes with a mild performance penalty as the type must be checked upon extracting any value.
 *
 * @tparam Value_ Type to return when values are extracted.
 * This is allowed to differ from the internal storage type. 
 */
template<typename Value_ = double>
class SomeNumericArray {
public:
    /**
     * @param[in] ptr Pointer to the array of interest, of run-time type specified by `type`.
     * The lifetime of the array should exceed that of the constructed `SomeNumericArray` and any of its copies.
     * @param number Length of the array pointed to by `ptr`.
     * @param type Type of the array. 
     */
    SomeNumericArray(void* ptr, size_t number, SomeNumericType type) : my_number(number), my_type(type) {
        switch(my_type) {
            case I8:
                my_i8 = static_cast<const int8_t*>(ptr);
                break;
            case U8:
                my_u8 = static_cast<const uint8_t*>(ptr);
                break;
            case I16:
                my_i16 = static_cast<const int16_t*>(ptr);
                break;
            case U16:
                my_u16 = static_cast<const uint16_t*>(ptr);
                break;
            case I32:
                my_i32 = static_cast<const int32_t*>(ptr);
                break;
            case U32:
                my_u32 = static_cast<const uint32_t*>(ptr);
                break;
            case I64:
                my_i64 = static_cast<const int64_t*>(ptr);
                break;
            case U64:
                my_u64 = static_cast<const uint64_t*>(ptr);
                break;
            case F32:
                my_f32 = static_cast<const float*>(ptr);
                break;
            case F64:
                my_f64 = static_cast<const double*>(ptr);
                break;
        }
    }

    /**
     * @param ptr Pointer to an existing array of `int8_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const int8_t* ptr, size_t number) : my_i8(ptr), my_number(number), my_type(I8) {}

    /**
     * @param ptr Pointer to an existing array of `uint8_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const uint8_t* ptr, size_t number) : my_u8(ptr), my_number(number), my_type(U8) {}

    /**
     * @param ptr Pointer to an existing array of `int16_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const int16_t* ptr, size_t number) : my_i16(ptr), my_number(number), my_type(I16) {}

    /**
     * @param ptr Pointer to an existing array of `uint16_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const uint16_t* ptr, size_t number) : my_u16(ptr), my_number(number), my_type(U16) {}
    
    /**
     * @param ptr Pointer to an existing array of `int32_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const int32_t* ptr, size_t number) : my_i32(ptr), my_number(number), my_type(I32) {}
    
    /**
     * @param ptr Pointer to an existing array of `uint32_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const uint32_t* ptr, size_t number) : my_u32(ptr), my_number(number), my_type(U32) {}
    
    /**
     * @param ptr Pointer to an existing array of `int64_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const int64_t* ptr, size_t number) : my_i64(ptr), my_number(number), my_type(I64) {}
    
    /**
     * @param ptr Pointer to an existing array of `uint64_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const uint64_t* ptr, size_t number) : my_u64(ptr), my_number(number), my_type(U64) {}
    
    /**
     * @param ptr Pointer to an existing array of `float`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const float* ptr, size_t number) : my_f32(ptr), my_number(number), my_type(F32) {}
    
    /**
     * @param ptr Pointer to an existing array of `double`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const double* ptr, size_t number) : my_f64(ptr), my_number(number), my_type(F64) {}

private:
    const int8_t* my_i8 = NULL;
    const uint8_t* my_u8 = NULL;
    const int16_t* my_i16 = NULL;
    const uint16_t* my_u16 = NULL;
    const int32_t* my_i32 = NULL;
    const uint32_t* my_u32 = NULL;
    const int64_t* my_i64 = NULL;
    const uint64_t* my_u64 = NULL;
    const float* my_f32 = NULL;
    const double* my_f64 = NULL;

    size_t my_number;
    SomeNumericType my_type;

public:
    /**
     * @param i Positional index on the array.
     * @return Value of the `i`-th element as a `Value_`.
     */
    Value_ operator[](size_t i) const {
        switch (my_type) {
            case I8:
                return my_i8[i];
            case U8:
                return my_u8[i];
            case I16:
                return my_i16[i];
            case U16:
                return my_u16[i];
            case I32:
                return my_i32[i];
            case U32:
                return my_u32[i];
            case I64:
                return my_i64[i];
            case U64:
                return my_u64[i];
            case F32:
                return my_f32[i];
            case F64:
                return my_f64[i];
        }
        return 0; // shouldn't reach here, but whatever.
    }

    /**
     * @return Length of the array, in terms of the number of elements of the specified type.
     */
    size_t size() const {
        return my_number;
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
        Iterator() : my_parent(NULL), my_index(0) {}

        /**
         * @param parent Pointer to the parental `SomeNumericArray` object.
         * @param index Index along the parental array, representing the current position of the iterator.
         *
         * Needless to say, we assume that the parental array outlives the iterator.
         */
        Iterator(const SomeNumericArray* parent, size_t index) : my_parent(parent), my_index(index) {}

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
        using value_type = Value_;

        /**
         * Pointer type, note the `const`.
         */
        using pointer = const Value_*;  

        /**
         * Reference type, note the `const`.
         */
        using reference = const Value_&; 

    private:
        const SomeNumericArray* my_parent;
        size_t my_index;

    public:
        /**
         * @return The value at the current position of the iterator on the parental `SomeNumericArray` object.
         */
        value_type operator*() const {
            return (*my_parent)[my_index];
        }

        /**
         * @param i The number of elements to add to the current position of the iterator to obtain a new position.
         * @return The value at the new position on the parental `SomeNumericArray` object.
         */
        value_type operator[](size_t i) const {
            return (*my_parent)[my_index + i];
        }

    public:
        /**
         * @param right Another `Iterator` object referencing the same parental array.
         * @return Whether the current iterator and `right` are pointing to the same position.
         */
        bool operator==(const Iterator& right) const {
            return my_index == right.my_index;
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
            return my_index < right.my_index;
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
            return my_index > right.my_index;
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
            my_index += n;
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
            my_index -= n;
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
            return Iterator(my_parent, my_index + n);
        }

        /**
         * @param n Number of elements to move back the iterator.
         * @return A new iterator is returned at the position of the current iterator minus `n`.
         */
        Iterator operator-(size_t n) const {
            return Iterator(my_parent, my_index - n);
        }

        /**
         * @param n Number of elements to advance the iterator.
         * @param it An existing `Iterator`.
         * @return A new iterator is returned at the position of `it` plus `n`.
         */
        friend Iterator operator+(size_t n, const Iterator& it) {
            return Iterator(it.my_parent, it.my_index + n);
        }

        /**
         * @param right Another `Iterator` object referencing the same parental array.
         * @return The difference in positions of the two iterators.
         */
        std::ptrdiff_t operator-(const Iterator& right) const {
            std::ptrdiff_t out;
            if (right.my_index > my_index) {
                out = right.my_index - my_index;
                out *= -1;
            } else {
                out = my_index - right.my_index;
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
