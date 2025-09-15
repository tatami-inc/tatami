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
enum SomeNumericType {
#ifdef INT8_MAX
    I8,
#endif
#ifdef UINT8_MAX
    U8,
#endif
#ifdef INT16_MAX
    I16,
#endif
#ifdef UINT16_MAX
    U16,
#endif
#ifdef INT32_MAX
    I32,
#endif
#ifdef UINT32_MAX
    U32,
#endif
#ifdef INT64_MAX
    I64,
#endif
#ifdef UINT64_MAX
    U64,
#endif
    F32,
    F64
};

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
template<typename Value_>
class SomeNumericArray {
public:
    /**
     * @param[in] ptr Pointer to the array of interest, of run-time type specified by `type`.
     * The lifetime of the array should exceed that of the constructed `SomeNumericArray` and any of its copies.
     * @param number Length of the array pointed to by `ptr`.
     * @param type Type of the array. 
     */
    SomeNumericArray(void* ptr, const std::size_t number, const SomeNumericType type) : my_number(number), my_type(type) {
        switch(my_type) {
#ifdef INT8_MAX
            case I8:
                my_i8 = static_cast<const std::int8_t*>(ptr);
                break;
#endif
#ifdef UINT8_MAX
            case U8:
                my_u8 = static_cast<const std::uint8_t*>(ptr);
                break;
#endif
#ifdef INT16_MAX
            case I16:
                my_i16 = static_cast<const std::int16_t*>(ptr);
                break;
#endif
#ifdef UINT16_MAX
            case U16:
                my_u16 = static_cast<const std::uint16_t*>(ptr);
                break;
#endif
#ifdef INT32_MAX
            case I32:
                my_i32 = static_cast<const std::int32_t*>(ptr);
                break;
#endif
#ifdef UINT32_MAX
            case U32:
                my_u32 = static_cast<const std::uint32_t*>(ptr);
                break;
#endif
#ifdef INT64_MAX
            case I64:
                my_i64 = static_cast<const std::int64_t*>(ptr);
                break;
#endif
#ifdef UINT64_MAX
            case U64:
                my_u64 = static_cast<const std::uint64_t*>(ptr);
                break;
#endif
            case F32:
                my_f32 = static_cast<const float*>(ptr);
                break;
            case F64:
                my_f64 = static_cast<const double*>(ptr);
                break;
        }
    }

#ifdef INT8_MAX
    /**
     * @param ptr Pointer to an existing array of `int8_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const std::int8_t* const ptr, const std::size_t number) : my_i8(ptr), my_number(number), my_type(I8) {}
#endif

#ifdef UINT8_MAX
    /**
     * @param ptr Pointer to an existing array of `uint8_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const std::uint8_t* const ptr, const std::size_t number) : my_u8(ptr), my_number(number), my_type(U8) {}
#endif

#ifdef INT16_MAX
    /**
     * @param ptr Pointer to an existing array of `int16_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const std::int16_t* const ptr, const std::size_t number) : my_i16(ptr), my_number(number), my_type(I16) {}
#endif

#ifdef UINT16_MAX
    /**
     * @param ptr Pointer to an existing array of `uint16_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const std::uint16_t* const ptr, const std::size_t number) : my_u16(ptr), my_number(number), my_type(U16) {}
#endif

#ifdef INT32_MAX
    /**
     * @param ptr Pointer to an existing array of `int32_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const std::int32_t* const ptr, const std::size_t number) : my_i32(ptr), my_number(number), my_type(I32) {}
#endif

#ifdef UINT32_MAX
    /**
     * @param ptr Pointer to an existing array of `uint32_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const std::uint32_t* const ptr, const std::size_t number) : my_u32(ptr), my_number(number), my_type(U32) {}
#endif

#ifdef INT64_MAX
    /**
     * @param ptr Pointer to an existing array of `int64_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const std::int64_t* const ptr, const std::size_t number) : my_i64(ptr), my_number(number), my_type(I64) {}
#endif

#ifdef UINT64_MAX
    /**
     * @param ptr Pointer to an existing array of `uint64_t`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const std::uint64_t* const ptr, const std::size_t number) : my_u64(ptr), my_number(number), my_type(U64) {}
#endif

    /**
     * @param ptr Pointer to an existing array of `float`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const float* const ptr, const std::size_t number) : my_f32(ptr), my_number(number), my_type(F32) {}

    /**
     * @param ptr Pointer to an existing array of `double`s.
     * @param number Length of the array pointed to by `ptr`.
     */
    SomeNumericArray(const double* const ptr, const std::size_t number) : my_f64(ptr), my_number(number), my_type(F64) {}

private:
#ifdef INT8_MAX
    const std::int8_t* my_i8 = NULL;
#endif
#ifdef UINT8_MAX
    const std::uint8_t* my_u8 = NULL;
#endif
#ifdef INT16_MAX
    const std::int16_t* my_i16 = NULL;
#endif
#ifdef UINT16_MAX
    const std::uint16_t* my_u16 = NULL;
#endif
#ifdef INT32_MAX
    const std::int32_t* my_i32 = NULL;
#endif
#ifdef UINT32_MAX
    const std::uint32_t* my_u32 = NULL;
#endif
#ifdef INT64_MAX
    const std::int64_t* my_i64 = NULL;
#endif
#ifdef UINT64_MAX
    const std::uint64_t* my_u64 = NULL;
#endif
    const float* my_f32 = NULL;
    const double* my_f64 = NULL;

    std::size_t my_number;
    SomeNumericType my_type;

public:
    /**
     * @param i Positional index on the array.
     * @return Value of the `i`-th element as a `Value_`.
     */
    Value_ operator[](const std::size_t i) const {
        switch (my_type) {
#ifdef INT8_MAX
            case I8:
                return my_i8[i];
#endif
#ifdef UINT8_MAX
            case U8:
                return my_u8[i];
#endif
#ifdef INT16_MAX
            case I16:
                return my_i16[i];
#endif
#ifdef UINT16_MAX
            case U16:
                return my_u16[i];
#endif
#ifdef INT32_MAX
            case I32:
                return my_i32[i];
#endif
#ifdef UINT32_MAX
            case U32:
                return my_u32[i];
#endif
#ifdef INT64_MAX
            case I64:
                return my_i64[i];
#endif
#ifdef UINT64_MAX
            case U64:
                return my_u64[i];
#endif
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
    std::size_t size() const {
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
        Iterator(const SomeNumericArray* const parent, const std::size_t index) : my_parent(parent), my_index(index) {}

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
        std::size_t my_index;

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
        value_type operator[](const std::size_t i) const {
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
        Iterator& operator+=(const std::size_t n) {
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
        Iterator operator++(const int) {
            auto copy = *this;
            ++(*this);
            return copy;
        }

        /**
         * @param n Number of elements to move back the current iterator.
         * @return The iterator's position is moved backward by `n`, and a reference to the iterator is returned.
         */
        Iterator& operator-=(const std::size_t n) {
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
        Iterator operator--(const int) {
            auto copy = *this;
            --(*this);
            return copy;
        }

    public:
        /**
         * @param n Number of elements to advance the iterator.
         * @return A new iterator is returned at the position of the current iterator plus `n`.
         */
        Iterator operator+(const std::size_t n) const {
            return Iterator(my_parent, my_index + n);
        }

        /**
         * @param n Number of elements to move back the iterator.
         * @return A new iterator is returned at the position of the current iterator minus `n`.
         */
        Iterator operator-(const std::size_t n) const {
            return Iterator(my_parent, my_index - n);
        }

        /**
         * @param n Number of elements to advance the iterator.
         * @param it An existing `Iterator`.
         * @return A new iterator is returned at the position of `it` plus `n`.
         */
        friend Iterator operator+(const std::size_t n, const Iterator& it) {
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
