#ifndef VECTOR_H
#define VECTOR_H

#include <string>
#include <memory>
#include <type_traits>
#include <cstdint>

namespace bioc {

struct vector {
    virtual ~vector() {}

    virtual size_t length() const = 0;

    virtual std::shared_ptr<vector> subset(size_t, size_t) const = 0;

    virtual std::shared_ptr<vector> subset(size_t, const size_t*) const = 0;

    virtual std::shared_ptr<vector> clone() const = 0;

    enum element_type { _double, _int32_t, _string };
    
    virtual element_type type() const = 0;
};

template<typename T>
struct typed_vector : public vector {
    ~typed_vector() {}

    virtual T get(size_t) const = 0;

    virtual void set(size_t, const T&) = 0;

    virtual const T* get_range(size_t, size_t, T*) const = 0;

    virtual void set_range(size_t, size_t, const T*) = 0;

    element_type type() const {
        if constexpr(std::is_same<T, double>::value) {
            return _double;
        } else if constexpr(std::is_same<T, int32_t>::value) {
            return _int32_t;
        } else if constexpr(std::is_same<T, std::string>::value) {
            return _string;
        }
    }
};

using numeric_vector = typed_vector<double>;

}

#endif
