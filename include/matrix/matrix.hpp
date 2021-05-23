#ifndef MATRIX_H
#define MATRIX_H

namespace bioc {

class matrix {
public:
    virtual ~matrix() {}

    virtual size_t nrow() const = 0;

    virtual size_t ncol() const = 0;
};

template <typename T>
class typed_matrix : public matrix {
public:
    ~typed_matrix() {}

    virtual const T* get_row(size_t, T*, size_t=0, size_t=-1) const = 0;

    virtual const T* get_column(size_t, T*, size_t=0, size_t=-1) const = 0;

    virtual void set_row(size_t, const T*, size_t=0, size_t=-1) = 0;

    virtual void set_column(size_t, const T*, size_t=0, size_t=-1) = 0;

    enum element_type { _double, _int32_t, _string };

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

using numeric_matrix = typed_matrix<double>;

}

#endif
