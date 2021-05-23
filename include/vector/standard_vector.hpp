#ifndef STANDARD_VECTOR_H
#define STANDARD_VECTOR_H

#include "vector.hpp"

#include <vector>
#include <deque>
#include <algorithm>

namespace bioc {

template<typename T, class V = std::vector<T> >
struct standard_vector : public typed_vector<T> {
    standard_vector(size_t n=0) : values(n) {}

    standard_vector(V&& source) : values(source) {}

    ~standard_vector() {}

    size_t length() const {
        return values.size();
    }

    // Cloners.
    std::shared_ptr<vector> subset(size_t n, const size_t* idx) const {
        standard_vector<T, V>* output = new standard_vector<T, V>(n);
        auto& replacement = output->values;

        for (size_t i = 0; i < n; ++i) {
            replacement[i] = values[idx[i]];
        }

        return std::shared_ptr<vector>(output);
    }

    std::shared_ptr<vector> subset(size_t start, size_t end) const {
        standard_vector<T, V>* output = new standard_vector<T, V>(end - start);
        auto& replacement = output->values;
        std::copy(values.begin() + start, values.begin() + end, replacement.begin());
        return std::shared_ptr<vector>(output);
    }

    std::shared_ptr<vector> clone() const {
        standard_vector<T, V>* output = new standard_vector<T, V>(*this);
        return std::shared_ptr<vector>(output);
    }

    // Getter/setters.
    T get(size_t idx) const {
        return values[idx];
    }

    void set(size_t idx, const T& val) {
        values[idx] = val;
        return;
    }

    const T* get_range(size_t start, size_t end, T* vals) const {
        auto it = values.begin() + start;
        if constexpr(std::is_same<V, std::vector<T> >::value) {
            // We can do a sneaky thing. The philosophy is that the return value
            // will always contain something. 'vals' may or may not be filled, 
            // and you'll have to check against the return value to see if it was used.
            return it; 
        } else {
            std::copy(it, values.begin() + end, vals);
            return vals;
        }
    }

    void set_range(size_t start, size_t end, const T* vals) {
        std::copy(vals, vals + end - start, values.begin() + start);
        return;
    }
private:
    V values;
};

}

#endif
