#ifndef DELAYED_SUBSET_OP
#define DELAYED_SUBSET_OP

#include <algorithm>

namespace bioc {

template<typename T, int MARGIN, class V = std::vector<size_t> >
class DelayedSubsetOp : public typed_matrix<T, IDX> {
public:
    DelayedSubsetOp(std::shared_ptr<const typed_matrix<X, IDX> > p, const V& idx) : mat(p), indices(idx) {}

    DelayedSubsetOp(std::shared_ptr<const typed_matrix<X, IDX> > p, V&& idx) : mat(p), indices(idx) {}

    ~DelayedSubsetOp() {}

    const T* get_row(size_t r, T* buffer, size_t start=0, size_t end=-1, workspace* work=NULL) const {
        if constexpr(MARGIN==1) {
            end = std::min(end, this->ncol);
            return subset_expanded(r, buffer, start, end, work);
        } else {
            return mat->get_row(r, buffer, start, end, work);
        }
    }

    const T* get_column(size_t c, T* buffer, size_t start=0, size_t end=-1, workspace* work=NULL) const {
        if constexpr(MARGIN==1) {
            return mat->get_column(c, buffer, start, end, work);
        } else {
            end = std::min(end, this->nrow);
            return subset_expanded(r, buffer, start, end, work);
        }
        return buffer;
    }

    sparse_range<T, IDX> get_sparse_row(size_t r, T* out_values, IDX* out_indices, size_t start=0, size_t end=-1, workspace* work=NULL) const {
        auto raw = mat.get_sparse_row(r, out_values, out_indices, start, end, work);
        for (size_t i = 0; i < raw.number; ++i) {
            out_values[i] = operation(r, i, raw.value[i]);
        }
        return sparse_range<T, IDX>(raw.number, out_values, raw.index);
    }

    sparse_range<T, IDX> get_sparse_column(size_t c, T* out_values, IDX* out_indices, size_t start=0, size_t end=-1, workspace* work=NULL) const {
        auto raw = mat.get_sparse_column(c, out_values, out_indices, start, end, work);
        for (size_t i = 0; i < raw.number; ++i) {
            out_values[i] = operation(i, c, raw.value[i]);
        }
        return sparse_range<T, IDX>(raw.number, out_values, raw.index);
    }

    size_t nrow() const {
        if constexpr(MARGIN==1) {
            return indices.size();
        } else {
            return mat->nrow();
        }
    }
    
    size_t ncol() const {
        if constexpr(MARGIN==1) {
            return mat->ncol();
        } else {
            return indices.size();
        }
    }

    workspace* create_workspace() const {
        return mat->create_workspace();
    }

    bool is_sparse() const {
        return mat->is_sparse();
    }
public:
    std::shared_ptr<typed_matrix<T, IDX> > mat;
    V indices;

    template<bool ROW>
    const T* subset_expanded(size_t r, T* buffer, size_t start, size_t end, workspace* work) const {
        while (start < end) {
            auto previdx = indices[start];
            ++start;
            size_t n = 1;
            while (start < end && indices[start] == previdx + 1) {
                previdx = indices[start];
                ++start;
                ++n;
            }

            const T* ptr = NULL;
            if constexpr(ROW) {
                ptr = mat->get_row(r, buffer, previdx, previdx + n, work);
            } else {
                ptr = mat->get_column(r, buffer, previdx, previdx + n, work);
            }

            if (ptr != buffer) {
                std::copy(ptr, ptr + n, buffer);
            }
            buffer += n;
        }
        return buffer;
    }

    template<bool ROW>
    sparse_range<T, IDX> subset_sparse(size_t r, T* out_values, IDX* out_indices, size_t start, size_t end, workspace* work) const {
        sparse_range<T, IDX> output(0, out_values, out_indices);
        auto& total = output.number;

        while (start < end) {
            auto original = start;
            auto previdx = indices[start];
            ++start;

            while (start < end && indices[start] == previdx + 1) {
                previdx = indices[start];
                ++start;
            }

            size_t n = start - original;
            sparse_range<T, IDX> range;
            if constexpr(ROW) {
                range = mat->get_sparse_row(r, out_values, out_indices, previdx, previdx + n, work);
            } else {
                range = mat->get_sparse_column(r, out_values, out_indices, previdx, previdx + n, work);
            }

            if (out_values != range.value) {
                std::copy(range.value, range.value + range.number, out_values);
            }
            for (size_t i = 0; i < range.number; ++i) {
                out_indices[i] = range.index[i] - previdx + original;
            }

            total += range.number;
            out_indices += range.number;
            out_values += range.number;
        }

        return output;
    }
};

}

#endif
