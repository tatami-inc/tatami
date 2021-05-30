#ifndef DELAYED_SUBSET_OP
#define DELAYED_SUBSET_OP

#include "../matrix/typed_matrix.hpp"
#include <algorithm>
#include <memory>

namespace bioc {

template<typename T, int MARGIN, class V = std::vector<size_t>, typename IDX = int>
class DelayedSubsetOp : public typed_matrix<T, IDX> {
public:
    DelayedSubsetOp(std::shared_ptr<const typed_matrix<T, IDX> > p, const V& idx) : mat(p), indices(idx) {}

    DelayedSubsetOp(std::shared_ptr<const typed_matrix<T, IDX> > p, V&& idx) : mat(p), indices(idx) {}

    ~DelayedSubsetOp() {}

public:
    const T* get_row(size_t r, T* buffer, size_t start=0, size_t end=-1, workspace* work=NULL) const {
        if constexpr(MARGIN==1) {
            end = std::min(end, this->ncol());
            subset_expanded<true>(r, buffer, start, end, work);
            return buffer;
        } else {
            return mat->get_row(indices[r], buffer, start, end, work);
        }
    }

    const T* get_column(size_t c, T* buffer, size_t start=0, size_t end=-1, workspace* work=NULL) const {
        if constexpr(MARGIN==1) {
            return mat->get_column(indices[c], buffer, start, end, work);
        } else {
            end = std::min(end, this->nrow());
            subset_expanded<false>(c, buffer, start, end, work);
            return buffer;
        }
    }

    using typed_matrix<T, IDX>::get_column;

    using typed_matrix<T, IDX>::get_row;

public:
    sparse_range<T, IDX> get_sparse_row(size_t r, T* out_values, IDX* out_indices, size_t start=0, size_t end=-1, workspace* work=NULL) const {
        if constexpr(MARGIN==1) {
            end = std::min(end, this->ncol());
            auto total = subset_sparse<true>(r, out_values, out_indices, start, end, work);
            return sparse_range<T, IDX>(total, out_values, out_indices);
        } else {
            return mat->get_sparse_row(indices[r], out_values, out_indices, start, end, work);
        }
    }

    sparse_range<T, IDX> get_sparse_column(size_t c, T* out_values, IDX* out_indices, size_t start=0, size_t end=-1, workspace* work=NULL) const {
        if constexpr(MARGIN==1) {
            return mat->get_sparse_column(indices[c], out_values, out_indices, start, end, work);
        } else {
            end = std::min(end, this->nrow());
            auto total = subset_sparse<false>(c, out_values, out_indices, start, end, work);
            return sparse_range<T, IDX>(total, out_values, out_indices);
        }
    }

    using typed_matrix<T, IDX>::get_sparse_column;

    using typed_matrix<T, IDX>::get_sparse_row;

public:
    size_t nrow() const {
        if constexpr(MARGIN==0) {
            return indices.size();
        } else {
            return mat->nrow();
        }
    }
    
    size_t ncol() const {
        if constexpr(MARGIN==0) {
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

private:
    std::shared_ptr<const typed_matrix<T, IDX> > mat;
    V indices;

    template<bool ROW>
    void subset_expanded(size_t r, T* buffer, size_t start, size_t end, workspace* work) const {
        while (start < end) {
            auto original = start;
            auto previdx = indices[start];
            ++start;
            while (start < end && indices[start] == previdx + 1) {
                previdx = indices[start];
                ++start;
            }

            const T* ptr = NULL;
            size_t n = start - original;
            previdx = indices[original];
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
        return;
    }

    template<bool ROW>
    size_t subset_sparse(size_t r, T* out_values, IDX* out_indices, size_t start, size_t end, workspace* work) const {
        size_t total = 0;
        while (start < end) {
            auto original = start;
            auto previdx = indices[start];
            ++start;
            while (start < end && indices[start] == previdx + 1) {
                previdx = indices[start];
                ++start;
            }

            size_t n = start - original;
            previdx = indices[original];
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

        return total;
    }
};

}

#endif
