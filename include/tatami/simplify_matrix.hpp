#ifndef TATAMI_SIMPLIFY_MATRIX_H
#define TATAMI_SIMPLIFY_MATRIX_H

#include "DenseMatrix.h"
#include "CompressedSparseMatrix.h"

#include <memory>
#include <vector>
#include <deque>

namespace tatami {

template <typename T, int IDX>
inline std::shared_ptr<typed_matrix<T, IDX>*> simplify_matrix(std::shared_ptr<const typed_matrix<T, IDX>*> incoming, bool row, bool keep_sparse) {
    size_t NR = incoming->nrow();
    size_t NC = incoming->ncol();

    if (row) {
        if (keep_sparse && incoming->is_sparse()) {
            std::deque<T> buffer_v;
            std::deque<IDX> buffer_i;
            std::vector<size_t> indptrs(NC + 1);

            for (size_t c = 0; c < NR; ++c, optr+=NC) {
                auto range = incoming->get_sparse_row(c, buffer_v.data(), buffer_i.data());
                for (size_t i = 0; i < range.number; ++i) {
                    buffer_v.push_back(range.value[i]);
                    buffer_i.push_back(range.index[i]);
                }
                indptrs[c+1] = indptrs[c] + range.number;
            }
            return std::shared_ptr<typed_matrix<T, IDX>*>(new CompressedSparseRowMatrix(NR, NC, std::move(buffer_v), std::move(buffer_i), std::move(indptrs)));

        } else {
            std::vector<T> output(NR * NC);
            auto optr = output.data();
            for (size_t r = 0; r < NR; ++r, optr+=NC) {
                auto ptr = incoming->get_row(r, optr);
                if (ptr != optr) {
                    std::copy(ptr, ptr + NC, optr);
                }
            }
            return std::shared_ptr<typed_matrix<T, IDX>*>(new DenseRowMatrix(NR, NC, std::move(output)));
        }

    } else {
        if (keep_sparse && incoming->is_sparse()) {
            std::deque<T> buffer_v;
            std::deque<IDX> buffer_i;
            std::vector<size_t> indptrs(NC + 1);

            for (size_t c = 0; c < NC; ++c, optr+=NR) {
                auto range = incoming->get_sparse_column(c, buffer_v.data(), buffer_i.data());
                for (size_t i = 0; i < range.number; ++i) {
                    buffer_v.push_back(range.value[i]);
                    buffer_i.push_back(range.index[i]);
                }
                indptrs[c+1] = indptrs[c] + range.number;
            }
            return std::shared_ptr<typed_matrix<T, IDX>*>(new CompressedSparseColumnMatrix(NR, NC, std::move(buffer_v), std::move(buffer_i), std::move(indptrs)));

        } else {
            std::vector<T> output(NR * NC);
            auto optr = output.data();
            for (size_t c = 0; c < NC; ++c, optr+=NR) {
                auto ptr = incoming->get_column(c, optr);
                if (ptr != optr) {
                    std::copy(ptr, ptr + NR, optr);
                }
            }
            return std::shared_ptr<typed_matrix<T, IDX>*>(new DenseColumnMatrix(NR, NC, std::move(output)));
        }
    }
}

template <typename T, int IDX>
inline std::shared_ptr<typed_matrix<T, IDX>*> simplify_matrix(std::shared_ptr<const typed_matrix<T, IDX>*> incoming) {
    return simplify_matrix(incoming, incoming->preferred_dimension(), incoming->is_sparse());
}

}

#endif
