#ifndef LOAD_SPARSE_H
#define LOAD_SPARSE_H

#include <vector>
#include <memory>

#include "tatami/CompressedSparseMatrix.hpp"

template<typename V>
inline std::unique_ptr<tatami::typed_matrix<double, int> > load_matrix_as_sparse_row_matrix(size_t nr, size_t nc, const V& source) {
    // Filling the sparse row matrix.
    std::vector<double> values;
    std::vector<int> indices;
    std::vector<size_t> indptr;

    auto sIt = source.begin();
    indptr.resize(nr+1);
    for (size_t i = 0; i < nr; ++i) {
        indptr[i+1] = indptr[i];
        for (size_t j = 0; j < nc; ++j, ++sIt) {
            if (*sIt != 0) {
                values.push_back(*sIt);
                indices.push_back(j);
                ++indptr[i+1];
            }
        }
    }
    return std::unique_ptr<tatami::typed_matrix<double, int> >(new tatami::CompressedSparseRowMatrix<double, int>(nr, nc, values, indices, indptr));
}

template<typename V>
inline std::unique_ptr<tatami::typed_matrix<double, int> > load_matrix_as_sparse_column_matrix(size_t nr, size_t nc, const V& source) {
    std::vector<double> values;
    std::vector<int> indices;
    std::vector<size_t> indptr;

    indptr.resize(nc+1);
    for (size_t i = 0; i < nc; ++i) {
        indptr[i+1] = indptr[i];
        auto sIt = source.begin() + i;
        for (size_t j = 0; j < nr; ++j, sIt+=nc) {
            if (*sIt != 0) {
                values.push_back(*sIt);
                indices.push_back(j);
                ++indptr[i+1];
            }
        }
    }

    return std::unique_ptr<tatami::typed_matrix<double, int> >(new tatami::CompressedSparseColumnMatrix<double, int>(nr, nc, values, indices, indptr));
}

#endif
