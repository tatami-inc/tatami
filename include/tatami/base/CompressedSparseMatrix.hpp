#ifndef TATAMI_COMPRESSED_SPARSE_MATRIX_H
#define TATAMI_COMPRESSED_SPARSE_MATRIX_H

#include "typed_matrix.hpp"
#include "has_data.hpp"
#include "sparse_range.hpp"

#include <vector>
#include <algorithm>

/**
 * @file CompressedSparseMatrix.hpp
 *
 * Compressed sparse matrix representation, with `typedef`s for the usual row and column formats.
 */

namespace tatami {

/**
 * @brief Compressed sparse matrix representation.
 *
 * @tparam ROW Whether this is a compressed sparse row representation.
 * If `false`, a compressed sparse column representation is used instead.
 * @tparam T Type of the matrix values.
 * @tparam IDX Type of the row/column indices.
 * @tparam U Vector class used to store the matrix values internally.
 * This does not necessarily have to contain `T`, as long as the type is convertible to `T`.
 * Methods should be available for `size()`, `begin()` and `end()`.
 * If a method is available for `data()` that returns a `const T*`, it will also be used.
 * @tparam V Vector class used to store the row/column indices internally.
 * This does not necessarily have to contain `IDX`, as long as the type is convertible to `IDX`.
 * Methods should be available for `size()`, `begin()` and `end()`.
 * If a method is available for `data()` that returns a `const IDX*`, it will also be used.
 * @tparam W Vector class used to store the column/row index pointers.
 * Methods should be available for `size()`, `begin()` and `end()`.
 */
template<bool ROW, typename T, typename IDX = int, class U = std::vector<T>, class V = std::vector<IDX>, class W = std::vector<size_t> >
class CompressedSparseMatrix : public typed_matrix<T, IDX> {
public:
    /**
     * @param nr Number of rows.
     * @param nc Number of columns.
     * @param vals Vector of non-zero elements.
     * @param idx Vector of row indices (if `ROW=false`) or column indices (if `ROW=true`) for the non-zero elements.
     * @param ptr Vector of index pointers.
     * @param check Should the input vectors be checked for validity?
     *
     * If `check=true`, the constructor will check that `vals` and `idx` have the same length;
     * `ptr` is ordered with first and last values set to 0 and the number of non-zero elements, respectively;
     * and `idx` is ordered within each interval defined by successive elements of `ptr`.
     */
    CompressedSparseMatrix(size_t nr, size_t nc, const U& vals, const V& idx, const W& ptr, bool check=true) : nrows(nr), ncols(nc), values(vals), indices(idx), indptrs(ptr) {
        check_values(check); 
        return;
    }

    /**
     * @param nr Number of rows.
     * @param nc Number of columns.
     * @param vals Vector of non-zero elements.
     * @param idx Vector of row indices (if `ROW=false`) or column indices (if `ROW=true`) for the non-zero elements.
     * @param ptr Vector of index pointers.
     * @param check Should the input vectors be checked for validity?
     *
     * If `check=true`, the constructor will check that `vals` and `idx` have the same length;
     * `ptr` is ordered with first and last values set to 0 and the number of non-zero elements, respectively;
     * and `idx` is ordered within each interval defined by successive elements of `ptr`.
     */
    CompressedSparseMatrix(size_t nr, size_t nc, U&& vals, V&& idx, W&& ptr, bool check=true) : nrows(nr), ncols(nc), values(vals), indices(idx), indptrs(ptr) {
        check_values(check); 
        return;
    }

    ~CompressedSparseMatrix() {}

public:
    size_t nrow() const { return nrows; }

    size_t ncol() const { return ncols; }

    /**
     * @return If `row == ROW`, a null pointer as no workspace is required for extraction along the preferred dimension.
     * Otherwise, a shared pointer to a `workspace` object is returned.
     *
     * @param row Should a workspace be created for row-wise extraction?
     */
    workspace_ptr new_workspace (bool row) const {
        if (row == ROW) {
            return nullptr;
        } else {
            return std::shared_ptr<workspace>(new compressed_sparse_workspace(indices, indptrs));
        }
    }

    /**
     * @return `true`.
     */
    bool sparse() const { return true; }

    /**
     * @return `true` if `ROW = true` (for `CompressedSparseRowMatrix` objects), otherwise returns `false` (for `CompressedSparseColumnMatrix` objects).
     */
    bool prefer_rows() const { return ROW; }

public:
    const T* row(size_t r, T* buffer, size_t first, size_t last, const workspace_ptr& work=nullptr) const {
        if constexpr(ROW) {
            primary_dimension_expanded(r, first, last, this->ncols, buffer, 0);
        } else {
            secondary_dimension_expanded(r, first, last, work, buffer, 0);
        }
        return buffer;
    }

    const T* column(size_t c, T* buffer, size_t first, size_t last, const workspace_ptr& work=nullptr) const {
        if constexpr(ROW) {
            secondary_dimension_expanded(c, first, last, work, buffer, 0);
        } else {
            primary_dimension_expanded(c, first, last, this->nrows, buffer, 0);
        }
        return buffer;
    }

    using typed_matrix<T, IDX>::row;

    using typed_matrix<T, IDX>::column;

public:
    /**
     * @copydoc typed_matrix::sparse_row()
     */
    sparse_range<T, IDX> sparse_row(size_t r, T* vbuffer, IDX* ibuffer, size_t first, size_t last, const workspace_ptr& work=nullptr) const {
        if constexpr(ROW) {
            return primary_dimension_raw(r, first, last, this->ncols, vbuffer, ibuffer);
        } else {
            return secondary_dimension_raw(r, first, last, work, vbuffer, ibuffer); 
        }
    }

    /**
     * @copydoc typed_matrix::sparse_column()
     */
    sparse_range<T, IDX> sparse_column(size_t c, T* vbuffer, IDX* ibuffer, size_t first, size_t last, const workspace_ptr& work=nullptr) const {
        if constexpr(ROW) {
            return secondary_dimension_raw(c, first, last, work, vbuffer, ibuffer); 
        } else {
            return primary_dimension_raw(c, first, last, this->nrows, vbuffer, ibuffer);
        }
    }

    using typed_matrix<T, IDX>::sparse_row;

    using typed_matrix<T, IDX>::sparse_column;

public:
    struct compressed_sparse_workspace : public workspace {
    public:
        compressed_sparse_workspace(const V& idx, const W& idp) : indices(idx), 
                                                                  indptrs(idp), 
                                                                  curptrs(idp), 
                                                                  previous(indptrs.size() - 1) {} 

        void update_indices(IDX i, size_t first, size_t last) {
            for (size_t current = first; current < last; ++current) {
                auto& prev_i = previous[current];
                if (i == prev_i) {
                    continue;
                }

                auto& curdex = curptrs[current];
                if (i == prev_i + 1) {
                    if (curdex != indptrs[current+1] && indices[curdex] < i) { 
                        ++curdex;
                    }
                } else if (i + 1 == prev_i) {
                    if (curdex != indptrs[current] && indices[curdex-1] >= i) {
                        --curdex;
                    }
                } else if (i > prev_i) {
                    curdex = std::lower_bound(indices.begin() + curdex, indices.begin() + indptrs[current+1], i) - indices.begin();
                } else if (i < prev_i) { 
                    curdex = std::lower_bound(indices.begin() + indptrs[current], indices.begin() + curdex, i) - indices.begin();
                }

                prev_i = i;
            }

            return;
        }

        const W& offsets() const { 
            return curptrs;
        }
    private:
        const V& indices;
        const W& indptrs; 
        W curptrs;
        std::vector<size_t> previous;
    };

private:
    size_t nrows, ncols;
    U values;
    V indices;
    W indptrs;

    void check_values(bool check) {
        if (!check) {
            return;
        }

        if (values.size() != indices.size()) {
            throw std::runtime_error("'values' and 'indices' should be of the same length");
        }

        if (ROW) {
            if (indptrs.size() != nrows + 1){
                throw std::runtime_error("length of 'indptrs' should be equal to 'nrows + 1'");
            }
        } else {
            if (indptrs.size() != ncols + 1){
                throw std::runtime_error("length of 'indptrs' should be equal to 'ncols + 1'");
            }
        }

        if (indptrs[0] != 0) {
            throw std::runtime_error("first element of 'indptrs' should be zero");
        }
        if (indptrs[indptrs.size() - 1] != indices.size()) {
            throw std::runtime_error("last element of 'indptrs' should be equal to length of 'indices'");
        }

        size_t counter = 0;
        for (size_t i = 1; i < indptrs.size(); ++i) {
            if (indptrs[i] < indptrs[i-1]) {
                throw std::runtime_error("'indptrs' should be in increasing order");
            }

            if (counter < indices.size()) {
                auto previous = indices[counter];
                ++counter;
                while (counter < indptrs[i]) {
                    if (previous >= indices[counter]) {
                        if (ROW) {
                            throw std::runtime_error("'indices' should be strictly increasing within each row");
                        } else {
                            throw std::runtime_error("'indices' should be strictly increasing within each column");
                        }
                    }
                    ++counter;
                }
            }
        }

        return;
    }

    std::pair<size_t, size_t> primary_dimension(size_t i, size_t first, size_t last, size_t otherdim) const {
        const auto pstart = indptrs[i];
        auto iIt = indices.begin() + pstart, 
             eIt = indices.begin() + indptrs[i + 1]; 
        auto xIt = values.begin() + pstart;

        if (first) { // Jumping ahead if non-zero.
            iIt = std::lower_bound(iIt, eIt, first);
        } 

        if (last != otherdim) { // Jumping to last element.
            eIt = std::lower_bound(iIt, eIt, last);
        }

        return std::make_pair(iIt - indices.begin(), eIt - iIt);
    }

    sparse_range<T, IDX> primary_dimension_raw(size_t i, size_t first, size_t last, size_t otherdim, T* out_values, IDX* out_indices) const {
        auto obtained = primary_dimension(i, first, last, otherdim);
        sparse_range<T, IDX> output(obtained.second);

        if constexpr(has_data<T, U>::value) {
            output.value = values.data() + obtained.first;
        } else {
            auto vIt = values.begin() + obtained.first;
            std::copy(vIt, vIt + obtained.second, out_values);
            output.value = out_values;
        }

        if constexpr(has_data<IDX, V>::value) {
            output.index = indices.data() + obtained.first;
        } else {
            auto iIt = indices.begin() + obtained.first;
            std::copy(iIt, iIt + obtained.second, out_indices);
            output.index = out_indices;
        }

        return output;
    }

    void primary_dimension_expanded(size_t i, size_t first, size_t last, size_t otherdim, T* out_values, T empty) const {
        std::fill(out_values, out_values + (last - first), empty);
        auto obtained = primary_dimension(i, first, last, otherdim);
        auto vIt = values.begin() + obtained.first;
        auto iIt = indices.begin() + obtained.first;
        for (size_t x = 0; x < obtained.second; ++x, ++vIt, ++iIt) {
            out_values[*iIt - first] = *vIt;            
        }
        return;
    }

    template<class STORE>
    void secondary_dimension(IDX i, size_t first, size_t last, const workspace_ptr& work, STORE& output) const {
        if (work == nullptr) {
            for (size_t c = first; c < last; ++c) { 
                auto start = indices.begin() + indptrs[c];
                auto end = indices.begin() + indptrs[c+1];
                auto iIt = std::lower_bound(start, end, i);
                if (iIt != end && *iIt == i) { 
                    auto offset = iIt - indices.begin();
                    output.add(c, values[offset]);
                }
            }
        } else {
            compressed_sparse_workspace& worker = *(dynamic_cast<compressed_sparse_workspace*>(work.get()));
            worker.update_indices(i, first, last);
            auto& new_indptrs = worker.offsets();

            auto pIt = indptrs.begin() + first + 1; // Points to first-past-the-end for each 'c'.
            for (size_t c = first; c < last; ++c, ++pIt) { 
                const int idex = new_indptrs[c];
                if (idex != *pIt && indices[idex] == i) { 
                    output.add(c, values[idex]);
                }
            }
        }
        return;
    }

    struct raw_store {
        T* out_values;
        IDX* out_indices;
        size_t n = 0;
        void add(IDX i, T val) {
            ++n;
            *out_indices = i;
            *out_values = val;
            ++out_values;
            ++out_indices;
            return;
        }
    };

    sparse_range<T, IDX> secondary_dimension_raw(IDX i, size_t first, size_t last, const workspace_ptr& work, T* out_values, IDX* out_indices) const {
        raw_store store;
        store.out_values = out_values;
        store.out_indices = out_indices;
        secondary_dimension(i, first, last, work, store);
        return sparse_range<T, IDX>(store.n, out_values, out_indices);
    }

    struct expanded_store {
        T* out_values;
        size_t first;
        void add(size_t i, T val) {
            out_values[i - first] = val;
            return;
        }
    };

    void secondary_dimension_expanded(IDX i, size_t first, size_t last, const workspace_ptr& work, T* out_values, T empty) const {
        std::fill(out_values, out_values + (last - first), empty);
        expanded_store store;
        store.out_values = out_values;
        store.first = first;
        secondary_dimension(i, first, last, work, store);
        return;
    }
};

/**
 * Compressed sparse column matrix.
 * See `tatami::CompressedSparseMatrix` for details on the template parameters.
 */
template<typename T, typename IDX, class U = std::vector<T>, class V = std::vector<IDX>, class W = std::vector<size_t> >
using CompressedSparseColumnMatrix = CompressedSparseMatrix<false, T, IDX, U, V, W>;

/**
 * Compressed sparse row matrix.
 * See `tatami::CompressedSparseMatrix` for details on the template parameters.
 */
template<typename T, typename IDX, class U = std::vector<T>, class V = std::vector<IDX>, class W = std::vector<size_t> >
using CompressedSparseRowMatrix = CompressedSparseMatrix<true, T, IDX, U, V, W>;

}

#endif
