#ifndef TATAMI_COMPRESSED_SPARSE_MATRIX_H
#define TATAMI_COMPRESSED_SPARSE_MATRIX_H

#include "Matrix.hpp"
#include "has_data.hpp"
#include "SparseRange.hpp"

#include <vector>
#include <algorithm>
#include <memory>
#include <utility>

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
 * Methods should be available for `size()`, `begin()`, `end()` and `[]`.
 * If a method is available for `data()` that returns a `const T*`, it will also be used.
 * @tparam V Vector class used to store the row/column indices internally.
 * This does not necessarily have to contain `IDX`, as long as the type is convertible to `IDX`.
 * Methods should be available for `size()`, `begin()`, `end()` and `[]`.
 * If a method is available for `data()` that returns a `const IDX*`, it will also be used.
 * @tparam W Vector class used to store the column/row index pointers.
 * Methods should be available for `size()`, `begin()`, `end()` and `[]`.
 */
template<bool ROW, typename T, typename IDX = int, class U = std::vector<T>, class V = std::vector<IDX>, class W = std::vector<size_t> >
class CompressedSparseMatrix : public Matrix<T, IDX> {
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
    CompressedSparseMatrix(size_t nr, size_t nc, U vals, V idx, W ptr, bool check=true) : nrows(nr), ncols(nc), values(std::move(vals)), indices(std::move(idx)), indptrs(std::move(ptr)) {
        check_values(check); 
        return;
    }

public:
    size_t nrow() const { return nrows; }

    size_t ncol() const { return ncols; }

    /**
     * @return `true`.
     */
    bool sparse() const { return true; }

    /**
     * @return `true` if `ROW = true` (for `CompressedSparseRowMatrix` objects), otherwise returns `false` (for `CompressedSparseColumnMatrix` objects).
     */
    bool prefer_rows() const { return ROW; }

public:
    const T* row(size_t r, T* buffer, size_t first, size_t last, Workspace* work=nullptr) const {
        if constexpr(ROW) {
            primary_dimension_expanded(r, first, last, this->ncols, buffer, 0);
        } else {
            secondary_dimension_expanded(r, first, last, work, buffer, 0);
        }
        return buffer;
    }

    const T* column(size_t c, T* buffer, size_t first, size_t last, Workspace* work=nullptr) const {
        if constexpr(ROW) {
            secondary_dimension_expanded(c, first, last, work, buffer, 0);
        } else {
            primary_dimension_expanded(c, first, last, this->nrows, buffer, 0);
        }
        return buffer;
    }

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::column;

public:
    /**
     * @copydoc Matrix::sparse_row()
     */
    SparseRange<T, IDX> sparse_row(size_t r, T* vbuffer, IDX* ibuffer, size_t first, size_t last, Workspace* work=nullptr, bool sorted=true) const {
        // It's always sorted anyway, no need to pass along 'sorted'.
        if constexpr(ROW) {
            return primary_dimension_raw(r, first, last, this->ncols, vbuffer, ibuffer);
        } else {
            return secondary_dimension_raw(r, first, last, work, vbuffer, ibuffer); 
        }
    }

    /**
     * @copydoc Matrix::sparse_column()
     */
    SparseRange<T, IDX> sparse_column(size_t c, T* vbuffer, IDX* ibuffer, size_t first, size_t last, Workspace* work=nullptr, bool sorted=true) const {
        // It's always sorted anyway, no need to pass along 'sorted'.
        if constexpr(ROW) {
            return secondary_dimension_raw(c, first, last, work, vbuffer, ibuffer); 
        } else {
            return primary_dimension_raw(c, first, last, this->nrows, vbuffer, ibuffer);
        }
    }

    using Matrix<T, IDX>::sparse_row;

    using Matrix<T, IDX>::sparse_column;

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

private:
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

    SparseRange<T, IDX> primary_dimension_raw(size_t i, size_t first, size_t last, size_t otherdim, T* out_values, IDX* out_indices) const {
        auto obtained = primary_dimension(i, first, last, otherdim);
        SparseRange<T, IDX> output(obtained.second);

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

public:
    /**
     * @param row Should a workspace be created for row-wise extraction?
     *
     * @return If `row == ROW`, a null pointer as no workspace is required for extraction along the preferred dimension.
     * Otherwise, a shared pointer to a `Workspace` object is returned.
     *
     * Extraction with a workspace is most efficient for accessing consecutive increasing indices,
     * due to the use of a caching mechanism for the index pointers and indices.
     * Access to increasing non-consecutive indices is the next-most efficient,
     * then decreasing consecutive indices,
     * then decreasing non-consecutive indices,
     * and finally random indices, which is probably about the same as not using a workspace at all.
     */
    std::shared_ptr<Workspace> new_workspace (bool row) const {
        if (row == ROW) {
            return nullptr;
        } else {
            return std::shared_ptr<Workspace>(new CompressedSparseWorkspace(max_secondary_index(), indices, indptrs));
        }
    }

    struct CompressedSparseWorkspace : public Workspace {
        CompressedSparseWorkspace(size_t max_index, const V& idx, const W& idp) : 
            previous_request(idp.size() - 1),
            current_indptrs(idp.begin(), idp.begin() + idp.size() - 1), // all but the last.
            current_indices(idp.size() - 1)
        {
            /* Here, the general idea is to store a local copy of the actual
             * row indices (for CSC matrices; column indices, for CSR matrices)
             * so that we don't have to keep on doing cache-unfriendly look-ups
             * for the indices based on the pointers that we do have. This assumes
             * that the density is so low that updates to the local indices are
             * rare relative to the number of comparisons to those same indices.
             * Check out the `secondary_dimension()` function for how this is used.
             */
            for (size_t i = 0; i < idp.size() - 1; ++i) {
                current_indices[i] = (idp[i] < idp[i+1] ? idx[idp[i]] : max_index);
            }
            return;
        } 

        std::vector<size_t> previous_request; // the last request for each column.
        std::vector<typename std::remove_reference<decltype(std::declval<W>()[0])>::type> current_indptrs; // the current position of the pointer
        std::vector<typename std::remove_reference<decltype(std::declval<V>()[0])>::type> current_indices; // the current index being pointed to
    };

private:
    size_t max_secondary_index() const {
        if constexpr(ROW) {
            return ncols;
        } else {
            return nrows;
        }
    }

    template<class STORE>
    void secondary_dimension(IDX i, size_t first, size_t last, Workspace* work, STORE& output) const {
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
            CompressedSparseWorkspace& worker = *(dynamic_cast<CompressedSparseWorkspace*>(work));
            auto max_index = max_secondary_index();

            for (size_t current = first; current < last; ++current) {
                auto& prev_i = worker.previous_request[current];
                auto& curptr = worker.current_indptrs[current];
                auto& curdex = worker.current_indices[current];

                if (i == prev_i + 1) {
                    if (curdex < i) { // if true, this implies that curptr < indptrs[current + 1], provided i < max_index.
                        ++curptr;
                        curdex = (curptr < indptrs[current + 1] ? indices[curptr] : max_index);
                    }
                } else if (i + 1 == prev_i) {
                    if (curptr != indptrs[current] && indices[curptr - 1] >= i) {
                        --curptr;
                        curdex = indices[curptr];
                    }
                } else if (i > prev_i) {
                    if (curdex < i) { // same implication as above.
                        curptr = std::lower_bound(indices.begin() + curptr, indices.begin() + indptrs[current + 1], i) - indices.begin();
                        curdex = (curptr < indptrs[current + 1] ? indices[curptr] : max_index);
                    }
                } else if (i < prev_i) { 
                    if (curptr != indptrs[current]) {
                        curptr = std::lower_bound(indices.begin() + indptrs[current], indices.begin() + curptr, i) - indices.begin();
                        curdex = indices[curptr];
                    }
                }

                prev_i = i;
                if (curdex == i) { // assuming i < max_index, of course.
                    output.add(current, values[curptr]);
                }
            }
        }
        return;
    }

private:
    struct raw_store {
        T* out_values;
        IDX* out_indices;
        size_t n = 0;
        void add(IDX i, T val) {
            ++n;
            *out_indices = i;
            ++out_indices;
            *out_values = val;
            ++out_values;
            return;
        }
    };

    SparseRange<T, IDX> secondary_dimension_raw(IDX i, size_t first, size_t last, Workspace* work, T* out_values, IDX* out_indices) const {
        raw_store store;
        store.out_values = out_values;
        store.out_indices = out_indices;
        secondary_dimension(i, first, last, work, store);
        return SparseRange<T, IDX>(store.n, out_values, out_indices);
    }

    struct expanded_store {
        T* out_values;
        size_t first;
        void add(size_t i, T val) {
            out_values[i - first] = val;
            return;
        }
    };

    void secondary_dimension_expanded(IDX i, size_t first, size_t last, Workspace* work, T* out_values, T empty) const {
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
