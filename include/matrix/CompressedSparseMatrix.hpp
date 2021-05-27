#ifndef COMPRESSED_SPARSE_MATRIX_H
#define COMPRESSED_SPARSE_MATRIX_H

#include "typed_matrix.hpp"
#include "../utils/is_contiguous.hpp"
#include "../utils/sparse_range.hpp"

#include <vector>
#include <algorithm>

namespace bioc {

template<bool ROW, typename T, typename IDX, class U = std::vector<T>, class V = std::vector<IDX>, class W = std::vector<size_t> >
class CompressedSparseMatrix : public typed_matrix<T, IDX> {
public: 
    CompressedSparseMatrix(size_t nr, size_t nc, const U& vals, const V& idx, const W& ptr, bool check=true) : nrows(nr), ncols(nc), values(vals), indices(idx), indptrs(ptr) {
        check_values(check); 
        return;
    }

    CompressedSparseMatrix(size_t nr, size_t nc, U&& vals, V&& idx, W&& ptr, bool check=true) : nrows(nr), ncols(nc), values(vals), indices(idx), indptrs(ptr) {
        check_values(check); 
        return;
    }

    ~CompressedSparseMatrix() {}

    size_t nrow() const { return nrows; }

    size_t ncol() const { return ncols; }

    const T* get_row(size_t i, T* out_values, size_t first=0, size_t last=-1, workspace* work=NULL) const {
        last = std::min(last, this->ncols);
        if constexpr(ROW) {
            get_primary_dimension_expanded(i, first, last, this->ncols, out_values, 0);
        } else {
            get_secondary_dimension_expanded(i, first, last, work, out_values, 0);
        }
        return out_values;
    }

    const T* get_column(size_t i, T* out_values, size_t first=0, size_t last=-1, workspace* work=NULL) const {
        last = std::min(last, this->nrows);
        if constexpr(ROW) {
            get_secondary_dimension_expanded(i, first, last, work, out_values, 0);
        } else {
            get_primary_dimension_expanded(i, first, last, this->nrows, out_values, 0);
        }
        return out_values;
    }
    
    sparse_range<T, IDX> get_sparse_row(size_t i, T* out_values, IDX* out_indices, size_t first=0, size_t last=-1, workspace* work=NULL) const {
        last = std::min(last, this->ncols);
        if constexpr(ROW) {
            return get_primary_dimension_raw(i, first, last, this->ncols, out_values, out_indices);
        } else {
            return get_secondary_dimension_raw(i, first, last, work, out_values, out_indices); 
        }
    }

    sparse_range<T, IDX> get_sparse_column(size_t i, T* out_values, IDX* out_indices, size_t first=0, size_t last=-1, workspace* work=NULL) const {
        last = std::min(last, this->nrows);
        if constexpr(ROW) {
            return get_secondary_dimension_raw(i, first, last, work, out_values, out_indices); 
        } else {
            return get_primary_dimension_raw(i, first, last, this->nrows, out_values, out_indices);
        }
    }

    struct compressed_sparse_workspace : public workspace {
    public:
        compressed_sparse_workspace(const V& idx, const W& idp) : indices(idx), indptrs(idp), 
            offsets(idp), prev_i(0), prev_first(0), prev_last(indptrs.size() - 1) {}

        void update_indices(IDX i, size_t first, size_t last) {
            /* If left/right slice are not equal to what is stored, we reset the indices,
             * so that the code below will know to recompute them. It's too much effort
             * to try to figure out exactly which columns need recomputing; just do them all.
             */
            if (first != prev_first || last != prev_last) {
                std::copy(indptrs.begin(), indptrs.end(), offsets.begin());
                prev_i = 0;
                prev_first = 0;
                prev_last = indptrs.size() - 1;
            }

            // No change necessary.
            if (i == prev_i) { 
                return; 
            } 

            auto pIt = indptrs.begin() + first;
            if (i == prev_i + 1) {
                ++pIt; // points to the first-past-the-end element, at any given 'c'.
                for (size_t c = first; c < last; ++c, ++pIt) {
                    auto& curdex = offsets[c];
                    if (curdex != *pIt && indices[curdex] < i) { 
                        ++curdex;
                    }
                }
            } else if (i + 1 == prev_i) {
                for (size_t c = first; c < last; ++c, ++pIt) {
                    auto& curdex = offsets[c];
                    if (curdex != *pIt && indices[curdex-1] >= i) {
                        --curdex;
                    }
                }

            } else { 
                if (i > prev_i) {
                    ++pIt; // points to the first-past-the-end element, at any given 'c'.
                    for (size_t c = first; c < last; ++c, ++pIt) { 
                        auto& curdex = offsets[c];
                        curdex = std::lower_bound(indices.begin() + curdex, indices.begin() + *pIt, i) - indices.begin();
                    }
                } else { 
                    for (size_t c = first; c < last; ++c, ++pIt) {
                        auto& curdex = offsets[c];
                        curdex = std::lower_bound(indices.begin() + *pIt, indices.begin() + curdex, i) - indices.begin();
                    }
                }
            }

            prev_i = i;
            prev_first = first;
            prev_last = last;
            return;
        }

        const W& get_offsets() const { 
            return offsets;
        }
    private:
        const V& indices;
        const W& indptrs; 
        W offsets;
        size_t prev_i, prev_first, prev_last;
    };

    workspace* create_workspace () const {
        workspace* output = new compressed_sparse_workspace(indices, indptrs);
        return output;
    }

    bool is_sparse() { return true; }

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

    std::pair<size_t, size_t> get_primary_dimension(size_t i, size_t first, size_t last, size_t otherdim) const {
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

    sparse_range<T, IDX> get_primary_dimension_raw(size_t i, size_t first, size_t last, size_t otherdim, T* out_values, IDX* out_indices) const {
        auto obtained = get_primary_dimension(i, first, last, otherdim);
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

    void get_primary_dimension_expanded(size_t i, size_t first, size_t last, size_t otherdim, T* out_values, T empty) const {
        std::fill(out_values, out_values + (last - first), empty);
        auto obtained = get_primary_dimension(i, first, last, otherdim);
        auto vIt = values.begin() + obtained.first;
        auto iIt = indices.begin() + obtained.first;
        for (size_t x = 0; x < obtained.second; ++x, ++vIt, ++iIt) {
            out_values[*iIt - first] = *vIt;            
        }
        return;
    }

    template<class STORE>
    void get_secondary_dimension(IDX i, size_t first, size_t last, workspace* work, STORE& output) const {
        if (work == NULL) {
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
            compressed_sparse_workspace& worker = *(dynamic_cast<compressed_sparse_workspace*>(work));
            worker.update_indices(i, first, last);
            auto& new_indptrs = worker.get_offsets();

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

    sparse_range<T, IDX> get_secondary_dimension_raw(IDX i, size_t first, size_t last, workspace* work, T* out_values, IDX* out_indices) const {
        raw_store store;
        store.out_values = out_values;
        store.out_indices = out_indices;
        get_secondary_dimension(i, first, last, work, store);
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

    void get_secondary_dimension_expanded(IDX i, size_t first, size_t last, workspace* work, T* out_values, T empty) const {
        std::fill(out_values, out_values + (last - first), empty);
        expanded_store store;
        store.out_values = out_values;
        store.first = first;
        get_secondary_dimension(i, first, last, work, store);
        return;
    }
};

template<typename T, typename IDX, class U = std::vector<T>, class V = std::vector<IDX>, class W = std::vector<size_t> >
using CompressedSparseColumnMatrix = CompressedSparseMatrix<false, T, IDX, U, V, W>;

template<typename T, typename IDX, class U = std::vector<T>, class V = std::vector<IDX>, class W = std::vector<size_t> >
using CompressedSparseRowMatrix = CompressedSparseMatrix<true, T, IDX, U, V, W>;

}

#endif
