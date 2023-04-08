#ifndef TATAMI_COMPRESSED_SPARSE_MATRIX_H
#define TATAMI_COMPRESSED_SPARSE_MATRIX_H

#include "Matrix.hpp"
#include "has_data.hpp"
#include "SparseRange.hpp"

#include <vector>
#include <algorithm>
#include <memory>
#include <utility>
#include <stdexcept>

/**
 * @file CompressedSparseMatrix.hpp
 *
 * @brief Compressed sparse matrix representation. 
 *
 * `typedef`s are provided for the usual row and column formats. 
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

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::column;

    using Matrix<T, IDX>::dense_row_workspace;

    using Matrix<T, IDX>::dense_column_workspace;

    using Matrix<T, IDX>::sparse_row_workspace;

    using Matrix<T, IDX>::sparse_column_workspace;

private:
    struct SecondaryWorkspaceBase {
        typedef typename std::remove_reference<decltype(std::declval<W>()[0])>::type indptr_type;
        typedef typename std::remove_reference<decltype(std::declval<V>()[0])>::type index_type;

        SecondaryWorkspaceBase(size_t max_index, const V& idx, const W& idp, size_t start, size_t length) : 
            current_indptrs(idp.begin() + start, idp.begin() + start + length), current_indices(length) 
        {
            /* Here, the general idea is to store a local copy of the actual
             * row indices (for CSC matrices; column indices, for CSR matrices)
             * so that we don't have to keep on doing cache-unfriendly look-ups
             * for the indices based on the pointers that we do have. This assumes
             * that the density is so low that updates to the local indices are
             * rare relative to the number of comparisons to those same indices.
             * Check out the `secondary_dimension()` function for how this is used.
             */
            auto idpIt = idp.begin() + start;
            for (size_t i = 0; i < length; ++i, ++idpIt) {
                current_indices[i] = (*idpIt < *(idpIt + 1) ? idx[*idpIt] : max_index);
            }
            return;
        } 

        SecondaryWorkspaceBase(size_t max_index, const V& idx, const W& idp) : SecondaryWorkspaceBase(max_index, idx, idp, 0, idp.size() - 1) {}

        SecondaryWorkspaceBase(size_t max_index, const V& idx, const W& idp, const std::vector<IDX>& subset) : current_indptrs(subset.size()), current_indices(subset.size()) {
            for (size_t i0 = 0, end = subset.size(); i0 < end; ++i0) {
                auto i = subset[i0];
                current_indptrs[i0] = idp[i];
                current_indices[i0] = (idp[i] < idp[i + 1] ? idx[idp[i]] : max_index);
            }
            return;
        } 

        size_t previous_request = 0;
        std::vector<indptr_type> current_indptrs; // the current position of the pointer

        // The current index being pointed to, i.e., current_indices[0] <= indices[current_ptrs[0]]. 
        // If current_ptrs[0] is out of range, the current index is instead set to the maximum index
        // (e.g., max rows for CSC matrices).
        std::vector<index_type> current_indices;

        // Efficiency boosts for reverse consecutive searches. These aren't populated until they're required,
        // given that forward consecutive searches are far more common.
        // 'below_indices' holds 1-past-the-index of the non zero element immediately below 'previous_request'.
        std::vector<index_type> below_indices;
        bool store_below = false;
    };

public:
    /**
     * @cond
     */
    struct CompressedSparsePrimaryDenseWorkspace : public DenseWorkspace<ROW> {
        CompressedSparsePrimaryDenseWorkspace() = default;
    };

    struct CompressedSparseSecondaryDenseWorkspace : public DenseWorkspace<!ROW> {
        CompressedSparseSecondaryDenseWorkspace(size_t max_index, const V& idx, const W& idp) : core(max_index, idx, idp) {}
        SecondaryWorkspaceBase core;
    };

    static SecondaryWorkspaceBase& extract_secondary_core(DenseWorkspace<!ROW>* wrk) {
        return static_cast<CompressedSparseSecondaryDenseWorkspace*>(wrk)->core;
    }
    /**
     * @endcond
     */

    std::shared_ptr<DenseRowWorkspace> dense_row_workspace(const WorkspaceOptions&) const {
        if constexpr(ROW) {
            return std::shared_ptr<DenseRowWorkspace>(new CompressedSparsePrimaryDenseWorkspace());
        } else {
            return std::shared_ptr<DenseRowWorkspace>(new CompressedSparseSecondaryDenseWorkspace(nrows, indices, indptrs));
        }
    }

    std::shared_ptr<DenseColumnWorkspace> dense_column_workspace(const WorkspaceOptions&) const {
        if constexpr(ROW) {
            return std::shared_ptr<DenseColumnWorkspace>(new CompressedSparseSecondaryDenseWorkspace(ncols, indices, indptrs));
        } else {
            return std::shared_ptr<DenseColumnWorkspace>(new CompressedSparsePrimaryDenseWorkspace());
        }
    }

    const T* row(size_t r, T* buffer, DenseRowWorkspace* work) const {
        if constexpr(ROW) {
            primary_dimension_expanded(r, 0, ncols, ncols, nullptr, buffer);
        } else {
            secondary_dimension_expanded(r, 0, ncols, extract_secondary_core(work), buffer);
        }
        return buffer;
    }

    const T* column(size_t c, T* buffer, DenseColumnWorkspace* work) const {
        if constexpr(ROW) {
            secondary_dimension_expanded(c, 0, nrows, extract_secondary_core(work), buffer);
        } else {
            primary_dimension_expanded(c, 0, nrows, nrows, nullptr, buffer);
        }
        return buffer;
    }

public:
    /**
     * @cond
     */
    struct CompressedSparsePrimarySparseWorkspace : public SparseWorkspace<ROW> {
        CompressedSparsePrimarySparseWorkspace(const WorkspaceOptions& opt) : SparseWorkspace<ROW>(opt.mode, opt.sorted) {}
    };

    struct CompressedSparseSecondarySparseWorkspace : public SparseWorkspace<!ROW> {
        CompressedSparseSecondarySparseWorkspace(const WorkspaceOptions& opt, size_t max_index, const V& idx, const W& idp) : 
            SparseWorkspace<!ROW>(opt.mode, opt.sorted), core(max_index, idx, idp) {}

        SecondaryWorkspaceBase core;
    };

    static SecondaryWorkspaceBase& extract_secondary_core(SparseWorkspace<!ROW>* wrk) {
        return static_cast<CompressedSparseSecondarySparseWorkspace*>(wrk)->core;
    }
    /**
     * @endcond
     */

    std::shared_ptr<SparseRowWorkspace> sparse_row_workspace(const WorkspaceOptions& opt) const {
        if constexpr(ROW) {
            return std::shared_ptr<SparseRowWorkspace>(new CompressedSparsePrimarySparseWorkspace(opt));
        } else {
            return std::shared_ptr<SparseRowWorkspace>(new CompressedSparseSecondarySparseWorkspace(opt, nrows, indices, indptrs));
        }
    }

    std::shared_ptr<SparseColumnWorkspace> sparse_column_workspace(const WorkspaceOptions& opt) const {
        if constexpr(ROW) {
            return std::shared_ptr<SparseColumnWorkspace>(new CompressedSparseSecondarySparseWorkspace(opt, ncols, indices, indptrs));
        } else {
            return std::shared_ptr<SparseColumnWorkspace>(new CompressedSparsePrimarySparseWorkspace(opt));
        }
    }

    SparseRange<T, IDX> row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowWorkspace* work) const {
        nullify_sparse_extract_pointers(work->mode, vbuffer, ibuffer);

        // It's always sorted anyway, no need to consider along 'sorted'.
        if constexpr(ROW) {
            return primary_dimension_raw(r, 0, ncols, ncols, nullptr, vbuffer, ibuffer);
        } else {
            return secondary_dimension_raw(r, 0, ncols, extract_secondary_core(work), vbuffer, ibuffer); 
        }
    }

    SparseRange<T, IDX> column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnWorkspace* work) const {
        nullify_sparse_extract_pointers(work->mode, vbuffer, ibuffer);

        // It's always sorted anyway, no need to pass along 'sorted'.
        if constexpr(ROW) {
            return secondary_dimension_raw(c, 0, nrows, extract_secondary_core(work), vbuffer, ibuffer); 
        } else {
            return primary_dimension_raw(c, 0, nrows, nrows, nullptr, vbuffer, ibuffer);
        }
    }

public:
    /**
     * @cond
     */
    struct PrimaryBlockWorkspaceBase {
        PrimaryBlockWorkspaceBase(size_t cache_size) : cached(cache_size, std::make_pair(-1, 0)) {}
        std::vector<std::pair<size_t, size_t> > cached;
    };

    struct CompressedSparsePrimaryDenseBlockWorkspace : public DenseBlockWorkspace<ROW> {
        CompressedSparsePrimaryDenseBlockWorkspace(size_t s, size_t l, const WorkspaceOptions& opt, size_t maxdim) : 
            DenseBlockWorkspace<ROW>(s, l), core(opt.cache ? maxdim : 0) {}

        PrimaryBlockWorkspaceBase core;
    };

    struct CompressedSparseSecondaryDenseBlockWorkspace : public DenseBlockWorkspace<!ROW> {
        CompressedSparseSecondaryDenseBlockWorkspace(size_t s, size_t l, size_t max_index, const V& idx, const W& idp) : 
            DenseBlockWorkspace<!ROW>(s, l), core(max_index, idx, idp, s, l) {}

        SecondaryWorkspaceBase core;
    };

    static PrimaryBlockWorkspaceBase& extract_primary_core(DenseBlockWorkspace<ROW>* wrk) {
        return static_cast<CompressedSparsePrimaryDenseBlockWorkspace*>(wrk)->core;
    }

    static SecondaryWorkspaceBase& extract_secondary_core(DenseBlockWorkspace<!ROW>* wrk) {
        return static_cast<CompressedSparseSecondaryDenseBlockWorkspace*>(wrk)->core;
    }
    /**
     * @endcond
     */

    std::shared_ptr<DenseRowBlockWorkspace> dense_row_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const { 
        if constexpr(ROW) {
            return std::shared_ptr<DenseRowBlockWorkspace>(new CompressedSparsePrimaryDenseBlockWorkspace(start, length, opt, nrows));
        } else {
            return std::shared_ptr<DenseRowBlockWorkspace>(new CompressedSparseSecondaryDenseBlockWorkspace(start, length, nrows, indices, indptrs));
        }
    }

    std::shared_ptr<DenseColumnBlockWorkspace> dense_column_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const { 
        if constexpr(ROW) {
            return std::shared_ptr<DenseColumnBlockWorkspace>(new CompressedSparseSecondaryDenseBlockWorkspace(start, length, ncols, indices, indptrs));
        } else {
            return std::shared_ptr<DenseColumnBlockWorkspace>(new CompressedSparsePrimaryDenseBlockWorkspace(start, length, opt, ncols));
        }
    }

    const T* row(size_t r, T* buffer, DenseRowBlockWorkspace* work) const {
        size_t start = work->start, len = work->length;
        if constexpr(ROW) {
            primary_dimension_expanded(r, start, len, ncols, &(extract_primary_core(work)), buffer);
        } else {
            secondary_dimension_expanded(r, start, len, extract_secondary_core(work), buffer);
        }
        return buffer;
    }

    const T* column(size_t c, T* buffer, DenseColumnBlockWorkspace* work) const {
        size_t start = work->start, len = work->length;
        if constexpr(ROW) {
            secondary_dimension_expanded(c, start, len, extract_secondary_core(work), buffer);
        } else {
            primary_dimension_expanded(c, start, len, nrows, &(extract_primary_core(work)), buffer);
        }
        return buffer;
    }

public:
    /**
     * @cond
     */
    struct CompressedSparsePrimarySparseBlockWorkspace : public SparseBlockWorkspace<ROW> {
        CompressedSparsePrimarySparseBlockWorkspace(size_t s, size_t l, const WorkspaceOptions& opt, size_t maxdim) : 
            SparseBlockWorkspace<ROW>(s, l, opt.mode, opt.sorted), core(opt.cache ? maxdim : 0) {}

        PrimaryBlockWorkspaceBase core;
    };

    struct CompressedSparseSecondarySparseBlockWorkspace : public SparseBlockWorkspace<!ROW> {
        CompressedSparseSecondarySparseBlockWorkspace(size_t s, size_t l, const WorkspaceOptions& opt, size_t max_index, const V& idx, const W& idp) : 
            SparseBlockWorkspace<!ROW>(s, l, opt.mode, opt.sorted), core(max_index, idx, idp, s, l) {}

        SecondaryWorkspaceBase core;
    };

    static PrimaryBlockWorkspaceBase& extract_primary_core(SparseBlockWorkspace<ROW>* wrk) {
        return static_cast<CompressedSparsePrimarySparseBlockWorkspace*>(wrk)->core;
    }

    static SecondaryWorkspaceBase& extract_secondary_core(SparseBlockWorkspace<!ROW>* wrk) {
        return static_cast<CompressedSparseSecondarySparseBlockWorkspace*>(wrk)->core;
    }
    /**
     * @endcond
     */

    std::shared_ptr<SparseRowBlockWorkspace> sparse_row_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const { 
        if constexpr(ROW) {
            return std::shared_ptr<SparseRowBlockWorkspace>(new CompressedSparsePrimarySparseBlockWorkspace(start, length, opt, nrows));
        } else {
            return std::shared_ptr<SparseRowBlockWorkspace>(new CompressedSparseSecondarySparseBlockWorkspace(start, length, opt, nrows, indices, indptrs));
        }
    }

    std::shared_ptr<SparseColumnBlockWorkspace> sparse_column_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const { 
        if constexpr(ROW) {
            return std::shared_ptr<SparseColumnBlockWorkspace>(new CompressedSparseSecondarySparseBlockWorkspace(start, length, opt, ncols, indices, indptrs));
        } else {
            return std::shared_ptr<SparseColumnBlockWorkspace>(new CompressedSparsePrimarySparseBlockWorkspace(start, length, opt, ncols));
        }
    }

    SparseRange<T, IDX> row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowBlockWorkspace* work) const {
        size_t start = work->start, len = work->length;
        nullify_sparse_extract_pointers(work->mode, vbuffer, ibuffer);

        if constexpr(ROW) {
            return primary_dimension_raw(r, start, len, ncols, &(extract_primary_core(work)), vbuffer, ibuffer);
        } else {
            return secondary_dimension_raw(r, start, len, extract_secondary_core(work), vbuffer, ibuffer); 
        }
    }

    SparseRange<T, IDX> column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnBlockWorkspace* work) const {
        size_t start = work->start, len = work->length;
        nullify_sparse_extract_pointers(work->mode, vbuffer, ibuffer);

        if constexpr(ROW) {
            return secondary_dimension_raw(c, start, len, extract_secondary_core(work), vbuffer, ibuffer);
        } else {
            return primary_dimension_raw(c, start, len, nrows, &(extract_primary_core(work)), vbuffer, ibuffer);
        }
    }

private:
    std::pair<size_t, size_t> primary_dimension(size_t i, size_t start, size_t length, size_t otherdim, PrimaryBlockWorkspaceBase* work) const {
        bool do_cache = work && !work->cached.empty();
        if (do_cache) {
            auto val = work->cached[i];
            if (val.first != -1) {
                return val;
            }
        }

        auto iIt = indices.begin() + indptrs[i], eIt = indices.begin() + indptrs[i + 1]; 

        if (start) { // Jumping ahead if non-zero.
            iIt = std::lower_bound(iIt, eIt, start);
        } 

        auto last = start + length;
        if (last != otherdim) { // Jumping to last element.
            eIt = std::lower_bound(iIt, eIt, last);
        }

        size_t outstart = iIt - indices.begin();
        size_t outlength = eIt - iIt;
        if (do_cache) {
            work->cached[i].first = outstart;
            work->cached[i].second = outlength;
        }

        return std::make_pair(outstart, outlength);
    }

    SparseRange<T, IDX> primary_dimension_raw(size_t i, size_t start, size_t length, size_t otherdim, PrimaryBlockWorkspaceBase* work, T* out_values, IDX* out_indices) const {
        auto obtained = primary_dimension(i, start, length, otherdim, work);
        SparseRange<T, IDX> output(obtained.second);

        if (out_values) {
            if constexpr(has_data<T, U>::value) {
                output.value = values.data() + obtained.first;
            } else {
                auto vIt = values.begin() + obtained.first;
                std::copy(vIt, vIt + obtained.second, out_values);
                output.value = out_values;
            }
        } else {
            output.value = NULL;
        }

        if (out_indices) {
            if constexpr(has_data<IDX, V>::value) {
                output.index = indices.data() + obtained.first;
            } else {
                auto iIt = indices.begin() + obtained.first;
                std::copy(iIt, iIt + obtained.second, out_indices);
                output.index = out_indices;
            }
        } else {
            output.index = NULL;
        }

        return output;
    }

    void primary_dimension_expanded(size_t i, size_t start, size_t length, size_t otherdim, PrimaryBlockWorkspaceBase* work, T* out_values) const {
        std::fill(out_values, out_values + length, static_cast<T>(0));
        auto obtained = primary_dimension(i, start, length, otherdim, work);
        auto vIt = values.begin() + obtained.first;
        auto iIt = indices.begin() + obtained.first;
        for (size_t x = 0; x < obtained.second; ++x, ++vIt, ++iIt) {
            out_values[*iIt - start] = *vIt;
        }
        return;
    }

private:
    size_t max_secondary_index() const {
        if constexpr(ROW) {
            return ncols;
        } else {
            return nrows;
        }
    }

    auto define_below_index(size_t pos, size_t curptr) const {
        // We set it to zero if there is no other non-zero element before
        // 'curptr' in 'pos'.  This aligns with the idea of 'below_indices'
        // holding 1-past-the-immediately-below-index; if it is zero,
        // nothing remaings below 'curptr'.
        return (curptr > indptrs[pos] ? indices[curptr - 1] + 1 : 0);
    }

    template<class Store>
    void add_to_store(size_t primary, size_t curptr, Store& output) const {
        // If the value isn't actually needed, we don't want to do a look-up on 'values'.
        if constexpr(Store::can_ignore_values) {
            if (output.ignore_values()) {
                output.add(primary);
                return;
            }
        }
        output.add(primary, values[curptr]);
    }

    template<class Store>
    void secondary_dimension_above(IDX secondary, size_t primary, size_t index_primary, SecondaryWorkspaceBase& work, Store& output) const {
        auto& curdex = work.current_indices[index_primary];
        auto& curptr = work.current_indptrs[index_primary];

        // Remember that we already store the index corresponding to the current indptr.
        // So, we only need to do more work if the request is greater than that index.
        if (secondary > curdex) {
            auto max_index = max_secondary_index();
            auto limit = indptrs[primary + 1];

            // Having a peek at the index of the next non-zero
            // element; maybe we're lucky enough that the requested
            // index is below this, as would be the case for
            // consecutive or near-consecutive accesses.
            ++curptr;
            if (curptr < limit) {
                auto candidate = indices[curptr];
                if (candidate >= secondary) {
                    curdex = candidate;

                } else if (secondary + 1 == max_index) {
                    // Special case if the requested index is at the end of the matrix,
                    // in which case we can just jump there directly rather than 
                    // doing an unnecessary binary search.
                    curptr = limit - 1;
                    if (indices[curptr] == secondary) {
                        curdex = indices[curptr];
                    } else {
                        ++curptr;
                        curdex = max_index;
                    }

                } else {
                    // Otherwise we need to search.
                    curptr = std::lower_bound(indices.begin() + curptr, indices.begin() + limit, secondary) - indices.begin();
                    curdex = (curptr < limit ? indices[curptr] : max_index);
                }
            } else {
                curdex = max_index;
            }

            if (work.store_below) {
                work.below_indices[index_primary] = define_below_index(primary, curptr);
            }
        }

        if (secondary == curdex) { // assuming secondary < max_index, of course.
            add_to_store(primary, curptr, output);
        } else {
            output.skip(primary);
        }
    }

    template<class Store>
    void secondary_dimension_below(IDX secondary, size_t primary, size_t index_primary, SecondaryWorkspaceBase& work, Store& output) const {
        auto& curdex = work.current_indices[index_primary];
        auto& curptr = work.current_indptrs[index_primary];

        // Remember that below_indices stores 'indices[curptr - 1] + 1'.
        // So, we only need to do more work if the request is less than this;
        // otherwise the curptr remains unchanged. We also skip if curptr is
        // equal to 'indptrs[index_primary]' because decreasing is impossible.
        auto& prevdex = work.below_indices[index_primary];
        auto limit = indptrs[primary];
        if (secondary < prevdex && curptr > limit) {

            // Having a peek at the index of the non-zero element below us.
            // Maybe we're lucky and the requested index is equal to this,
            // in which case we can use it directly. This would be the case
            // for consecutive accesses in the opposite direction.
            auto candidate = indices[curptr - 1];
            if (candidate == secondary) {
                --curptr;
                curdex = candidate;

            } else if (secondary == 0) {
                // Special case if the requested index is at the start of the matrix,
                // in which case we can just jump there directly rather than 
                // doing an unnecessary binary search. 
                curdex = indices[limit];
                curptr = limit;

            } else {
                // Otherwise we need to search. Note that there's no need to
                // handle candidate < secondary, as that would imply that
                // 'secondary >= prevdex', which isn't possible.
                curptr = std::lower_bound(indices.begin() + limit, indices.begin() + curptr - 1, secondary) - indices.begin();
                curdex = indices[curptr]; // guaranteed to be valid as curptr - 1 >= limit.
            }

            prevdex = define_below_index(primary, curptr);
        }

        if (secondary == curdex) { // assuming i < max_index, of course.
            add_to_store(primary, curptr, output);
        } else {
            output.skip(primary);
        }
    }

    template<class STORE>
    void secondary_dimension(IDX i, size_t start, size_t length, SecondaryWorkspaceBase& work, STORE& output) const {
        IDX prev_i = work.previous_request;

        if (i >= prev_i) {
            for (size_t current = 0; current < length; ++current) {
                secondary_dimension_above(i, current + start, current, work, output);
            }
        } else {
            // Backfilling everything so that all uses of 'below_indices' are valid.
            if (!work.store_below) {
                work.store_below = true;
                work.below_indices.resize(length);
                for (size_t j = 0; j < length; ++j) {
                    work.below_indices[j] = define_below_index(j + start, work.current_indptrs[j]);
                }
            }

            for (size_t current = 0; current < length; ++current) {
                secondary_dimension_below(i, current + start, current, work, output);
            }
        }

        work.previous_request = i;
        return;
    }

    // Overview of the Store contract:
    struct raw_store {
        T* out_values;
        IDX* out_indices;
        size_t n = 0;

        // Whether or not 'values' might be ignored, typically
        // because a NULL pointer is supplied in the 'vbuffer'
        static constexpr bool can_ignore_values = true;

        // If 'can_ignore_values = true', we need this method
        // to define whether the values are ignored for this instance.
        bool ignore_values() const {
            return !out_values;
        }

        // If 'can_ignore_values = true', we need this method to handle
        // hits, but _without_ a costly fetch of the value from memory.
        void add(IDX i) {
            ++n;
            if (out_indices) {
                *out_indices = i;
                ++out_indices;
            }
        }

        // Storing a hit in terms of its index and its value. If 'can_ignore_values = true',
        // it can be assumed that 'skip_values()' is false if this method is called. Note 
        // that add() and skip() is always called on indices in increasing order.
        void add(IDX i, T val) {
            add(i);
            *out_values = val;
            ++out_values;
            return;
        }

        // Skip an index if it isn't found in the data. 
        void skip(IDX) {} 
    };

    SparseRange<T, IDX> secondary_dimension_raw(IDX i, size_t start, size_t length, SecondaryWorkspaceBase& work, T* out_values, IDX* out_indices) const {
        raw_store store;
        store.out_values = out_values;
        store.out_indices = out_indices;
        secondary_dimension(i, start, length, work, store);
        return SparseRange<T, IDX>(store.n, out_values, out_indices);
    }

    struct expanded_store {
        T* out_values;
        IDX first;

        static constexpr bool can_ignore_values = false;

        void add(IDX i, T val) {
            out_values[i - first] = val;
            return;
        }

        void skip(IDX) {} 
    };

    void secondary_dimension_expanded(IDX i, size_t start, size_t length, SecondaryWorkspaceBase& work, T* out_values) const {
        std::fill(out_values, out_values + length, static_cast<T>(0));
        expanded_store store;
        store.out_values = out_values;
        store.first = start;
        secondary_dimension(i, start, length, work, store);
        return;
    }

public:
    /**
     * @cond
     */
    struct PrimaryIndexedWorkspaceBase {
        PrimaryIndexedWorkspaceBase(size_t cache_size) : starts(cache_size, -1) {}
        std::vector<size_t> starts;
    };

    struct CompressedSparsePrimaryDenseIndexWorkspace : public DenseIndexWorkspace<IDX, ROW> {
        CompressedSparsePrimaryDenseIndexWorkspace(std::vector<IDX> subset, const WorkspaceOptions& opt, size_t maxdim) : 
            DenseIndexWorkspace<IDX, ROW>(subset.size()), core(opt.cache ? maxdim : 0), indices_(std::move(subset)) {}

        PrimaryIndexedWorkspaceBase core;

        std::vector<IDX> indices_; 
        const std::vector<IDX>& indices() const { return indices_; }
    };

    struct CompressedSparseSecondaryDenseIndexWorkspace : public DenseIndexWorkspace<IDX, !ROW> {
        CompressedSparseSecondaryDenseIndexWorkspace(std::vector<IDX> subset, size_t max_index, const V& idx, const W& idp) : 
            DenseIndexWorkspace<IDX, !ROW>(subset.size()),  core(max_index, idx, idp, subset), indices_(std::move(subset)) {}

        SecondaryWorkspaceBase core;

        std::vector<IDX> indices_; // must be after 'core', as we're moving 'subset' after 'core' is done with it.
        const std::vector<IDX>& indices() const { return indices_; }
    };

    static PrimaryIndexedWorkspaceBase& extract_primary_core(DenseIndexWorkspace<IDX, ROW>* wrk) {
        return static_cast<CompressedSparsePrimaryDenseIndexWorkspace*>(wrk)->core;
    }

    static SecondaryWorkspaceBase& extract_secondary_core(DenseIndexWorkspace<IDX, !ROW>* wrk) {
        return static_cast<CompressedSparseSecondaryDenseIndexWorkspace*>(wrk)->core;
    }
    /**
     * @endcond
     */

    std::shared_ptr<DenseRowIndexWorkspace<IDX> > dense_row_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const { 
        if constexpr(ROW) {
            return std::shared_ptr<DenseRowIndexWorkspace<IDX> >(new CompressedSparsePrimaryDenseIndexWorkspace(std::move(subset), opt, nrows));
        } else {
            return std::shared_ptr<DenseRowIndexWorkspace<IDX> >(new CompressedSparseSecondaryDenseIndexWorkspace(std::move(subset), nrows, indices, indptrs));
        }
    }

    std::shared_ptr<DenseColumnIndexWorkspace<IDX> > dense_column_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const { 
        if constexpr(ROW) {
            return std::shared_ptr<DenseColumnIndexWorkspace<IDX> >(new CompressedSparseSecondaryDenseIndexWorkspace(std::move(subset), ncols, indices, indptrs));
        } else {
            return std::shared_ptr<DenseColumnIndexWorkspace<IDX> >(new CompressedSparsePrimaryDenseIndexWorkspace(std::move(subset), opt, ncols));
        }
    }

    const T* row(size_t r, T* buffer, DenseRowIndexWorkspace<IDX>* work) const {
        if constexpr(ROW) {
            primary_dimension_expanded(r, work->indices(), ncols, extract_primary_core(work), buffer);
        } else {
            secondary_dimension_expanded(r, work->indices(), extract_secondary_core(work), buffer);
        }
        return buffer;
    }

    const T* column(size_t c, T* buffer, DenseColumnIndexWorkspace<IDX>* work) const {
        if constexpr(ROW) {
            secondary_dimension_expanded(c, work->indices(), extract_secondary_core(work), buffer);
        } else {
            primary_dimension_expanded(c, work->indices(), nrows, extract_primary_core(work), buffer);
        }
        return buffer;
    }

public:
    /**
     * @cond
     */
    struct CompressedSparsePrimarySparseIndexWorkspace : public SparseIndexWorkspace<IDX, ROW> {
        CompressedSparsePrimarySparseIndexWorkspace(std::vector<IDX> subset, const WorkspaceOptions& opt, size_t maxdim) : 
            SparseIndexWorkspace<IDX, ROW>(subset.size(), opt.mode, opt.sorted), core(opt.cache ? maxdim : 0), indices_(std::move(subset)) {}

        PrimaryIndexedWorkspaceBase core;

        std::vector<IDX> indices_; 
        const std::vector<IDX>& indices() const { return indices_; }
    };

    struct CompressedSparseSecondarySparseIndexWorkspace : public SparseIndexWorkspace<IDX, !ROW> {
        CompressedSparseSecondarySparseIndexWorkspace(std::vector<IDX> subset, const WorkspaceOptions& opt, size_t max_index, const V& idx, const W& idp) : 
            SparseIndexWorkspace<IDX, !ROW>(subset.size(), opt.mode, opt.sorted), core(max_index, idx, idp, subset), indices_(std::move(subset)) {}

        SecondaryWorkspaceBase core;

        std::vector<IDX> indices_; // must be after 'core', as we're moving 'subset' after 'core' is done with it.
        const std::vector<IDX>& indices() const { return indices_; }
    };

    static PrimaryIndexedWorkspaceBase& extract_primary_core(SparseIndexWorkspace<IDX, ROW>* wrk) {
        return static_cast<CompressedSparsePrimarySparseIndexWorkspace*>(wrk)->core;
    }

    static SecondaryWorkspaceBase& extract_secondary_core(SparseIndexWorkspace<IDX, !ROW>* wrk) {
        return static_cast<CompressedSparseSecondarySparseIndexWorkspace*>(wrk)->core;
    }
    /**
     * @endcond
     */

    std::shared_ptr<SparseRowIndexWorkspace<IDX> > sparse_row_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const { 
        if constexpr(ROW) {
            return std::shared_ptr<SparseRowIndexWorkspace<IDX> >(new CompressedSparsePrimarySparseIndexWorkspace(std::move(subset), opt, nrows));
        } else {
            return std::shared_ptr<SparseRowIndexWorkspace<IDX> >(new CompressedSparseSecondarySparseIndexWorkspace(std::move(subset), opt, nrows, indices, indptrs));
        }
    }

    std::shared_ptr<SparseColumnIndexWorkspace<IDX> > sparse_column_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const { 
        if constexpr(ROW) {
            return std::shared_ptr<SparseColumnIndexWorkspace<IDX> >(new CompressedSparseSecondarySparseIndexWorkspace(std::move(subset), opt, ncols, indices, indptrs));
        } else {
            return std::shared_ptr<SparseColumnIndexWorkspace<IDX> >(new CompressedSparsePrimarySparseIndexWorkspace(std::move(subset), opt, ncols));
        }
    }

    SparseRange<T, IDX> row(size_t r, T* vbuffer, IDX* ibuffer, SparseRowIndexWorkspace<IDX>* work) const {
        nullify_sparse_extract_pointers(work->mode, vbuffer, ibuffer);

        if constexpr(ROW) {
            return primary_dimension_raw(r, work->indices(), ncols, extract_primary_core(work), vbuffer, ibuffer);
        } else {
            return secondary_dimension_raw(r, work->indices(), extract_secondary_core(work), vbuffer, ibuffer);
        }
    }

    SparseRange<T, IDX> column(size_t c, T* vbuffer, IDX* ibuffer, SparseColumnIndexWorkspace<IDX>* work) const {
        nullify_sparse_extract_pointers(work->mode, vbuffer, ibuffer);

        if constexpr(ROW) {
            return secondary_dimension_raw(c, work->indices(), extract_secondary_core(work), vbuffer, ibuffer);
        } else {
            return primary_dimension_raw(c, work->indices(), nrows, extract_primary_core(work), vbuffer, ibuffer);
        }
    }

private:
    template<class Store>
    void primary_dimension(size_t i, const std::vector<IDX>& subset, PrimaryIndexedWorkspaceBase& work, Store& store) const {
        auto iIt = indices.begin() + indptrs[i], eIt = indices.begin() + indptrs[i + 1]; 

        if (!subset.empty() && subset[0]) { // Jumping ahead if non-zero.
            bool do_cache = !work.starts.empty();
            if (do_cache) {
                if (work.starts[i] != -1) { // retrieving the jump from cache, if we came here before.
                    iIt += work.starts[i];
                } else {
                    auto iIt2 = std::lower_bound(iIt, eIt, subset[0]);
                    work.starts[i] = iIt2 - iIt;
                    iIt = iIt2;
                }
            } else {
                iIt = std::lower_bound(iIt, eIt, subset[0]);
            }
        } 

        if (iIt == eIt) {
            return;
        }

        size_t counter = 0, end = subset.size();
        while (counter < end) {
            auto current = subset[counter];

            while (iIt != eIt && current > *iIt) {
                ++iIt;
            }
            if (iIt == eIt) {
                break;
            }

            if (current == *iIt) {
                add_to_store(current, iIt - indices.begin(), store);
            } else {
                store.skip(current);
            }
            ++counter;
        }

        return;
    }

    SparseRange<T, IDX> primary_dimension_raw(size_t i, const std::vector<IDX>& subset, size_t otherdim, PrimaryIndexedWorkspaceBase& work, T* out_values, IDX* out_indices) const {
        raw_store store;
        store.out_values = out_values;
        store.out_indices = out_indices;
        primary_dimension(i, subset, work, store);
        return SparseRange<T, IDX>(store.n, out_values, out_indices);
    }

    struct expanded_store2 {
        T* out_values;

        static constexpr bool can_ignore_values = false;

        void add(IDX, T val) {
            *out_values = val;
            ++out_values;
            return;
        }

        void skip(IDX) {
            ++out_values;
        }
    };

    void primary_dimension_expanded(size_t i, const std::vector<IDX>& subset, size_t otherdim, PrimaryIndexedWorkspaceBase& work, T* out_values) const {
        std::fill(out_values, out_values + subset.size(), static_cast<T>(0));
        expanded_store2 store;
        store.out_values = out_values;
        primary_dimension(i, subset, work, store);
        return;
    }

private:
    template<class STORE>
    void secondary_dimension(IDX i, const std::vector<IDX>& subset, SecondaryWorkspaceBase& work, STORE& output) const {
        IDX prev_i = work.previous_request;
        size_t length = subset.size();

        if (i >= prev_i) {
            for (size_t current = 0; current < length; ++current) {
                secondary_dimension_above(i, subset[current], current, work, output);
            }
        } else {
            // Backfilling everything so that all uses of 'below_indices' are valid.
            if (!work.store_below) {
                work.store_below = true;
                work.below_indices.resize(length);
                for (size_t j = 0; j < length; ++j) {
                    work.below_indices[j] = define_below_index(subset[j], work.current_indptrs[j]);
                }
            }

            for (size_t current = 0; current < length; ++current) {
                secondary_dimension_below(i, subset[current], current, work, output);
            }
        }

        work.previous_request = i;
        return;
    }

    SparseRange<T, IDX> secondary_dimension_raw(IDX i, const std::vector<IDX>& subset, SecondaryWorkspaceBase& work, T* out_values, IDX* out_indices) const {
        raw_store store;
        store.out_values = out_values;
        store.out_indices = out_indices;
        secondary_dimension(i, subset, work, store);
        return SparseRange<T, IDX>(store.n, out_values, out_indices);
    }

    void secondary_dimension_expanded(IDX i, const std::vector<IDX>& subset, SecondaryWorkspaceBase& work, T* out_values) const {
        std::fill(out_values, out_values + subset.size(), static_cast<T>(0));
        expanded_store2 store;
        store.out_values = out_values;
        secondary_dimension(i, subset, work, store);
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
