#ifndef TATAMI_COMPRESSED_SPARSE_MATRIX_H
#define TATAMI_COMPRESSED_SPARSE_MATRIX_H

#include "../base/Matrix.hpp"
#include "primary_extraction.hpp"
#include "secondary_extraction.hpp"

#include "../utils/ElementType.hpp"
#include "../utils/has_data.hpp"
#include "../utils/PseudoOracularExtractor.hpp"

#include <vector>
#include <algorithm>
#include <memory>
#include <utility>
#include <stdexcept>

/**
 * @file CompressedSparseMatrix.hpp
 *
 * @brief Compressed sparse matrix representation. 
 */

namespace tatami {

/**
 * @cond
 */
namespace CompressedSparseMatrix_internal {

/********************
 *** Primary full ***
 ********************/

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
class PrimaryMyopicFullDense : public MyopicDenseExtractor<Value_, Index_> {
public:
    PrimaryMyopicFullDense(const ValueStorage_& values, const IndexStorage_& indices, const PointerStorage_& pointers, Index_ secondary) :
        my_values(values),
        my_indices(indices), 
        my_pointers(pointers),
        my_secondary(secondary) 
    {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto offset = my_pointers[i];
        size_t delta = my_pointers[i+1] - my_pointers[i];
        std::fill_n(buffer, my_secondary, static_cast<Value_>(0));
        for (size_t x = 0; x < delta; ++x) {
            auto cur_offset = offset + x;
            buffer[my_indices[cur_offset]] = my_values[cur_offset];
        }
        return buffer;
    }

private:
    const ValueStorage_& my_values;
    const IndexStorage_& my_indices;
    const PointerStorage_& my_pointers;
    Index_ my_secondary;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
class PrimaryMyopicFullSparse : public MyopicSparseExtractor<Value_, Index_> {
public:
    PrimaryMyopicFullSparse(const ValueStorage_& values, const IndexStorage_& indices, const PointerStorage_& pointers, Index_ secondary, const Options& opt) :
        my_values(values),
        my_indices(indices),
        my_pointers(pointers),
        my_secondary(secondary),
        my_needs_value(opt.sparse_extract_value),
        my_needs_index(opt.sparse_extract_index) 
    {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        auto offset = my_pointers[i];
        auto delta = my_pointers[i+1] - my_pointers[i];

        SparseRange<Value_, Index_> output(delta, NULL, NULL);
        if (my_needs_value) {
            output.value = sparse_utils::extract_primary_vector(my_values, offset, delta, value_buffer);
        }
        if (my_needs_index) {
            output.index = sparse_utils::extract_primary_vector(my_indices, offset, delta, index_buffer);
        }
        return output;
    }

private:
    const ValueStorage_& my_values;
    const IndexStorage_& my_indices;
    const PointerStorage_& my_pointers;
    Index_ my_secondary;
    bool my_needs_value, my_needs_index;
};

/*********************
 *** Primary block ***
 *********************/

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
class PrimaryMyopicBlockDense : public MyopicDenseExtractor<Value_, Index_> {
public:
    PrimaryMyopicBlockDense(const ValueStorage_& values, const IndexStorage_& indices, const PointerStorage_& pointers, Index_ secondary, Index_ block_start, Index_ block_length) :
        my_values(values), 
        my_indices(indices), 
        my_pointers(pointers), 
        my_secondary(secondary), 
        my_block_start(block_start), 
        my_block_length(block_length) 
    {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto iStart = my_indices.begin() + my_pointers[i];
        auto iEnd = my_indices.begin() + my_pointers[i + 1];
        sparse_utils::refine_primary_block_limits(iStart, iEnd, my_secondary, my_block_start, my_block_length);
        size_t offset = (iStart - my_indices.begin());
        size_t number = iEnd - iStart;

        std::fill_n(buffer, my_block_length, static_cast<Value_>(0));
        for (size_t i = 0; i < number; ++i) {
            auto cur_offset = offset + i;
            buffer[my_indices[cur_offset] - my_block_start] = my_values[cur_offset];
        }
        return buffer;
    }

private:
    const ValueStorage_& my_values;
    const IndexStorage_& my_indices;
    const PointerStorage_& my_pointers;
    Index_ my_secondary;
    Index_ my_block_start, my_block_length;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
class PrimaryMyopicBlockSparse : public MyopicSparseExtractor<Value_, Index_> {
public:
    PrimaryMyopicBlockSparse(const ValueStorage_& values, const IndexStorage_& indices, const PointerStorage_& pointers, Index_ secondary, Index_ block_start, Index_ block_length, const Options& opt) :
        my_values(values),
        my_indices(indices),
        my_pointers(pointers),
        my_secondary(secondary),
        my_block_start(block_start),
        my_block_length(block_length),
        my_needs_value(opt.sparse_extract_value),
        my_needs_index(opt.sparse_extract_index) 
    {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        auto iStart = my_indices.begin() + my_pointers[i];
        auto iEnd = my_indices.begin() + my_pointers[i + 1];
        sparse_utils::refine_primary_block_limits(iStart, iEnd, my_secondary, my_block_start, my_block_length);
        size_t offset = iStart - my_indices.begin();
        size_t delta = iEnd - iStart;

        SparseRange<Value_, Index_> output(delta, NULL, NULL);
        if (my_needs_value) {
            output.value = sparse_utils::extract_primary_vector(my_values, offset, delta, value_buffer);
        }
        if (my_needs_index) {
            output.index = sparse_utils::extract_primary_vector(my_indices, offset, delta, index_buffer);
        }
        return output;
    }

private:
    const ValueStorage_& my_values;
    const IndexStorage_& my_indices;
    const PointerStorage_& my_pointers;
    Index_ my_secondary;
    Index_ my_block_start, my_block_length;
    bool my_needs_value, my_needs_index;
};

/***********************
 *** Primary indexed ***
 ***********************/

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
class PrimaryMyopicIndexDense : public MyopicDenseExtractor<Value_, Index_> {
public:
    PrimaryMyopicIndexDense(const ValueStorage_& values, const IndexStorage_& indices, const PointerStorage_& pointers, Index_ secondary, const VectorPtr<Index_>& indices_ptr) :
        my_values(values),
        my_indices(indices),
        my_pointers(pointers),
        my_retriever(*indices_ptr, secondary),
        my_num_indices(indices_ptr->size()) 
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        std::fill_n(buffer, my_num_indices, static_cast<Value_>(0));
        auto vIt = my_values.begin() + my_pointers[i];
        my_retriever.populate(
            my_indices.begin() + my_pointers[i], 
            my_indices.begin() + my_pointers[i+1],
            [&](size_t s, size_t offset) {
                buffer[s] = *(vIt + offset);
            }
        );
        return buffer;
    }

private:
    const ValueStorage_& my_values;
    const IndexStorage_& my_indices;
    const PointerStorage_& my_pointers;
    sparse_utils::RetrievePrimarySubsetDense<Index_> my_retriever;
    size_t my_num_indices;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
class PrimaryMyopicIndexSparse : public MyopicSparseExtractor<Value_, Index_> {
public:
    PrimaryMyopicIndexSparse(const ValueStorage_& values, const IndexStorage_& indices, const PointerStorage_& pointers, Index_ secondary, const VectorPtr<Index_>& indices_ptr, const Options& opt) :
        my_values(values),
        my_indices(indices), 
        my_pointers(pointers),
        my_retriever(*indices_ptr, secondary),
        my_needs_value(opt.sparse_extract_value),
        my_needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        Index_ count = 0;
        auto vcopy = value_buffer;
        auto icopy = index_buffer;

        auto vIt = my_values.begin() + my_pointers[i];
        my_retriever.populate(
            my_indices.begin() + my_pointers[i], 
            my_indices.begin() + my_pointers[i+1],
            [&](size_t offset, Index_ ix) {
                ++count;
                if (my_needs_value) {
                    *vcopy = *(vIt + offset);
                    ++vcopy;
                }
                if (my_needs_index) {
                    *icopy = ix;
                    ++icopy;
                }
            }
        );

        return SparseRange<Value_, Index_>(count, my_needs_value ? value_buffer : NULL, my_needs_index ? index_buffer : NULL);
    }

private:
    const ValueStorage_& my_values;
    const IndexStorage_& my_indices;
    const PointerStorage_& my_pointers;
    sparse_utils::RetrievePrimarySubsetSparse<Index_> my_retriever;
    bool my_needs_value, my_needs_index;
};

/**********************
 *** Secondary full ***
 **********************/

template<typename Index_, class IndexStorage_, class PointerStorage_>
class ServeIndices {
public:
    ServeIndices(const IndexStorage_& i, const PointerStorage_& p) : my_indices(i), my_pointers(p) {}

private:
    const IndexStorage_& my_indices;
    const PointerStorage_& my_pointers;

public:
    typedef ElementType<PointerStorage_> pointer_type;

    pointer_type start_offset(Index_ primary) const {
        return my_pointers[primary];
    }

    pointer_type end_offset(Index_ primary) const {
        return my_pointers[primary + 1];
    }

    auto raw(Index_) const {
        return my_indices.begin();
    }
};

template<typename Index_, class IndexStorage_, class PointerStorage_>
auto make_ServeIndices(const IndexStorage_& i, const PointerStorage_& p) {
    return ServeIndices<Index_, IndexStorage_, PointerStorage_>(i, p);
}

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
class SecondaryMyopicFullDense : public MyopicDenseExtractor<Value_, Index_> {
public:
    SecondaryMyopicFullDense(const ValueStorage_& values, const IndexStorage_& indices, const PointerStorage_& pointers, Index_ secondary) :
        my_values(values),
        my_cache(make_ServeIndices<Index_>(indices, pointers), secondary, pointers.size() - 1) 
    {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        std::fill_n(buffer, my_cache.size(), static_cast<Value_>(0));
        my_cache.search(i, [&](Index_, Index_ index_primary, ElementType<PointerStorage_> ptr) {
            buffer[index_primary] = my_values[ptr];
        });
        return buffer;
    }

private:
    const ValueStorage_& my_values;
    sparse_utils::FullSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexStorage_, PointerStorage_> > my_cache;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
class SecondaryMyopicFullSparse : public MyopicSparseExtractor<Value_, Index_> {
public:
    SecondaryMyopicFullSparse(const ValueStorage_& values, const IndexStorage_& indices, const PointerStorage_& pointers, Index_ secondary, const Options& opt) :
        my_values(values),
        my_cache(make_ServeIndices<Index_>(indices, pointers), secondary, pointers.size() - 1),
        my_needs_value(opt.sparse_extract_value),
        my_needs_index(opt.sparse_extract_index) 
    {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        Index_ count = 0;
        my_cache.search(i, [&](Index_ primary, Index_, ElementType<PointerStorage_> ptr) {
            if (my_needs_value) {
                value_buffer[count] = my_values[ptr];
            }
            if (my_needs_index) {
                index_buffer[count] = primary;
            }
            ++count;
        });
        return SparseRange<Value_, Index_>(count, my_needs_value ? value_buffer : NULL, my_needs_index ? index_buffer : NULL);
    }

private:
    const ValueStorage_& my_values;
    sparse_utils::FullSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexStorage_, PointerStorage_> > my_cache;
    bool my_needs_value, my_needs_index;
};

/***********************
 *** Secondary block ***
 ***********************/

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
class SecondaryMyopicBlockDense : public MyopicDenseExtractor<Value_, Index_> {
public:
    SecondaryMyopicBlockDense(const ValueStorage_& values, const IndexStorage_& indices, const PointerStorage_& pointers, Index_ secondary, Index_ block_start, Index_ block_length) :
        my_values(values),
        my_cache(make_ServeIndices<Index_>(indices, pointers), secondary, block_start, block_length) 
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        std::fill_n(buffer, my_cache.size(), static_cast<Value_>(0));
        my_cache.search(i, [&](Index_, Index_ index_primary, ElementType<PointerStorage_> ptr) {
            buffer[index_primary] = my_values[ptr];
        });
        return buffer;
    }

private:
    const ValueStorage_& my_values;
    sparse_utils::BlockSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexStorage_, PointerStorage_> > my_cache;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
class SecondaryMyopicBlockSparse : public MyopicSparseExtractor<Value_, Index_> {
public:
    SecondaryMyopicBlockSparse(const ValueStorage_& values, const IndexStorage_& indices, const PointerStorage_& pointers, Index_ secondary, Index_ block_start, Index_ block_length, const Options& opt) :
        my_values(values),
        my_cache(make_ServeIndices<Index_>(indices, pointers), secondary, block_start, block_length), 
        my_needs_value(opt.sparse_extract_value),
        my_needs_index(opt.sparse_extract_index) 
    {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        Index_ count = 0;
        my_cache.search(i, [&](Index_ primary, Index_, ElementType<PointerStorage_> ptr) {
            if (my_needs_value) {
                value_buffer[count] = my_values[ptr];
            }
            if (my_needs_index) {
                index_buffer[count] = primary;
            }
            ++count;
        });
        return SparseRange<Value_, Index_>(count, my_needs_value ? value_buffer : NULL, my_needs_index ? index_buffer : NULL);
    }

private:
    const ValueStorage_& my_values;
    sparse_utils::BlockSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexStorage_, PointerStorage_> > my_cache;
    bool my_needs_value, my_needs_index;
};

/***********************
 *** Secondary index ***
 ***********************/

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
class SecondaryMyopicIndexDense : public MyopicDenseExtractor<Value_, Index_> {
public:
    SecondaryMyopicIndexDense(const ValueStorage_& values, const IndexStorage_& indices, const PointerStorage_& pointers, Index_ secondary, VectorPtr<Index_> sub_ptr) :
        my_values(values), my_cache(make_ServeIndices<Index_>(indices, pointers), secondary, std::move(sub_ptr)) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        std::fill_n(buffer, my_cache.size(), static_cast<Value_>(0));
        my_cache.search(i, [&](Index_, Index_ index_primary, ElementType<PointerStorage_> ptr) {
            buffer[index_primary] = my_values[ptr];
        });
        return buffer;
    }

private:
    const ValueStorage_& my_values;
    sparse_utils::IndexSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexStorage_, PointerStorage_> > my_cache;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
class SecondaryMyopicIndexSparse : public MyopicSparseExtractor<Value_, Index_> {
public:
    SecondaryMyopicIndexSparse(const ValueStorage_& values, const IndexStorage_& indices, const PointerStorage_& pointers, Index_ secondary, VectorPtr<Index_> sub_ptr, const Options& opt) :
        my_values(values),
        my_cache(make_ServeIndices<Index_>(indices, pointers), secondary, std::move(sub_ptr)),
        my_needs_value(opt.sparse_extract_value),
        my_needs_index(opt.sparse_extract_index) 
    {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        Index_ count = 0;
        my_cache.search(i, [&](Index_ primary, Index_, ElementType<PointerStorage_> ptr) {
            if (my_needs_value) {
                value_buffer[count] = my_values[ptr];
            }
            if (my_needs_index) {
                index_buffer[count] = primary;
            }
            ++count;
        });
        return SparseRange<Value_, Index_>(count, my_needs_value ? value_buffer : NULL, my_needs_index ? index_buffer : NULL);
    }

private:
    const ValueStorage_& my_values;
    sparse_utils::IndexSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexStorage_, PointerStorage_> > my_cache;
    bool my_needs_value, my_needs_index;
};

}
/**
 * @endcond
 */

/**
 * @brief Compressed sparse matrix representation.
 *
 * @tparam Value_ Type of the matrix values.
 * @tparam Index_ Type of the row/column indices.
 * @tparam ValueStorage_ Vector class used to store the matrix values internally.
 * This does not necessarily have to contain `Value_`, as long as the type is convertible to `Value_`.
 * Methods should be available for `size()`, `begin()`, `end()` and `[]`.
 * If a method is available for `data()` that returns a `const Value_*`, it will also be used.
 * @tparam IndexStorage_ Vector class used to store the row/column indices internally.
 * This does not necessarily have to contain `Index_`, as long as the type is convertible to `Index_`.
 * Methods should be available for `size()`, `begin()`, `end()` and `[]`.
 * If a method is available for `data()` that returns a `const Index_*`, it will also be used.
 * @tparam PointerStorage_ Vector class used to store the column/row index pointers.
 * Methods should be available for `size()`, `begin()`, `end()` and `[]`.
 */
template<
    typename Value_, 
    typename Index_,
    class ValueStorage_ = std::vector<Value_>, 
    class IndexStorage_ = std::vector<Index_>, 
    class PointerStorage_ = std::vector<size_t> 
>
class CompressedSparseMatrix : public Matrix<Value_, Index_> {
public:
    /**
     * @param nrow Number of rows.
     * @param ncol Number of columns.
     * @param values Vector of non-zero elements.
     * @param indices Vector of row indices (if `csr = false`) or column indices (if `csr = true`) for the non-zero elements.
     * @param pointers Vector of index pointers.
     * @param csr Whether this is a compressed sparse row representation.
     * @param check Should the input vectors be checked for validity?
     *
     * If `check=true`, the constructor will check that `values` and `indices` have the same length, equal to the number of structural non-zero elements;
     * `pointers` has length equal to the number of rows (if `csr = true`) or columns (otherwise) plus one;
     * `pointers` is non-decreasing with first and last values set to 0 and the number of structural non-zeroes, respectively;
     * `indices` is strictly increasing within each interval defined by successive elements of `pointers`;
     * and all values of `indices` are non-negative and less than the number of columns (if `csr = true`) or rows (otherwise).
     */
    CompressedSparseMatrix(Index_ nrow, Index_ ncol, ValueStorage_ values, IndexStorage_ indices, PointerStorage_ pointers, bool csr, bool check = true) : 
        my_nrow(nrow), my_ncols(ncol), my_values(std::move(values)), my_indices(std::move(indices)), my_pointers(std::move(pointers)), my_csr(csr)
    {
        if (check) {
            size_t nnzero = my_values.size();
            if (nnzero != my_indices.size()) {
                throw std::runtime_error("'my_values' and 'my_indices' should be of the same length");
            }

            size_t npointers = my_pointers.size(); 
            if (my_csr) {
                if (npointers != static_cast<size_t>(my_nrow) + 1){
                    throw std::runtime_error("length of 'pointers' should be equal to 'nrow + 1'");
                }
            } else {
                if (npointers != static_cast<size_t>(my_ncols) + 1){
                    throw std::runtime_error("length of 'pointers' should be equal to 'ncols + 1'");
                }
            }

            if (my_pointers[0] != 0) {
                throw std::runtime_error("first element of 'pointers' should be zero");
            }

            auto last = my_pointers[npointers - 1]; // don't use back() as this is not guaranteed to be available for arbitrary PointerStorage_.
            if (static_cast<size_t>(last) != nnzero) {
                throw std::runtime_error("last element of 'pointers' should be equal to length of 'indices'");
            }

            ElementType<IndexStorage_> max_index = (my_csr ? my_ncols : my_nrow);
            for (size_t i = 1; i < npointers; ++i) {
                auto start = my_pointers[i- 1], end = my_pointers[i];
                if (end < start || end > last) {
                    throw std::runtime_error("'pointers' should be in non-decreasing order");
                }

                for (auto x = start; x < end; ++x) {
                    if (my_indices[x] < 0 || my_indices[x] >= max_index) {
                        throw std::runtime_error("'indices' should contain non-negative integers less than the number of " + (my_csr ? std::string("columns") : std::string("rows")));
                    }
                }

                for (decltype(start) j = start + 1; j < end; ++j) {
                    if (my_indices[j] <= my_indices[j - 1]) {
                        throw std::runtime_error("'indices' should be strictly increasing within each " + (my_csr ? std::string("row") : std::string("column")));
                    }
                }
            }
        }
    }

private:
    Index_ my_nrow, my_ncols;
    ValueStorage_ my_values;
    IndexStorage_ my_indices;
    PointerStorage_ my_pointers;
    bool my_csr;

public:
    Index_ nrow() const { return my_nrow; }

    Index_ ncol() const { return my_ncols; }

    bool is_sparse() const { return true; }

    double is_sparse_proportion() const { return 1; }

    bool prefer_rows() const { return my_csr; }

    double prefer_rows_proportion() const { return static_cast<double>(my_csr); }

    bool uses_oracle(bool) const { return false; }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

private:
    Index_ secondary() const {
        if (my_csr) {
            return my_ncols;
        } else {
            return my_nrow;
        }
    }

    /*****************************
     ******* Dense myopic ********
     *****************************/
public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, const Options&) const {
        if (my_csr == row) {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicFullDense<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(
                my_values, my_indices, my_pointers, secondary()
            );
        } else {
            return std::make_unique<CompressedSparseMatrix_internal::SecondaryMyopicFullDense<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(
                my_values, my_indices, my_pointers, secondary()
            ); 
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, Index_ block_start, Index_ block_end, const Options&) const {
        if (my_csr == row) {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicBlockDense<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(
                my_values, my_indices, my_pointers, secondary(), block_start, block_end
            );
        } else {
            return std::make_unique<CompressedSparseMatrix_internal::SecondaryMyopicBlockDense<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(
                my_values, my_indices, my_pointers, secondary(), block_start, block_end
            );
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, VectorPtr<Index_> indices_ptr, const Options&) const {
        if (my_csr == row) {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicIndexDense<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(
                my_values, my_indices, my_pointers, secondary(), std::move(indices_ptr)
            );
        } else {
            return std::make_unique<CompressedSparseMatrix_internal::SecondaryMyopicIndexDense<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(
                my_values, my_indices, my_pointers, secondary(), std::move(indices_ptr)
            );
        }
    }

    /******************************
     ******* Sparse myopic ********
     ******************************/
public:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const Options& opt) const {
        if (my_csr == row) {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicFullSparse<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(
                my_values, my_indices, my_pointers, secondary(), opt
            );
        } else {
            return std::make_unique<CompressedSparseMatrix_internal::SecondaryMyopicFullSparse<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(
                my_values, my_indices, my_pointers, secondary(), opt
            );
        }
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_end, const Options& opt) const {
        if (my_csr == row) {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicBlockSparse<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(
                my_values, my_indices, my_pointers, secondary(), block_start, block_end, opt
            );
        } else {
            return std::make_unique<CompressedSparseMatrix_internal::SecondaryMyopicBlockSparse<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(
                my_values, my_indices, my_pointers, secondary(), block_start, block_end, opt
            );
        }
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        if (my_csr == row) {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicIndexSparse<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(
                my_values, my_indices, my_pointers, secondary(), std::move(indices_ptr), opt
            );
        } else {
            return std::make_unique<CompressedSparseMatrix_internal::SecondaryMyopicIndexSparse<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(
                my_values, my_indices, my_pointers, secondary(), std::move(indices_ptr), opt
            );
        }
    }

    /*******************************
     ******* Dense oracular ********
     *******************************/
public:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense(row, opt));
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_end, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense(row, block_start, block_end, opt));
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> my_indices_ptr, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense(row, std::move(my_indices_ptr), opt));
    }

    /********************************
     ******* Sparse oracular ********
     ********************************/
public:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse(row, opt));
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_end, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse(row, block_start, block_end, opt));
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> my_indices_ptr, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse(row, std::move(my_indices_ptr), opt));
    }
};

/**
 * @brief Compressed sparse column matrix.
 *
 * See `tatami::CompressedSparseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_, class ValueStorage_ = std::vector<Value_>, class IndexStorage_ = std::vector<Index_>, class PointerStorage_ = std::vector<size_t> >
class CompressedSparseColumnMatrix : public CompressedSparseMatrix<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> {
public:
    /**
     * @param nrow Number of rows.
     * @param ncol Number of columns.
     * @param values Vector of non-zero elements.
     * @param indices Vector of row indices for the non-zero elements.
     * @param pointers Vector of index pointers, of length equal to the number of columns plus 1.
     * @param check Should the input vectors be checked for validity?
     */
    CompressedSparseColumnMatrix(Index_ nrow, Index_ ncol, ValueStorage_ values, IndexStorage_ indices, PointerStorage_ pointers, bool check = true) :
        CompressedSparseMatrix<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_>(nrow, ncol, std::move(values), std::move(indices), std::move(pointers), false, check) {}
};

/**
 * @brief Compressed sparse row matrix.
 *
 * See `tatami::CompressedSparseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_, class ValueStorage_ = std::vector<Value_>, class IndexStorage_ = std::vector<Index_>, class PointerStorage_ = std::vector<size_t> >
class CompressedSparseRowMatrix : public CompressedSparseMatrix<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> {
public:
    /**
     * @param nrow Number of rows.
     * @param ncol Number of columns.
     * @param values Vector of non-zero elements.
     * @param indices Vector of row indices for the non-zero elements.
     * @param pointers Vector of index pointers, of length equal to the number of columns plus 1.
     * @param check Should the input vectors be checked for validity?
     */
    CompressedSparseRowMatrix(Index_ nrow, Index_ ncol, ValueStorage_ values, IndexStorage_ indices, PointerStorage_ pointers, bool check = true) :
        CompressedSparseMatrix<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_>(nrow, ncol, std::move(values), std::move(indices), std::move(pointers), true, check) {}
};

}

#endif
