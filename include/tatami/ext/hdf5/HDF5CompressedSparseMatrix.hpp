#ifndef TATAMI_HDF5_SPARSE_MATRIX_HPP
#define TATAMI_HDF5_SPARSE_MATRIX_HPP

#include "H5Cpp.h"

#include <string>
#include <cstdint>
#include <type_traits>
#include <cmath>
#include <list>

#include "../../base/Matrix.hpp"
#include "../../utils/Oracles.hpp"
#include "utils.hpp"

/**
 * @file HDF5CompressedSparseMatrix.hpp
 *
 * @brief Defines a class for a HDF5-backed compressed sparse matrix.
 */

namespace tatami {

/**
 * @brief Compressed sparse matrix in a HDF5 file.
 *
 * This class retrieves sparse data from the HDF5 file on demand rather than loading it all in at the start.
 * This allows us to handle very large datasets in limited memory at the cost of speed.
 *
 * We manually handle the chunk caching to speed up access for consecutive rows or columns (for compressed sparse row and column matrices, respectively).
 * The policy is to minimize the number of calls to the HDF5 library by requesting large contiguous slices where possible.
 * The size of the slice is determined by the cache limit in the constructor.
 *
 * Callers should follow the `prefer_rows()` suggestion when extracting data,
 * as this tries to minimize the number of chunks that need to be read per access request.
 * This recommendation is even stronger than for the `HDF5DenseMatrix`,
 * as the access pattern on disk for the non-preferred dimension is very suboptimal.
 *
 * As the HDF5 library is not generally thread-safe, the HDF5-related operations should only be run in a single thread.
 * For OpenMP, this is handled automatically by putting all HDF5 operations in a critical region.
 * For other parallelization schemes, callers should define the `TATAMI_HDF5_PARALLEL_LOCK` macro;
 * this should be a function that accepts and executes a no-argument lambda within an appropriate serial region (e.g., based on a global mutex).
 *
 * @tparam row_ Whether the matrix is stored in compressed sparse row format.
 * @tparam Value_ Type of the matrix values.
 * @tparam Index_ Type of the row/column indices.
 */
template<bool row_, typename Value_, typename Index_>
class HDF5CompressedSparseMatrix : public tatami::Matrix<Value_, Index_> {
    Index_ nrows, ncols;
    std::string file_name;
    std::string data_name, index_name;
    std::vector<hsize_t> pointers;

    size_t cache_size_limit;
    Index_ max_non_zeros;

public:
    /**
     * @param nr Number of rows in the matrix.
     * @param nc Number of columns in the matrix.
     * @param file Path to the file.
     * @param vals Name of the 1D dataset inside `file` containing the non-zero elements.
     * @param idx Name of the 1D dataset inside `file` containing the indices of the non-zero elements.
     * If `row_ = true`, this should contain column indices sorted within each row, otherwise it should contain row indices sorted within each column.
     * @param ptr Name of the 1D dataset inside `file` containing the index pointers for the start and end of each row (if `row_ = true`) or column (otherwise).
     * This should have length equal to the number of rows (if `row_ = true`) or columns (otherwise).
     * @param cache_limit Limit to the size of the cache, in bytes.
     *
     * The cache is created by extracting multiple columns (for CSC matrices) or rows (CSR) on every call to the HDF5 library.
     * These are held in memory in the `Extractor` while the relevant column/row is returned to the user by `row()` or `column()`.
     * The aim is to minimize the number of calls to the HDF5 library - and thus expensive file reads - for consecutive accesses.
     */
    HDF5CompressedSparseMatrix(Index_ nr, Index_ nc, std::string file, std::string vals, std::string idx, std::string ptr, size_t cache_limit = 100000000) :
        nrows(nr),
        ncols(nc),
        file_name(file),
        data_name(std::move(vals)),
        index_name(std::move(idx)),
        pointers(static_cast<size_t>(row_ ? nr : nc) + 1)
    {
#ifndef TATAMI_HDF5_PARALLEL_LOCK
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        H5::H5File file_handle(file_name, H5F_ACC_RDONLY);

        auto dhandle = HDF5::open_and_check_dataset<false>(file_handle, data_name);
        const Index_ nonzeros = HDF5::get_array_dimensions<1>(dhandle, "vals")[0];

        auto ihandle = HDF5::open_and_check_dataset<true>(file_handle, index_name);
        if (HDF5::get_array_dimensions<1>(ihandle, "idx")[0] != nonzeros) {
            throw std::runtime_error("number of non-zero elements is not consistent between 'data' and 'idx'");
        }

        auto phandle = HDF5::open_and_check_dataset<true>(file_handle, ptr);
        const Index_ ptr_size = HDF5::get_array_dimensions<1>(phandle, "ptr")[0];
        if (ptr_size != pointers.size()) {
            throw std::runtime_error("'ptr' dataset should have length equal to the number of " + (row_ ? std::string("rows") : std::string("columns")) + " plus 1");
        }

        // Checking the contents of the index pointers.
        phandle.read(pointers.data(), H5::PredType::NATIVE_HSIZE);
        if (pointers[0] != 0) {
            throw std::runtime_error("first index pointer should be zero");
        }
        if (pointers.back() != nonzeros) {
            throw std::runtime_error("last index pointer should be equal to the number of non-zero elements");
        }

#ifndef TATAMI_HDF5_PARALLEL_LOCK
        }
#else
        });
#endif

        cache_size_limit = cache_limit;
        max_non_zeros = 0;
        for (size_t i = 1; i < pointers.size(); ++i) {
            Index_ diff = pointers[i] - pointers[i-1];
            if (diff > max_non_zeros) {
                max_non_zeros = diff;
            }
        }
    }

public:
    Index_ nrow() const {
        return nrows;
    }

    Index_ ncol() const {
        return ncols;
    }

    /**
     * @return `true`.
     */
    bool sparse() const {
        return true;
    }

    /**
     * @return `true` if this is in compressed sparse row format.
     */
    bool prefer_rows() const {
        return row_;
    }

    bool uses_oracle(bool) const {
        return false; // placeholder for proper support.
    }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

    /********************************************
     ************ Primary extraction ************
     ********************************************/
private:
    struct OracleCache {
        struct Element {
            size_t data_offset;
            size_t mem_offset;
            Index_ length;
            bool bounded;
        };

        std::vector<Value_> cache_value;
        std::vector<Index_> cache_index;
        std::unordered_map<Index_, Index_> cache_exists, next_cache_exists;
        std::vector<Element> cache_data, next_cache_data;

        OracleStream<Index_> prediction_stream;
        std::vector<Index_> predictions_made;
        std::vector<Index_> needed;
        std::vector<Index_> present;
        size_t predictions_fulfilled = 0;

        size_t max_cache_elements = -1;
    };

    struct LruCache {
        struct Element {
            std::vector<Value_> value;
            std::vector<Index_> index;
            Index_ length;
            Index_ id;
            bool bounded;
        };

        std::list<Element> cache_data;
        std::unordered_map<Index_, typename std::list<Element>::iterator> cache_exists;
        size_t max_cache_number = -1;
    };

    struct PrimaryWorkspace {
        void fill(const HDF5CompressedSparseMatrix* parent, Index_ extraction_cache_size) {
            // TODO: set more suitable chunk cache values here, to avoid re-reading
            // chunks on the boundaries of the primary cache.
            file.openFile(parent->file_name, H5F_ACC_RDONLY);

            data = file.openDataSet(parent->data_name);
            index = file.openDataSet(parent->index_name);
            dataspace = data.getSpace();

            extraction_bounds.resize(extraction_cache_size, std::pair<size_t, size_t>(-1, 0));

            historian.reset(new LruCache);
        }

    public:
        H5::H5File file;
        H5::DataSet data, index;
        H5::DataSpace dataspace;
        H5::DataSpace memspace;

    public:
        // Cache with an oracle.
        std::unique_ptr<OracleCache> futurist;

        // Cache without an oracle.
        std::unique_ptr<LruCache> historian;

    public:
        // Cache for re-use.
        std::vector<std::pair<size_t, size_t> > extraction_bounds;
    };

private:
    struct Extracted {
        Extracted() = default;

        Extracted(const typename LruCache::Element& cache) {
            value = cache.value.data();
            index = cache.index.data();
            length = cache.length;
            bounded = cache.bounded;
        }

        Extracted(const OracleCache& cache, Index_ i, bool needs_value) {
            const auto& element = cache.cache_data[i];
            auto offset = element.mem_offset;
            if (needs_value) {
                value = cache.cache_value.data() + offset;
            }
            index = cache.cache_index.data() + offset;
            length = element.length;
            bounded = element.bounded;
        }

        const Value_* value;
        const Index_* index;
        Index_ length;
        bool bounded;
    };

    Extracted extract_primary_without_oracle(Index_ i, PrimaryWorkspace& work, bool needs_value) const {
        auto& historian = *(work.historian);
        auto it = historian.cache_exists.find(i);
        if (it != historian.cache_exists.end()) {
            auto chosen = it->second;
            historian.cache_data.splice(historian.cache_data.end(), historian.cache_data, chosen); // move to end.
            return Extracted(*chosen);
        }

        // Check if bounds already exist from the reusable cache. If so,
        // we can use them to reduce the amount of data I/O.
        hsize_t extraction_start = pointers[i];
        hsize_t extraction_len = pointers[i + 1] - pointers[i];
        bool bounded = false;

        if (work.extraction_bounds.size()) {
            const auto& current = work.extraction_bounds[i];
            if (current.first != -1) {
                bounded = true;
                extraction_start = current.first;
                extraction_len = current.second;
            }
        }

        // Need to use the max_cache_number as the capacity of each recycled
        // vector may be much larger than the reported size of the chunk; this
        // would cause us to overrun the cache_size_limit if we recycled the
        // vector enough times such that each cache element had capacity equal
        // to the maximum number of non-zero elements in any dimension element.
        //
        // Alternatives would be to create a new Element on every recycling
        // iteration, or to hope that shrink_to_fit() behaves. Both would allow
        // us to store more cache elements but would involve reallocations,
        // which degrades perf in the most common case where a dimension
        // element is accessed no more than once during iteration.
        if (historian.max_cache_number == -1) {
            historian.max_cache_number = cache_size_limit / (max_non_zeros * ((needs_value ? sizeof(Value_) : 0) + sizeof(Index_)));
            if (historian.max_cache_number == 0) {
                historian.max_cache_number = 1;
            }
        }

        // Adding a new last element, or recycling the front to the back.
        // We initialize each element with the maximum number of non-zeros to avoid reallocations.
        typename std::list<typename LruCache::Element>::iterator location;
        if (historian.cache_data.size() < historian.max_cache_number) {
            historian.cache_data.push_back(typename LruCache::Element());
            location = std::prev(historian.cache_data.end());
            location->index.resize(max_non_zeros);
            if (needs_value) {
                location->value.resize(max_non_zeros);
            }
        } else {
            location = historian.cache_data.begin();
            historian.cache_exists.erase(location->id);
            historian.cache_data.splice(historian.cache_data.end(), historian.cache_data, location); // move to end.
        }
        historian.cache_exists[i] = location;

        auto& current_cache = *location;
        current_cache.id = i;
        current_cache.length = extraction_len;
        current_cache.bounded = bounded;

#ifndef TATAMI_HDF5_PARALLEL_LOCK
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        work.dataspace.selectHyperslab(H5S_SELECT_SET, &extraction_len, &extraction_start);
        work.memspace.setExtentSimple(1, &extraction_len);
        work.memspace.selectAll();
        work.index.read(current_cache.index.data(), HDF5::define_mem_type<Index_>(), work.memspace, work.dataspace);

        if (needs_value) {
            work.data.read(current_cache.value.data(), HDF5::define_mem_type<Value_>(), work.memspace, work.dataspace);
        }

#ifndef TATAMI_HDF5_PARALLEL_LOCK
        }
#else
        });
#endif

        return Extracted(current_cache);
    }

    template<class Function_>
    static void sort_by_field(std::vector<int>& indices, Function_ field) {
        auto comp = [&field](size_t l, size_t r) -> bool {
            return field(l) < field(r);
        };

        bool sorted = true;
        for (size_t i = 1, end = indices.size(); i < end; ++i) {
            if (!comp(indices[i-1], indices[i])) {
                sorted = false;
                break;
            }
        }

        if (!sorted) {
            std::sort(indices.begin(), indices.end(), comp);
        }
    }

    Extracted extract_primary_with_oracle(PrimaryWorkspace& work, bool needs_value) const {
        auto& pred = *work.futurist;
        if (pred.predictions_made.size() > pred.predictions_fulfilled) {
            auto chosen = pred.predictions_made[pred.predictions_fulfilled++];
            return Extracted(pred, chosen, needs_value);
        }

        // Grow the number of predictions over time, until we get to a point
        // where we consistently fill the cache.
        size_t max_predictions = pred.predictions_made.size() * 2; 
        if (max_predictions < 100) {
            max_predictions = 100;
        } else {
            size_t upper = (row_ ? nrows : ncols);
            if (max_predictions > upper) {
                max_predictions = upper;
            }
        }

        pred.predictions_made.clear();
        pred.needed.clear();
        pred.present.clear();
        pred.next_cache_data.clear();
        pred.next_cache_exists.clear();

        // Here, we use a giant contiguous buffer to optimize for
        // near-consecutive iteration. This allows the HDF5 library to pull out
        // long strips of data from the file.  It also allows us to maximize
        // the use of the cache_size_limit by accounting for differences in the
        // non-zeros for each element, rather than conservatively assuming
        // they're all at max (as in the LRU case). The downside is that we
        // need to do some copying within the cache to make space for new
        // reads, but that works out to be no more than one extra copy per
        // fetch() call, which is tolerable. I suppose we could do better
        // by defragmenting within this buffer but that's probably overkill.
        if (pred.max_cache_elements == -1) {
            pred.max_cache_elements = cache_size_limit / ((needs_value ? sizeof(Value_) : 0) + sizeof(Index_));
            if (pred.max_cache_elements < max_non_zeros) {
                pred.max_cache_elements = max_non_zeros; // make sure we have enough space to store the largest possible primary dimension element.
            }
            pred.cache_index.resize(pred.max_cache_elements);
            if (needs_value) {
                pred.cache_value.resize(pred.max_cache_elements);
            }
        }
        size_t filled_elements = 0;

        for (size_t p = 0; p < max_predictions; ++p) {
            Index_ current;
            if (!pred.prediction_stream.next(current)) {
                break;
            }

            // Seeing if this element already exists somewhere.
            auto nit = pred.next_cache_exists.find(current);
            if (nit != pred.next_cache_exists.end()) {
                pred.predictions_made.push_back(nit->second);
                continue;
            }

            auto it = pred.cache_exists.find(current);
            if (it != pred.cache_exists.end()) {
                auto& candidate = pred.cache_data[it->second];
                filled_elements += candidate.length;
                if (filled_elements > pred.max_cache_elements) {
                    pred.prediction_stream.back();
                    break;
                }

                Index_ used = pred.next_cache_data.size();
                pred.predictions_made.push_back(used);
                pred.present.push_back(used);
                pred.next_cache_exists[current] = used;
                pred.next_cache_data.push_back(std::move(candidate));
                continue;
            }

            // Check if bounds already exist from the reusable cache. If so,
            // we can use them to reduce the amount of data I/O.
            hsize_t extraction_start = pointers[current];
            hsize_t extraction_len = pointers[current + 1] - pointers[current];
            bool bounded = false;

            if (work.extraction_bounds.size()) {
                const auto& bounds = work.extraction_bounds[current];
                if (bounds.first != -1) {
                    bounded = true;
                    extraction_start = bounds.first;
                    extraction_len = bounds.second;
                }
            }

            filled_elements += extraction_len;
            if (filled_elements > pred.max_cache_elements) {
                pred.prediction_stream.back();
                break;
            }

            Index_ used = pred.next_cache_data.size();
            pred.predictions_made.push_back(used);
            pred.needed.emplace_back(used);
            pred.next_cache_exists[current] = used;

            typename OracleCache::Element latest;
            latest.data_offset = extraction_start;
            latest.length = extraction_len;
            latest.bounded = bounded;
            pred.next_cache_data.push_back(std::move(latest));
        }

        if (pred.needed.size()) {
            size_t dest_offset = 0;

            if (pred.present.size()) {
                // Shuffling all re-used elements to the start of the buffer,
                // so that we can perform a contiguous extraction of the needed
                // elements in the rest of the buffer. This needs some sorting
                // to ensure that we're not clobbering one re-used element's
                // contents when shifting another element to the start.
                sort_by_field(pred.present, [&pred](size_t i) -> size_t { return pred.next_cache_data[i].mem_offset; });

                for (const auto& p : pred.present) {
                    auto& info = pred.next_cache_data[p];
                    auto isrc = pred.cache_index.begin() + info.mem_offset;

#ifdef DEBUG
                    if (info.mem_offset < dest_offset) {
                        throw std::runtime_error("detected clobbering of memory cache from overlapping offsets");
                    }
#endif

                    std::copy(isrc, isrc + info.length, pred.cache_index.begin() + dest_offset);
                    if (needs_value) {
                        auto vsrc = pred.cache_value.begin() + info.mem_offset;
                        std::copy(vsrc, vsrc + info.length, pred.cache_value.begin() + dest_offset); 
                    }
                    info.mem_offset = dest_offset;
                    dest_offset += info.length;
                }
            }

            // Sorting so that we get consecutive accesses in the hyperslab construction.
            // This should improve re-use of partially read chunks inside the HDF5 call.
            sort_by_field(pred.needed, [&pred](size_t i) -> size_t { return pred.next_cache_data[i].data_offset; });

#ifndef TATAMI_HDF5_PARALLEL_LOCK
            #pragma omp critical
            {
#else
            TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

            size_t sofar = 0;
            hsize_t combined_len = 0;
            work.dataspace.selectNone();

            while (sofar < pred.needed.size()) {
                auto& first = pred.next_cache_data[pred.needed[sofar]];
                first.mem_offset = dest_offset + combined_len;
                hsize_t src_offset = first.data_offset;
                hsize_t len = first.length;
                ++sofar;

                // Finding the stretch of consecutive extractions, and bundling them into a single hyperslab.
                for (; sofar < pred.needed.size(); ++sofar) {
                    auto& next = pred.next_cache_data[pred.needed[sofar]];
                    if (src_offset + len < next.data_offset) {
                        break;
                    }
                    next.mem_offset = first.mem_offset + len;
                    len += next.length;
                }

                work.dataspace.selectHyperslab(H5S_SELECT_OR, &len, &src_offset);
                combined_len += len;
            }

            work.memspace.setExtentSimple(1, &combined_len);
            work.memspace.selectAll();
            work.index.read(pred.cache_index.data() + dest_offset, HDF5::define_mem_type<Index_>(), work.memspace, work.dataspace);
            if (needs_value) {
                work.data.read(pred.cache_value.data() + dest_offset, HDF5::define_mem_type<Value_>(), work.memspace, work.dataspace);
            }

#ifndef TATAMI_HDF5_PARALLEL_LOCK
            }
#else
            });
#endif
        }

        pred.cache_data.swap(pred.next_cache_data);
        pred.cache_exists.swap(pred.next_cache_exists);
        pred.predictions_fulfilled = 1; // using the first one now.
        return Extracted(pred, pred.predictions_made.front(), needs_value);
    }

    /********************************************
     ************ Primary extraction ************
     ********************************************/
private:
    template<class Function_>
    void extract_primary_raw(size_t i, Function_ fill, Index_ start, PrimaryWorkspace& work, bool needs_value) const {
        Extracted details;
        if (work.futurist) {
            details = extract_primary_with_oracle(work, needs_value);
        } else {
            details = extract_primary_without_oracle(i, work, needs_value);
        }

        auto istart = details.index;
        auto iend = details.index + details.length;
        size_t offset = 0;

        // If we used the extraction_bounds during extraction, there's no need
        // to do another search. Similarly, if we didn't use the extraction_bounds
        // (e.g., it was already cached) but we have extraction_bounds available,
        // we can again skip the binary search.
        if (!details.bounded && start) {
            bool hit = false;
            if (work.extraction_bounds.size()) {
                auto& target = work.extraction_bounds[i];
                if (target.first != -1) {
                    hit = true;
                    offset = target.first - pointers[i];
                    istart += offset;
                }
            } 

            if (!hit) {
                istart = std::lower_bound(details.index, iend, start);
                offset = istart - details.index;
            }
        }

        size_t iterated = fill(istart, iend, (needs_value ? details.value + offset : NULL));

        if (work.extraction_bounds.size()) {
            auto& target = work.extraction_bounds[i];
            if (target.first == -1) {
                target.first = pointers[i] + offset;
                target.second = iterated;
            }
        }

        return;
    }

    const Value_* extract_primary(size_t i, Value_* buffer, Index_ start, Index_ length, PrimaryWorkspace& work) const {
        std::fill(buffer, buffer + length, 0);

        if (length) {
            extract_primary_raw(i, 

                [&](const Index_* is, const Index_* ie, const Value_* vs) -> size_t {
                    auto ioriginal = is;
                    Index_ end = start + length;
                    for (; is != ie && *is < end; ++is, ++vs) {
                        buffer[*is - start] = *vs;
                    }
                    return is - ioriginal;
                },

                start, 
                work,
                true
            );
        }

        return buffer;
    }

    SparseRange<Value_, Index_> extract_primary(size_t i, Value_* dbuffer, Index_* ibuffer, Index_ start, Index_ length, PrimaryWorkspace& work, bool needs_value, bool needs_index) const {
        Index_ counter = 0;

        if (length) {
            extract_primary_raw(i, 

                [&](const Index_* is, const Index_* ie, const Value_* vs) -> Index_ {
                    auto ioriginal = is;
                    Index_ end = start + length;
                    for (; is != ie && *is < end; ++is) {
                        ++counter;                        
                    }

                    if (needs_index) {
                        std::copy(ioriginal, ioriginal + counter, ibuffer);
                    }
                    if (needs_value) {
                        std::copy(vs, vs + counter, dbuffer);
                    }

                    return counter;
                },

                start, 
                work,
                needs_value
            );
        }

        if (!needs_value) {
            dbuffer = NULL;
        }
        if (!needs_index) {
            ibuffer = NULL;
        }

        return SparseRange<Value_, Index_>(counter, dbuffer, ibuffer);
    }

private:
    template<class Fill_, class Skip_>
    static size_t indexed_extraction(const Index_* istart, const Index_* iend, const Value_* vstart, bool needs_value, const std::vector<Index_>& indices, Fill_ fill, Skip_ skip) {
        auto ioriginal = istart;
        if (needs_value) {
            for (auto idx : indices) {
                while (istart != iend && *istart < idx) {
                    ++istart;
                    ++vstart;
                }
                if (istart == iend) {
                    break;
                }
                if (*istart == idx) {
                    fill(idx, *vstart);
                    ++istart;
                    ++vstart;
                } else {
                    skip();
                }
            }
        } else {
            for (auto idx : indices) {
                while (istart != iend && *istart < idx) {
                    ++istart;
                }
                if (istart == iend) {
                    break;
                }
                if (*istart == idx) {
                    fill(idx, 0);
                    ++istart;
                } else {
                    skip();
                }
            }
        }

        return istart - ioriginal;
    }

    const Value_* extract_primary(size_t i, Value_* buffer, const std::vector<Index_>& indices, PrimaryWorkspace& work) const {
        std::fill(buffer, buffer + indices.size(), 0);
        auto original = buffer;

        if (indices.size()) {
            extract_primary_raw(i, 

                [&](const Index_* is, const Index_* ie, const Value_* vs) -> size_t {
                    return indexed_extraction(is, ie, vs, true, indices, 
                        [&](Index_, Value_ value) -> void {
                            *buffer = value;
                            ++buffer;
                        },
                        [&]() -> void {
                            ++buffer;
                        }
                    );
                },

                indices.front(),
                work,
                true
            );
        }

        return original;
    }

    SparseRange<Value_, Index_> extract_primary(size_t i, Value_* dbuffer, Index_* ibuffer, const std::vector<Index_>& indices, PrimaryWorkspace& work, bool needs_value, bool needs_index) const {
        Index_ counter = 0;

        if (indices.size()) {
            extract_primary_raw(i, 

                [&](const Index_* is, const Index_* ie, const Value_* vs) -> size_t {
                    return indexed_extraction(is, ie, vs, needs_value, indices,
                        [&](Index_ pos, Value_ value) -> void {
                            if (needs_value) {
                                dbuffer[counter] = value;
                            }
                            if (needs_index) {
                                ibuffer[counter] = pos;
                            }
                            ++counter;
                        },
                        []() -> void {}
                    );
                },

                indices.front(),
                work,
                needs_value
            );
        }

        if (!needs_index) {
            ibuffer = NULL;
        }
        if (!needs_value) {
            dbuffer = NULL;
        }

        return SparseRange<Value_, Index_>(counter, dbuffer, ibuffer);
    }

    /**********************************************
     ************ Secondary extraction ************
     **********************************************/
private:
    // This could be improved by extracting multiple rows at any given call and
    // caching them for subsequent requests. However, even then, we'd require
    // multiple re-reads from file when we exceed the cache. So, any caching
    // would be just turning an extremely bad access pattern into a very bad
    // pattern, when users shouldn't even be calling this at all... 
    struct SecondaryWorkspace {
        H5::H5File file;
        H5::DataSet data, index;
        H5::DataSpace dataspace;
        H5::DataSpace memspace;

        void fill(const HDF5CompressedSparseMatrix* parent) {
            // TODO: set more suitable chunk cache values here, to avoid re-reading
            // chunks on the boundaries of the primary cache.
            file.openFile(parent->file_name, H5F_ACC_RDONLY);

            data = file.openDataSet(parent->data_name);
            index = file.openDataSet(parent->index_name);
            dataspace = data.getSpace();
        }

        std::vector<Index_> index_cache;
    };

    template<class Function_>
    bool extract_secondary_raw(Index_ primary, Index_ secondary, Function_& fill, SecondaryWorkspace& core, bool needs_value) const {
        hsize_t left = pointers[primary], right = pointers[primary + 1];
        core.index_cache.resize(right - left);

        // Serial locks should be applied by the callers.
        hsize_t offset = left;
        hsize_t count = core.index_cache.size();
        core.dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
        core.memspace.setExtentSimple(1, &count);
        core.memspace.selectAll();
        core.index.read(core.index_cache.data(), HDF5::define_mem_type<Index_>(), core.memspace, core.dataspace);

        auto it = std::lower_bound(core.index_cache.begin(), core.index_cache.end(), secondary);
        if (it != core.index_cache.end() && *it == secondary) {
            if (needs_value) {
                offset = left + (it - core.index_cache.begin());
                count = 1;
                core.dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
                core.memspace.setExtentSimple(1, &count);
                core.memspace.selectAll();

                Value_ dest;
                core.data.read(&dest, HDF5::define_mem_type<Value_>(), core.memspace, core.dataspace);
                fill(primary, dest);
            } else {
                fill(primary, 0);
            }
            return true;
        } else {
            return false;
        }
    }

    template<class Function_>
    void extract_secondary_raw_loop(size_t i, Function_ fill, Index_ start, Index_ length, SecondaryWorkspace& core, bool needs_value) const {
#ifndef TATAMI_HDF5_PARALLEL_LOCK
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        Index_ end = start + length;
        for (size_t j = start; j < end; ++j) {
            extract_secondary_raw(j, i, fill, core, needs_value);
        }

#ifndef TATAMI_HDF5_PARALLEL_LOCK 
        }
#else
        });
#endif
    }

    const Value_* extract_secondary(size_t i, Value_* buffer, Index_ start, Index_ length, SecondaryWorkspace& core) const {
        std::fill(buffer, buffer + length, 0);

        extract_secondary_raw_loop(i, 
            [&](Index_ pos, Value_ value) -> void {
                buffer[pos - start] = value;
            }, 
            start, 
            length, 
            core,
            true
        );

        return buffer;
    }

    SparseRange<Value_, Index_> extract_secondary(size_t i, Value_* dbuffer, Index_* ibuffer, Index_ start, Index_ length, SecondaryWorkspace& core, bool needs_value, bool needs_index) const {
        Index_ counter = 0;

        extract_secondary_raw_loop(i, 
            [&](Index_ pos, Value_ value) -> void {
                if (needs_value) {
                    dbuffer[counter] = value;
                }
                if (needs_index) {
                    ibuffer[counter] = pos;
                }
                ++counter;
            }, 
            start, 
            length, 
            core,
            needs_value
        );

        if (!needs_value) {
            dbuffer = NULL;
        }
        if (!needs_index) {
            ibuffer = NULL;
        }

        return SparseRange<Value_, Index_>(counter, dbuffer, ibuffer);
    }

    template<class Function_, class Skip_>
    void extract_secondary_raw_loop(size_t i, Function_ fill, Skip_ skip, const std::vector<Index_>& indices, SecondaryWorkspace& core, bool needs_value) const {
#ifndef TATAMI_HDF5_PARALLEL_LOCK
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        for (auto j : indices) {
            if (!extract_secondary_raw(j, i, fill, core, needs_value)) {
                skip();
            }
        }

#ifndef TATAMI_HDF5_PARALLEL_LOCK 
        }
#else
        });
#endif
    }

    const Value_* extract_secondary(size_t i, Value_* buffer, const std::vector<Index_>& indices, SecondaryWorkspace& core) const {
        std::fill(buffer, buffer + indices.size(), 0);
        auto original = buffer;
        extract_secondary_raw_loop(i, 
            [&](Index_ pos, Value_ value) -> void {
                *buffer = value;
                ++buffer;
            }, 
            [&]() -> void {
                ++buffer;
            },
            indices, 
            core,
            true
        );
        return original;
    }

    SparseRange<Value_, Index_> extract_secondary(size_t i, Value_* dbuffer, Index_* ibuffer, const std::vector<Index_>& indices, SecondaryWorkspace& core, bool needs_value, bool needs_index) const {
        Index_ counter = 0;

        extract_secondary_raw_loop(i, 
            [&](Index_ pos, Value_ value) -> void {
                if (needs_value) {
                    dbuffer[counter] = value;
                }
                if (needs_index) {
                    ibuffer[counter] = pos;
                }
                ++counter;
            }, 
            []() -> void {},
            indices, 
            core,
            needs_value
        );

        if (!needs_value) {
            dbuffer = NULL;
        }
        if (!needs_index) {
            ibuffer = NULL;
        }

        return SparseRange<Value_, Index_>(counter, dbuffer, ibuffer);
    }

    /******************************************
     ************ Public overrides ************
     ******************************************/
private:
    template<bool accrow_, DimensionSelectionType selection_, bool sparse_>
    struct Hdf5SparseExtractor : public Extractor<selection_, sparse_, Value_, Index_> {
        Hdf5SparseExtractor(const HDF5CompressedSparseMatrix* p, const Options& opt) : parent(p) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = (accrow_ ? parent->ncols : parent->nrows);
            }

            if constexpr(row_ == accrow_) {
                if (opt.cache_for_reuse) {
                    core.fill(parent, accrow_ ? parent->nrows : parent->ncols);
                } else {
                    core.fill(parent, 0);
                }
            } else {
                core.fill(parent);
            }
        }

        Hdf5SparseExtractor(const HDF5CompressedSparseMatrix* p, const Options& opt, Index_ bs, Index_ bl) : Hdf5SparseExtractor(p, opt) {
            if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = bs;
                this->block_length = bl;
            }
        }

        Hdf5SparseExtractor(const HDF5CompressedSparseMatrix* p, const Options& opt, std::vector<Index_> idx) : Hdf5SparseExtractor(p, opt) {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                this->index_length = idx.size();
                indices = std::move(idx);
            }
        }

    protected:
        const HDF5CompressedSparseMatrix* parent;
        typename std::conditional<row_ == accrow_, PrimaryWorkspace, SecondaryWorkspace>::type core;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type indices;

    public:
        const Index_* index_start() const {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                return indices.data();
            } else {
                return NULL;
            }
        }

    public:
        void set_oracle(std::unique_ptr<Oracle<Index_> > o) {
            if constexpr(row_ == accrow_) {
                core.futurist.reset(new OracleCache);
                core.futurist->prediction_stream.set(std::move(o));
                core.historian.reset();
            }
        }
    };

    template<bool accrow_, DimensionSelectionType selection_>
    struct DenseHdf5SparseExtractor : public Hdf5SparseExtractor<accrow_, selection_, false> {
        template<typename... Args_>
        DenseHdf5SparseExtractor(const HDF5CompressedSparseMatrix* p, const Options& opt, Args_... args) : 
            Hdf5SparseExtractor<accrow_, selection_, false>(p, opt, std::move(args)...) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                if constexpr(row_ == accrow_) {
                    return this->parent->extract_primary(i, buffer, 0, this->full_length, this->core);
                } else {
                    return this->parent->extract_secondary(i, buffer, 0, this->full_length, this->core);
                }
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                if constexpr(row_ == accrow_) {
                    return this->parent->extract_primary(i, buffer, this->block_start, this->block_length, this->core);
                } else {
                    return this->parent->extract_secondary(i, buffer, this->block_start, this->block_length, this->core);
                }
            } else {
                if constexpr(row_ == accrow_) {
                    return this->parent->extract_primary(i, buffer, this->indices, this->core);
                } else {
                    return this->parent->extract_secondary(i, buffer, this->indices, this->core);
                }
            }
        }
    };

    template<bool accrow_, DimensionSelectionType selection_>
    struct SparseHdf5SparseExtractor : public Hdf5SparseExtractor<accrow_, selection_, true> {
        template<typename... Args_>
        SparseHdf5SparseExtractor(const HDF5CompressedSparseMatrix* p, const Options& opt, Args_... args) : 
            Hdf5SparseExtractor<accrow_, selection_, true>(p, opt, std::move(args)...), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                if constexpr(row_ == accrow_) {
                    if (needs_index || needs_value) {
                        return this->parent->extract_primary(i, vbuffer, ibuffer, 0, this->full_length, this->core, needs_value, needs_index);
                    } else {
                        // Quick return is possible if we don't need any indices or values.
                        return SparseRange<Value_, Index_>(this->parent->pointers[i+1] - this->parent->pointers[i], NULL, NULL);
                    }
                } else {
                    return this->parent->extract_secondary(i, vbuffer, ibuffer, 0, this->full_length, this->core, needs_value, needs_index);
                }
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                if constexpr(row_ == accrow_) {
                    return this->parent->extract_primary(i, vbuffer, ibuffer, this->block_start, this->block_length, this->core, needs_value, needs_index);
                } else {
                    return this->parent->extract_secondary(i, vbuffer, ibuffer, this->block_start, this->block_length, this->core, needs_value, needs_index);
                }
            } else {
                if constexpr(row_ == accrow_) {
                    return this->parent->extract_primary(i, vbuffer, ibuffer, this->indices, this->core, needs_value, needs_index);
                } else {
                    return this->parent->extract_secondary(i, vbuffer, ibuffer, this->indices, this->core, needs_value, needs_index);
                }
            }
        }

    protected:
        bool needs_value;
        bool needs_index;
    };

    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > populate(const Options& opt, Args_... args) const {
        std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > output;

#ifndef TATAMI_HDF5_PARALLEL_LOCK
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        if constexpr(sparse_) {
            output.reset(new SparseHdf5SparseExtractor<accrow_, selection_>(this, opt, std::move(args)...));
        } else {
            output.reset(new DenseHdf5SparseExtractor<accrow_, selection_>(this, opt, std::move(args)...));
        }

#ifndef TATAMI_HDF5_PARALLEL_LOCK
        }
#else
        });
#endif

        return output;
    }

public:
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        return populate<true, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<true, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        return populate<true, DimensionSelectionType::INDEX, false>(opt, std::move(indices));
    }

    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        return populate<false, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<false, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        return populate<false, DimensionSelectionType::INDEX, false>(opt, std::move(indices));
    }

public:
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_row(const Options& opt) const {
        return populate<true, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<true, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& opt) const {
        return populate<true, DimensionSelectionType::INDEX, true>(opt, std::move(indices));
    }

    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_column(const Options& opt) const {
        return populate<false, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<false, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options& opt) const {
        return populate<false, DimensionSelectionType::INDEX, true>(opt, std::move(indices));
    }
};

}

#endif
