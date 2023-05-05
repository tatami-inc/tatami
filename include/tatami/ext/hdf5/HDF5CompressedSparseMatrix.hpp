#ifndef TATAMI_HDF5_SPARSE_MATRIX_HPP
#define TATAMI_HDF5_SPARSE_MATRIX_HPP

#include "H5Cpp.h"

#include <string>
#include <cstdint>
#include <type_traits>
#include <cmath>

#include "../../base/Matrix.hpp"
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

    std::vector<Index_> primary_cache_id;
    std::vector<std::pair<Index_, Index_> > primary_cache_limits;

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
        pointers(row_ ? nr + 1 : nc + 1)
    {
        size_t total_element_size;

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

        // Setting up the cache parameters.
        Index_ primary_dim = pointers.size() - 1;
        primary_cache_id.resize(primary_dim);
        if (primary_dim) {
            primary_cache_limits.emplace_back(0, pointers[1]);
        }

        size_t effective_cache_limit = cache_limit / (sizeof(Value_) + sizeof(Index_));
        Index_ counter = 0, start = 0;
        for (size_t i = 1; i < primary_dim; ++i) {
            Index_ end = pointers[i + 1];
            if (end - start <= effective_cache_limit) {
                primary_cache_limits.back().second = end;
            } else {
                ++counter;
                start = pointers[i];
                primary_cache_limits.emplace_back(start, end);
            }
            primary_cache_id[i] = counter;
        }

        return;
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
    struct CacheElement {
        std::vector<Value_> value;
        std::vector<Index_> index;
        size_t index_offset = 0;
        size_t valid_length = 0;
        size_t byte_size = 0;
        bool init = false;
    };
    
    template<bool accrow_>
    struct OracleCache {
        std::unordered_map<Index_, Index_> cache_exists, next_cache_exists;
        std::vector<CacheElement> cache_data, next_cache_data;
        std::vector<std::pair<Index_, Index_> > chunks_in_need; 

        OracleStream<Index_> prediction_stream;
        std::vector<Index_> predictions_made;
        size_t predictions_fulfilled = 0;
    };

    struct LruCache {
        typedef std::pair<CacheElement, Index_> Element;
        std::list<Element> cache_data;
        std::unordered_map<Index_, typename std::list<Element>::iterator> cache_exists;
        size_t current_cache_size;
    };

    struct Workspace {
        void fill(const HDF5DenseMatrix* parent, Index_ cache_size) {
            // TODO: set more suitable chunk cache values here, to avoid re-reading
            // chunks on the boundaries of the primary cache.
            file.openFile(parent->file_name, H5F_ACC_RDONLY);

            data = core.file.openDataSet(parent->data_name);
            index = core.file.openDataSet(parent->index_name);
            dataspace = core.data.getSpace();

            extraction_bounds.resize(cache_size, std::pair<size_t, size_t>(-1, 0));

            historian.reset(new LruCache);
        }

    public:
        H5::H5File file;
        H5::DataSet data, index;
        H5::DataSpace dataspace;
        H5::DataSpace memspace;

    public:
        // Caching members.
        size_t cache_size_limit;

        // Cache with an oracle.
        std::unique_ptr<OracleCache<accrow_> > futurist;

        // Cache without an oracle.
        std::unique_ptr<LruCache> historian;

    public:
        // Cache for re-use.
        std::vector<std::pair<size_t, size_t> > extraction_bounds;
    };

private:
    void extract_base(Index_ i, CacheElement& cache, Workspace& work, Index_ start_index, Index_ end_index, bool needs_value) const {
        auto start_ptr = pointers[i], len = pointers[i+1] - start_ptr;
        cache.index.resize(len); // this should be full length, to avoid discrepancies when cache_for_reuse = true.

        bool excache_hit = false;
        if (work.extraction_bounds.size()) {
            const auto& current = work.extraction_bounds[i];
            if (current.first != -1) {
                excache_hit = true;
                len = current.second;
                start_ptr += cache.index_offset;
            }
        }

        dataspace.selectHyperslab(H5S_SELECT_SET, &len, &start_ptr);
        memspace.setExtentSimple(1, &len);
        memspace.selectAll();
        index.read(cache.index.data(), HDF5::define_mem_type<Index_>(), core.memspace, core.dataspace);

        cache.byte_size = len * sizeof(Index_); // for cache usage calculations.
        cache.index_offset = 0; // to avoid another binary search in callers.

        // Now we decide whether or not to constrict the retrieval zone for values (and for future indices).
        if (!excache_hit) {
            if (start_index) {
                cache.index_offset = std::lower_bound(cache.index.data(), cache.index.data() + len, start_index) - cache.index.data();
                start_ptr += cache.index_offset;
            }
            if (end_index < (row_ ? ncols : nrows)) {
                auto start = cache.index.data(), + cache.index_offset;
                len = std::lower_bound(start, cache.index.data() + len, end_index) - start;
            }

            if (needs_value) {
                dataspace.selectHyperslab(H5S_SELECT_SET, &len, &start_ptr);
                memspace.setExtentSimple(1, &len);
                memspace.selectAll();
            }
            if (work.extraction_bounds.size()) {
                work.extraction_bounds[i] = std::pair<size_t, size_t>(cache.index_offset, len);
            }
        }

        // Storing the final number of non-zeros.
        cache.valid_length = len;

        if (needs_value) {
            cache.value.resize(len);
            data.read(cache.value.data(), HDF5::define_mem_type<Value_>(), core.memspace, core.dataspace);
            cache.byte_size += len * sizeof(Value_);
        }
    }

private:
    CacheElement& populate_cache_without_oracle(Index_ i, Workspace& work, Index_ start_index, Index_ end_index, bool needs_value) const {
        auto& historian = *(work.historian);
        auto it = historian.cache_exists.find(chunk);
        if (it != historian.cache_exists.end()) {
            auto chosen = it->second;
            historian.cache_data.splice(historian.cache_data.end(), historian.cache_data, chosen); // move to end.
            return chosen->first;
        }

        // Estimating the size of the next cache element.
        size_t next_element_size = (pointers[i + 1] - pointers[i]) * (sizeof(Value_) + sizeof(Index_));

        // Adding a new last element, or recycling the front to the back.
        typename std::list<typename LruCache::Element>::iterator location;
        if (historian.current_cache_size + next_element_size <= work.cache_size_limit) {
            historian.cache_data.emplace_back(std::vector<Value_>(work.per_cache_size), chunk);
            location = std::prev(historian.cache_data.end());
        } else {
            location = historian.cache_data.begin();
            historian.cache_exists.erase(location->second);
            work.current_cache_size -= location->first.byte_size;
            location->second = chunk;
            historian.cache_data.splice(historian.cache_data.end(), historian.cache_data, location); // move to end.
        }
        historian.cache_exists[chunk] = location;

#ifndef TATAMI_HDF5_PARALLEL_LOCK
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

            extract_base(i, location->first, work, start_index, end_index, needs_value);

#ifndef TATAMI_HDF5_PARALLEL_LOCK
        }
#else
        });
#endif

        historian.current_cache_size += (location->first).byte_size;
        return location->first;
    }

    const CacheElement& extract_with_oracle(Workspace& work, Index_ start_index, Index_ end_index, bool needs_value)
        auto& pred = *work.futurist;
        if (pred.predictions_made.size() > pred.predictions_fulfilled) {
            auto chosen = pred.predictions_made[pred.predictions_fulfilled++];
            return pred.cache_data[chosen];
        }

        // Same logic as in the dense case.
        Index_ used = 0;

        bool starting = pred.cache_data.empty();
        if (starting) {
            pred.cache_data.resize(num_chunks);
            pred.next_cache_data.resize(num_chunks); 
            for (auto& x : pred.next_cache_data) {
                x.init = true;
            }
            pred.chunks_in_need.reserve(num_chunks);
        } else {
            pred.next_cache_exists.clear();
            pred.chunks_in_need.clear();
        }

        size_t max_predictions = static_cast<size_t>(num_chunks) * chunk_mydim * 2; // double the cache size, basically.
        pred.predictions_made.clear();
        pred.predictions_made.reserve(max_predictions);
        size_t anticipated_size = 0;

        for (size_t p = 0; p < max_predictions; ++p) {
            Index_ current;
            if (!pred.prediction_stream.next(current)) {
                break;
            }
            pred.predictions_made.push_back(current);

            auto it = pred.next_cache_exists.find(current);
            if (it == pred.next_cache_exists.end()) {
                auto it2 = pred.cache_exists.find(curchunk);
                if (it2 != pred.cache_exists.end()) {
                    pred.next_cache_data[used].swap(pred.cache_data[it2->second]);
                    anticipated_size += pred.next_cache_data[used].byte_size;
                    if (antitipated_size > work.cache_size_limit) {
                        break;
                    }
                } else {
                    size_t next_element_size = (pointers[current + 1] - pointers[current]) * (sizeof(Value_) + sizeof(Index_));
                    anticipated_size += next_element_size;
                    if (antitipated_size > work.cache_size_limit) {
                        break;
                    }
                    pred.chunks_in_need.emplace_back(current, used);
                }
                pred.next_cache_exists[current] = used;

            } else {
                pred.predictions_made.emplace_back(current);
            }
        }

        // See corresponding logic in the HDF5DenseMatrix.
        size_t search = 0;
        for (const auto& c : pred.chunks_in_need) {
            if (pred.next_cache_data[c.second].init) { 
                while (pred.cache_data[search].init) { 
                    ++search;
                }

#ifdef DEBUG
                if (search >= pred.cache_data.size()) {
                    throw std::runtime_error("internal cache management error for HDF5DenseMatrix with oracles");
                }
#endif

                pred.next_cache_data[c.second].swap(pred.cache_data[search]);
                ++search;
            }
        }


#ifndef TATAMI_HDF5_PARALLEL_LOCK
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        for (const auto& c : pred.chunks_in_need) {
            auto& cache_target = pred.next_cache_data[c.second];
            auto actual_dim = extract_chunk<accrow_>(c.first, mydim, chunk_mydim, cache_target.data(), extract_value, extract_length, work);
            if constexpr(accrow_ == transpose_) {
                pred.cache_transpose_info.emplace_back(c.second, actual_dim);
            }
        }

#ifndef TATAMI_HDF5_PARALLEL_LOCK
        }
#else
        });
#endif

        // Applying transpositions to all cached buffers for easier retrieval, but only once the lock is released.
        if constexpr(accrow_ == transpose_) {
            for (const auto& x : pred.cache_transpose_info) {
                transpose(pred.next_cache_data[x.first], work.transposition_buffer, x.second, extract_length);
            }
        }

        pred.cache_data.swap(pred.next_cache_data);
        pred.cache_exists.swap(pred.next_cache_exists);

        pred.predictions_fulfilled = 1; // well, because we just used one.
        const auto& chosen = pred.predictions_made.front(); // assuming at least one prediction was made, otherwise, why was fetch() even called?
        return pred.cache_data[chosen.first].data() + extract_length * chosen.second;
    }

    /********************************************
     ************ Primary extraction ************
     ********************************************/
private:
    void populate_primary_cache(size_t i, PrimaryH5Core& core, bool needs_value) const {
        if (primary_cache_id[i] == core.current_cache_id && core.init) {
            return;
        }

        core.init = true;
        core.current_cache_id = primary_cache_id[i];

        // Pulling out the entire chunk of primary dimensions containing
        // 'i'. We have to do this for all indices, regardless of the
        // slicing/indexing. 
        const auto& limits = primary_cache_limits[core.current_cache_id];
        hsize_t offset = limits.first;
        hsize_t count = limits.second - limits.first;

#ifndef TATAMI_HDF5_PARALLEL_LOCK
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        core.dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
        core.memspace.setExtentSimple(1, &count);
        core.memspace.selectAll();

        core.index_cache.resize(count);
        core.index.read(core.index_cache.data(), HDF5::define_mem_type<Index_>(), core.memspace, core.dataspace);

        if (needs_value) {
            // In theory, we could avoid extracting data for the entire column when
            // populating the cache for sliced or indexed queries. In practice,
            // each chunk is often larger than the number of non-zeroes in a
            // column, so we end up having to pull out the entire column anyway.
            // Also, I want to get out of this critical region ASAP and do my
            // indexing elsewhere. 

            core.data_cache.resize(count);
            core.data.read(core.data_cache.data(), HDF5::define_mem_type<Value_>(), core.memspace, core.dataspace);
        }

#ifndef TATAMI_HDF5_PARALLEL_LOCK 
        }
#else
        });
#endif

    };

    template<class Function_>
    void extract_primary_raw(size_t i, Function_ fill, Index_ start, Index_ length, PrimaryH5Core& core, bool needs_value) const {
        if (length == 0) {
            return;
        }

        Index_ primary_length = pointers[i + 1] - pointers[i];
        if (primary_length ==0) {
            return;
        }

        populate_primary_cache(i, core, needs_value);

        const auto& limits = primary_cache_limits[core.current_cache_id];
        Index_ offset = pointers[i] - limits.first;
        auto istart = core.index_cache.begin() + offset;

        Index_ request_start = 0;
        if (start > *istart) { // 'istart' guaranteed to be valid from length check above.
            bool do_cache = !core.starts.empty();
            if (do_cache && core.starts[i] != -1) {
                request_start = core.starts[i];
            } else {
                request_start = std::lower_bound(istart, istart + primary_length, start) - istart;
                if (do_cache) {
                    core.starts[i] = request_start;
                }
            }
        }

        istart += request_start;
        Index_ end = start + length;
        if (needs_value) {
            auto dstart = core.data_cache.begin() + offset + request_start;
            for (size_t i = request_start; i < primary_length && *istart < end; ++i, ++istart, ++dstart) {
                fill(*istart, *dstart);
            }
        } else {
            for (size_t i = request_start; i < primary_length && *istart < end; ++i, ++istart) {
                fill(*istart, 0);
            }
        }

        return;
    }

    const Value_* extract_primary(size_t i, Value_* dbuffer, Index_ start, Index_ length, PrimaryH5Core& core) const {
        std::fill(dbuffer, dbuffer + length, 0);

        extract_primary_raw(i, 
            [&](Index_ pos, Value_ value) -> void {
                dbuffer[pos - start] = value;
            }, 
            start, 
            length, 
            core, 
            true
        );

        return dbuffer;
    }

    SparseRange<Value_, Index_> extract_primary(size_t i, Value_* dbuffer, Index_* ibuffer, Index_ start, Index_ length, PrimaryH5Core& core, bool needs_value, bool needs_index) const {
        Index_ counter = 0;

        extract_primary_raw(i, 
            [&](Index_ pos, Value_ value) -> void {
                if (needs_index) {
                    ibuffer[counter] = pos;
                }
                if (needs_value) {
                    dbuffer[counter] = value;
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
    void extract_primary_raw(size_t i, Function_ fill, Skip_ skip, const std::vector<Index_>& indices, PrimaryH5Core& core, bool needs_value) const {
        if (indices.empty()) {
            return;
        }

        Index_ primary_length = pointers[i + 1] - pointers[i];
        if (primary_length ==0) {
            return;
        }

        populate_primary_cache(i, core, needs_value);

        const auto& limits = primary_cache_limits[core.current_cache_id];
        Index_ offset = pointers[i] - limits.first;
        auto istart = core.index_cache.begin() + offset;
        auto iend = istart + primary_length;

        Index_ quick_shift = 0;
        if (indices[0] > *istart) { // Both are guaranteed to valid, from length check above.
            bool do_cache = !core.starts.empty();
            if (do_cache && core.starts[i] != -1) {
                quick_shift = core.starts[i];
            } else {
                quick_shift = std::lower_bound(istart, istart + primary_length, indices[0]) - istart;
                if (do_cache) {
                    core.starts[i] = quick_shift;
                }
            }
        }

        istart += quick_shift;
        auto dstart = core.data_cache.begin();
        if (needs_value) {
            dstart += offset + quick_shift;
        }

        for (auto idx : indices) {
            while (istart != iend && *istart < idx) {
                ++istart;
                ++dstart;
            }
            if (istart == iend) {
                break;
            }
            if (*istart == idx) {
                if (needs_value) {
                    fill(idx, *dstart);
                } else {
                    fill(idx, 0);
                }
            } else {
                skip();
            }
        }
    }

    const Value_* extract_primary(size_t i, Value_* buffer, const std::vector<Index_>& indices, PrimaryH5Core& core) const {
        std::fill(buffer, buffer + indices.size(), 0);
        auto original = buffer;

        extract_primary_raw(i, 
            [&](Index_, Value_ value) -> void {
                *buffer = value;
                ++buffer;
            }, 
            [&]() -> void{
                ++buffer;
            },
            indices, 
            core,
            true
        );

        return original;
    }

    SparseRange<Value_, Index_> extract_primary(size_t i, Value_* dbuffer, Index_* ibuffer, const std::vector<Index_>& indices, PrimaryH5Core& core, bool needs_value, bool needs_index) const {
        Index_ counter = 0;

        extract_primary_raw(i, 
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
    struct H5Core {
        H5::H5File file;
        H5::DataSet data, index;
        H5::DataSpace dataspace;
        H5::DataSpace memspace;

        std::vector<Index_> index_cache;
    };

    void fill_core(H5Core& core) const {
        // TODO: set more suitable chunk cache values here, to avoid re-reading
        // chunks on the boundaries of the primary cache.
        core.file.openFile(file_name, H5F_ACC_RDONLY);

        core.data = core.file.openDataSet(data_name);
        core.index = core.file.openDataSet(index_name);
        core.dataspace = core.data.getSpace();
    }

    template<class Function_>
    bool extract_secondary_raw(Index_ primary, Index_ secondary, Function_& fill, H5Core& core, bool needs_value) const {
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
    void extract_secondary_raw_loop(size_t i, Function_ fill, Index_ start, Index_ length, H5Core& core, bool needs_value) const {
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

    const Value_* extract_secondary(size_t i, Value_* buffer, Index_ start, Index_ length, H5Core& core) const {
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

    SparseRange<Value_, Index_> extract_secondary(size_t i, Value_* dbuffer, Index_* ibuffer, Index_ start, Index_ length, H5Core& core, bool needs_value, bool needs_index) const {
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
    void extract_secondary_raw_loop(size_t i, Function_ fill, Skip_ skip, const std::vector<Index_>& indices, H5Core& core, bool needs_value) const {
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

    const Value_* extract_secondary(size_t i, Value_* buffer, const std::vector<Index_>& indices, H5Core& core) const {
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

    SparseRange<Value_, Index_> extract_secondary(size_t i, Value_* dbuffer, Index_* ibuffer, const std::vector<Index_>& indices, H5Core& core, bool needs_value, bool needs_index) const {
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
                    parent->fill_core(core, accrow_ ? parent->nrows : parent->ncols);
                } else {
                    parent->fill_core(core, 0);
                }
            } else {
                parent->fill_core(core);
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
        ConditionalH5Core<accrow_> core;
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
        void set_oracle(std::unique_ptr<Oracle<Index_> >) {
            return; // TODO: add proper support for oracle handling.
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
