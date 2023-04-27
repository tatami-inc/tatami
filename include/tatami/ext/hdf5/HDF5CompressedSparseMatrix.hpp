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

private:
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

    struct PrimaryH5Core : public H5Core {
        std::vector<Value_> data_cache;
        Index_ current_cache_id = 0;
        bool init = false;

        std::vector<size_t> starts;
    };

    void fill_core(PrimaryH5Core& core, Index_ cache_size) const {
        fill_core(core);
        core.starts.resize(cache_size, -1);
    }

    template<bool accrow_>
    using ConditionalH5Core = typename std::conditional<accrow_ == row_, PrimaryH5Core, H5Core>::type;

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
        void set_oracle(std::unique_ptr<SequenceOracle<Index_> >) {
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
