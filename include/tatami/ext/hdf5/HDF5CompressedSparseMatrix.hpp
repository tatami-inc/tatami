#ifndef TATAMI_HDF5_SPARSE_MATRIX_HPP
#define TATAMI_HDF5_SPARSE_MATRIX_HPP

#include "H5Cpp.h"

#include <string>
#include <cstdint>
#include <type_traits>
#include <cmath>

#include "../../base/Matrix.hpp"
#include "../../base/CompressedSparseMatrix.hpp"
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
 * @tparam ROW Whether the matrix is stored in compressed sparse row format.
 * @tparam T Type of the matrix values.
 * @tparam IDX Type of the row/column indices.
 */
template<bool ROW, typename T, typename IDX>
class HDF5CompressedSparseMatrix : public tatami::Matrix<T, IDX> {
    size_t nrows, ncols;
    std::string file_name;
    std::string data_name, index_name;
    std::vector<hsize_t> pointers;

    std::vector<size_t> primary_cache_id;
    std::vector<std::pair<size_t, size_t> > primary_cache_limits;

public:
    /**
     * @param nr Number of rows in the matrix.
     * @param nc Number of columns in the matrix.
     * @param file Path to the file.
     * @param vals Name of the 1D dataset inside `file` containing the non-zero elements.
     * @param idx Name of the 1D dataset inside `file` containing the indices of the non-zero elements.
     * If `ROW = true`, this should contain column indices sorted within each row, otherwise it should contain row indices sorted within each column.
     * @param ptr Name of the 1D dataset inside `file` containing the index pointers for the start and end of each row (if `ROW = true`) or column (otherwise).
     * This should have length equal to the number of rows (if `ROW = true`) or columns (otherwise).
     * @param cache_limit Limit to the size of the cache, in bytes.
     *
     * The cache is created by extracting multiple columns (for CSC matrices) or rows (CSR) on every call to the HDF5 library.
     * These are held in memory in the `Workspace` while the relevant column/row is returned to the user by `row()` or `column()`.
     * The aim is to minimize the number of calls to the HDF5 library - and thus expensive file reads - for consecutive accesses.
     */
    HDF5CompressedSparseMatrix(size_t nr, size_t nc, std::string file, std::string vals, std::string idx, std::string ptr, size_t cache_limit = 100000000) :
        nrows(nr),
        ncols(nc),
        file_name(file),
        data_name(std::move(vals)),
        index_name(std::move(idx)),
        pointers(ROW ? nr + 1 : nc + 1)
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
        const size_t nonzeros = HDF5::get_array_dimensions<1>(dhandle, "vals")[0];

        auto ihandle = HDF5::open_and_check_dataset<true>(file_handle, index_name);
        if (HDF5::get_array_dimensions<1>(ihandle, "idx")[0] != nonzeros) {
            throw std::runtime_error("number of non-zero elements is not consistent between 'data' and 'idx'");
        }

        auto phandle = HDF5::open_and_check_dataset<true>(file_handle, ptr);
        const size_t ptr_size = HDF5::get_array_dimensions<1>(phandle, "ptr")[0];
        if (ptr_size != pointers.size()) {
            throw std::runtime_error("'ptr' dataset should have length equal to the number of " + (ROW ? std::string("rows") : std::string("columns")) + " plus 1");
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
        size_t primary_dim = pointers.size() - 1;
        primary_cache_id.resize(primary_dim);
        if (primary_dim) {
            primary_cache_limits.emplace_back(0, pointers[1]);
        }

        size_t effective_cache_limit = cache_limit / (sizeof(T) + sizeof(IDX));
        size_t counter = 0, start = 0;
        for (size_t i = 1; i < primary_dim; ++i) {
            size_t end = pointers[i + 1];
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
    size_t nrow() const {
        return nrows;
    }

    size_t ncol() const {
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
        return ROW;
    }

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::column;

    using Matrix<T, IDX>::dense_row_workspace;

    using Matrix<T, IDX>::dense_column_workspace;

    using Matrix<T, IDX>::sparse_row_workspace;

    using Matrix<T, IDX>::sparse_column_workspace;

private:
    struct SparseBase {
        SparseBase(const WorkspaceOptions& opt) : needs_index(sparse_extract_index(opt.sparse_extract_mode)), needs_value(sparse_extract_value(opt.sparse_extract_mode)) {}
        bool needs_index;
        bool needs_value;
    };

    struct DenseBase {
        DenseBase(const WorkspaceOptions& opt) {}
    };

    template<bool SPARSE>
    using ConditionalBase = typename std::conditional<SPARSE, SparseBase, DenseBase>::type;

private:
    struct H5Core {
        H5::H5File file;
        H5::DataSet data, index;
        H5::DataSpace dataspace;
        H5::DataSpace memspace;

        std::vector<IDX> index_cache;
    };

    void fill_core(H5Core& core) const {
#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        // TODO: set more suitable chunk cache values here, to avoid re-reading
        // chunks on the boundaries of the primary cache.
        core.file.openFile(file_name, H5F_ACC_RDONLY);

        core.data = core.file.openDataSet(data_name);
        core.index = core.file.openDataSet(index_name);
        core.dataspace = core.data.getSpace();

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        }
#else
        });
#endif
    }

    struct PrimaryH5Core : public H5Core {
        std::vector<T> data_cache;
        size_t current_cache_id = 0;
        bool init = false;

        std::vector<size_t> starts;
    };

    void fill_core(PrimaryH5Core& core, size_t cache_size) const {
        fill_core(core);
        core.starts.resize(cache_size, -1);
    }

    template<bool WORKROW>
    using ConditionalH5Core = typename std::conditional<WORKROW == ROW, PrimaryH5Core, H5Core>::type;

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
        core.index.read(core.index_cache.data(), HDF5::define_mem_type<IDX>(), core.memspace, core.dataspace);

        if (needs_value) {
            // In theory, we could avoid extracting data for the entire column when
            // populating the cache for sliced or indexed queries. In practice,
            // each chunk is often larger than the number of non-zeroes in a
            // column, so we end up having to pull out the entire column anyway.
            // Also, I want to get out of this critical region ASAP and do my
            // indexing elsewhere. 

            core.data_cache.resize(count);
            core.data.read(core.data_cache.data(), HDF5::define_mem_type<T>(), core.memspace, core.dataspace);
        }

#ifndef TATAMI_HDF5_PARALLEL_LOCK 
        }
#else
        });
#endif

    };

    template<class Function>
    void extract_primary_raw(size_t i, Function fill, size_t start, size_t length, PrimaryH5Core& core, bool needs_value) const {
        if (length == 0) {
            return;
        }

        size_t primary_length = pointers[i + 1] - pointers[i];
        if (primary_length ==0) {
            return;
        }

        populate_primary_cache(i, core, needs_value);

        const auto& limits = primary_cache_limits[core.current_cache_id];
        size_t offset = pointers[i] - limits.first;
        auto istart = core.index_cache.begin() + offset;

        size_t request_start = 0;
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
        size_t end = start + length;
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

    const T* extract_primary(size_t i, T* dbuffer, size_t start, size_t length, PrimaryH5Core& core) const {
        std::fill(dbuffer, dbuffer + length, 0);

        extract_primary_raw(i, 
            [&](IDX pos, T value) -> void {
                dbuffer[pos - start] = value;
            }, 
            start, 
            length, 
            core, 
            true
        );

        return dbuffer;
    }

    SparseRange<T, IDX> extract_primary(size_t i, T* dbuffer, IDX* ibuffer, size_t start, size_t length, PrimaryH5Core& core, bool needs_index, bool needs_value) const {
        size_t counter = 0;

        extract_primary_raw(i, 
            [&](IDX pos, T value) -> void {
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

        return SparseRange<T, IDX>(counter, dbuffer, ibuffer);
    }

private:
    // This could be improved by extracting multiple rows at any given call and
    // caching them for subsequent requests. However, even then, we'd require
    // multiple re-reads from file when we exceed the cache. So, any caching
    // would be just turning an extremely bad access pattern into a very bad
    // pattern, when users shouldn't even be calling this at all... 

    template<class Function>
    bool extract_secondary_raw(IDX primary, size_t secondary, Function& fill, H5Core& core, bool needs_value) const {
        size_t left = pointers[primary], right = pointers[primary + 1];
        core.index_cache.resize(right - left);

        // Serial locks should be applied by the callers.
        hsize_t offset = left;
        hsize_t count = core.index_cache.size();
        core.dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
        core.memspace.setExtentSimple(1, &count);
        core.memspace.selectAll();
        core.index.read(core.index_cache.data(), HDF5::define_mem_type<IDX>(), core.memspace, core.dataspace);

        auto it = std::lower_bound(core.index_cache.begin(), core.index_cache.end(), secondary);
        if (it != core.index_cache.end() && *it == secondary) {
            if (needs_value) {
                offset = left + (it - core.index_cache.begin());
                count = 1;
                core.dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
                core.memspace.setExtentSimple(1, &count);
                core.memspace.selectAll();

                T dest;
                core.data.read(&dest, HDF5::define_mem_type<T>(), core.memspace, core.dataspace);
                fill(primary, dest);
            } else {
                fill(primary, 0);
            }
            return true;
        } else {
            return false;
        }
    }

    template<class Function>
    void extract_secondary_raw_loop(size_t i, Function fill, size_t start, size_t length, H5Core& core, bool needs_value) const {
#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        size_t end = start + length;
        for (size_t j = start; j < end; ++j) {
            extract_secondary_raw(j, i, fill, core, needs_value);
        }

#ifndef TATAMI_HDF5_PARALLEL_LOCK 
        }
#else
        });
#endif
    }

    const T* extract_secondary(size_t i, T* dbuffer, size_t start, size_t length, H5Core& core) const {
        std::fill(dbuffer, dbuffer + length, 0);

        extract_secondary_raw_loop(i, 
            [&](IDX pos, T value) -> void {
                dbuffer[pos - start] = value;
            }, 
            start, 
            length, 
            core,
            true
        );

        return dbuffer;
    }

    SparseRange<T, IDX> extract_secondary(size_t i, T* dbuffer, IDX* ibuffer, size_t start, size_t length, H5Core& core, bool needs_index, bool needs_value) const {
        size_t counter = 0;

        extract_secondary_raw_loop(i, 
            [&](IDX pos, T value) -> void {
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

        if (needs_value) {
            dbuffer = NULL;
        }
        if (needs_index) {
            ibuffer = NULL;
        }

        return SparseRange<T, IDX>(counter, dbuffer, ibuffer);
    }

    /*****************************************
     ************ Full extraction ************
     *****************************************/
private:
    template<bool WORKROW, bool SPARSE>
    struct Hdf5SparseWorkspace : public Workspace<WORKROW, SPARSE>, public ConditionalBase<SPARSE> {
        Hdf5SparseWorkspace(const WorkspaceOptions& opt) : ConditionalBase<SPARSE>(opt) {}
        ConditionalH5Core<WORKROW> core;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<Workspace<WORKROW, SPARSE> > create_new_workspace(const WorkspaceOptions& opt) const {
        auto ptr = new Hdf5SparseWorkspace<WORKROW, SPARSE>(opt);
        std::shared_ptr<Workspace<WORKROW, SPARSE> > output(ptr);
        if constexpr(WORKROW == ROW) {
            fill_core(ptr->core, 0);
        } else {
            fill_core(ptr->core);
        }
        return output;
    }

    template<bool WORKROW>
    const T* get_dense(size_t i, T* buffer, DenseWorkspace<WORKROW>* work) const {
        size_t n = (WORKROW ? ncols : nrows);
        auto wptr = static_cast<Hdf5SparseWorkspace<WORKROW, false>*>(work);
        if constexpr(ROW == WORKROW) {
            return extract_primary(i, buffer, 0, n, wptr->core);
        } else {
            return extract_secondary(i, buffer, 0, n, wptr->core);
        }
    }

    template<bool WORKROW>
    SparseRange<T, IDX> get_sparse(size_t i, T* vbuffer, IDX* ibuffer, SparseWorkspace<WORKROW>* work) const {
        size_t n = (WORKROW ? ncols : nrows);
        auto wptr = static_cast<Hdf5SparseWorkspace<WORKROW, true>*>(work);
        if constexpr(ROW == WORKROW) {
            if (wptr->needs_index || wptr->needs_value) {
                return extract_primary(i, vbuffer, ibuffer, 0, n, wptr->core, wptr->needs_index, wptr->needs_value);
            } else {
                // Quick return is possible if we don't need any indices or values.                        
                return SparseRange<T, IDX>(pointers[i+1] - pointers[i], NULL, NULL);
            }
        } else {
            return extract_secondary(i, vbuffer, ibuffer, 0, n, wptr->core, wptr->needs_index, wptr->needs_value);
        }
    }

public:
    std::shared_ptr<DenseRowWorkspace> dense_row_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<true, false>(opt);
    }

    std::shared_ptr<DenseColumnWorkspace> dense_column_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<false, false>(opt);
    }

    const T* row(size_t r, T* buffer, DenseRowWorkspace* work) const {
        return get_dense(r, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnWorkspace* work) const {
        return get_dense(c, buffer, work);
    }

    std::shared_ptr<SparseRowWorkspace> sparse_row_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<true, true>(opt);
    }

    std::shared_ptr<SparseColumnWorkspace> sparse_column_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<false, true>(opt);
    }

    SparseRange<T, IDX> row(size_t r, T* dbuffer, IDX* ibuffer, SparseRowWorkspace* work) const {
        return get_sparse(r, dbuffer, ibuffer, work);
    }

    SparseRange<T, IDX> column(size_t c, T* dbuffer, IDX* ibuffer, SparseColumnWorkspace* work) const {
        return get_sparse(c, dbuffer, ibuffer, work);
    }

    /*****************************************
     *********** Block extraction ************
     *****************************************/
private:
    template<bool WORKROW, bool SPARSE>
    struct Hdf5SparseBlockWorkspace : public BlockWorkspace<WORKROW, SPARSE>, public ConditionalBase<SPARSE> {
        Hdf5SparseBlockWorkspace(size_t s, size_t l, const WorkspaceOptions& opt) : BlockWorkspace<WORKROW, SPARSE>(s, l), ConditionalBase<SPARSE>(opt) {}
        ConditionalH5Core<WORKROW> core;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > create_new_workspace(size_t s, size_t l, const WorkspaceOptions& opt) const {
        auto ptr = new Hdf5SparseBlockWorkspace<WORKROW, SPARSE>(s, l, opt);
        std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > output(ptr);

        if constexpr(WORKROW == ROW) {
            if constexpr(WORKROW) {
                fill_core(ptr->core, opt.cache_for_reuse ? nrows : 0);
            } else {
                fill_core(ptr->core, opt.cache_for_reuse ? ncols : 0);
            }
        } else {
            fill_core(ptr->core);
        }

        return output;
    }

    template<bool WORKROW>
    const T* get_dense(size_t i, T* buffer, DenseBlockWorkspace<WORKROW>* work) const {
        auto wptr = static_cast<Hdf5SparseBlockWorkspace<WORKROW, false>*>(work);
        if constexpr(ROW == WORKROW) {
            return extract_primary(i, buffer, wptr->start, wptr->length, wptr->core);
        } else {
            return extract_secondary(i, buffer, wptr->start, wptr->length, wptr->core);
        }
    }

    template<bool WORKROW>
    SparseRange<T, IDX> get_sparse(size_t i, T* vbuffer, IDX* ibuffer, SparseBlockWorkspace<WORKROW>* work) const {
        auto wptr = static_cast<Hdf5SparseBlockWorkspace<WORKROW, true>*>(work);
        if constexpr(ROW == WORKROW) {
            return extract_primary(i, vbuffer, ibuffer, wptr->start, wptr->length, wptr->core, wptr->needs_index, wptr->needs_value);
        } else {
            return extract_secondary(i, vbuffer, ibuffer, wptr->start, wptr->length, wptr->core, wptr->needs_index, wptr->needs_value);
        }
    }

public:
    std::shared_ptr<DenseRowBlockWorkspace> dense_row_workspace(size_t s, size_t l, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, false>(s, l, opt);
    }

    std::shared_ptr<DenseColumnBlockWorkspace> dense_column_workspace(size_t s, size_t l, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, false>(s, l, opt);
    }

    const T* row(size_t r, T* buffer, DenseRowBlockWorkspace* work) const {
        return get_dense(r, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnBlockWorkspace* work) const {
        return get_dense(c, buffer, work);
    }

    std::shared_ptr<SparseRowBlockWorkspace> sparse_row_workspace(size_t s, size_t l, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, true>(s, l, opt);
    }

    std::shared_ptr<SparseColumnBlockWorkspace> sparse_column_workspace(size_t s, size_t l, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, true>(s, l, opt);
    }

    SparseRange<T, IDX> row(size_t r, T* dbuffer, IDX* ibuffer, SparseRowBlockWorkspace* work) const {
        return get_sparse(r, dbuffer, ibuffer, work);
    }

    SparseRange<T, IDX> column(size_t c, T* dbuffer, IDX* ibuffer, SparseColumnBlockWorkspace* work) const {
        return get_sparse(c, dbuffer, ibuffer, work);
    }

    /*****************************************
     *********** Index extraction ************
     *****************************************/
private:
    template<class Function, class Skip>
    void extract_primary_raw(size_t i, Function fill, Skip skip, const std::vector<IDX>& indices, PrimaryH5Core& core, bool needs_value) const {
        if (indices.empty()) {
            return;
        }

        size_t primary_length = pointers[i + 1] - pointers[i];
        if (primary_length ==0) {
            return;
        }

        populate_primary_cache(i, core, needs_value);

        const auto& limits = primary_cache_limits[core.current_cache_id];
        size_t offset = pointers[i] - limits.first;
        auto istart = core.index_cache.begin() + offset;
        auto iend = istart + primary_length;

        size_t quick_shift = 0;
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

    const T* extract_primary(size_t i, T* dbuffer, const std::vector<IDX>& indices, PrimaryH5Core& core) const {
        std::fill(dbuffer, dbuffer + indices.size(), 0);
        auto original = dbuffer;

        extract_primary_raw(i, 
            [&](IDX, T value) -> void {
                *dbuffer = value;
                ++dbuffer;
            }, 
            [&]() -> void{
                ++dbuffer;
            },
            indices, 
            core,
            true
        );

        return original;
    }

    SparseRange<T, IDX> extract_primary(size_t i, T* dbuffer, IDX* ibuffer, const std::vector<IDX>& indices, PrimaryH5Core& core, bool needs_index, bool needs_value) const {
        size_t counter = 0;

        extract_primary_raw(i, 
            [&](IDX pos, T value) -> void {
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

        if (needs_index) {
            ibuffer = NULL;
        }
        if (needs_value) {
            dbuffer = NULL;
        }

        return SparseRange<T, IDX>(counter, dbuffer, ibuffer);
    }

    template<class Function, class Skip>
    size_t extract_secondary_raw_loop(size_t i, Function fill, Skip skip, const std::vector<IDX>& indices, H5Core& core, bool needs_value) const {
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

    const T* extract_secondary(size_t i, T* dbuffer, const std::vector<IDX>& indices, H5Core& core) const {
        std::fill(dbuffer, dbuffer + indices.size(), 0);
        auto original = dbuffer;
        extract_secondary_raw_loop(i, 
            [&](IDX pos, T value) -> void {
                *dbuffer = value;
                ++dbuffer;
            }, 
            [&]() -> void {
                ++dbuffer;
            },
            indices, 
            core,
            true
        );
        return original;
    }

    SparseRange<T, IDX> extract_secondary(size_t i, T* dbuffer, IDX* ibuffer, const std::vector<IDX>& indices, H5Core& core, bool needs_index, bool needs_value) const {
        size_t counter = 0;

        extract_secondary_raw_loop(i, 
            [&](IDX pos, T value) -> void {
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

        if (needs_value) {
            dbuffer = NULL;
        }
        if (needs_index) {
            ibuffer = NULL;
        }

        return SparseRange<T, IDX>(counter, dbuffer, ibuffer);
    }

private:
    template<bool WORKROW, bool SPARSE>
    struct Hdf5SparseIndexWorkspace : public IndexWorkspace<IDX, WORKROW, SPARSE>, public ConditionalBase<SPARSE> {
        Hdf5SparseIndexWorkspace(std::vector<IDX> i, const WorkspaceOptions& opt) : IndexWorkspace<IDX, WORKROW, SPARSE>(i.size()), ConditionalBase<SPARSE>(opt), indices_(std::move(i)) {}

        ConditionalH5Core<WORKROW> core;

        std::vector<IDX> indices_;
        const std::vector<IDX>& indices() const { return indices_; }
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > create_new_workspace(std::vector<IDX> i, const WorkspaceOptions& opt) const {
        auto ptr = new Hdf5SparseIndexWorkspace<WORKROW, SPARSE>(std::move(i), opt);
        std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > output(ptr);

        if constexpr(WORKROW == ROW) {
            if constexpr(WORKROW) {
                fill_core(ptr->core, opt.cache_for_reuse ? nrows : 0);
            } else {
                fill_core(ptr->core, opt.cache_for_reuse ? ncols : 0);
            }
        } else {
            fill_core(ptr->core);
        }

        return output;
    }

    template<bool WORKROW>
    const T* get_dense(size_t i, T* buffer, DenseIndexWorkspace<IDX, WORKROW>* work) const {
        auto wptr = static_cast<Hdf5SparseIndexWorkspace<WORKROW, false>*>(work);
        if constexpr(ROW == WORKROW) {
            return extract_primary(i, buffer, wptr->indices_, wptr->core);
        } else {
            return extract_secondary(i, buffer, wptr->indices_, wptr->core);
        }
    }

    template<bool WORKROW>
    SparseRange<T, IDX> get_sparse(size_t i, T* vbuffer, IDX* ibuffer, SparseIndexWorkspace<IDX, WORKROW>* work) const {
        auto wptr = static_cast<Hdf5SparseIndexWorkspace<WORKROW, true>*>(work);
        if constexpr(ROW == WORKROW) {
            return extract_primary(i, vbuffer, ibuffer, wptr->indices_, wptr->core, wptr->needs_index, wptr->needs_value);
        } else {
            return extract_secondary(i, vbuffer, ibuffer, wptr->indices_, wptr->core, wptr->needs_index, wptr->needs_value);
        }
    }

public:
    std::shared_ptr<DenseRowIndexWorkspace<IDX> > dense_row_workspace(std::vector<IDX> i, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, false>(std::move(i), opt);
    }

    std::shared_ptr<DenseColumnIndexWorkspace<IDX> > dense_column_workspace(std::vector<IDX> i, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, false>(std::move(i), opt);
    }

    const T* row(size_t r, T* buffer, DenseRowIndexWorkspace<IDX>* work) const {
        return get_dense<true>(r, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnIndexWorkspace<IDX>* work) const {
        return get_dense<false>(c, buffer, work);
    }

    std::shared_ptr<SparseRowIndexWorkspace<IDX> > sparse_row_workspace(std::vector<IDX> i, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, true>(std::move(i), opt);
    }

    std::shared_ptr<SparseColumnIndexWorkspace<IDX> > sparse_column_workspace(std::vector<IDX> i, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, true>(std::move(i), opt);
    }

    SparseRange<T, IDX> row(size_t r, T* dbuffer, IDX* ibuffer, SparseRowIndexWorkspace<IDX>* work) const {
        return get_sparse<true>(r, dbuffer, ibuffer, work);
    }

    SparseRange<T, IDX> column(size_t c, T* dbuffer, IDX* ibuffer, SparseColumnIndexWorkspace<IDX>* work) const {
        return get_sparse<false>(c, dbuffer, ibuffer, work);
    }
};

}

#endif
