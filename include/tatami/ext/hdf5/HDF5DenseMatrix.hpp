#ifndef TATAMI_HDF5_DENSE_MATRIX_HPP
#define TATAMI_HDF5_DENSE_MATRIX_HPP

#include "H5Cpp.h"

#include <string>
#include <cstdint>
#include <type_traits>
#include <cmath>

#include "../../base/Matrix.hpp"
#include "utils.hpp"

/**
 * @file HDF5DenseMatrix.hpp
 *
 * @brief Defines a class for a HDF5-backed dense matrix.
 */

namespace tatami {

/**
 * @brief Dense matrix backed by a DataSet in a HDF5 file.
 *
 * This class retrieves data from the HDF5 file on demand rather than loading it all in at the start.
 * This allows us to handle very large datasets in limited memory at the cost of some speed.
 *
 * We manually handle the chunk caching to speed up access for consecutive rows and columns.
 * The policy is to minimize the number of calls to the HDF5 library by requesting large contiguous slices where possible.
 * The size of the slice is determined by the cache limit in the constructor.
 *
 * Callers should follow the `prefer_rows()` suggestion when extracting data,
 * as this tries to minimize the number of chunks that need to be read per access request.
 * If they do not, the access pattern on disk may be slightly to highly suboptimal, depending on the chunk dimensions.
 *
 * As the HDF5 library is not generally thread-safe, the HDF5-related operations should only be run in a single thread.
 * For OpenMP, this is handled automatically by putting all HDF5 operations in a critical region.
 * For other parallelization schemes, callers should define the `TATAMI_HDF5_PARALLEL_LOCK` macro;
 * this should be a function that accepts and executes a no-argument lambda within an appropriate serial region (e.g., based on a global mutex).
 *
 * @tparam T Type of the matrix values.
 * @tparam IDX Type of the row/column indices.
 * @tparam transpose Whether the dataset is transposed in its storage order, i.e., rows in HDF5 are columns in this matrix.
 */
template<typename T, typename IDX, bool transpose = false>
class HDF5DenseMatrix : public tatami::Matrix<T, IDX> {
    size_t firstdim, seconddim;
    std::string file_name, dataset_name;

    hsize_t cache_firstdim, cache_seconddim;
    bool prefer_firstdim;

public:
    /**
     * @param file Path to the file.
     * @param name Path to the dataset inside the file.
     * @param cache_limit Limit to the size of the chunk cache, in bytes.
     *
     * The cache size should be large enough to fit all chunks spanned by a row or column, for (near-)consecutive row and column access respectively.
     * Otherwise, performance will degrade as the same chunks may need to be repeatedly read back into memory.
     */
    HDF5DenseMatrix(std::string file, std::string name, size_t cache_limit = 100000000) : 
        file_name(std::move(file)), 
        dataset_name(std::move(name))
    {
#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        H5::H5File fhandle(file_name, H5F_ACC_RDONLY);
        auto dhandle = HDF5::open_and_check_dataset<false>(fhandle, dataset_name);
        auto dims = HDF5::get_array_dimensions<2>(dhandle, dataset_name);
        firstdim = dims[0];
        seconddim = dims[1];

        hsize_t chunk_firstdim, chunk_seconddim;
        auto dparms = dhandle.getCreatePlist();
        if (dparms.getLayout() != H5D_CHUNKED) {
            // If contiguous, each firstdim is treated as a chunk.
            chunk_firstdim = 1;
            chunk_seconddim = seconddim;
        } else {
            hsize_t chunk_dims[2];
            dparms.getChunk(2, chunk_dims);
            chunk_firstdim = chunk_dims[0];
            chunk_seconddim = chunk_dims[1];
        }

        auto limit_cache = [](size_t chunk_dim, size_t mydim, size_t otherdim, size_t limit) -> size_t {
            double size_along_otherdim = sizeof(T) * otherdim;
            size_t nelements_along_chunkdim = std::min(mydim, static_cast<size_t>(limit / size_along_otherdim));
            size_t nchunks = nelements_along_chunkdim / chunk_dim;
            return std::min(nelements_along_chunkdim, nchunks * chunk_dim); 
        };

        cache_firstdim = limit_cache(chunk_firstdim, firstdim, seconddim, cache_limit);
        cache_seconddim = limit_cache(chunk_seconddim, seconddim, firstdim, cache_limit);

        // Favoring extraction along the first dimension if it's not crippled
        // by the chunking. This favors pulling of contiguous sections along
        // the last (i.e., second) dimension, which is how HDF5 stores and
        // returns 2D arrays, so we avoid needing some manual transposition.
        if (cache_firstdim >= chunk_firstdim) {
            prefer_firstdim = true;
        } else if (cache_seconddim >= chunk_seconddim) {
            // Of course, if the first dimension is chunk-compromised, then we
            // would rather do the transpositions than repeat a read from file.
            prefer_firstdim = false;
        } else {
            // If both are crippled, then we might as well extract along the first dimension.
            prefer_firstdim = true;
        }

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        }
#else
        });
#endif

        return;
    }

public:
    size_t nrow() const {
        if constexpr(transpose) {
            return seconddim;
        } else {
            return firstdim;
        }
    }

    size_t ncol() const {
        if constexpr(transpose) {
            return firstdim;
        } else {
            return seconddim;
        }
    }

    /**
     * @return Boolean indicating whether to prefer row extraction.
     *
     * We favor extraction on the first dimension (rows by default, columns when `transpose = true`) as this matches the HDF5 storage order.
     * However, for some chunking scheme and `cache_limit`, this might require repeated reads from file;
     * in such cases, we switch to extraction on the second dimension.
     */
    bool prefer_rows() const {
        if constexpr(transpose) {
            return !prefer_firstdim;
        } else {
            return prefer_firstdim;
        }
    }

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::column;

private:
    struct Hdf5WorkspaceBase {
        H5::H5File file;
        H5::DataSet dataset;
        H5::DataSpace dataspace;
        H5::DataSpace memspace;

        std::vector<T> cache;
        size_t cached_chunk = 0;
        bool init = true;

        std::vector<T> buffer; // buffer is provided for transpositions for easier extraction.
    };

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct Hdf5Workspace : public Workspace<ROW> {
        Hdf5WorkspaceBase base;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowWorkspace> new_row_workspace(bool = false) const {
        auto ptr = new Hdf5Workspace<true>;
        std::shared_ptr<RowWorkspace> output(ptr);
        fill_base(ptr->base);
        return output;
    }

    std::shared_ptr<ColumnWorkspace> new_column_workspace(bool = false) const {
        auto ptr = new Hdf5Workspace<false>;
        std::shared_ptr<ColumnWorkspace> output(ptr);
        fill_base(ptr->base);
        return output;
    }

    const T* row(size_t r, T* buffer, RowWorkspace* work) const {
        return extract<true>(r, buffer, 0, ncol(), static_cast<Hdf5Workspace<true>*>(work)->base);
    }

    const T* column(size_t c, T* buffer, ColumnWorkspace* work) const {
        return extract<false>(c, buffer, 0, nrow(), static_cast<Hdf5Workspace<false>*>(work)->base);
    }

private:
    void fill_base(Hdf5WorkspaceBase& base) const  {
#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        // Turn off HDF5's caching, as we'll be handling that.
        H5::FileAccPropList fapl(H5::FileAccPropList::DEFAULT.getId());
        fapl.setCache(0, 0, 0, 0);

        base.file.openFile(file_name, H5F_ACC_RDONLY, fapl);
        base.dataset = base.file.openDataSet(dataset_name);
        base.dataspace = base.dataset.getSpace();

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        }
#else
        });
#endif
    }

    template<bool ROW, typename ExtractType>
    const T* extract(size_t primary_start, size_t primary_length, T* target, const ExtractType& extract_value, size_t extract_length, Hdf5WorkspaceBase& work) const {
        hsize_t offset[2];
        hsize_t count[2];

        constexpr int dimdex = (ROW != transpose);
        offset[1-dimdex] = primary_start;
        count[1-dimdex] = primary_length;

        constexpr bool indexed = std::is_same<ExtractType, std::vector<IDX> >::value;

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        if constexpr(indexed) {
            // Take slices across the current chunk for each index. This should be okay if consecutive,
            // but hopefully they've fixed the problem with non-consecutive slices in:
            // https://forum.hdfgroup.org/t/union-of-non-consecutive-hyperslabs-is-very-slow/5062
            count[dimdex] = 1;
            work.dataspace.selectNone();
            for (auto idx : extract_value) {
                offset[dimdex] = idx;
                work.dataspace.selectHyperslab(H5S_SELECT_OR, count, offset);
            }
            count[dimdex] = extract_length; // for the memspace setter.
        } else {
            offset[dimdex] = extract_value;
            count[dimdex] = extract_length;
            work.dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);
        }

        // HDF5 is a lot faster when the memspace and dataspace match in dimensionality.
        // Presumably there is some shuffling that happens inside when dimensions don't match.
        work.memspace.setExtentSimple(2, count);
        work.memspace.selectAll();

        work.dataset.read(target, HDF5::define_mem_type<T>(), work.memspace, work.dataspace);

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        }
#else
        });
#endif
    }

    template<bool ROW, typename ExtractType>
    const T* extract(size_t i, T* buffer, const ExtractType& extract_value, size_t extract_length, Hdf5WorkspaceBase& work) const {
        // Figuring out which chunk the request belongs to.
        hsize_t cache_mydim, mydim, otherdim;
        if constexpr(ROW != transpose) {
            cache_mydim = cache_firstdim;
            mydim = firstdim;
            otherdim = seconddim;
        } else {
            cache_mydim = cache_seconddim;
            mydim = seconddim;
            otherdim = firstdim;
        }

        // No caching can be done here, so we just extract directly.
        if (cache_mydim == 0) {
            extract<ROW>(i, 1, buffer, extract_value, extract_length, work);
            return buffer;
        }

        size_t chunk = i / cache_mydim;
        if (chunk != work.cached_chunk || work.init) {
            work.init = false;

            hsize_t cache_mydim_start = chunk * cache_mydim;
            hsize_t cache_mydim_end = std::min(mydim, cache_mydim_start + cache_mydim);
            hsize_t cache_mydim_actual = cache_mydim_end - cache_mydim_start;
            size_t new_cache_size = extract_length * cache_mydim_actual;

            T* destination;
            work.cache.resize(new_cache_size);
            if constexpr(ROW != transpose) {
                destination = work.cache.data();
            } else {
                work.buffer.resize(new_cache_size);
                destination = work.buffer.data();
            }

            extract<ROW>(cache_mydim_start, cache_mydim_actual, destination, extract_value, extract_length, work);

            if constexpr(ROW == transpose) {
                auto output = work.cache.begin();
                for (hsize_t x = 0; x < cache_mydim_actual; ++x, output += extract_length) {
                    auto in = work.buffer.begin() + x;
                    for (hsize_t y = 0; y < extract_length; ++y, in += cache_mydim_actual) {
                        *(output + y) = *in;
                    }
                }
            }

            work.cached_chunk = chunk;
        }

        size_t index = i % cache_mydim;
        auto wIt = work.cache.begin() + index * extract_length;
        std::copy(wIt, wIt + extract_length, buffer);
        return buffer;
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct Hdf5BlockWorkspace : public BlockWorkspace<ROW> {
        Hdf5BlockWorkspace(size_t s, size_t l) : details(s, l) {}

        std::pair<size_t, size_t> details;
        const std::pair<size_t, size_t>& block() const { return details; }

        Hdf5WorkspaceBase base;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowBlockWorkspace> new_row_workspace(size_t s, size_t l, bool = false) const {
        auto ptr = new Hdf5BlockWorkspace<true>(s, l);
        std::shared_ptr<RowBlockWorkspace> output(ptr);
        fill_base(ptr->base);
        return output;
    }

    std::shared_ptr<ColumnBlockWorkspace> new_column_workspace(size_t s, size_t l, bool = false) const {
        auto ptr = new Hdf5BlockWorkspace<false>(s, l);
        std::shared_ptr<ColumnBlockWorkspace> output(ptr);
        fill_base(ptr->base);
        return output;
    }

    const T* row(size_t r, T* buffer, RowBlockWorkspace* work) const {
        auto ptr = static_cast<Hdf5BlockWorkspace<true>*>(work);
        return extract<true>(r, buffer, ptr->details.first, ptr->details.second, ptr->base);
    }

    const T* column(size_t c, T* buffer, ColumnBlockWorkspace* work) const {
        auto ptr = static_cast<Hdf5BlockWorkspace<false>*>(work);
        return extract<false>(c, buffer, ptr->details.first, ptr->details.second, ptr->base);
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct Hdf5IndexWorkspace : public IndexWorkspace<IDX, ROW> {
        Hdf5IndexWorkspace(std::vector<IDX> i) : indices_(std::move(i)) {}

        std::vector<IDX> indices_;
        const std::vector<IDX>& indices() const { return indices_; }

        Hdf5WorkspaceBase base;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(std::vector<IDX> i, bool = false) const { 
        auto ptr = new Hdf5IndexWorkspace<true>(std::move(i));
        std::shared_ptr<RowIndexWorkspace<IDX> > output(ptr);
        fill_base(ptr->base);
        return output;
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX> > new_column_workspace(std::vector<IDX> i, bool = false) const { 
        auto ptr = new Hdf5IndexWorkspace<false>(std::move(i));
        std::shared_ptr<ColumnIndexWorkspace<IDX> > output(ptr);
        fill_base(ptr->base);
        return output;
    }

    const T* row(size_t r, T* buffer, RowIndexWorkspace<IDX>* work) const {
        auto ptr = static_cast<Hdf5IndexWorkspace<true>*>(work);
        return extract<true>(r, buffer, ptr->indices_, ptr->indices_.size(), ptr->base);
    }

    const T* column(size_t c, T* buffer, ColumnIndexWorkspace<IDX>* work) const {
        auto ptr = static_cast<Hdf5IndexWorkspace<false>*>(work);
        return extract<false>(c, buffer, ptr->indices_, ptr->indices_.size(), ptr->base);
    }
};

}

#endif
