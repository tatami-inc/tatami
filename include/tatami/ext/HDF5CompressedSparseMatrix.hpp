#ifndef TATAMI_HDF5_SPARSE_MATRIX_HPP
#define TATAMI_HDF5_SPARSE_MATRIX_HPP

#include "H5Cpp.h"

#include <string>
#include <cstdint>
#include <type_traits>
#include <cmath>

#include "../base/Matrix.hpp"
#include "HDF5DenseMatrix.hpp"

/**
 * @file HDF5CompressedSparseMatrix.hpp
 *
 * @brief Defines a class for a HDF5-backed compressed sparse matrix.
 */

namespace tatami {

/**
 * @brief Compressed sparse matrix in a HDF5 file.
 *
 * This class retrieves data from the HDF5 file on demand rather than loading it all in at the start.
 * This allows us to handle very large datasets in limited memory at the cost of speed.
 * It is strongly advised to follow the `prefer_rows()` suggestion when extracting data,
 * otherwise the access pattern on disk will be highly suboptimal.
 *
 * @tparam T Type of the matrix values.
 * @tparam IDX Type of the row/column indices.
 * @tparam transpose Whether the dataset is transposed in its storage order, i.e., rows in HDF5 are columns in this matrix.
 */
template<typename T, typename IDX, bool row>
class HDF5CompressedSparseMatrix : public tatami::Matrix<T, IDX> {
    size_t nrows, ncols;
    std::string file_name, data_name, index_name;
    std::vector<hsize_t> pointers;

public:
    /**
     * @param file Path to the file.
     * @param name Path to the dataset inside the file.
     * @param cache_limit Limit to the size of the cache, in bytes.
     *
     * The actual cache size is automatically chosen to optimize for row or column extraction,
     * as long as this choice is less than `cache_limit`.
     * Otherwise, access is likely to be pretty slow.
     */
    HDF5DenseMatrix(size_t nr, size_t nc, std::string file, std::string vals, std::string idx, std::string ptr) :
        nrows(nr),
        ncols(nc),
        file_name(std::move(file)), 
        data_name(std::move(vals)),
        index_name(std::move(idx)),
        pointers(row ? nr + 1 : nc + 1)
    {
        H5::H5File fhandle(file_name, H5F_ACC_RDONLY);

        auto check_size = [&](H5::DataSpace& dhandle) -> size_t {
            auto dtype = dhandle.getTypeClass();
            if (dtype != H5T_INTEGER && dtype != H5T_FLOAT) { 
                throw std::runtime_error(std::string("expected numeric values in the '") + name + "' dataset");
            }
            auto space = dhandle.getSpace();
            
            int ndim = space.getSimpleExtentNdims();
            if (ndim == 1) {
                hsize_t dims_out[1];
                space.getSimpleExtentDims(dims_out, NULL);
                return dims_out[0];

            } else if (ndim == 2) {
                hsize_t dims_out[2];
                space.getSimpleExtentDims(dims_out, NULL);
                if (dims_out[0] == 1) {
                    return dims_out[1];
                } else if (dims_out[1] == 1) {
                    return dims_out[0];
                }
            }

            throw std::runtime_error(std::string("'") + name + "' should be a one-dimensional array");
            return 0;
        };

        auto dhandle = fhandle.openDataSet(data_name);
        const size_t nonzeros = check_size(dhandle, "vals")'

        auto ihandle = fhandle.openDataSet(index_name, "idx");
        if (check_size(ihandle) != nonzeros) {
            throw std::runtime_error("number of non-zero elements is not consistent between 'data' and 'idx'");
        }

        auto phandle = fhandle.openDataSet(ptr);
        const size_t ptr_size = check_size(phandle, "ptr");
        if (ptr_size != pointers.size()) {
            throw std::runtime_error("'ptr' is not of the appropriate length");
        }

        // Checking the contents of the index pointers.
        phandle.read(pointers.data(), H5::PredType::NATIVE_HSIZE);
        if (pointers[0] != 0) {
            throw std::runtime_error("first index pointer should be zero");
        }
        if (pointers.back() != nonzeros) {
            throw std::runtime_error("last index pointer should be equal to the number of non-zero elements");
        }
        for (size_t i = 1; i < pointers.size(); ++i) {
            if (pointers[i] < pointers[i-1]) {
                throw std::runtime_error("index pointers should be sorted");
            }
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
     * @return `true` if this is in compressed sparse row format.
     */
    bool prefer_rows() const {
        return row;
    }

public:
    /**
     * @cond
     */
    struct HDF5DenseWorkspace : public Workspace {
        HDF5DenseWorkspace(const std::string& fpath, const std::string& dpath, size_t nslots, size_t cache_size, size_t dim) : 
            current_size(dim),
            memspace(1, &current_size)
        { 
            // Configuring the chunk cache.
            if (nslots) {
                H5::FileAccPropList plist(H5::FileAccPropList::DEFAULT.getId());

                /* The first argument is ignored, according to https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html.
                 * Setting w0 to 0 to evict the last used chunk; no need to worry about full vs partial reads here.
                 */
                plist.setCache(10000, nslots, cache_size, 0);

                file.openFile(fpath, H5F_ACC_RDONLY, plist);
            } else {
                file.openFile(fpath, H5F_ACC_RDONLY);
            }

            dataset = file.openDataSet(dpath);
            dataspace = dataset.getSpace();
            dataspace.selectNone();
            memspace.selectAll();
            return;
        }

        H5::H5File file;
        H5::DataSet dataset;

        H5::DataSpace dataspace;
        hsize_t offset[2];
        hsize_t count[2];

        hsize_t current_size;
        H5::DataSpace memspace;
    };
    /**
     * @endcond
     */

    /**
     * @param row Should a workspace be created for row-wise extraction?
     * @return A shared pointer to a `Workspace` object is returned.
     */
    std::shared_ptr<Workspace> new_workspace(bool row) const {
        return std::shared_ptr<Workspace>(
            new HDF5DenseWorkspace(
                file_name, 
                dataset_name, 
                nslots, 
                row != transpose ? cache_size_firstdim : cache_size_seconddim,
                row ? ncol() : nrow()
            )
        );
    }

private:
    template<bool row>
    const T* extract(size_t i, T* buffer, size_t first, size_t last, HDF5DenseWorkspace& work) const {
        if constexpr(row != transpose) {
            work.offset[0] = i;
            work.offset[1] = first;
            work.count[0] = 1;
            work.count[1] = last - first;
        } else {
            work.offset[0] = first;
            work.offset[1] = i;
            work.count[0] = last - first;
            work.count[1] = 1;
        }

        work.dataspace.selectHyperslab(H5S_SELECT_SET, work.count, work.offset);

        if (static_cast<size_t>(work.current_size) != last - first) {
            work.current_size = last - first;
            work.memspace.setExtentSimple(1, &(work.current_size));
            work.memspace.selectAll();
        }

        work.dataset.read(buffer, HDF5::define_mem_type<T>(), work.memspace, work.dataspace);
        return buffer;
    }

    template<bool row>
    const T* extract(size_t i, T* buffer, size_t first, size_t last, Workspace* work) const {
        if (work) {
            auto wptr = dynamic_cast<HDF5DenseWorkspace*>(work);
            return extract<row>(i, buffer, first, last, *wptr);
        } else {
            auto full = new_workspace(row);
            auto wptr = dynamic_cast<HDF5DenseWorkspace*>(full.get());
            return extract<row>(i, buffer, first, last, *wptr);
        }
    }

public:
    const T* row(size_t r, T* buffer, size_t first, size_t last, Workspace* work=nullptr) const {
        return extract<true>(r, buffer, first, last, work);
    }

    const T* column(size_t c, T* buffer, size_t first, size_t last, Workspace* work=nullptr) const {
        return extract<false>(c, buffer, first, last, work);
    }

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::column;
};

}

#endif
