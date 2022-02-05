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
 */
template<bool ROW, typename T, typename IDX>
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
    HDF5CompressedSparseMatrix(size_t nr, size_t nc, std::string file, std::string vals, std::string idx, std::string ptr) :
        nrows(nr),
        ncols(nc),
        file_name(std::move(file)), 
        data_name(std::move(vals)),
        index_name(std::move(idx)),
        pointers(ROW ? nr + 1 : nc + 1)
    {
        H5::H5File fhandle(file_name, H5F_ACC_RDONLY);

        auto check_size = [&](const H5::DataSet& dhandle, const std::string& name) -> size_t {
            auto dtype = dhandle.getTypeClass();
            if (dtype != H5T_INTEGER && dtype != H5T_FLOAT) { 
                throw std::runtime_error(std::string("expected numeric values in the '") + name + "' dataset");
            }
            auto space = dhandle.getSpace();
            
            int ndim = space.getSimpleExtentNdims();
            if (ndim != 1) {
                throw std::runtime_error(std::string("'") + name + "' should be a one-dimensional array");
            }

            hsize_t dims_out[1];
            space.getSimpleExtentDims(dims_out, NULL);
            return dims_out[0];
        };

        auto dhandle = fhandle.openDataSet(data_name);
        const size_t nonzeros = check_size(dhandle, "vals");

        auto ihandle = fhandle.openDataSet(index_name);
        if (check_size(ihandle, "idx") != nonzeros) {
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

public:
    /**
     * @cond
     */
    struct HDF5SparseWorkspace : public Workspace {
        HDF5SparseWorkspace(const std::string& fpath, const std::string& dpath, const std::string& ipath) :
            file(fpath, H5F_ACC_RDONLY),
            data(file.openDataSet(dpath)),
            index(file.openDataSet(ipath)),
            dataspace(data.getSpace())
        {}

        H5::H5File file;
        H5::DataSet data, index;

        H5::DataSpace dataspace;
        H5::DataSpace memspace;
        std::vector<T> dbuffer;
        std::vector<IDX> ibuffer;
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
            new HDF5SparseWorkspace(
                file_name, 
                data_name, 
                index_name
            )
        );
    }

private:
    size_t extract_primary(size_t i, T* dbuffer, IDX* ibuffer, size_t first, size_t last, HDF5SparseWorkspace& work) const {
        hsize_t offset = pointers[i];
        hsize_t count = pointers[i+1] - pointers[i];
        if (!count) {
            return 0;
        }

        work.dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
        work.memspace.setExtentSimple(1, &count);
        work.memspace.selectAll();

        if (dbuffer && ibuffer && last - first >= count) {
            // Read indices directly to the output space, avoid unnecessary copying.
            work.index.read(ibuffer, HDF5::define_mem_type<IDX>(), work.memspace, work.dataspace);

            size_t start = 0;
            size_t end = std::lower_bound(ibuffer, ibuffer + count, last) - ibuffer;
            if (ibuffer[0] < first) {
                start = std::lower_bound(ibuffer, ibuffer + end, first) - ibuffer;
                std::copy(ibuffer + start, ibuffer + end, ibuffer);
            }

            // Seeing if we can narrow the extraction space for the data.
            hsize_t delta = end - start;
            if (delta) {
                if (delta != count) {
                    offset += start;
                    work.dataspace.selectHyperslab(H5S_SELECT_SET, &delta, &offset);
                    work.memspace.setExtentSimple(1, &delta);
                    work.memspace.selectAll();
                }
                work.data.read(dbuffer, HDF5::define_mem_type<T>(), work.memspace, work.dataspace);
            }

            return delta;

        } else {
            // Copying is required here as the supplied space is not enough to hold the run contents.
            work.ibuffer.resize(count);
            work.index.read(work.ibuffer.data(), HDF5::define_mem_type<IDX>(), work.memspace, work.dataspace);

            size_t start = std::lower_bound(work.ibuffer.begin(), work.ibuffer.end(), first) - work.ibuffer.begin();
            size_t end = std::lower_bound(work.ibuffer.begin() + start, work.ibuffer.end(), last) - work.ibuffer.begin();

            if (ibuffer) {
                std::copy(work.ibuffer.begin() + start, work.ibuffer.begin() + end, ibuffer);
            } else {
                std::copy(work.ibuffer.begin() + start, work.ibuffer.begin() + end, work.ibuffer.begin());
            }

            hsize_t delta = end - start;
            if (delta) {
                if (delta != count) {
                    offset += start;
                    work.dataspace.selectHyperslab(H5S_SELECT_SET, &delta, &offset);
                    work.memspace.setExtentSimple(1, &delta);
                    work.memspace.selectAll();
                }

                // At this point, the supplied space _must_ be large enough to hold the run contents, so whatever.
                if (dbuffer) {
                    work.data.read(dbuffer, HDF5::define_mem_type<T>(), work.memspace, work.dataspace);
                } else {
                    work.dbuffer.resize(delta);
                    work.data.read(work.dbuffer.data(), HDF5::define_mem_type<T>(), work.memspace, work.dataspace);
                }
            }

            return delta;
        }
    }

    SparseRange<T, IDX> extract_primary(size_t i, T* dbuffer, IDX* ibuffer, size_t first, size_t last, Workspace* work) const {
        if (work) {
            auto wptr = dynamic_cast<HDF5SparseWorkspace*>(work);
            return SparseRange<T, IDX>(extract_primary(i, dbuffer, ibuffer, first, last, *wptr), dbuffer, ibuffer);
        } else {
            auto full = new_workspace(ROW); // extract along the primary dimension.
            auto wptr = dynamic_cast<HDF5SparseWorkspace*>(full.get());
            return SparseRange<T, IDX>(extract_primary(i, dbuffer, ibuffer, first, last, *wptr), dbuffer, ibuffer);
        }
    }

private:
    size_t extract_secondary(size_t i, T* dbuffer, IDX* ibuffer, size_t first, size_t last, HDF5SparseWorkspace& work) const {
        auto& dhandle = work.data;
        auto& ihandle = work.index;

        auto& dataspace = work.dataspace;
        auto& memspace = work.memspace;

        std::vector<IDX>& index = work.ibuffer;
        size_t counter = 0;

        // Looping over all secondary dimensions.
        for (size_t j = first; j < last; ++j) {
            size_t left = pointers[j], right = pointers[j + 1];
            index.resize(right - left);

            hsize_t offset = left;
            hsize_t count = index.size();
            dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
            memspace.setExtentSimple(1, &count);
            memspace.selectAll();
            ihandle.read(index.data(), HDF5::define_mem_type<IDX>(), memspace, dataspace);

            auto it = std::lower_bound(index.begin(), index.end(), i);
            if (it != index.end() && *it == i) {
                offset = left + (it - index.begin());
                count = 1;
                dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
                memspace.setExtentSimple(1, &count);
                memspace.selectAll();
                dhandle.read(dbuffer + counter, HDF5::define_mem_type<T>(), memspace, dataspace);
                ibuffer[counter] = j;
                ++counter;
            }
        }

        return counter;
    }

    SparseRange<T, IDX> extract_secondary(size_t i, T* dbuffer, IDX* ibuffer, size_t first, size_t last, Workspace* work) const {
        if (work) {
            auto wptr = dynamic_cast<HDF5SparseWorkspace*>(work);
            return SparseRange<T, IDX>(extract_secondary(i, dbuffer, ibuffer, first, last, *wptr), dbuffer, ibuffer);
        } else {
            auto full = new_workspace(!ROW); // extract along the secondary dimension.
            auto wptr = dynamic_cast<HDF5SparseWorkspace*>(full.get());
            return SparseRange<T, IDX>(extract_secondary(i, dbuffer, ibuffer, first, last, *wptr), dbuffer, ibuffer);
        }
    }

public:
    SparseRange<T, IDX> sparse_row(size_t r, T* dbuffer, IDX* ibuffer, size_t first, size_t last, Workspace* work=nullptr, bool sorted=true) const {
        if constexpr(ROW) {
            return extract_primary(r, dbuffer, ibuffer, first, last, work);
        } else {
            return extract_secondary(r, dbuffer, ibuffer, first, last, work);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* dbuffer, IDX* ibuffer, size_t first, size_t last, Workspace* work=nullptr, bool sorted=true) const {
        if constexpr(ROW) {
            return extract_secondary(c, dbuffer, ibuffer, first, last, work);
        } else {
            return extract_primary(c, dbuffer, ibuffer, first, last, work);
        }
    }

    using Matrix<T, IDX>::sparse_row;

    using Matrix<T, IDX>::sparse_column;

private:
    const T* expand(size_t number, const T* value, const IDX* index, size_t first, size_t last, T* buffer) const {
        std::fill(buffer, buffer + last - first, 0);
        for (size_t i = 0; i < number; ++i) {
            buffer[index[i] - first] = value[i];
        }
        return buffer;
    }

    const T* expand_primary(size_t i, T* buffer, size_t first, size_t last, Workspace* work) const {
        T* dummy_d = 0;
        IDX* dummy_i = 0;
        if (work) {
            auto wptr = dynamic_cast<HDF5SparseWorkspace*>(work);
            size_t n = extract_primary(i, dummy_d, dummy_i, first, last, *wptr);
            return expand(n, wptr->dbuffer.data(), wptr->ibuffer.data(), first, last, buffer);
        } else {
            auto full = new_workspace(ROW);
            auto wptr = dynamic_cast<HDF5SparseWorkspace*>(full.get());
            size_t n = extract_primary(i, dummy_d, dummy_i, first, last, *wptr); 
            return expand(n, wptr->dbuffer.data(), wptr->ibuffer.data(), first, last, buffer);
        }
    }

    const T* expand_secondary(size_t i, T* buffer, size_t first, size_t last, Workspace* work) const {
        std::vector<T> dbuffer(last - first);
        std::vector<IDX> ibuffer(last - first);
        auto dptr = dbuffer.data();
        auto iptr = ibuffer.data();

        if (work) {
            auto wptr = dynamic_cast<HDF5SparseWorkspace*>(work);
            size_t n = extract_secondary(i, dptr, iptr, first, last, *wptr);
            return expand(n, dptr, iptr, first, last, buffer);
        } else {
            auto full = new_workspace(ROW);
            auto wptr = dynamic_cast<HDF5SparseWorkspace*>(full.get());
            size_t n = extract_secondary(i, dptr, iptr, first, last, *wptr);
            return expand(n, dptr, iptr, first, last, buffer);
        }
    }

public:
    const T* row(size_t r, T* buffer, size_t first, size_t last, Workspace* work=nullptr) const {
        if constexpr(ROW) {
            return expand_primary(r, buffer, first, last, work);
        } else {
            return expand_secondary(r, buffer, first, last, work);
        }
    }

    const T* column(size_t c, T* buffer, size_t first, size_t last, Workspace* work=nullptr) const {
        if constexpr(ROW) {
            return expand_secondary(c, buffer, first, last, work);
        } else {
            return expand_primary(c, buffer, first, last, work);
        }
    }

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::column;
};

}

#endif
