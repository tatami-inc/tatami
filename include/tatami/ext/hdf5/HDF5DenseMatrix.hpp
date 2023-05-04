#ifndef TATAMI_HDF5_DENSE_MATRIX_HPP
#define TATAMI_HDF5_DENSE_MATRIX_HPP

#include "H5Cpp.h"

#include <string>
#include <cstdint>
#include <type_traits>
#include <cmath>
#include <list>
#include <vector>

#include "../../base/dense/VirtualDenseMatrix.hpp"
#include "../../utils/Oracles.hpp"
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
 * @tparam Value_ Type of the matrix values.
 * @tparam Index_ Type of the row/column indices.
 * @tparam transpose_ Whether the dataset is transposed in its storage order, i.e., rows in HDF5 are columns in this matrix.
 */
template<typename Value_, typename Index_, bool transpose_ = false>
class HDF5DenseMatrix : public VirtualDenseMatrix<Value_, Index_> {
    Index_ firstdim, seconddim;
    std::string file_name, dataset_name;

    Index_ chunk_firstdim, chunk_seconddim;
    size_t total_cache_size;
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

        // Favoring extraction on the dimension that involves pulling out fewer chunks per dimension element.
        double nchunks_firstdim = static_cast<double>(firstdim)/static_cast<double>(chunk_firstdim);
        double nchunks_seconddim = static_cast<double>(seconddim)/static_cast<double>(chunk_seconddim);
        prefer_firstdim = (nchunks_firstdim > nchunks_seconddim);

        total_cache_size = static_cast<double>(cache_limit) / sizeof(Value_);

#ifndef TATAMI_HDF5_PARALLEL_LOCK
        }
#else
        });
#endif

        return;
    }

public:
    Index_ nrow() const {
        if constexpr(transpose_) {
            return seconddim;
        } else {
            return firstdim;
        }
    }

    Index_ ncol() const {
        if constexpr(transpose_) {
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
        if constexpr(transpose_) {
            return !prefer_firstdim;
        } else {
            return prefer_firstdim;
        }
    }

    bool uses_oracle(bool) const {
        return true;
    }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

private:
    template<bool accrow_>
    struct OracleCache {
        std::unordered_map<Index_, Index_> cache_exists, next_cache_exists;
        std::vector<std::vector<Value_> > cache_data, next_cache_data;
        std::vector<std::pair<Index_, Index_> > chunks_in_need; 
        typename std::conditional<accrow_ == transpose_, std::vector<std::pair<Index_, Index_> >, bool>::type cache_transpose_info;

        OracleStream<Index_> prediction_stream;
        std::vector<std::pair<Index_, Index_> > predictions_made;
        size_t predictions_fulfilled = 0;
    };

    struct LruCache {
        typedef std::pair<std::vector<Value_>, Index_> Element;
        std::list<Element> cache_data;
        std::unordered_map<Index_, typename std::list<Element>::iterator> cache_exists;
    };

    template<bool accrow_>
    struct Workspace {
        void fill(const HDF5DenseMatrix* parent, Index_ other_dim) {
            // Turn off HDF5's caching, as we'll be handling that. This allows us
            // to parallelize extractions without locking when the data has already
            // been loaded into memory; if we just used HDF5's cache, we would have
            // to lock on every extraction, given the lack of thread safety.
            H5::FileAccPropList fapl(H5::FileAccPropList::DEFAULT.getId());
            fapl.setCache(0, 0, 0, 0);

            file.openFile(parent->file_name, H5F_ACC_RDONLY, fapl);
            dataset = file.openDataSet(parent->dataset_name);
            dataspace = dataset.getSpace();

            auto chunk_dim = (accrow_ != transpose_ ? parent->chunk_firstdim : parent->chunk_seconddim);
            per_cache_size = static_cast<size_t>(chunk_dim) * static_cast<size_t>(other_dim);
            num_chunks = static_cast<double>(parent->total_cache_size) / per_cache_size;

            historian.reset(new LruCache);
        }

    public:
        // HDF5 members.
        H5::H5File file;
        H5::DataSet dataset;
        H5::DataSpace dataspace;
        H5::DataSpace memspace;

    public:
        // Caching members.
        size_t per_cache_size;
        Index_ num_chunks;
        typename std::conditional<accrow_ == transpose_, std::vector<Value_>, bool>::type transposition_buffer;

        // Cache with an oracle.
        std::unique_ptr<OracleCache<accrow_> > futurist;

        // Cache without an oracle.
        std::unique_ptr<LruCache> historian;
    };

private:
    template<bool accrow_, typename ExtractType_>
    static void extract_base(Index_ primary_start, Index_ primary_length, Value_* target, const ExtractType_& extract_value, Index_ extract_length, Workspace<accrow_>& work) {
        hsize_t offset[2];
        hsize_t count[2];

        constexpr int dimdex = (accrow_ != transpose_);
        offset[1-dimdex] = primary_start;
        count[1-dimdex] = primary_length;

        constexpr bool indexed = std::is_same<ExtractType_, std::vector<Index_> >::value;

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

        work.dataset.read(target, HDF5::define_mem_type<Value_>(), work.memspace, work.dataspace);
    }

    template<bool accrow_, typename ExtractType_>
    static Index_ extract_chunk(Index_ chunk_id, Index_ dim, Index_ chunk_dim, Value_* target, const ExtractType_& extract_value, Index_ extract_length, Workspace<accrow_>& work) {
        Index_ chunk_start = chunk_id * chunk_dim;
        Index_ chunk_end = std::min(dim, chunk_start + chunk_dim);
        Index_ chunk_actual = chunk_end - chunk_start;
        extract_base<accrow_>(chunk_start, chunk_actual, target, extract_value, extract_length, work);
        return chunk_actual;
    }

    static void transpose(std::vector<Value_>& cache, std::vector<Value_>& buffer, Index_ actual_dim, Index_ extract_length) {
        buffer.resize(cache.size());
        auto output = buffer.begin();
        for (Index_ x = 0; x < actual_dim; ++x, output += extract_length) {
            auto in = cache.begin() + x;
            for (Index_ y = 0; y < extract_length; ++y, in += actual_dim) {
                *(output + y) = *in;
            }
        }
        cache.swap(buffer);
        return;
    }

private:
    template<bool accrow_, typename ExtractType_>
    const Value_* extract_without_cache(Index_ i, Value_* buffer, const ExtractType_& extract_value, Index_ extract_length, Workspace<accrow_>& work) const {
#ifndef TATAMI_HDF5_PARALLEL_LOCK
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

            extract_base<accrow_>(i, 1, buffer, extract_value, extract_length, work);

#ifndef TATAMI_HDF5_PARALLEL_LOCK
        }
#else
        });
#endif

        return buffer;
    }

    template<bool accrow_, typename ExtractType_>
    const Value_* extract_with_oracle(Index_ mydim, Index_ chunk_mydim, const ExtractType_& extract_value, Index_ extract_length, Workspace<accrow_>& work) const {
        auto& pred = *work.futurist;
        if (pred.predictions_made.size() > pred.predictions_fulfilled) {
            const auto& chosen = pred.predictions_made[pred.predictions_fulfilled++];
            return pred.cache_data[chosen.first].data() + extract_length * chosen.second;
        }

        /* We use the oracle to batch together multiple HDF5 library
         * calls so that we don't have to reenter the critical section
         * for every fetch() call. Extracted data are cached for 
         * easy parallel retrieval without involving the HDF5 library.
         *
         * We don't bother trying to bundle everything into a single
         * HDF5 call, as this is a pain to pull data out into the
         * individual cache buffers afterwards.
         */
        Index_ num_chunks = work.num_chunks;
        Index_ used = 0;

        bool starting = pred.cache_data.empty();
        if (starting) {
            pred.cache_data.resize(num_chunks);
            pred.next_cache_data.resize(num_chunks, std::vector<Value_>(1)); // using a placeholder value for availability checks.
            pred.chunks_in_need.reserve(num_chunks);
        } else {
            pred.next_cache_exists.clear();
            pred.chunks_in_need.clear();
        }

        size_t max_predictions = static_cast<size_t>(num_chunks) * chunk_mydim * 2; // double the cache size, basically.
        pred.predictions_made.clear();
        pred.predictions_made.reserve(max_predictions);

        for (size_t p = 0; p < max_predictions; ++p) {
            Index_ current;
            if (!pred.prediction_stream.next(current)) {
                break;
            }

            auto curchunk = current / chunk_mydim;
            auto curindex = current % chunk_mydim;

            auto it = pred.next_cache_exists.find(curchunk);
            if (it == pred.next_cache_exists.end()) {
                pred.next_cache_exists[curchunk] = used;
                pred.predictions_made.emplace_back(used, curindex);

                auto it2 = pred.cache_exists.find(curchunk);
                if (it2 != pred.cache_exists.end()) {
                    pred.next_cache_data[used].swap(pred.cache_data[it2->second]);
                } else {
                    pred.chunks_in_need.emplace_back(curchunk, used);
                }

                ++used;
                if (used == num_chunks) {
                    break;
                }
            } else {
                pred.predictions_made.emplace_back(it->second, curindex);
            }
        }

        /**
         * Doing a linear scan across chunks to find the currently unused cache
         * elements. This is the simplest and safest approach; trying to keep
         * track of the unused caches would require a scan to fill a map/set
         * anyway, and deleting an iterator from cache_exists only works if
         * cache_exists actually contains all cache elements, which it might
         * not be if we reached max_predictions without filling up the cache.
         *
         * In any case, a linear scan should be pretty good for consecutive
         * access; the first cache elements would be the oldest, so there
         * wouldn't be any wasted iterations to find available cache elements.
         */
        size_t search = 0;
        for (const auto& c : pred.chunks_in_need) {
            // If we're lucky enough to already be sitting on an allocation
            // (e.g., if the last call to this function didn't consume all
            // cache elements), we just use it directly.
            if (pred.next_cache_data[c.second].empty()) { 
                while (pred.cache_data[search].empty()) { 
                    ++search;
                }

#ifdef DEBUG
                // This should never reach the end, as there should always be
                // 'num_chunks' non-empty elements across cache_data and
                // next_cache_data. Nonetheless, we add an assertion here.
                if (search >= pred.cache_data.size()) {
                    throw std::runtime_error("internal cache management error for HDF5DenseMatrix with oracles");
                }
#endif

                pred.next_cache_data[c.second].swap(pred.cache_data[search]);
                ++search;
            }

            pred.next_cache_data[c.second].resize(work.per_cache_size); // eventually no-op when all available caches are of the right size.
        }

        if constexpr(accrow_ == transpose_) {
            pred.cache_transpose_info.clear();
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

    template<bool accrow_, typename ExtractType_>
    const Value_* extract_without_oracle(Index_ i, Index_ mydim, Index_ chunk_mydim, const ExtractType_& extract_value, Index_ extract_length, Workspace<accrow_>& work) const {
        auto chunk = i / chunk_mydim;
        auto index = i % chunk_mydim;

        auto& historian = *(work.historian);
        auto it = historian.cache_exists.find(chunk);
        if (it != historian.cache_exists.end()) {
            auto chosen = it->second;
            historian.cache_data.splice(historian.cache_data.end(), historian.cache_data, chosen); // move to end.
            return chosen->first.data() + index * extract_length;
        }

        // Adding a new last element, or recycling the front to the back.
        typename std::list<typename LruCache::Element>::iterator location;
        if (historian.cache_data.size() < static_cast<size_t>(work.num_chunks)) {
            historian.cache_data.emplace_back(std::vector<Value_>(work.per_cache_size), chunk);
            location = std::prev(historian.cache_data.end());
        } else {
            location = historian.cache_data.begin();
            historian.cache_exists.erase(location->second);
            location->second = chunk;
            historian.cache_data.splice(historian.cache_data.end(), historian.cache_data, location); // move to end.
        }
        historian.cache_exists[chunk] = location;

        auto& cache_target = location->first;
        Index_ actual_dim;

#ifndef TATAMI_HDF5_PARALLEL_LOCK
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        actual_dim = extract_chunk<accrow_>(chunk, mydim, chunk_mydim, cache_target.data(), extract_value, extract_length, work);

#ifndef TATAMI_HDF5_PARALLEL_LOCK
        }
#else
        });
#endif

        // Applying a transposition for easier retrieval, but only once the lock is released.
        if constexpr(accrow_ == transpose_) {
            transpose(cache_target, work.transposition_buffer, actual_dim, extract_length);
        }

        return cache_target.data() + index * extract_length;
    }

    template<bool accrow_, typename ExtractType_>
    const Value_* extract(Index_ i, Value_* buffer, const ExtractType_& extract_value, Index_ extract_length, Workspace<accrow_>& work) const {
        // If there isn't any space for caching, we just extract directly.
        if (work.num_chunks == 0) {
            return extract_without_cache(i, buffer, extract_value, extract_length, work);
        }

        Index_ chunk_mydim, mydim;
        if constexpr(accrow_ != transpose_) {
            chunk_mydim = chunk_firstdim;
            mydim = firstdim;
        } else {
            chunk_mydim = chunk_seconddim;
            mydim = seconddim;
        }

        const Value_* cache;
        if (work.futurist) {
            cache = extract_with_oracle(mydim, chunk_mydim, extract_value, extract_length, work);
        } else {
            cache = extract_without_oracle(i, mydim, chunk_mydim, extract_value, extract_length, work);;
        }

        std::copy(cache, cache + extract_length, buffer);
        return buffer;
    }

private:
    template<bool accrow_, DimensionSelectionType selection_>
    struct Hdf5Extractor : public Extractor<selection_, false, Value_, Index_> {
        Hdf5Extractor(const HDF5DenseMatrix* p) : parent(p) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = (accrow_ ? parent->ncol() : parent->nrow());
                base.fill(parent, this->full_length); 
            }
        }

        Hdf5Extractor(const HDF5DenseMatrix* p, Index_ start, Index_ length) : parent(p) {
            if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = start;
                this->block_length = length;
                base.fill(parent, this->block_length); 
            }
        }

        Hdf5Extractor(const HDF5DenseMatrix* p, std::vector<Index_> idx) : parent(p) {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                this->index_length = idx.size();
                indices = std::move(idx);
                base.fill(parent, this->index_length); 
            }
        }

    protected:
        const HDF5DenseMatrix* parent;
        Workspace<accrow_> base;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type indices;

    public:
        const Index_* index_start() const {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                return indices.data();
            } else {
                return NULL;
            }
        }

        const Value_* fetch(Index_ i, Value_* buffer) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                return parent->extract<accrow_>(i, buffer, 0, this->full_length, this->base);
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                return parent->extract<accrow_>(i, buffer, this->block_start, this->block_length, this->base);
            } else {
                return parent->extract<accrow_>(i, buffer, this->indices, this->index_length, this->base);
            }
        }

        void set_oracle(std::unique_ptr<Oracle<Index_> > o) {
            base.futurist.reset(new OracleCache<accrow_>);
            base.futurist->prediction_stream.set(std::move(o));
            base.historian.reset();
        }
    };

    template<bool accrow_, DimensionSelectionType selection_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, false, Value_, Index_> > populate(const Options& opt, Args_... args) const {
        std::unique_ptr<Extractor<selection_, false, Value_, Index_> > output;

#ifndef TATAMI_HDF5_PARALLEL_LOCK
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        output.reset(new Hdf5Extractor<accrow_, selection_>(this, args...));

#ifndef TATAMI_HDF5_PARALLEL_LOCK
        }
#else
        });
#endif

        return output;
    }

public:
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        return populate<true, DimensionSelectionType::FULL>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<true, DimensionSelectionType::BLOCK>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        return populate<true, DimensionSelectionType::INDEX>(opt, std::move(indices));
    }

    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        return populate<false, DimensionSelectionType::FULL>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<false, DimensionSelectionType::BLOCK>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        return populate<false, DimensionSelectionType::INDEX>(opt, std::move(indices));
    }
};

}

#endif
