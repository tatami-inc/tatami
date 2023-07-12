#ifndef TATAMI_CUSTOM_CHUNKED_MATRIX_HPP
#define TATAMI_CUSTOM_CHUNKED_MATRIX_HPP

#include "../base/Matrix.hpp"
#include "CustomChunkManager.hpp"
#include "utils.hpp"

namespace tatami {

/**
 * @brief Options for custom chunk extraction.
 */
struct CustomChunkedOptions : public TypicalChunkCacheOptions {};

/**
 * @cond
 */
template<typename Index_, bool sparse_, class Chunk_>
class CustomChunkedMatrixMethods {
protected:
    CustomChunkedMatrixMethods(Index_ mat_nr, Index_ mat_nc, Index_ chunk_nr, Index_ chunk_nc, std::vector<Chunk_> chunks, bool rm) :
        mat_nrow(mat_nr),
        mat_ncol(mat_nc),
        chunk_nrow(chunk_nr), 
        chunk_ncol(chunk_nc),
        num_chunks_per_row(integer_ceil(mat_ncol, chunk_ncol)),
        num_chunks_per_column(integer_ceil(mat_nrow, chunk_nrow)),
        chunk_array(std::move(chunks)),
        row_major(rm)
    {}

protected:
    Index_ mat_nrow, mat_ncol;
    Index_ chunk_nrow, chunk_ncol;
    Index_ num_chunks_per_row, num_chunks_per_column;
    std::vector<Chunk_> chunk_array;
    bool row_major;

    bool prefer_rows_internal() const {
        // Prefer rows if we have to extract fewer chunks per row.
        return num_chunks_per_column > num_chunks_per_row; 
    }

protected:
    template<bool accrow_>
    Index_ get_primary_dim() const {
        if constexpr(accrow_) {
            return mat_nrow;
        } else {
            return mat_ncol;
        }
    }

    template<bool accrow_>
    Index_ get_secondary_dim() const {
        if constexpr(accrow_) {
            return mat_ncol;
        } else {
            return mat_nrow;
        }
    }

    template<bool accrow_>
    Index_ get_primary_chunkdim() const {
        if constexpr(accrow_) {
            return chunk_nrow;
        } else {
            return chunk_ncol;
        }
    }

    template<bool accrow_>
    Index_ get_secondary_chunkdim() const {
        if constexpr(accrow_) {
            return chunk_ncol;
        } else {
            return chunk_nrow;
        }
    }

    template<bool accrow_>
    Index_ get_secondary_num_chunks() const {
        if constexpr(accrow_) {
            return num_chunks_per_row;
        } else {
            return num_chunks_per_column;
        }
    }

protected:
    static Index_ integer_ceil(Index_ left, Index_ right) {
        return left / right + (left % right > 0); // avoids overflow.
    }

protected:
    typedef std::vector<typename Chunk_::value_type> DenseSlab;

    struct SparseSlab {
        SparseSlab() = default;
        SparseSlab(size_t primary_dim) : indices(primary_dim), values(primary_dim) {}
        std::vector<std::vector<typename Chunk_::value_type> > values;
        std::vector<std::vector<typename Chunk_::index_type> > indices;

        bool empty() const {
            return values.empty();
        }

        void swap(SparseSlab& right) {
            values.swap(right.values);
            indices.swap(right.indices);
        }

        void resize(size_t primary_dim) {
            indices.resize(primary_dim);
            values.resize(primary_dim);
        }
    };

    typedef typename std::conditional<sparse_, SparseSlab, DenseSlab>::type Slab;

protected:
    template<bool accrow_, DimensionSelectionType selection_, bool exact_, class Extractor_>
    void extract(Index_ chunk_id, Index_ chunk_offset, Slab& slab, Extractor_* ext) const {
        auto primary_chunkdim = get_primary_chunkdim<accrow_>();
        auto secondary_chunkdim = get_secondary_chunkdim<accrow_>();

        size_t offset, increment; // use size_t to avoid overflow.
        if (row_major) {
            if constexpr(accrow_) {
                offset = chunk_id * num_chunks_per_row;
                increment = 1;
            } else {
                offset = chunk_id;
                increment = num_chunks_per_row;
            }
        } else {
            if constexpr(accrow_) {
                offset = chunk_id;
                increment = num_chunks_per_column;
            } else {
                offset = chunk_id * num_chunks_per_column;
                increment = 1;
            }
        }

        Index_ primary_dim = get_primary_dim<accrow_>();
        Index_ secondary_dim = get_secondary_dim<accrow_>();

        Index_ primary_start_pos, primary_len;
        if constexpr(exact_) {
            primary_start_pos = chunk_offset;
            primary_len = 1;
        } else {
            primary_start_pos = 0;
            primary_len = std::min(primary_chunkdim, primary_dim - chunk_id * primary_chunkdim); // avoid running off the end.
        }

        auto slab_ptr = [&]{
            if constexpr(!sparse_) {
                return slab.data();
            } else {
                for (auto& x : slab.indices) {
                    x.clear();
                }
                for (auto& x : slab.values) {
                    x.clear();
                }
                return false;
            }
        }();

        if constexpr(selection_ == DimensionSelectionType::INDEX) {
            if (!ext->indices.empty()) {
                auto secondary_num_chunks = get_secondary_num_chunks<accrow_>();
                Index_ start_chunk_index = ext->indices.front() / secondary_chunkdim; // 'indices' is guaranteed to be non-empty at this point.
                offset += static_cast<size_t>(start_chunk_index) * increment; // use size_t to avoid integer overflow.
                Index_ secondary_start_pos = start_chunk_index * secondary_chunkdim;

                auto iIt = ext->indices.begin();
                auto iEnd = ext->indices.end();
                for (Index_ c = start_chunk_index; c < secondary_num_chunks; ++c) {
                    const auto& chunk = chunk_array[offset];

                    Index_ secondary_end_pos = std::min(secondary_dim - secondary_start_pos, secondary_chunkdim) + secondary_start_pos; // avoid overflow.
                    ext->chunk_indices.clear();
                    while (iIt != iEnd && *iIt < secondary_end_pos) {
                        ext->chunk_indices.push_back(*iIt - secondary_start_pos);
                        ++iIt;
                    }

                    if (!ext->chunk_indices.empty()) {
                        if constexpr(sparse_) {
                            chunk.template extract<accrow_>(primary_start_pos, primary_len, ext->chunk_indices, ext->chunk_workspace, slab.values, slab.indices, secondary_start_pos);
                        } else {
                            chunk.template extract<accrow_>(primary_start_pos, primary_len, ext->chunk_indices, ext->chunk_workspace, slab_ptr, ext->indices.size());
                        }
                    }

                    secondary_start_pos += secondary_chunkdim; 
                    offset += increment;
                    if constexpr(!sparse_) {
                        slab_ptr += ext->chunk_indices.size();
                    }
                }
            }

        } else {
            Index_ start, len;
            if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                start = ext->block_start;
                len = ext->block_length;
            } else {
                start = 0;
                len = ext->full_length;
            }
            Index_ end = start + len;

            Index_ start_chunk_index = start / secondary_chunkdim;
            Index_ secondary_start_pos = start_chunk_index * secondary_chunkdim;
            Index_ end_chunk_index = integer_ceil(end, secondary_chunkdim);
            offset += increment * static_cast<size_t>(start_chunk_index); // size_t to avoid integer overflow.

            for (Index_ c = start_chunk_index; c < end_chunk_index; ++c) {
                const auto& chunk = chunk_array[offset];
                Index_ from = (c == start_chunk_index ? start - secondary_start_pos : 0);
                Index_ to = (c + 1 == end_chunk_index ? end - secondary_start_pos : secondary_chunkdim);

                // No need to protect against a zero length, as it should be impossible
                // here (otherwise, start_chunk_index == end_chunk_index and we'd never iterate).
                if constexpr(sparse_) {
                    chunk.template extract<accrow_>(primary_start_pos, primary_len, from, to - from, ext->chunk_workspace, slab.values, slab.indices, secondary_start_pos);
                } else {
                    chunk.template extract<accrow_>(primary_start_pos, primary_len, from, to - from, ext->chunk_workspace, slab_ptr, len);
                }

                secondary_start_pos += secondary_chunkdim;
                offset += increment;
                if constexpr(!sparse_) {
                    slab_ptr += (to - from);
                }
            }
        }
    }

public:
    template<bool accrow_, DimensionSelectionType selection_, class Extractor_>
    std::pair<const Slab*, Index_> fetch_cache(Index_ i, Extractor_* ext) const {
        Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
        size_t alloc = sparse_ ? primary_chunkdim : primary_chunkdim * static_cast<size_t>(extracted_length<selection_, Index_>(*ext)); // use size_t to avoid integer overflow.
        auto& cache_workspace = ext->cache_workspace;

        if (cache_workspace.oracle_cache) {
            return cache_workspace.oracle_cache->next_chunk(
                /* identify = */ [&](Index_ i) -> std::pair<Index_, Index_> {
                    return std::make_pair(i / primary_chunkdim, i % primary_chunkdim);
                },
                /* swap = */ [&](Slab& left, Slab& right) -> void {
                    left.swap(right);
                },
                /* ready = */ [&](const Slab& slab) -> bool {
                    return !slab.empty();
                },
                /* allocate = */ [&](Slab& slab) -> void {
                    slab.resize(alloc);
                },
                /* populate =*/ [&](const std::vector<std::pair<Index_, Index_> >& chunks_in_need, std::vector<Slab>& chunk_data) -> void {
                    for (const auto& p : chunks_in_need) {
                        extract<accrow_, selection_, false>(p.first, /* no-op */ 0, chunk_data[p.second], ext);
                    }
                }
            );

        } else {
            auto chunk_id = i / primary_chunkdim;
            auto chunk_offset = i % primary_chunkdim;

            if (cache_workspace.num_chunk_sets_in_cache == 0) {
                extract<accrow_, selection_, true>(chunk_id, chunk_offset, ext->solo, ext);
                return std::make_pair(&(ext->solo), static_cast<Index_>(0));

            } else {
                auto& cache = cache_workspace.lru_cache->find_chunk(
                    chunk_id,
                    /* create = */ [&]() -> Slab {
                        return Slab(alloc);
                    },
                    /* populate = */ [&](Index_ id, Slab& slab) -> void {
                        extract<accrow_, selection_, false>(id, /* no-op */ 0, slab, ext);
                    }
                );
                return std::make_pair(&cache, chunk_offset);
            }
        }
    }
};
/**
 * @endcond
 */

/**
 * @brief Matrix of custom dense chunks.
 *
 * @tparam Value_ Numeric type for the matrix value.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Chunk_ Class of the chunk.
 *
 * Implements a `Matrix` subclass where data is contained in dense rectangular chunks.
 * These chunks are typically compressed in some manner to reduce memory usage;
 * on access, each chunk is decompressed and the desired values are extracted.
 * The `Chunk_` class should provide the following:
 *
 * - A `typedef value_type` specifying the type of the decompressed chunk data.
 * - A nested `Workspace` class that allocates memory for decompression.
 *   This should be default-constructible and may be re-used for decompressing multiple chunks.
 * - A `void extract(Index_ primary_start, Index_ primary_length, Index_ secondary_start, Index_ secondary_length, Workspace& work, Output_* output, size_t stride) const` method,
 *   which extracts contiguous ranges from the chunk and stores them in `output`.
 *   See `SimpleDenseChunkWrapper::extract()` for an example that describes the expected layout of the output.
 * - A `void extract(Index_ primary_start, Index_ primary_length, const std::vector<Index_>& secondary_indices, Workspace& work, Output_* output, size_t stride) const` method.
 *   which extracts an indexed subset from the chunk and stores them in `output`.
 *   See `SimpleDenseChunkWrapper::extract()` for an example that describes the expected layout of the output.
 *
 * All chunks should have the same dimensions, i.e., covering the same shape/area of the matrix.
 * The matrix should be partitioned at regular intervals starting from zero -
 * the first chunk should start at (0, 0), the next chunk should be immediately adjacent in one of the dimensions, and so on.
 * The exception is for chunks at the non-zero boundaries of the matrix dimensions, which may be truncated.
 */
template<typename Value_, typename Index_, typename Chunk_>
class CustomChunkedDenseMatrix : public VirtualDenseMatrix<Value_, Index_>, public CustomChunkedMatrixMethods<Index_, false, Chunk_> {
public:
    /**
     * @param mat_nrow Number of rows in the matrix.
     * @param mat_ncol Number of columns in the matrix.
     * @param chunk_nrow Number of rows in each chunk.
     * @param chunk_ncol Number of columns in each chunk.
     * @param chunks Vector containing a 2D array of chunks that cover the entire matrix.
     * This should have length equal to the product of the number of chunks along the rows and columns of the matrix, i.e., `ceil(mat_nrow / chunk_nrow) * ceil(mat_ncol / chunk_ncol)`.
     * @param row_major Whether `chunks` is in row-major format.
     * @param opt Further options for chunked extraction.
     */
    CustomChunkedDenseMatrix(Index_ mat_nrow, Index_ mat_ncol, Index_ chunk_nrow, Index_ chunk_ncol, std::vector<Chunk_> chunks, bool row_major, const CustomChunkedOptions& opt) : 
        CustomChunkedMatrixMethods<Index_, false, Chunk_>(mat_nrow, mat_ncol, chunk_nrow, chunk_ncol, std::move(chunks), row_major),
        cache_size_in_elements(opt.maximum_cache_size / sizeof(typename Chunk_::value_type)),
        require_minimum_cache(opt.require_minimum_cache)
    {
        if (static_cast<size_t>(this->num_chunks_per_column) * static_cast<size_t>(this->num_chunks_per_row) != this->chunk_array.size()) {
            throw std::runtime_error("length of 'chunks' should be equal to the product of the number of chunks along each row and column");
        }
    }

private:
    size_t cache_size_in_elements;
    bool require_minimum_cache;

public:
    Index_ nrow() const { return this->mat_nrow; }

    Index_ ncol() const { return this->mat_ncol; }

    bool prefer_rows() const { 
        return this->prefer_rows_internal();
    }

    bool uses_oracle(bool) const { 
        return true; 
    }

    double prefer_rows_proportion() const { 
        return static_cast<double>(this->prefer_rows_internal());
    }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

private:
    template<bool accrow_, DimensionSelectionType selection_>
    struct CustomExtractor : public Extractor<selection_, false, Value_, Index_> {
        CustomExtractor(const CustomChunkedDenseMatrix* p) : parent(p) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = parent->template get_secondary_dim<accrow_>();
            }
            initialize_cache();
        }

        CustomExtractor(const CustomChunkedDenseMatrix* p, Index_ bs, Index_ bl) : parent(p) {
            if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = bs;
                this->block_length = bl;
            }
            initialize_cache();
        }

        CustomExtractor(const CustomChunkedDenseMatrix* p, std::vector<Index_> idx) : parent(p) {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                indices = std::move(idx);
                this->index_length = indices.size();
            }
            initialize_cache();
        }

    private:
        const CustomChunkedDenseMatrix* parent;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type indices;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type chunk_indices;

        typedef typename CustomChunkedMatrixMethods<Index_, false, Chunk_>::Slab Slab;
        typename Chunk_::Workspace chunk_workspace;
        TypicalChunkCacheWorkspace<Index_, Slab> cache_workspace;
        Slab solo;

        void initialize_cache() {
            auto len = extracted_length<selection_, Index_>(*this);

            cache_workspace = TypicalChunkCacheWorkspace<Index_, Slab>(
                accrow_ ? parent->chunk_nrow : parent->chunk_ncol,
                len,
                parent->cache_size_in_elements,
                parent->require_minimum_cache
            );

            if (cache_workspace.num_chunk_sets_in_cache == 0) {
                solo.resize(len);
            }
        }

    public:
        const Index_* index_start() const {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                return indices.data();
            } else {
                return NULL;
            }
        }

        void set_oracle(std::unique_ptr<Oracle<Index_> > o) {
            cache_workspace.set_oracle(std::move(o));
            return;
        }

    public:
        friend class CustomChunkedMatrixMethods<Index_, false, Chunk_>;

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto fetched = parent->template fetch_cache<accrow_, selection_>(i, this);
            size_t len = extracted_length<selection_, Index_>(*this); // size_t to avoid overflow.
            auto ptr = fetched.first->data() + fetched.second * len;
            std::copy(ptr, ptr + len, buffer);
            return buffer;
        }
    };

public:
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        auto ptr = new CustomExtractor<true, DimensionSelectionType::FULL>(this);
        return std::unique_ptr<FullDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        auto ptr = new CustomExtractor<true, DimensionSelectionType::BLOCK>(this, block_start, block_length);
        return std::unique_ptr<BlockDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        auto ptr = new CustomExtractor<true, DimensionSelectionType::INDEX>(this, std::move(indices));
        return std::unique_ptr<IndexDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        auto ptr = new CustomExtractor<false, DimensionSelectionType::FULL>(this);
        return std::unique_ptr<FullDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        auto ptr = new CustomExtractor<false, DimensionSelectionType::BLOCK>(this, block_start, block_length);
        return std::unique_ptr<BlockDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        auto ptr = new CustomExtractor<false, DimensionSelectionType::INDEX>(this, std::move(indices));
        return std::unique_ptr<IndexDenseExtractor<Value_, Index_> >(ptr);
    }
};

/**
 * @brief Matrix of custom sparse chunks.
 *
 * @tparam Value_ Numeric type for the matrix value.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Chunk_ Class of the chunk.
 *
 * Implements a `Matrix` subclass where data is contained in sparse rectangular chunks.
 * These chunks are typically compressed in some manner to reduce memory usage;
 * on access, each chunk is decompressed and the desired values are extracted.
 * The `Chunk_` class should provide the following:
 *
 * - A `typedef value_type` specifying the type of the decompressed chunk data.
 * - A `typedef index_type` specifying the type of the indices in the decompressed chunk.
 * - A nested `Workspace` class that allocates memory for decompression.
 *   This should be default-constructible and may be re-used for decompressing multiple chunks.
 * - A `void extract(Index_ primary_start, Index_ primary_length, Index_ secondary_start, Index_ secondary_length, Workspace& work, 
 *   std::vector<std::vector<value_type> >& output_values, std::vector<std::vector<index_type> >& output_indices) const` method,
 *   which extracts contiguous ranges from the chunk and stores them in `output_values` and `output_indices`.
 *   See `SimpleSparseChunkWrapper::extract()` for an example that describes the expected layout of the output.
 * - A `void extract(Index_ primary_start, Index_ primary_length, Index_ secondary_start, Index_ secondary_length, Workspace& work, 
 *   std::vector<std::vector<value_type> >& output_values, std::vector<std::vector<index_type> >& output_indices) const` method,
 *   which extracts an indexed subset from the chunk and stores them in `output_values` and `output_indices`.
 *   See `SimpleSparseChunkWrapper::extract()` for an example that describes the expected layout of the output.
 *
 * All chunks should have the same dimensions, i.e., covering the same shape/area of the matrix.
 * The matrix should be partitioned at regular intervals starting from zero -
 * the first chunk should start at (0, 0), the next chunk should be immediately adjacent in one of the dimensions, and so on.
 * The exception is for chunks at the non-zero boundaries of the matrix dimensions, which may be truncated.
 */
template<typename Value_, typename Index_, typename Chunk_>
class CustomChunkedSparseMatrix : public VirtualDenseMatrix<Value_, Index_>, public CustomChunkedMatrixMethods<Index_, true, Chunk_> {
public:
    /**
     * @param mat_nrow Number of rows in the matrix.
     * @param mat_ncol Number of columns in the matrix.
     * @param chunk_nrow Number of rows in each chunk.
     * @param chunk_ncol Number of columns in each chunk.
     * @param chunks Vector containing a 2D array of chunks that cover the entire matrix.
     * This should have length equal to the product of the number of chunks along the rows and columns of the matrix, i.e., `ceil(mat_nrow / chunk_nrow) * ceil(mat_ncol / chunk_ncol)`.
     * @param row_major Whether `chunks` is in row-major format.
     * @param opt Further options for chunked extraction.
     */
    CustomChunkedSparseMatrix(Index_ mat_nrow, Index_ mat_ncol, Index_ chunk_nrow, Index_ chunk_ncol, std::vector<Chunk_> chunks, bool row_major, const CustomChunkedOptions& opt) : 
        CustomChunkedMatrixMethods<Index_, true, Chunk_>(mat_nrow, mat_ncol, chunk_nrow, chunk_ncol, std::move(chunks), row_major),
        cache_size_in_elements(opt.maximum_cache_size / (sizeof(typename Chunk_::value_type) + sizeof(typename Chunk_::index_type))),
        require_minimum_cache(opt.require_minimum_cache)
    {
        if (static_cast<size_t>(this->num_chunks_per_column) * static_cast<size_t>(this->num_chunks_per_row) != this->chunk_array.size()) {
            throw std::runtime_error("length of 'chunks' should be equal to the product of the number of chunks along each row and column");
        }
    }

private:
    size_t cache_size_in_elements;
    bool require_minimum_cache;

public:
    Index_ nrow() const { return this->mat_nrow; }

    Index_ ncol() const { return this->mat_ncol; }

    bool prefer_rows() const { 
        return this->prefer_rows_internal();
    }

    bool uses_oracle(bool) const { 
        return true; 
    }

    double prefer_rows_proportion() const { 
        return static_cast<double>(this->prefer_rows_internal());
    }

    bool sparse() const {
        return true;
    }

    double sparse_proportion() const {
        return 1;
    }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

private:
    template<bool accrow_, bool sparse_, DimensionSelectionType selection_>
    struct CustomExtractorBase : public Extractor<selection_, sparse_, Value_, Index_> {
        CustomExtractorBase(const CustomChunkedSparseMatrix* p) : parent(p) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = parent->template get_secondary_dim<accrow_>();
            }
            initialize_cache();
        }

        CustomExtractorBase(const CustomChunkedSparseMatrix* p, Index_ bs, Index_ bl) : parent(p) {
            if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = bs;
                this->block_length = bl;
            }
            initialize_cache();
        }

        CustomExtractorBase(const CustomChunkedSparseMatrix* p, std::vector<Index_> idx) : parent(p) {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                indices = std::move(idx);
                this->index_length = indices.size();
            }
            initialize_cache();
        }

    protected:
        const CustomChunkedSparseMatrix* parent;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type indices;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type chunk_indices;

        typedef typename CustomChunkedMatrixMethods<Index_, true, Chunk_>::Slab Slab;
        typename Chunk_::Workspace chunk_workspace;
        TypicalChunkCacheWorkspace<Index_, Slab> cache_workspace;
        Slab solo;

        void initialize_cache() {
            cache_workspace = TypicalChunkCacheWorkspace<Index_, Slab>(
                accrow_ ? parent->chunk_nrow : parent->chunk_ncol,
                extracted_length<selection_, Index_>(*this),
                parent->cache_size_in_elements,
                parent->require_minimum_cache
            );

            if (cache_workspace.num_chunk_sets_in_cache == 0) {
                solo.resize(1);
            }
        }

    public:
        const Index_* index_start() const {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                return indices.data();
            } else {
                return NULL;
            }
        }

        void set_oracle(std::unique_ptr<Oracle<Index_> > o) {
            cache_workspace.set_oracle(std::move(o));
            return;
        }

        friend class CustomChunkedMatrixMethods<Index_, true, Chunk_>;
    };

private:
    template<bool accrow_, DimensionSelectionType selection_>
    struct CustomDenseExtractor : public CustomExtractorBase<accrow_, false, selection_> {
    public:
        template<typename ... Args_>
        CustomDenseExtractor(const CustomChunkedSparseMatrix* p, Args_&& ... args) : CustomExtractorBase<accrow_, false, selection_>(p, std::forward<Args_>(args)...) {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                const auto& indices = this->indices;
                if (!indices.empty()) {
                    remap.resize(indices.back() + 1);
                    for (Index_ i = 0, end = indices.size(); i < end; ++i) {
                        remap[indices[i]] = i;
                    }
                }
            }
        }

    public:
        const Value_* fetch(Index_ i, Value_* buffer) {
            auto contents = this->parent->template fetch_cache<accrow_, selection_>(i, this);
            const auto& values = contents.first->values[contents.second];
            const auto& indices = contents.first->indices[contents.second];

            auto len = extracted_length<selection_, Index_>(*this);
            std::fill(buffer, buffer + len, 0);

            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                for (size_t i = 0, end = indices.size(); i < end; ++i) {
                    buffer[remap[indices[i]]] = values[i];
                }
            } else {
                Index_ start = 0;
                if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                    start = this->block_start;
                }

                for (size_t i = 0, end = indices.size(); i < end; ++i) {
                    buffer[indices[i] - start] = values[i];
                }
            }

            return buffer;
        }

    private:
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type remap;
    };

public:
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        auto ptr = new CustomDenseExtractor<true, DimensionSelectionType::FULL>(this);
        return std::unique_ptr<FullDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        auto ptr = new CustomDenseExtractor<true, DimensionSelectionType::BLOCK>(this, block_start, block_length);
        return std::unique_ptr<BlockDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        auto ptr = new CustomDenseExtractor<true, DimensionSelectionType::INDEX>(this, std::move(indices));
        return std::unique_ptr<IndexDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        auto ptr = new CustomDenseExtractor<false, DimensionSelectionType::FULL>(this);
        return std::unique_ptr<FullDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        auto ptr = new CustomDenseExtractor<false, DimensionSelectionType::BLOCK>(this, block_start, block_length);
        return std::unique_ptr<BlockDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        auto ptr = new CustomDenseExtractor<false, DimensionSelectionType::INDEX>(this, std::move(indices));
        return std::unique_ptr<IndexDenseExtractor<Value_, Index_> >(ptr);
    }

private:
    template<bool accrow_, DimensionSelectionType selection_>
    struct CustomSparseExtractor : public CustomExtractorBase<accrow_, true, selection_> {
    public:
        template<typename ... Args_>
        CustomSparseExtractor(const CustomChunkedSparseMatrix* p, bool ev, bool ei, Args_&& ... args) : 
            CustomExtractorBase<accrow_, true, selection_>(p, std::forward<Args_>(args)...), 
            extract_value(ev), 
            extract_index(ei)
        {}

    public:
        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            auto contents = this->parent->template fetch_cache<accrow_, selection_>(i, this);
            const auto& values = contents.first->values[contents.second];
            const auto& indices = contents.first->indices[contents.second];

            if (extract_value) {
                std::copy(values.begin(), values.end(), vbuffer);
            } else {
                vbuffer = NULL;
            }

            if (extract_index) {
                std::copy(indices.begin(), indices.end(), ibuffer);
            } else {
                ibuffer = NULL;
            }

            return SparseRange<Value_, Index_>(values.size(), vbuffer, ibuffer);
        }

    private:
        bool extract_value, extract_index;
    };

public:
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_row(const Options& opt) const {
        auto ptr = new CustomSparseExtractor<true, DimensionSelectionType::FULL>(this, opt.sparse_extract_value, opt.sparse_extract_index);
        return std::unique_ptr<FullSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        auto ptr = new CustomSparseExtractor<true, DimensionSelectionType::BLOCK>(this, opt.sparse_extract_value, opt.sparse_extract_index, block_start, block_length);
        return std::unique_ptr<BlockSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& opt) const {
        auto ptr = new CustomSparseExtractor<true, DimensionSelectionType::INDEX>(this, opt.sparse_extract_value, opt.sparse_extract_index, std::move(indices));
        return std::unique_ptr<IndexSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_column(const Options& opt) const {
        auto ptr = new CustomSparseExtractor<false, DimensionSelectionType::FULL>(this, opt.sparse_extract_value, opt.sparse_extract_index);
        return std::unique_ptr<FullSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        auto ptr = new CustomSparseExtractor<false, DimensionSelectionType::BLOCK>(this, opt.sparse_extract_value, opt.sparse_extract_index, block_start, block_length);
        return std::unique_ptr<BlockSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options& opt) const {
        auto ptr = new CustomSparseExtractor<false, DimensionSelectionType::INDEX>(this, opt.sparse_extract_value, opt.sparse_extract_index, std::move(indices));
        return std::unique_ptr<IndexSparseExtractor<Value_, Index_> >(ptr);
    }

};

}

#endif
