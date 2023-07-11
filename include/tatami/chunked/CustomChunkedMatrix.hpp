#ifndef TATAMI_CUSTOM_CHUNKED_MATRIX_HPP
#define TATAMI_CUSTOM_CHUNKED_MATRIX_HPP

#include "../base/Matrix.hpp"
#include "CustomChunkManager.hpp"
#include "utils.hpp"

namespace tatami {

/**
 * @brief Options for custom chunk extraction.
 */
struct CustomChunkedOptions : public TypicalBlockCacheOptions {};

template<typename Value_, typename Index_, typename Chunk_>
class CustomChunkedDenseMatrix : public VirtualDenseMatrix<Value_, Index_> {
public:
    CustomChunkedMatrix(Index_ mat_nrow, Index_ mat_ncol, Index_ chunk_nrow, Index_ chunk_ncol, Index_ nrow_in_chunks, Index_ ncol_in_chunks, std::vector<Chunk_> chunks, bool row_major, const CustomChunkedOptions& opt) : 
        manager(std::move(chunks)), 
        nrows(mat_nrow), 
        ncols(mat_ncol), 
        cache_size_in_elements(opt.maximum_cache_size / sizeof(Chunk_::value_type)),
        require_minimum_cache(opt.require_minimum_cache)
    {
        manager.chunk_nrow = chunk_nrow;
        manager.chunk_ncol = chunk_ncol;
        manager.num_chunks_per_column = nrow_in_chunks;
        manager.num_chunks_per_row = ncol_in_chunks;
        manager.row_major = row_major;
    }

private:
    CustomChunkManager<Chunk_> manager;
    Index_ nrows, ncols;
    size_t cache_size_in_elements;
    bool require_minimum_cache;

    static_assert(!Chunk_::sparse); // Should be dense, obviously.


public:
    Index_ nrow() const { return nrows; }

    Index_ ncol() const { return ncols; }

    bool prefer_rows() const { 
        // Prefer rows if we have to extract fewer chunks per row.
        return manager.num_chunks_per_column > manager.num_chunks_per_row; 
    }

    bool uses_oracle(bool) const { 
        return true; 
    }

    double prefer_rows_proportion() const { 
        return static_cast<double>(row_); 
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
                this->full_length = (accrow_ ? p->ncols : p->nrows);
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

        typedef typename CustomChunkManager::Slab Slab;
        typename CustomChunkManager::Workspace buffer;
        TypicalChunkCacheWorkspace<Index_, Slab> workspace;
        std::unique_ptr<Slab> solo;

        void initialize_cache() {
            workspace = TypicalChunkCacheWorkspace<Index_, Slab>(
                accrow_ ? parent->manager->chunk_nrow, parent->manager->chunk_ncol,
                extracted_length(this)
                parent->cache_size_in_elements,
                parent->require_minimum_cache
            );

            if (workspace.num_chunk_sets_in_cache == 0) {
                solo.reset(new Slab(parent->manager.create_slab<accrow_, true>(len)));
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
            workspace.set_oracle(std::move(o));
            return;
        }

    private:
        Index_ get_primary_dim() const {
            if constexpr(accrow_) {
                return parent->nrows;
            } else {
                return parent->ncols;
            }
        }

        Index_ get_primary_chunkdim() const {
            if constexpr(accrow_) {
                return parent->manager.chunk_nrow;
            } else {
                return parent->manager.chunk_ncol;
            }
        }

        template<bool exact_>
        void extract(Index_ chunk_id, Index_ chunk_offset, Slab& slab) const {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                manager.extract<accrow_, exact_>(chunk_id, chunk_offset, get_primary_dim(), static_cast<Index_>(0), this->full_length, slab, buffer);
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                manager.extract<accrow_, exact_>(chunk_id, chunk_offset, get_primary_dim(), this->block_start, this->block_length, slab, buffer);
            } else {
                manager.extract<accrow_, exact_>(chunk_id, chunk_offset, get_primary_dim(), indices, slab, buffer);
            }
        }

    public:
        const Value_* fetch(Index_ i, Value_* buffer) {
            const typename Slab::value_type* ptr;
            size_t len = extracted_length(this);

            if (workspace.oracle_cache) {
                auto predicted = oracle_cache.next_chunk(
                    /* identify = */ [&](Index_ i) -> std::pair<Index_, Index_> {
                        return std::make_pair(i / get_primary_chunkdim(), i % get_primary_chunkdim());
                    },
                    /* swap = */ [&](Slab& left, Slab& right) -> void {
                        left.swap(right);
                    },
                    /* ready = */ [&](const Slab& slab) -> bool {
                        return !cache.values.empty();
                    },
                    /* allocate = */ [&](const Slab& slab) -> void {
                        cache.values.resize(len * get_primary_chunkdim());
                    },
                    /* populate =*/ [&](const std::vector<std::pair<Index_, Index_> >& chunks_in_need, std::vector<Slab>& chunk_data) -> void {
                        for (const auto& p : chunks_in_need) {
                            extract<false>(p.first, chunk_offset, chunk_data[p.second]);
                        }
                    }
                );
                ptr = predicted.first.values.data() + len * predicted.second;

            } else {
                auto chunk_id = i / get_primary_chunkdim();
                auto chunk_offset = i % get_primary_chunkdim();

                if (workspace.num_chunk_sets_in_cache) {
                    extract<true>(chunk_id, chunk_offset, *solo);
                    ptr = solo->values.data();

                } else {
                    auto& cache = lru_cache.find_chunk(
                        chunk_id,
                        /* create = */ [&]() -> Slab {
                            return Slab(workspace.chunk_set_size_in_elements);
                        },
                        /* populate = */ [&](Index_ id, Slab& slab) -> void {
                            extract<false>(id, chunk_offset, slab);
                        }
                    );
                    ptr = cache.values.data() + len * chunk_offset;
                }
            }

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

template<typename Value_, typename Index_, typename Chunk_>
class CustomChunkedSparseMatrix : public Matrix<Value_, Index_> {
public:
    CustomChunkedMatrix(Index_ mat_nrow, Index_ mat_ncol, Index_ chunk_nrow, Index_ chunk_ncol, Index_ nrow_in_chunks, Index_ ncol_in_chunks, std::vector<Chunk_> chunks, bool row_major, const CustomChunkedOptions& opt) : 
        manager(std::move(chunks)), 
        nrows(mat_nrow), 
        ncols(mat_ncol), 
        cache_size_in_elements(opt.maximum_cache_size / sizeof(Chunk_::value_type)),
        require_minimum_cache(opt.require_minimum_cache)
    {
        manager.chunk_nrow = chunk_nrow;
        manager.chunk_ncol = chunk_ncol;
        manager.num_chunks_per_column = nrow_in_chunks;
        manager.num_chunks_per_row = ncol_in_chunks;
        manager.row_major = row_major;
    }

private:
    CustomChunkManager<Chunk_> manager;
    Index_ nrows, ncols;
    size_t cache_size_in_elements;
    bool require_minimum_cache;

    static_assert(Chunk_::sparse); // Should be sparse, obviously.


public:
    Index_ nrow() const { return nrows; }

    Index_ ncol() const { return ncols; }

    bool prefer_rows() const { 
        // Prefer rows if we have to extract fewer chunks per row.
        return manager.num_chunks_per_column > manager.num_chunks_per_row; 
    }

    bool uses_oracle(bool) const { 
        return true; 
    }

    double prefer_rows_proportion() const { 
        return static_cast<double>(row_); 
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
                this->full_length = (accrow_ ? p->ncols : p->nrows);
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

        typedef typename CustomChunkManager::Cache Slab;
        typename CustomChunkManager::Workspace buffer;
        TypicalChunkCacheWorkspace<Index_, Slab> workspace;
        std::unique_ptr<Slab> solo;

        void initialize_cache() {
            workspace = TypicalChunkCacheWorkspace<Index_, Slab>(
                accrow_ ? parent->manager->chunk_nrow, parent->manager->chunk_ncol,
                extracted_length(this)
                parent->cache_size_in_elements,
                parent->require_minimum_cache
            );

            if (workspace.num_chunk_sets_in_cache == 0) {
                solo.reset(new Slab(parent->manager.create_slab<accrow_, true>(len)));
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
            workspace.set_oracle(std::move(o));
            return;
        }

    private:
        Index_ get_primary_dim() const {
            if constexpr(accrow_) {
                return parent->nrows;
            } else {
                return parent->ncols;
            }
        }

        Index_ get_primary_chunkdim() const {
            if constexpr(accrow_) {
                return parent->manager.chunk_nrow;
            } else {
                return parent->manager.chunk_ncol;
            }
        }

        template<bool exact_>
        void extract(Index_ chunk_id, Index_ chunk_offset, Slab& slab) const {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                manager.extract<accrow_, exact_>(chunk_id, chunk_offset, get_primary_dim(), static_cast<Index_>(0), this->full_length, slab, buffer);
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                manager.extract<accrow_, exact_>(chunk_id, chunk_offset, get_primary_dim(), this->block_start, this->block_length, slab, buffer);
            } else {
                manager.extract<accrow_, exact_>(chunk_id, chunk_offset, get_primary_dim(), indices, slab, buffer);
            }
        }

    public:
        std::pair<const Slab*, Index_> fetch_cache(Index_ i) {
            if (workspace.oracle_cache) {
                auto predicted = oracle_cache.next_chunk(
                    /* identify = */ [&](Index_ i) -> std::pair<Index_, Index_> {
                        return std::make_pair(i / get_primary_chunkdim(), i % get_primary_chunkdim());
                    },
                    /* swap = */ [&](Slab& left, Slab& right) -> void {
                        left.swap(right);
                    },
                    /* ready = */ [&](const Slab& slab) -> bool {
                        return !cache.cache.empty();
                    },
                    /* allocate = */ [&](const Slab& slab) -> void {
                        cache.cache.resize(len * get_primary_chunkdim());
                    },
                    /* populate =*/ [&](const std::vector<std::pair<Index_, Index_> >& chunks_in_need, std::vector<Slab>& chunk_data) -> void {
                        for (const auto& p : chunks_in_need) {
                            extract<false>(p.first, chunk_offset, chunk_data[p.second]);
                        }
                    }
                );
                return predicted;

            } else {
                auto chunk_id = i / get_primary_chunkdim();
                auto chunk_offset = i % get_primary_chunkdim();

                if (workspace.num_chunk_sets_in_cache) {
                    extract<true>(chunk_id, chunk_offset, *solo);
                    return std::make_pair(solo.get(), 0);

                } else {
                    auto& cache = lru_cache.find_chunk(
                        chunk_id,
                        /* create = */ [&]() -> Slab {
                            return Slab(workspace.chunk_set_size_in_elements);
                        },
                        /* populate = */ [&](Index_ id, Slab& slab) -> void {
                            extract<false>(id, chunk_offset, slab);
                        }
                    );
                    return std::make_pair(&cache, chunk_offset);
                }
            }
        }
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
            auto contents = this->fetch_cache(i);
            const auto& values = contents.first->values[contents.second];
            const auto& indices = contents.first->indices[contents.second];

            auto len = extracted_length(this);
            std::fill(buffer, buffer + len, 0);

            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                for (size_t i = 0; end = indices.size(); ++i) {
                    buffer[remap[indices[i]]] = values[i];
                }
            } else {
                Index_ start = (selection_ == DimensionSelectionType::BLOCK ? this->block_start : 0);
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
        CustomDenseExtractor(const CustomChunkedSparseMatrix* p, bool ev, bool ei, Args_&& ... args) : 
            CustomExtractorBase<accrow_, true, selection_>(p, std::forward<Args_>(args)...), 
            extract_value(ev), 
            extract_index(ei)
        {}

    public:
        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            auto contents = this->fetch_cache(i);
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

            return SparseRange<Value, Index_>(values.size(), vbuffer, ibuffer);
        }

    private:
        bool extract_value, extract_index;
    };

public:
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        auto ptr = new CustomSparseExtractor<true, DimensionSelectionType::FULL>(this);
        return std::unique_ptr<FullSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        auto ptr = new CustomSparseExtractor<true, DimensionSelectionType::BLOCK>(this, block_start, block_length);
        return std::unique_ptr<BlockSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        auto ptr = new CustomSparseExtractor<true, DimensionSelectionType::INDEX>(this, std::move(indices));
        return std::unique_ptr<IndexSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<FullSparseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        auto ptr = new CustomSparseExtractor<false, DimensionSelectionType::FULL>(this);
        return std::unique_ptr<FullSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        auto ptr = new CustomSparseExtractor<false, DimensionSelectionType::BLOCK>(this, block_start, block_length);
        return std::unique_ptr<BlockSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        auto ptr = new CustomSparseExtractor<false, DimensionSelectionType::INDEX>(this, std::move(indices));
        return std::unique_ptr<IndexSparseExtractor<Value_, Index_> >(ptr);
    }

};


}

#endif
