#ifndef TATAMI_TEST_ORACLE_ACCESS_HPP
#define TATAMI_TEST_ORACLE_ACCESS_HPP

#include <gtest/gtest.h>

#include "../tatami/base/Matrix.hpp"
#include "../tatami/utils/Oracles.hpp"
#include <random>
#include <deque>

namespace tatami_test {

template<typename Value_, typename Index_> 
class CrankyMatrix : public tatami::Matrix<Value_, Index_> {
public:
    CrankyMatrix(std::shared_ptr<tatami::Matrix<Value_, Index_> > x, int pt) : matrix(std::move(x)), predict_to(pt) {}

private:
    std::shared_ptr<tatami::Matrix<Value_, Index_> > matrix;
    int predict_to;

public:
    Index_ nrow() const { return matrix->nrow(); }

    Index_ ncol() const { return matrix->ncol(); }

    bool sparse() const { return matrix->sparse(); }

    double sparse_proportion() const { return matrix->sparse_proportion(); }

    bool prefer_rows() const { return matrix->prefer_rows(); }

    double prefer_rows_proportion() const { return matrix->prefer_rows_proportion(); }

    bool uses_oracle(bool) const { return true; }

public:
    template<tatami::DimensionSelectionType selection_, bool sparse_>
    struct CrankyExtractor : public tatami::Extractor<selection_, sparse_, Value_, Index_> {
        CrankyExtractor(std::unique_ptr<tatami::Extractor<selection_, sparse_, Value_, Index_> > inner, bool u, int pt) : internal(std::move(inner)), unused(u), predict_to(pt) {
            if constexpr(selection_ == tatami::DimensionSelectionType::FULL) {
                this->full_length = internal->full_length;
            } else if constexpr(selection_ == tatami::DimensionSelectionType::BLOCK) {
                this->block_start = internal->block_start;
                this->block_length = internal->block_length;
            } else {
                this->index_length = internal->index_length;
            }
        }

    protected:
        std::unique_ptr<tatami::Extractor<selection_, sparse_, Value_, Index_> > internal;
        std::unique_ptr<tatami::Oracle<Index_> > source;
        std::deque<Index_> filled;
        std::vector<Index_> buffer;
        bool unused;
        int predict_to;

    public:
        const Index_* index_start() const {
            return internal->index_start();
        }

    private:
        struct CrankyOracle : public tatami::Oracle<Index_> {
            CrankyOracle(CrankyExtractor* p) : parent(p) {}

            size_t predict(Index_* buffer, size_t number) {
                size_t n = parent->source->predict(buffer, number);
                parent->filled.insert(parent->filled.end(), buffer, buffer + n);
                return n;
            }
        private:
            CrankyExtractor* parent;
        };

    protected:
        void check(Index_ i) {
            if (unused && filled.empty()) {
                buffer.resize(predict_to);
                size_t n = source->predict(buffer.data(), buffer.size());
                filled.insert(filled.end(), buffer.begin(), buffer.begin() + n);
            }

            if (filled.empty() || i != filled.front()) { // not happy if the predictions don't match up!
                throw std::runtime_error("mismatch in the predictions versus the iteration order");
            } else {
                filled.pop_front();
            }
        }

    public:
        void set_oracle(std::unique_ptr<tatami::Oracle<Index_> > o) {
            source = std::move(o);
            internal->set_oracle(std::make_unique<CrankyOracle>(this));
        }
    };

    template<tatami::DimensionSelectionType selection_>
    struct DenseCrankyExtractor : public CrankyExtractor<selection_, false> {
        DenseCrankyExtractor(std::unique_ptr<tatami::Extractor<selection_, false, Value_, Index_> > inner, bool unused, int predict_to) : 
            CrankyExtractor<selection_, false>(std::move(inner), unused, predict_to) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto out = this->internal->fetch(i, buffer);
            this->check(i);
            return out;
        }
    };

    template<tatami::DimensionSelectionType selection_>
    struct SparseCrankyExtractor : public CrankyExtractor<selection_, true> {
        SparseCrankyExtractor(std::unique_ptr<tatami::Extractor<selection_, true, Value_, Index_> > inner, bool unused, int predict_to) : 
            CrankyExtractor<selection_, true>(std::move(inner), unused, predict_to) {}

        tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            auto out = this->internal->fetch(i, vbuffer, ibuffer);
            this->check(i);
            return out;
        }
    };

public:
    std::unique_ptr<tatami::FullDenseExtractor<Value_, Index_> > dense_row(const tatami::Options& opt) const {
        return std::make_unique<DenseCrankyExtractor<tatami::DimensionSelectionType::FULL> >(matrix->dense_row(opt), !matrix->uses_oracle(true), predict_to);
    }

    std::unique_ptr<tatami::BlockDenseExtractor<Value_, Index_> > dense_row(Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return std::make_unique<DenseCrankyExtractor<tatami::DimensionSelectionType::BLOCK> >(matrix->dense_row(bs, bl, opt), !matrix->uses_oracle(true),  predict_to);
    }

    std::unique_ptr<tatami::IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> i, const tatami::Options& opt) const {
        return std::make_unique<DenseCrankyExtractor<tatami::DimensionSelectionType::INDEX> >(matrix->dense_row(std::move(i), opt), !matrix->uses_oracle(true), predict_to);
    }

    std::unique_ptr<tatami::FullDenseExtractor<Value_, Index_> > dense_column(const tatami::Options& opt) const {
        return std::make_unique<DenseCrankyExtractor<tatami::DimensionSelectionType::FULL> >(matrix->dense_column(opt), !matrix->uses_oracle(false), predict_to);
    }

    std::unique_ptr<tatami::BlockDenseExtractor<Value_, Index_> > dense_column(Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return std::make_unique<DenseCrankyExtractor<tatami::DimensionSelectionType::BLOCK> >(matrix->dense_column(bs, bl, opt), !matrix->uses_oracle(false), predict_to);
    }

    std::unique_ptr<tatami::IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> i, const tatami::Options& opt) const {
        return std::make_unique<DenseCrankyExtractor<tatami::DimensionSelectionType::INDEX> >(matrix->dense_column(std::move(i), opt), !matrix->uses_oracle(false), predict_to);
    }

public:
    std::unique_ptr<tatami::FullSparseExtractor<Value_, Index_> > sparse_row(const tatami::Options& opt) const {
        return std::make_unique<SparseCrankyExtractor<tatami::DimensionSelectionType::FULL> >(matrix->sparse_row(opt), !matrix->uses_oracle(true), predict_to);
    }

    std::unique_ptr<tatami::BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return std::make_unique<SparseCrankyExtractor<tatami::DimensionSelectionType::BLOCK> >(matrix->sparse_row(bs, bl, opt), !matrix->uses_oracle(true), predict_to);
    }

    std::unique_ptr<tatami::IndexSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> i, const tatami::Options& opt) const {
        return std::make_unique<SparseCrankyExtractor<tatami::DimensionSelectionType::INDEX> >(matrix->sparse_row(std::move(i), opt), !matrix->uses_oracle(true), predict_to);
    }

    std::unique_ptr<tatami::FullSparseExtractor<Value_, Index_> > sparse_column(const tatami::Options& opt) const {
        return std::make_unique<SparseCrankyExtractor<tatami::DimensionSelectionType::FULL> >(matrix->sparse_column(opt), !matrix->uses_oracle(false), predict_to);
    }

    std::unique_ptr<tatami::BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return std::make_unique<SparseCrankyExtractor<tatami::DimensionSelectionType::BLOCK> >(matrix->sparse_column(bs, bl, opt), !matrix->uses_oracle(false), predict_to);
    }

    std::unique_ptr<tatami::IndexSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> i, const tatami::Options& opt) const {
        return std::make_unique<SparseCrankyExtractor<tatami::DimensionSelectionType::INDEX> >(matrix->sparse_column(std::move(i), opt), !matrix->uses_oracle(false), predict_to);
    }
};

template<typename Value_, typename Index_> 
std::shared_ptr<tatami::Matrix<Value_, Index_> > make_CrankyMatrix(std::shared_ptr<tatami::Matrix<Value_, Index_> > p, int predict_to = 10) {
    return std::shared_ptr<tatami::Matrix<Value_, Index_> >(new CrankyMatrix<Value_, Index_>(std::move(p), predict_to));
}

template<class Matrix, typename ... Args_>
void test_oracle_column_access(const Matrix* ptr, const Matrix* ref, bool randomized, Args_... args) {
    int NR = ptr->nrow();
    int NC = ptr->ncol();

    auto pwork = ref->dense_column(args...);
    auto swork = ref->sparse_column(args...);

    auto pwork_o = ptr->dense_column(args...);
    auto swork_o = ptr->sparse_column(args...);

    auto iterator = [&](int i) -> void {
        auto expected = pwork->fetch(i);
        auto observed = pwork_o->fetch(i);
        EXPECT_EQ(expected, observed);

        auto sexpected = swork->fetch(i);
        auto sobserved = swork_o->fetch(i);
        EXPECT_EQ(sexpected.index, sobserved.index);
        EXPECT_EQ(sexpected.value, sobserved.value);
    };

    typedef typename Matrix::index_type Index_;
    if (randomized) {
        std::mt19937_64 rng(NR + NC * 10); // making up an interesting seed.
        std::vector<Index_> fixed(NC * 2);
        for (auto& x : fixed) {
            x = rng() % NC;
        }

        pwork_o->set_oracle(std::make_unique<tatami::FixedOracle<Index_> >(fixed.data(), fixed.size()));
        swork_o->set_oracle(std::make_unique<tatami::FixedOracle<Index_> >(fixed.data(), fixed.size()));

        for (auto i : fixed) {
            iterator(i);
        }

    } else {
        pwork_o->set_oracle(std::make_unique<tatami::ConsecutiveOracle<Index_> >(0, NC));
        swork_o->set_oracle(std::make_unique<tatami::ConsecutiveOracle<Index_> >(0, NC));

        for (int i = 0; i < NC; ++i) {
            iterator(i);
        }
    }
}

template<class Matrix, typename ... Args_>
void test_oracle_row_access(const Matrix* ptr, const Matrix* ref, bool randomized, Args_... args) {
    int NR = ptr->nrow();
    int NC = ptr->ncol();

    auto pwork = ref->dense_row(args...);
    auto swork = ref->sparse_row(args...);

    auto pwork_o = ptr->dense_row(args...);
    auto swork_o = ptr->sparse_row(args...);

    auto iterator = [&](int i) -> void {
        auto expected = pwork->fetch(i);
        auto observed = pwork_o->fetch(i);
        EXPECT_EQ(expected, observed);

        auto sexpected = swork->fetch(i);
        auto sobserved = swork_o->fetch(i);
        EXPECT_EQ(sexpected.index, sobserved.index);
        EXPECT_EQ(sexpected.value, sobserved.value);
    };

    typedef typename Matrix::index_type Index_;
    if (randomized) {
        std::mt19937_64 rng(3 * NR + NC); // making an interesting seed again.
        std::vector<Index_> fixed(NR * 2);
        for (auto& x : fixed) {
            x = rng() % NR;
        }

        pwork_o->set_oracle(std::make_unique<tatami::FixedOracle<Index_> >(fixed.data(), fixed.size()));
        swork_o->set_oracle(std::make_unique<tatami::FixedOracle<Index_> >(fixed.data(), fixed.size()));

        for (auto i : fixed) {
            iterator(i);
        }

    } else {
        pwork_o->set_oracle(std::make_unique<tatami::ConsecutiveOracle<Index_> >(0, NR));
        swork_o->set_oracle(std::make_unique<tatami::ConsecutiveOracle<Index_> >(0, NR));

        for (int i = 0; i < NR; ++i) {
            iterator(i);
        }
    }
}

}

#endif
