#ifndef TATAMI_TEST_FORCED_ORACLE_HPP
#define TATAMI_TEST_FORCED_ORACLE_HPP

#include "../tatami/base/Matrix.hpp"
#include "../tatami/utils/copy.hpp"

#include <algorithm>
#include <memory>

namespace tatami_test {

template<typename Value_, typename Index_>
struct TestWrapper : public tatami::Matrix<Value_, Index_> {
    TestWrapper(std::shared_ptr<const tatami::Matrix<Value_, Index_> > p) : mat(std::move(p)) {}

protected:
    std::shared_ptr<const tatami::Matrix<Value_, Index_> > mat;

public:
    Index_ nrow() const {
        return mat->nrow();
    }

    Index_ ncol() const {
        return mat->ncol();
    }

    bool is_sparse() const {
        return mat->is_sparse();
    }

    double is_sparse_proportion() const {
        return mat->is_sparse_proportion();
    }

    bool prefer_rows() const {
        return mat->prefer_rows();
    }

    double prefer_rows_proportion() const {
        return mat->prefer_rows_proportion();
    }

    bool uses_oracle(bool) const {
        return false;
    }

public:
    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense(bool row, const tatami::Options& opt) const { 
        return mat->dense(row, opt); 
    }

    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense(bool row, Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return mat->dense(row, bs, bl, opt);
    }

    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense(bool row, tatami::VectorPtr<Index_> idx, const tatami::Options& opt) const {
        return mat->dense(row, std::move(idx), opt);
    }

public:
    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const tatami::Options& opt) const { 
        return mat->sparse(row, opt); 
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return mat->sparse(row, bs, bl, opt);
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, tatami::VectorPtr<Index_> idx, const tatami::Options& opt) const {
        return mat->sparse(row, std::move(idx), opt);
    }

public:
    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const tatami::Oracle<Index_> > ora, const tatami::Options& opt) const { 
        return mat->dense(row, std::move(ora), opt); 
    }

    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const tatami::Oracle<Index_> > ora, Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return mat->dense(row, std::move(ora), bs, bl, opt);
    }

    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const tatami::Oracle<Index_> > ora, tatami::VectorPtr<Index_> idx, const tatami::Options& opt) const {
        return mat->dense(row, std::move(ora), std::move(idx), opt);
    }

public:
    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const tatami::Oracle<Index_> > ora, const tatami::Options& opt) const { 
        return mat->sparse(row, std::move(ora), opt); 
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const tatami::Oracle<Index_> > ora, Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return mat->sparse(row, std::move(ora), bs, bl, opt);
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const tatami::Oracle<Index_> > ora, tatami::VectorPtr<Index_> idx, const tatami::Options& opt) const {
        return mat->sparse(row, std::move(ora), std::move(idx), opt);
    }
};

template<typename Value_, typename Index_>
struct ForcedOracleWrapper : public TestWrapper<Value_, Index_> {
    ForcedOracleWrapper(std::shared_ptr<const tatami::Matrix<Value_, Index_> > p) : TestWrapper<Value_, Index_>(std::move(p)) {}

    using tatami::Matrix<Value_, Index_>::uses_oracle;

    bool uses_oracle(bool) const {
        return true;
    }
};

template<typename Value_, typename Index_>
struct UnsortedWrapper : public TestWrapper<Value_, Index_> {
    UnsortedWrapper(std::shared_ptr<const tatami::Matrix<Value_, Index_> > p) : TestWrapper<Value_, Index_>(std::move(p)) {}

private:
    template<bool oracle_>
    struct UnsortedExtractor : public tatami::SparseExtractor<oracle_, Value_, Index_> {
        UnsortedExtractor(std::unique_ptr<tatami::SparseExtractor<oracle_, Value_, Index_> > h, bool srtd) : host(std::move(h)), must_sort(srtd) {}

    private:
        std::unique_ptr<tatami::SparseExtractor<oracle_, Value_, Index_> > host;
        bool must_sort;

    public:
        tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            auto range = host->fetch(i, vbuffer, ibuffer);
            if (!must_sort) {
                if (range.value) {
                    tatami::copy_n(range.value, range.number, vbuffer);
                    std::reverse(vbuffer, vbuffer + range.number);
                    range.value = vbuffer;
                }
                if (range.index) {
                    tatami::copy_n(range.index, range.number, ibuffer);
                    std::reverse(ibuffer, ibuffer + range.number);
                    range.index = ibuffer;
                }
            }
            return range;
        }
    };

public:
    using tatami::Matrix<Value_, Index_>::sparse;
    using TestWrapper<Value_, Index_>::mat;

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const tatami::Options& opt) const { 
        return std::make_unique<UnsortedExtractor<false> >(mat->sparse(row, opt), opt.sparse_ordered_index); 
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return std::make_unique<UnsortedExtractor<false> >(mat->sparse(row, bs, bl, opt), opt.sparse_ordered_index);
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, tatami::VectorPtr<Index_> idx, const tatami::Options& opt) const {
        return std::make_unique<UnsortedExtractor<false> >(mat->sparse(row, std::move(idx), opt), opt.sparse_ordered_index);
    }


public:
    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const tatami::Oracle<Index_> > ora, const tatami::Options& opt) const { 
        return std::make_unique<UnsortedExtractor<true> >(mat->sparse(row, std::move(ora), opt), opt.sparse_ordered_index);
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const tatami::Oracle<Index_> > ora, Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return std::make_unique<UnsortedExtractor<true> >(mat->sparse(row, std::move(ora), bs, bl, opt), opt.sparse_ordered_index);
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const tatami::Oracle<Index_> > ora, tatami::VectorPtr<Index_> idx, const tatami::Options& opt) const {
        return std::make_unique<UnsortedExtractor<true> >(mat->sparse(row, std::move(ora), std::move(idx), opt), opt.sparse_ordered_index);
    }
};

}

#endif
