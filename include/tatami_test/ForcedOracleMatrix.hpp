#ifndef TATAMI_TEST_FORCED_ORACLE_HPP
#define TATAMI_TEST_FORCED_ORACLE_HPP

#include "../tatami/base/Matrix.hpp"

namespace tatami_test {

enum TestAccessOracle {
    NO_ORACLE,
    DEFAULT_ORACLE,
    FORCED_ORACLE
};

template<typename Value_, typename Index_>
struct ForcedOracleMatrix : public tatami::Matrix<Value_, Index_> {
    ForcedOracleMatrix(std::shared_ptr<const tatami::Matrix<Value_, Index_> > p) : mat(std::move(p)) {}

private:
    std::shared_ptr<const tatami::Matrix<Value_, Index_> > mat;

public:
    Index_ nrow() const {
        return mat->nrow();
    }

    Index_ ncol() const {
        return mat->ncol();
    }

    bool sparse() const {
        return mat->sparse();
    }

    double sparse_proportion() const {
        return mat->sparse_proportion();
    }

    bool prefer_rows() const {
        return mat->prefer_rows();
    }

    double prefer_rows_proportion() const {
        return mat->prefer_rows_proportion();
    }

    bool uses_oracle(bool) const {
        return true;
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

}

#endif
