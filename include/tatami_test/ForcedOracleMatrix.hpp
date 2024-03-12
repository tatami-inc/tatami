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
    ForcedOracleMatrix(std::shared_ptr<tatami::Matrix<Value_, Index_> > p) : mat(std::move(p)) {}

private:
    std::shared_ptr<tatami::Matrix<Value_, Index_> > mat;

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
    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense_row(const tatami::Options& opt) const { 
        return mat->dense_row(opt); 
    }

    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense_row(Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return mat->dense_row(bs, bl, opt);
    }

    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> idx, const tatami::Options& opt) const {
        return mat->dense_row(std::move(idx), opt);
    }

    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense_column(const tatami::Options& opt) const { 
        return mat->dense_column(opt); 
    }

    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense_column(Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return mat->dense_column(bs, bl, opt);
    }

    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> idx, const tatami::Options& opt) const {
        return mat->dense_column(std::move(idx), opt);
    }

public:
    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse_row(const tatami::Options& opt) const { 
        return mat->sparse_row(opt); 
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse_row(Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return mat->sparse_row(bs, bl, opt);
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> idx, const tatami::Options& opt) const {
        return mat->sparse_row(std::move(idx), opt);
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse_column(const tatami::Options& opt) const { 
        return mat->sparse_column(opt); 
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse_column(Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return mat->sparse_column(bs, bl, opt);
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> idx, const tatami::Options& opt) const {
        return mat->sparse_column(std::move(idx), opt);
    }

public:
    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<tatami::Oracle<Index_> > ora, const tatami::Options& opt) const { 
        return mat->dense_row(std::move(ora), opt); 
    }

    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<tatami::Oracle<Index_> > ora, Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return mat->dense_row(std::move(ora), bs, bl, opt);
    }

    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<tatami::Oracle<Index_> > ora, std::vector<Index_> idx, const tatami::Options& opt) const {
        return mat->dense_row(std::move(ora), std::move(idx), opt);
    }

    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<tatami::Oracle<Index_> > ora, const tatami::Options& opt) const { 
        return mat->dense_column(std::move(ora), opt); 
    }

    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<tatami::Oracle<Index_> > ora, Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return mat->dense_column(std::move(ora), bs, bl, opt);
    }

    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<tatami::Oracle<Index_> > ora, std::vector<Index_> idx, const tatami::Options& opt) const {
        return mat->dense_column(std::move(ora), std::move(idx), opt);
    }

public:
    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<tatami::Oracle<Index_> > ora, const tatami::Options& opt) const { 
        return mat->sparse_row(std::move(ora), opt); 
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<tatami::Oracle<Index_> > ora, Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return mat->sparse_row(std::move(ora), bs, bl, opt);
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<tatami::Oracle<Index_> > ora, std::vector<Index_> idx, const tatami::Options& opt) const {
        return mat->sparse_row(std::move(ora), std::move(idx), opt);
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<tatami::Oracle<Index_> > ora, const tatami::Options& opt) const { 
        return mat->sparse_column(std::move(ora), opt); 
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<tatami::Oracle<Index_> > ora, Index_ bs, Index_ bl, const tatami::Options& opt) const {
        return mat->sparse_column(std::move(ora), bs, bl, opt);
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<tatami::Oracle<Index_> > ora, std::vector<Index_> idx, const tatami::Options& opt) const {
        return mat->sparse_column(std::move(ora), std::move(idx), opt);
    }
};

}

#endif
