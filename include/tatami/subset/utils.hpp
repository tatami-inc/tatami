#ifndef TATAMI_DELAYED_SUBSET_UTILS_HPP
#define TATAMI_DELAYED_SUBSET_UTILS_HPP

#include "../base/Matrix.hpp"
#include "../utils/new_extractor.hpp"

#include <vector>
#include <algorithm>
#include <type_traits>

namespace tatami {

namespace subset_utils {

template<typename Index_, class SubsetStorage_>
class SubsetOracle final : public Oracle<Index_> {
public:
    SubsetOracle(std::shared_ptr<const Oracle<Index_> > oracle, const SubsetStorage_& subset) : my_oracle(std::move(oracle)), my_subset(subset) {}

    Index_ get(PredictionIndex i) const {
        return my_subset[my_oracle->get(i)];
    }

    PredictionIndex total() const {
        return my_oracle->total();
    }

private:
    std::shared_ptr<const Oracle<Index_> > my_oracle;
    const SubsetStorage_& my_subset;
};

template<typename Value_, typename Index_, class SubsetStorage_>
class MyopicPerpendicularDense final : public MyopicDenseExtractor<Value_, Index_> {
public:
    template<typename ... Args_>
    MyopicPerpendicularDense(const Matrix<Value_, Index_>* matrix, const SubsetStorage_& subset, bool row, Args_&& ... args) : 
        my_subset(subset), my_ext(new_extractor<false, false>(matrix, row, false, std::forward<Args_>(args)...)) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        return my_ext->fetch(my_subset[i], buffer);
    }

protected:
    const SubsetStorage_& my_subset;
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > my_ext;
};

template<typename Value_, typename Index_, class SubsetStorage_>
class MyopicPerpendicularSparse final : public MyopicSparseExtractor<Value_, Index_> {
public:
    template<typename ... Args_>
    MyopicPerpendicularSparse(const Matrix<Value_, Index_>* matrix, const SubsetStorage_& subset, bool row, Args_&& ... args) : 
        my_subset(subset), my_ext(new_extractor<true, false>(matrix, row, false, std::forward<Args_>(args)...)) {}

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        return my_ext->fetch(my_subset[i], value_buffer, index_buffer);
    }

protected:
    const SubsetStorage_& my_subset;
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > my_ext;
};

template<typename Value_, typename Index_>
class OracularPerpendicularDense final : public OracularDenseExtractor<Value_, Index_> {
public:
    template<class SubsetStorage_, typename ... Args_>
    OracularPerpendicularDense(const Matrix<Value_, Index_>* matrix, const SubsetStorage_& subset, bool row, std::shared_ptr<const Oracle<Index_> > oracle, Args_&& ... args) :
        my_ext(new_extractor<false, true>(matrix, row, std::make_shared<SubsetOracle<Index_, SubsetStorage_> >(std::move(oracle), subset), std::forward<Args_>(args)...)) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        return my_ext->fetch(i, buffer);
    }

protected:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > my_ext;
};

template<typename Value_, typename Index_>
class OracularPerpendicularSparse final : public OracularSparseExtractor<Value_, Index_> {
public:
    template<class SubsetStorage_, typename ... Args_>
    OracularPerpendicularSparse(const Matrix<Value_, Index_>* matrix, const SubsetStorage_& subset, bool row, std::shared_ptr<const Oracle<Index_> > oracle, Args_&& ... args) :
        my_ext(new_extractor<true, true>(matrix, row, std::make_shared<SubsetOracle<Index_, SubsetStorage_> >(std::move(oracle), subset), std::forward<Args_>(args)...)) {}

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        return my_ext->fetch(i, value_buffer, index_buffer);
    }

protected:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > my_ext;
};

}

}

#endif
