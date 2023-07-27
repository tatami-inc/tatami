#include "../sparse/convert_to_compressed_sparse.hpp"

namespace tatami {

// Present only for back-compatibility.
template <bool row_, 
    typename Value_ = double,
    typename Index_ = int,
    typename StoredValue_ = Value_,
    typename StoredIndex_ = Index_,
    typename InputValue_,
    typename InputIndex_
>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_sparse(const Matrix<InputValue_, InputIndex_>* incoming, int threads = 1) {
    if (threads == 0) {
        threads = 1; // for back-compatibility with calls using the old 'reserve' argument.
    }
    return convert_to_compressed_sparse<row_, Value_, Index_, StoredValue_, StoredIndex_>(incoming, false, threads);
}

template <
    typename Value_ = double,
    typename Index_ = int,
    typename StoredValue_ = Value_,
    typename StoredIndex_ = Index_,
    typename InputValue_,
    typename InputIndex_
>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_sparse(const Matrix<InputValue_, InputIndex_>* incoming, int order, int threads = 1) {
    return convert_to_compressed_sparse<Value_, Index_, StoredValue_, StoredIndex_>(incoming, order, threads);
}

}
