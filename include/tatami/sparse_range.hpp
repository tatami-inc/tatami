#ifndef TATAMI_SPARSE_RANGE_H
#define TATAMI_SPARSE_RANGE_H

namespace tatami {

template <typename T, typename IDX>
struct sparse_range {
    sparse_range(size_t n, const T* v=NULL, const IDX* i=NULL) : number(n), value(v), index(i) {}
    sparse_range() {}
    size_t number = 0;
    const T* value = NULL;
    const IDX* index = NULL;
};

}

#endif
