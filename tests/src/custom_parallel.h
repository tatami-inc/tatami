#ifndef CUSTOM_PARALLEL_H
#define CUSTOM_PARALLEL_H

// Put this before any tatami imports.
#ifdef CUSTOM_PARALLEL_TEST
#include "subpar/subpar.hpp"
#define TATAMI_CUSTOM_PARALLEL(fun, ntasks, nthreads) ::subpar::test_parallelize_range(nthreads, ntasks, std::move(fun))
#endif

#endif
