#ifndef TATAMI_H
#define TATAMI_H

#include "base/DenseMatrix.hpp"
#include "base/CompressedSparseMatrix.hpp"
#include "base/DelayedIsometricOp.hpp"
#include "base/DelayedSubset.hpp"
#include "base/DelayedSubsetBlock.hpp"
#include "base/DelayedBind.hpp"
#include "base/DelayedTranspose.hpp"

#include "utils/compress_sparse_triplets.hpp"
#include "utils/convert_to_sparse.hpp"
#include "utils/convert_to_dense.hpp"
#include "utils/wrap_shared_ptr.hpp"
#include "utils/NakedArray.hpp"
#include "utils/bind_intersection.hpp"

#include "stats/sums.hpp"
#include "stats/variances.hpp"
#include "stats/medians.hpp"

#define TATAMI_VERSION_MAJOR 0
#define TATAMI_VERSION_MINOR 99
#define TATAMI_VERSION_PATCH 0

#endif
