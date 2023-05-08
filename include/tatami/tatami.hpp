#ifndef TATAMI_TATAMI_HPP
#define TATAMI_TATAMI_HPP

#include "base/dense/DenseMatrix.hpp"

#include "base/sparse/CompressedSparseMatrix.hpp"
#include "base/sparse/SemiCompressedSparseMatrix.hpp"

#include "base/isometric/DelayedIsometricOp.hpp"

#include "base/other/DelayedBind.hpp"
#include "base/other/DelayedCast.hpp"
#include "base/other/DelayedTranspose.hpp"

#include "base/subset/DelayedSubsetBlock.hpp"
#include "base/subset/make_DelayedSubset.hpp"

#include "utils/compress_sparse_triplets.hpp"
#include "utils/convert_to_sparse.hpp"
#include "utils/convert_to_dense.hpp"
#include "utils/wrap_shared_ptr.hpp"
#include "utils/ArrayView.hpp"
#include "utils/SomeNumericArray.hpp"
#include "utils/bind_intersection.hpp"
#include "utils/Oracles.hpp"

#include "stats/sums.hpp"
#include "stats/variances.hpp"
#include "stats/medians.hpp"

#define TATAMI_VERSION_MAJOR 1
#define TATAMI_VERSION_MINOR 0
#define TATAMI_VERSION_PATCH 0

/**
 * @file tatami.hpp
 * @brief Flexible representations of matrix data
 */

/**
 * @namespace tatami
 * @brief Flexible representations for matrix data
 */
namespace tatami {}

#endif
