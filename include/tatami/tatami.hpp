#ifndef TATAMI_TATAMI_HPP
#define TATAMI_TATAMI_HPP

#include "base/DenseMatrix.hpp"
#include "base/CompressedSparseMatrix.hpp"
#include "base/DelayedIsometricOp.hpp"
#include "base/DelayedSubset.hpp"
#include "base/DelayedSubsetBlock.hpp"
#include "base/DelayedBind.hpp"
#include "base/DelayedCast.hpp"
#include "base/DelayedTranspose.hpp"

#include "utils/compress_sparse_triplets.hpp"
#include "utils/convert_to_sparse.hpp"
#include "utils/convert_to_dense.hpp"
#include "utils/wrap_shared_ptr.hpp"
#include "utils/ArrayView.hpp"
#include "utils/SomeNumericArray.hpp"
#include "utils/bind_intersection.hpp"

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
