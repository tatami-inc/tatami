#ifndef TATAMI_TATAMI_HPP
#define TATAMI_TATAMI_HPP

#include "dense/DenseMatrix.hpp"

#include "sparse/CompressedSparseMatrix.hpp"
#include "sparse/SemiCompressedSparseMatrix.hpp"
#include "sparse/FragmentedSparseMatrix.hpp"

#include "isometric/unary/DelayedUnaryIsometricOp.hpp"
#include "isometric/binary/DelayedBinaryIsometricOp.hpp"

#include "other/DelayedBind.hpp"
#include "other/DelayedCast.hpp"
#include "other/DelayedTranspose.hpp"

#include "subset/DelayedSubsetBlock.hpp"
#include "subset/make_DelayedSubset.hpp"

#include "utils/compress_sparse_triplets.hpp"
#include "utils/convert_to_sparse.hpp"
#include "utils/convert_to_dense.hpp"
#include "utils/wrap_shared_ptr.hpp"
#include "utils/ArrayView.hpp"
#include "utils/SomeNumericArray.hpp"
#include "utils/bind_intersection.hpp"
#include "utils/Oracles.hpp"
#include "utils/process_consecutive_indices.hpp"

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
