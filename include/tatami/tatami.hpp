#ifndef TATAMI_TATAMI_HPP
#define TATAMI_TATAMI_HPP

#include "dense/DenseMatrix.hpp"
#include "dense/convert_to_dense.hpp"
#include "dense/transpose.hpp"
#include "dense/ForcedDense.hpp"

#include "sparse/CompressedSparseMatrix.hpp"
#include "sparse/FragmentedSparseMatrix.hpp"
#include "sparse/convert_to_compressed_sparse.hpp"
#include "sparse/convert_to_fragmented_sparse.hpp"
#include "sparse/compress_sparse_triplets.hpp"

#include "isometric/unary/DelayedUnaryIsometricOperation.hpp"
#include "isometric/unary/arithmetic_helpers.hpp"
#include "isometric/unary/math_helpers.hpp"
#include "isometric/unary/compare_helpers.hpp"
#include "isometric/unary/boolean_helpers.hpp"
#include "isometric/unary/substitute_helpers.hpp"
#include "isometric/unary/helper_interface.hpp"

#include "isometric/binary/DelayedBinaryIsometricOperation.hpp"
#include "isometric/binary/helper_interface.hpp"
#include "isometric/binary/arithmetic_helpers.hpp"
#include "isometric/binary/compare_helpers.hpp"
#include "isometric/binary/boolean_helpers.hpp"

#include "other/DelayedBind.hpp"
#include "other/DelayedCast.hpp"
#include "other/DelayedTranspose.hpp"
#include "other/ConstantMatrix.hpp"

#include "subset/DelayedSubsetBlock.hpp"
#include "subset/make_DelayedSubset.hpp"

#include "utils/wrap_shared_ptr.hpp"
#include "utils/ArrayView.hpp"
#include "utils/SomeNumericArray.hpp"
#include "utils/ConsecutiveOracle.hpp"
#include "utils/parallelize.hpp"
#include "utils/FixedOracle.hpp"
#include "utils/process_consecutive_indices.hpp"

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
