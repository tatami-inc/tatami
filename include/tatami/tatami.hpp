#ifndef TATAMI_TATAMI_HPP
#define TATAMI_TATAMI_HPP

#include "dense/DenseMatrix.hpp"
#include "dense/convert_to_dense.hpp"

#include "sparse/CompressedSparseMatrix.hpp"
#include "sparse/SemiCompressedSparseMatrix.hpp"
#include "sparse/FragmentedSparseMatrix.hpp"
#include "sparse/convert_to_compressed_sparse.hpp"
#include "sparse/convert_to_fragmented_sparse.hpp"
#include "sparse/compress_sparse_triplets.hpp"

#include "isometric/unary/DelayedUnaryIsometricOp.hpp"
#include "isometric/binary/DelayedBinaryIsometricOp.hpp"

#include "other/DelayedBind.hpp"
#include "other/DelayedCast.hpp"
#include "other/DelayedTranspose.hpp"

#include "subset/DelayedSubsetBlock.hpp"
#include "subset/make_DelayedSubset.hpp"

#include "utils/wrap_shared_ptr.hpp"
#include "utils/ArrayView.hpp"
#include "utils/SomeNumericArray.hpp"
#include "utils/bind_intersection.hpp"
#include "utils/Oracles.hpp"
#include "utils/process_consecutive_indices.hpp"

#include "stats/sums.hpp"
#include "stats/variances.hpp"
#include "stats/medians.hpp"
#include "stats/ranges.hpp"
#include "stats/counts.hpp"
#include "stats/utils.hpp"
#include "stats/grouped_medians.hpp"
#include "stats/grouped_sums.hpp"

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
