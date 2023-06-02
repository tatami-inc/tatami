<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.9.7">
  <compound kind="file">
    <name>DenseMatrix.hpp</name>
    <path>tatami/base/dense/</path>
    <filename>DenseMatrix_8hpp.html</filename>
    <class kind="class">tatami::DenseMatrix</class>
    <namespace>tatami</namespace>
    <member kind="typedef">
      <type>DenseMatrix&lt; false, Value_, Index_, Storage_ &gt;</type>
      <name>DenseColumnMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac47a769e00660eb7e9b5fcd543bcf2d3</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>DenseMatrix&lt; true, Value_, Index_, Storage_ &gt;</type>
      <name>DenseRowMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a51122d20490b377cd3f4609cc044f314</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>VirtualDenseMatrix.hpp</name>
    <path>tatami/base/dense/</path>
    <filename>VirtualDenseMatrix_8hpp.html</filename>
    <class kind="class">tatami::VirtualDenseMatrix</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>Extractor.hpp</name>
    <path>tatami/base/</path>
    <filename>Extractor_8hpp.html</filename>
    <class kind="struct">tatami::ExtractorBase</class>
    <class kind="struct">tatami::FullExtractor</class>
    <class kind="struct">tatami::BlockExtractor</class>
    <class kind="struct">tatami::IndexExtractor</class>
    <class kind="class">tatami::DenseExtractor</class>
    <class kind="class">tatami::SparseExtractor</class>
    <namespace>tatami</namespace>
    <member kind="typedef">
      <type>typename std::conditional&lt; selection_==DimensionSelectionType::FULL, FullExtractor&lt; Index_ &gt;, typename std::conditional&lt; selection_==DimensionSelectionType::BLOCK, BlockExtractor&lt; Index_ &gt;, IndexExtractor&lt; Index_ &gt; &gt;::type &gt;::type</type>
      <name>ConditionalSelectionExtractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a85ed61a4f772a2f7be4a12f739554e6e</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>typename std::conditional&lt; sparse_, SparseExtractor&lt; selection_, Value_, Index_ &gt;, DenseExtractor&lt; selection_, Value_, Index_ &gt; &gt;::type</type>
      <name>Extractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ae9f8db5316521603085577d977a6955f</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Extractor&lt; DimensionSelectionType::FULL, false, Value_, Index_ &gt;</type>
      <name>FullDenseExtractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a6008dbced6de41e5619156b5335f5762</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Extractor&lt; DimensionSelectionType::BLOCK, false, Value_, Index_ &gt;</type>
      <name>BlockDenseExtractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ae75de1fc78b7d361ea8b59a5379ea4da</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Extractor&lt; DimensionSelectionType::INDEX, false, Value_, Index_ &gt;</type>
      <name>IndexDenseExtractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a47ce406c32c3914c2ecce187e21b6ced</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Extractor&lt; DimensionSelectionType::FULL, true, Value_, Index_ &gt;</type>
      <name>FullSparseExtractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a0fbb0624c8e1913a87e8fb5c975400e1</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Extractor&lt; DimensionSelectionType::BLOCK, true, Value_, Index_ &gt;</type>
      <name>BlockSparseExtractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac8d0024399a66ce61f6315f5f46ebb63</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Extractor&lt; DimensionSelectionType::INDEX, true, Value_, Index_ &gt;</type>
      <name>IndexSparseExtractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a4b67b4d1b6c00cd0bd449703432a5f7b</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>extracted_length</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>af9d13ceaa112d2c091265510d741488d</anchor>
      <arglist>(const ConditionalSelectionExtractor&lt; selection_, Index_ &gt; &amp;ex)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Matrix.hpp</name>
    <path>tatami/base/</path>
    <filename>Matrix_8hpp.html</filename>
    <class kind="class">tatami::Matrix</class>
    <namespace>tatami</namespace>
    <member kind="typedef">
      <type>Matrix&lt; double, int &gt;</type>
      <name>NumericMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a35c670894994f1d620abb55953f98441</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Options.hpp</name>
    <path>tatami/base/</path>
    <filename>Options_8hpp.html</filename>
    <class kind="struct">tatami::Options</class>
    <class kind="struct">tatami::Oracle</class>
    <namespace>tatami</namespace>
    <member kind="enumeration">
      <type></type>
      <name>DimensionSelectionType</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a0a2ecaf58e2b69bb4a808e814aeb16a1</anchor>
      <arglist></arglist>
      <enumvalue file="namespacetatami.html" anchor="a0a2ecaf58e2b69bb4a808e814aeb16a1aba7de5bc6888294e5884b024a4c894f1">FULL</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="a0a2ecaf58e2b69bb4a808e814aeb16a1a4d34f53389ed7f28ca91fc31ea360a66">BLOCK</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="a0a2ecaf58e2b69bb4a808e814aeb16a1acb4ae3b37047fb4b2c0d16f8bf84f076">INDEX</enumvalue>
    </member>
  </compound>
  <compound kind="file">
    <name>DelayedBind.hpp</name>
    <path>tatami/base/other/</path>
    <filename>DelayedBind_8hpp.html</filename>
    <class kind="class">tatami::DelayedBind</class>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>make_DelayedBind</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a5f8e0f69139575707aa9314174b415b3</anchor>
      <arglist>(std::vector&lt; std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt; &gt; ps)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>DelayedCast.hpp</name>
    <path>tatami/base/other/</path>
    <filename>DelayedCast_8hpp.html</filename>
    <class kind="class">tatami::DelayedCast</class>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_out_, Index_out_ &gt; &gt;</type>
      <name>make_DelayedCast</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aca62c3bf751cdd06a08e8e503b0b591a</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_in_, Index_in_ &gt; &gt; p)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>DelayedTranspose.hpp</name>
    <path>tatami/base/other/</path>
    <filename>DelayedTranspose_8hpp.html</filename>
    <class kind="class">tatami::DelayedTranspose</class>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>make_DelayedTranspose</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>afa35d8e9fe286967f327ec0eb6bd5005</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>CompressedSparseMatrix.hpp</name>
    <path>tatami/base/sparse/</path>
    <filename>CompressedSparseMatrix_8hpp.html</filename>
    <class kind="class">tatami::CompressedSparseMatrix</class>
    <namespace>tatami</namespace>
    <member kind="typedef">
      <type>CompressedSparseMatrix&lt; false, Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_ &gt;</type>
      <name>CompressedSparseColumnMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a18cee3a5d9734f0092b03d023cfe4b6a</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>CompressedSparseMatrix&lt; true, Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_ &gt;</type>
      <name>CompressedSparseRowMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a0111adeeb583aeb7e24e9e1e25be4aa0</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>SemiCompressedSparseMatrix.hpp</name>
    <path>tatami/base/sparse/</path>
    <filename>SemiCompressedSparseMatrix_8hpp.html</filename>
    <class kind="class">tatami::SemiCompressedSparseMatrix</class>
    <namespace>tatami</namespace>
    <member kind="typedef">
      <type>SemiCompressedSparseMatrix&lt; false, Value_, Index_, IndexStorage_, PointerStorage_ &gt;</type>
      <name>SemiCompressedSparseColumnMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a6f03d0d880bc056e09c2cbb80eb2c2ec</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>SemiCompressedSparseMatrix&lt; true, Value_, Index_, IndexStorage_, PointerStorage_ &gt;</type>
      <name>SemiCompressedSparseRowMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a412cb6ee12f3ee81d404d6eb0e494e4d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>SparseRange.hpp</name>
    <path>tatami/base/</path>
    <filename>SparseRange_8hpp.html</filename>
    <class kind="struct">tatami::SparseRange</class>
    <class kind="struct">tatami::SparseRangeCopy</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>DelayedSubset.hpp</name>
    <path>tatami/base/subset/</path>
    <filename>DelayedSubset_8hpp.html</filename>
    <class kind="class">tatami::DelayedSubset</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>DelayedSubsetBlock.hpp</name>
    <path>tatami/base/subset/</path>
    <filename>DelayedSubsetBlock_8hpp.html</filename>
    <class kind="class">tatami::DelayedSubsetBlock</class>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>make_DelayedSubsetBlock</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a8f4aafc0a1fbdc0c31bc122d24122a63</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p, Index_ f, Index_ l)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>DelayedSubsetSorted.hpp</name>
    <path>tatami/base/subset/</path>
    <filename>DelayedSubsetSorted_8hpp.html</filename>
    <class kind="class">tatami::DelayedSubsetSorted</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>DelayedSubsetSortedUnique.hpp</name>
    <path>tatami/base/subset/</path>
    <filename>DelayedSubsetSortedUnique_8hpp.html</filename>
    <class kind="class">tatami::DelayedSubsetSortedUnique</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>DelayedSubsetUnique.hpp</name>
    <path>tatami/base/subset/</path>
    <filename>DelayedSubsetUnique_8hpp.html</filename>
    <class kind="class">tatami::DelayedSubsetUnique</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>make_DelayedSubset.hpp</name>
    <path>tatami/base/subset/</path>
    <filename>make__DelayedSubset_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>make_DelayedSubset</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>abfc63177b00e6e3e2fa47754b8d87704</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p, IndexStorage_ idx)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>utils.hpp</name>
    <path>tatami/stats/</path>
    <filename>stats_2utils_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>void</type>
      <name>parallelize</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a29ce7a2219ea60d45de1aa3d4de66063</anchor>
      <arglist>(Function_ fun, size_t tasks, size_t threads)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>consecutive_extractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a36c6ecf33bcb87e1ed33c0a7d744dd82</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *mat, Index_ iter_start, Index_ iter_length, Args_ &amp;&amp;... args)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>LruChunkCache.hpp</name>
    <path>tatami/chunked/</path>
    <filename>LruChunkCache_8hpp.html</filename>
    <class kind="class">tatami::LruChunkCache</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>OracleChunkCache.hpp</name>
    <path>tatami/chunked/</path>
    <filename>OracleChunkCache_8hpp.html</filename>
    <class kind="class">tatami::OracleChunkCache</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>convert_to_layered_sparse.hpp</name>
    <path>tatami/ext/layered/</path>
    <filename>convert__to__layered__sparse_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>LayeredMatrixData&lt; T, IDX &gt;</type>
      <name>convert_to_layered_sparse</name>
      <anchorfile>convert__to__layered__sparse_8hpp.html</anchorfile>
      <anchor>a8a10c9d0fdbf5076e8583d1068a7f6d6</anchor>
      <arglist>(const Matrix *incoming)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>LayeredMatrixData.hpp</name>
    <path>tatami/ext/layered/</path>
    <filename>LayeredMatrixData_8hpp.html</filename>
    <class kind="struct">tatami::LayeredMatrixData</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>MatrixMarket.hpp</name>
    <path>tatami/ext/</path>
    <filename>MatrixMarket_8hpp.html</filename>
  </compound>
  <compound kind="file">
    <name>layered.hpp</name>
    <path>tatami/ext/</path>
    <filename>layered_8hpp.html</filename>
  </compound>
  <compound kind="file">
    <name>layered.hpp</name>
    <path>tatami/ext/MatrixMarket/</path>
    <filename>MatrixMarket_2layered_8hpp.html</filename>
    <namespace>tatami</namespace>
    <namespace>tatami::MatrixMarket</namespace>
    <member kind="function">
      <type>LayeredMatrixData&lt; T, IDX &gt;</type>
      <name>load_layered_sparse_matrix_from_file</name>
      <anchorfile>namespacetatami_1_1MatrixMarket.html</anchorfile>
      <anchor>a1957591502ac286ad6b1bcae78ee75f9</anchor>
      <arglist>(const char *filepath, int compression=0, size_t bufsize=65536)</arglist>
    </member>
    <member kind="function">
      <type>LayeredMatrixData&lt; T, IDX &gt;</type>
      <name>load_layered_sparse_matrix_from_buffer</name>
      <anchorfile>namespacetatami_1_1MatrixMarket.html</anchorfile>
      <anchor>af7655515516f20bf22da4311b7e8dcc5</anchor>
      <arglist>(const unsigned char *buffer, size_t n, int compression=0, size_t bufsize=65536)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>simple.hpp</name>
    <path>tatami/ext/MatrixMarket/</path>
    <filename>simple_8hpp.html</filename>
    <class kind="struct">tatami::MatrixMarket::HeaderDetails</class>
    <namespace>tatami</namespace>
    <namespace>tatami::MatrixMarket</namespace>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; T, IDX &gt; &gt;</type>
      <name>load_sparse_matrix_from_file</name>
      <anchorfile>namespacetatami_1_1MatrixMarket.html</anchorfile>
      <anchor>acb8ab6a2dba7c9b9088245792d8dd4c5</anchor>
      <arglist>(const char *filepath, int compression=0, size_t bufsize=65536)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; T, IDX &gt; &gt;</type>
      <name>load_sparse_matrix_from_buffer</name>
      <anchorfile>namespacetatami_1_1MatrixMarket.html</anchorfile>
      <anchor>abca83cf0a5a51b93f78991719c8bd772</anchor>
      <arglist>(const unsigned char *buffer, size_t n, int compression=0, size_t bufsize=65536)</arglist>
    </member>
    <member kind="function">
      <type>HeaderDetails</type>
      <name>extract_header_from_file</name>
      <anchorfile>namespacetatami_1_1MatrixMarket.html</anchorfile>
      <anchor>afc195d253d64d0a67216cfd6d52e8639</anchor>
      <arglist>(const char *filepath, int compression=0, size_t bufsize=65536)</arglist>
    </member>
    <member kind="function">
      <type>HeaderDetails</type>
      <name>extract_header_from_buffer</name>
      <anchorfile>namespacetatami_1_1MatrixMarket.html</anchorfile>
      <anchor>ae0ced3181c20a3281a2ab4844b978313</anchor>
      <arglist>(const unsigned char *buffer, size_t n, int compression=0, size_t bufsize=65536)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>arith_scalar_helpers.hpp</name>
    <path>tatami/isometric_unary/</path>
    <filename>arith__scalar__helpers_8hpp.html</filename>
    <class kind="struct">tatami::DelayedAddScalarHelper</class>
    <class kind="struct">tatami::DelayedMultiplyScalarHelper</class>
    <class kind="struct">tatami::DelayedSubtractScalarHelper</class>
    <class kind="struct">tatami::DelayedDivideScalarHelper</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>arith_vector_helpers.hpp</name>
    <path>tatami/isometric_unary/</path>
    <filename>arith__vector__helpers_8hpp.html</filename>
    <class kind="struct">tatami::DelayedAddVectorHelper</class>
    <class kind="struct">tatami::DelayedSubtractVectorHelper</class>
    <class kind="struct">tatami::DelayedMultiplyVectorHelper</class>
    <class kind="struct">tatami::DelayedDivideVectorHelper</class>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>DelayedAddVectorHelper&lt; MARGIN, T, V &gt;</type>
      <name>make_DelayedAddVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aa82419090e54a4674d971429e93e6af5</anchor>
      <arglist>(V v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedSubtractVectorHelper&lt; RIGHT, MARGIN, T, V &gt;</type>
      <name>make_DelayedSubtractVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>abb68a8d27e821f213965edab7f029c48</anchor>
      <arglist>(V v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedMultiplyVectorHelper&lt; MARGIN, T, V &gt;</type>
      <name>make_DelayedMultiplyVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a316f4cf0ae5a9e83c2dad6ccbbbcfb43</anchor>
      <arglist>(V v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedDivideVectorHelper&lt; RIGHT, MARGIN, T, V &gt;</type>
      <name>make_DelayedDivideVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ab104e5029db4d09afe0b47a125f457bf</anchor>
      <arglist>(V v)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>DelayedUnaryIsometricOp.hpp</name>
    <path>tatami/isometric_unary/</path>
    <filename>DelayedUnaryIsometricOp_8hpp.html</filename>
    <class kind="class">tatami::DelayedUnaryIsometricOp</class>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>make_DelayedUnaryIsometricOp</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a07b072944c766d85dc4ce4dc734e0b75</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p, Operation_ op)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>math_helpers.hpp</name>
    <path>tatami/isometric_unary/</path>
    <filename>math__helpers_8hpp.html</filename>
    <class kind="struct">tatami::DelayedAbsHelper</class>
    <class kind="struct">tatami::DelayedLogHelper</class>
    <class kind="struct">tatami::DelayedSqrtHelper</class>
    <class kind="struct">tatami::DelayedLog1pHelper</class>
    <class kind="struct">tatami::DelayedRoundHelper</class>
    <class kind="struct">tatami::DelayedExpHelper</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>medians.hpp</name>
    <path>tatami/stats/</path>
    <filename>medians_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>column_medians</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a44202790861791b1ed9df6d480c69f6a</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>row_medians</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ab877f67694b3fa500f59469f044b1c01</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>ranges.hpp</name>
    <path>tatami/stats/</path>
    <filename>ranges_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>column_maxs</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a8102a2423efab5e329543f4235ac2290</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>row_maxs</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a47890bfd538f65971f9e5dba5e1ed785</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>column_mins</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a60e44772f82183100c7c25c181479b56</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>row_mins</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a3bcdcb0499b36b5b844ebd88fc30eefe</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; std::vector&lt; Output_ &gt;, std::vector&lt; Output_ &gt; &gt;</type>
      <name>column_ranges</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ab5eff79824930c458e730d23f640c1ab</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; std::vector&lt; Output_ &gt;, std::vector&lt; Output_ &gt; &gt;</type>
      <name>row_ranges</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a2b2475a166ad744b31b1275e26e562be</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>sums.hpp</name>
    <path>tatami/stats/</path>
    <filename>sums_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>column_sums</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a3aab6733d637b66abdf45f7300d2e1ba</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>row_sums</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a53a658059404691856bef57fb85d83d6</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>variances.hpp</name>
    <path>tatami/stats/</path>
    <filename>variances_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>std::pair&lt; Output_, Output_ &gt;</type>
      <name>compute_direct</name>
      <anchorfile>variances_8hpp.html</anchorfile>
      <anchor>a69156143bd931041e7a59f2ee469ff38</anchor>
      <arglist>(const Value_ *ptr, size_t n)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; Output_, Output_ &gt;</type>
      <name>compute_direct</name>
      <anchorfile>variances_8hpp.html</anchorfile>
      <anchor>aa55280f08a8d90b49d37c8d1e87869c0</anchor>
      <arglist>(const SparseRange&lt; Value_, Index_ &gt; &amp;range, size_t n)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute_running</name>
      <anchorfile>variances_8hpp.html</anchorfile>
      <anchor>a8fee3fa83a665089963ff712e0dfb6bb</anchor>
      <arglist>(const Value_ *ptr, size_t n, Output_ *means, Output_ *vars, int &amp;count)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute_running</name>
      <anchorfile>variances_8hpp.html</anchorfile>
      <anchor>a7441bd3d6a018a2ea941f778900dbd18</anchor>
      <arglist>(const SparseRange&lt; Value_, Index_ &gt; &amp;range, Output_ *means, Output_ *vars, Nonzero_ *nonzeros, int &amp;count, bool skip_zeros=true, Index_ subtract=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>finish_running</name>
      <anchorfile>variances_8hpp.html</anchorfile>
      <anchor>a3437fa820311dc3f221048e5ded97ef3</anchor>
      <arglist>(size_t n, Output_ *means, Output_ *vars, int count)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>finish_running</name>
      <anchorfile>variances_8hpp.html</anchorfile>
      <anchor>aaea56bd6999abb7c4979d07599c63824</anchor>
      <arglist>(size_t n, Output_ *means, Output_ *vars, const Nonzero_ *nonzeros, int count)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>column_variances</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac37af9a76d15f08a3634881696a27f65</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>row_variances</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ab806279616e19f6150376a7c07d5b64b</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>tatami.hpp</name>
    <path>tatami/</path>
    <filename>tatami_8hpp.html</filename>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>ArrayView.hpp</name>
    <path>tatami/utils/</path>
    <filename>ArrayView_8hpp.html</filename>
    <class kind="class">tatami::ArrayView</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>bind_intersection.hpp</name>
    <path>tatami/utils/</path>
    <filename>bind__intersection_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>std::pair&lt; std::shared_ptr&lt; Matrix &gt;, std::vector&lt; size_t &gt; &gt;</type>
      <name>bind_intersection</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a87c52cf8a974ad0025550b002d391d39</anchor>
      <arglist>(const std::vector&lt; std::shared_ptr&lt; Matrix &gt; &gt; &amp;inputs, const std::vector&lt; const Id * &gt; &amp;ids)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>compress_sparse_triplets.hpp</name>
    <path>tatami/utils/</path>
    <filename>compress__sparse__triplets_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>std::vector&lt; size_t &gt;</type>
      <name>compress_sparse_triplets</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ab17e92414b0bff60f7b7a6431ac8a330</anchor>
      <arglist>(size_t nr, size_t nc, U &amp;values, V &amp;rows, W &amp;cols)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>convert_to_dense.hpp</name>
    <path>tatami/utils/</path>
    <filename>convert__to__dense_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>void</type>
      <name>convert_to_dense</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>afaccb18d22f35a4ed6b6a7091342da06</anchor>
      <arglist>(const Matrix_ *incoming, StoredValue_ *store, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index &gt; &gt;</type>
      <name>convert_to_dense</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>af80143ac537339fe8dafd892632e96de</anchor>
      <arglist>(const Matrix_ *incoming, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_dense</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a64e8745f8aaaa5ad93c6ce65a0a6591e</anchor>
      <arglist>(const Matrix_ *incoming, int order, int threads=1)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>convert_to_sparse.hpp</name>
    <path>tatami/utils/</path>
    <filename>convert__to__sparse_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_sparse</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a220c47e470b3ade64ede7303e50ca92d</anchor>
      <arglist>(const InputMatrix_ *incoming, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_sparse</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a5847f3cd78ce89170ca43262d7f5f8e3</anchor>
      <arglist>(const InputMatrix_ *incoming, int order, int threads=1)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Oracles.hpp</name>
    <path>tatami/utils/</path>
    <filename>Oracles_8hpp.html</filename>
    <class kind="struct">tatami::FixedOracle</class>
    <class kind="struct">tatami::ConsecutiveOracle</class>
    <class kind="struct">tatami::OracleStream</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>process_consecutive_indices.hpp</name>
    <path>tatami/utils/</path>
    <filename>process__consecutive__indices_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>void</type>
      <name>process_consecutive_indices</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>af01c93a616eb99f6a17861a8b19f7ee0</anchor>
      <arglist>(const Index_ *indices, Index_ length, Function_ fun)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>SomeNumericArray.hpp</name>
    <path>tatami/utils/</path>
    <filename>SomeNumericArray_8hpp.html</filename>
    <class kind="struct">tatami::SomeNumericArray</class>
    <class kind="struct">tatami::SomeNumericArray::Iterator</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>wrap_shared_ptr.hpp</name>
    <path>tatami/utils/</path>
    <filename>wrap__shared__ptr_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>std::shared_ptr&lt; const Matrix&lt; T, IDX &gt; &gt;</type>
      <name>wrap_shared_ptr</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a7518f5e8e09a6f6d7d3955b8ea286689</anchor>
      <arglist>(const Matrix&lt; T, IDX &gt; *ptr)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::ArrayView</name>
    <filename>classtatami_1_1ArrayView.html</filename>
    <templarg>typename T</templarg>
    <member kind="function">
      <type></type>
      <name>ArrayView</name>
      <anchorfile>classtatami_1_1ArrayView.html</anchorfile>
      <anchor>a183bb27de4abe81292ecdd7be2633516</anchor>
      <arglist>(const T *p, size_t n)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>size</name>
      <anchorfile>classtatami_1_1ArrayView.html</anchorfile>
      <anchor>a7983c151067c852f6382b76661937dfd</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const T *</type>
      <name>data</name>
      <anchorfile>classtatami_1_1ArrayView.html</anchorfile>
      <anchor>a6a02f0144809d21ff20cc6018279f218</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const T *</type>
      <name>begin</name>
      <anchorfile>classtatami_1_1ArrayView.html</anchorfile>
      <anchor>a476ff30f6f88958e34401d254a2ff77d</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const T *</type>
      <name>end</name>
      <anchorfile>classtatami_1_1ArrayView.html</anchorfile>
      <anchor>a904e6d03b59a20f2cc7831fbb318cf33</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator[]</name>
      <anchorfile>classtatami_1_1ArrayView.html</anchorfile>
      <anchor>a329bcfebf2ebc7ca73d559836170722c</anchor>
      <arglist>(size_t i) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::BlockExtractor</name>
    <filename>structtatami_1_1BlockExtractor.html</filename>
    <templarg>typename Index_</templarg>
    <base>tatami::ExtractorBase</base>
    <member kind="variable">
      <type>Index_</type>
      <name>block_start</name>
      <anchorfile>structtatami_1_1BlockExtractor.html</anchorfile>
      <anchor>a6c013beccbd49cfa2cf5aae42098a527</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Index_</type>
      <name>block_length</name>
      <anchorfile>structtatami_1_1BlockExtractor.html</anchorfile>
      <anchor>a2421eff5ead7d4b29306374de30ca8a2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::CompressedSparseMatrix</name>
    <filename>classtatami_1_1CompressedSparseMatrix.html</filename>
    <templarg>bool row_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class ValueStorage_</templarg>
    <templarg>class IndexStorage_</templarg>
    <templarg>class PointerStorage_</templarg>
    <base>tatami::Matrix</base>
    <member kind="function">
      <type></type>
      <name>CompressedSparseMatrix</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>a537979311be374af1b03e5770a1dde60</anchor>
      <arglist>(Index_ nr, Index_ nc, ValueStorage_ vals, IndexStorage_ idx, PointerStorage_ ptr, bool check=true)</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>nrow</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>a3e67d6e54117433851d7c9453439f759</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>ncol</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>aa3869e4e92e47566b5d238c2ac87526c</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>a96666f952ad63ede2f7920ce72fdaf37</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sparse_proportion</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>a4f86d32bb4e4e74e3e20b60b02c45795</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>prefer_rows</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>a615e4caecc2e69d0adf16458908603e9</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>prefer_rows_proportion</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>acbe930f2e78252b1631d7657bbda2455</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>uses_oracle</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>a3cb93a9b9f549f26eb847e29d8fa61e2</anchor>
      <arglist>(bool) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>aa516aa869fbb00c8c71fd0742fdc8498</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>a4b6c0fd7cef115174287f339922d0796</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>ae350260bfe054c54eaf10e500bfb5f39</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>ac6ac5a086f25e46725f5fb0c574b9437</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>a0a4fcb350adb08def60565714524c43d</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>a30815fb2c7a4940e2f72e18e2ebf5dbf</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>a3573e9c46b77c06e88aa8142f611654c</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>aa82a4748a02e119301ce52baacc31c73</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>afe3e44e7ab9df3b45f0724951a662d17</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>ae56e0afe6873b2186a0945b9312519dc</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>ade534aca393bd02fe77ed58cc93454ad</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1CompressedSparseMatrix.html</anchorfile>
      <anchor>aa300f54d8f02e1fab49d6700c9564480</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::ConsecutiveOracle</name>
    <filename>structtatami_1_1ConsecutiveOracle.html</filename>
    <templarg>typename Index_</templarg>
    <base>tatami::Oracle</base>
    <member kind="function">
      <type></type>
      <name>ConsecutiveOracle</name>
      <anchorfile>structtatami_1_1ConsecutiveOracle.html</anchorfile>
      <anchor>a4e896a9ad264a7a88e28fd1c4bce76c3</anchor>
      <arglist>(Index_ s, Index_ l)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>predict</name>
      <anchorfile>structtatami_1_1ConsecutiveOracle.html</anchorfile>
      <anchor>a6ed7501a227b2ee9b96e78733f166f16</anchor>
      <arglist>(Index_ *buffer, size_t number)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedAbsHelper</name>
    <filename>structtatami_1_1DelayedAbsHelper.html</filename>
    <templarg>typename T</templarg>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structtatami_1_1DelayedAbsHelper.html</anchorfile>
      <anchor>ae1f88bd02c562b770648f23cc2711470</anchor>
      <arglist>(size_t r, size_t c, T val) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>sparse_</name>
      <anchorfile>structtatami_1_1DelayedAbsHelper.html</anchorfile>
      <anchor>aafc67de259524e824030550ab13011aa</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_row_</name>
      <anchorfile>structtatami_1_1DelayedAbsHelper.html</anchorfile>
      <anchor>aac72feb5d194fccb6809655e97db71e2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_column_</name>
      <anchorfile>structtatami_1_1DelayedAbsHelper.html</anchorfile>
      <anchor>a40fe4bb821193a0834f68b778c472656</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedAddScalarHelper</name>
    <filename>structtatami_1_1DelayedAddScalarHelper.html</filename>
    <templarg>typename T</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedAddScalarHelper</name>
      <anchorfile>structtatami_1_1DelayedAddScalarHelper.html</anchorfile>
      <anchor>a038479e8924a7acc01587a2ea0508a1b</anchor>
      <arglist>(T s)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structtatami_1_1DelayedAddScalarHelper.html</anchorfile>
      <anchor>ac285729c698438913524ced072fb6cbf</anchor>
      <arglist>(size_t r, size_t c, T val) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>sparse_</name>
      <anchorfile>structtatami_1_1DelayedAddScalarHelper.html</anchorfile>
      <anchor>a0e8ddacf3d91493b400f653adef66059</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_row_</name>
      <anchorfile>structtatami_1_1DelayedAddScalarHelper.html</anchorfile>
      <anchor>a15de20eb4e0babd35bd449c7c0356f78</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_column_</name>
      <anchorfile>structtatami_1_1DelayedAddScalarHelper.html</anchorfile>
      <anchor>ab1fff146fca9e14548a9f08188de0e35</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedAddVectorHelper</name>
    <filename>structtatami_1_1DelayedAddVectorHelper.html</filename>
    <templarg>int MARGIN</templarg>
    <templarg>typename T</templarg>
    <templarg>class V</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedAddVectorHelper</name>
      <anchorfile>structtatami_1_1DelayedAddVectorHelper.html</anchorfile>
      <anchor>ab4d8c15e312b6b3bafdcd7e9f8f2513e</anchor>
      <arglist>(V v)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structtatami_1_1DelayedAddVectorHelper.html</anchorfile>
      <anchor>ab4efab1cf8add8c7314b7f6da5bd3fa4</anchor>
      <arglist>(size_t r, size_t c, T val) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>sparse_</name>
      <anchorfile>structtatami_1_1DelayedAddVectorHelper.html</anchorfile>
      <anchor>a7f9ce0431bcad5ab9dc8a79b5e71e18b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_row_</name>
      <anchorfile>structtatami_1_1DelayedAddVectorHelper.html</anchorfile>
      <anchor>ae99620d5453f6f207da18771fd9f3d63</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_column_</name>
      <anchorfile>structtatami_1_1DelayedAddVectorHelper.html</anchorfile>
      <anchor>a935c85c09fb22b85db6b3c9149f9a2dd</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::DelayedBind</name>
    <filename>classtatami_1_1DelayedBind.html</filename>
    <templarg>int margin_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <base>Matrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>DelayedBind</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>a22ca07b4b4bea21b0bf2de420e780e5f</anchor>
      <arglist>(std::vector&lt; std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; &gt; ps)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DelayedBind</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>afb3555f25b08b90de7f3730698aef355</anchor>
      <arglist>(const std::vector&lt; std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt; &gt; &amp;ps)</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>nrow</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>ad4b4b14a53e79294bb3c51e2d03edb54</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>ncol</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>af3e218d269c0dd0bf61358a2df462def</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>a015cdac3de086438b5305195b84a2eee</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sparse_proportion</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>a34cf182a2f93e85e8d408c25e01a5ac6</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>prefer_rows</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>a29672009babe2d5bea7a95ab94a827df</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>prefer_rows_proportion</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>ab225a92ca17c09b1eed5d4086cb1e33d</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>uses_oracle</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>abd1a1ba10e6bd17707cd70852119e38c</anchor>
      <arglist>(bool row) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>ac56b0c60e1604760d061053d2a3f1cda</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>a9606ec7e5d1bfd116242aa8e8db5776f</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>aed6d8ffb63bb0b52c8535b458508214b</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>afa05d17a2ef5d844d4a79c13e42ed11c</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>a98d3d0040a8bd1442c7c48165dbc214c</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>af45d99f9425031c3cd3ab9a1e1e7adeb</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>a09843333a18e5e1917c96d6399078ce8</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>aca855b3323d167a3087e2426c97c990f</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>a905ee4de8c18087c1af4988b7d11fd05</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>a23015c47f274dd2d83f2e3b8ea90821c</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>a3cf86883e75357e64d89adbc48a70563</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedBind.html</anchorfile>
      <anchor>a1431ba4def7af7275f0a2fc2ae14ed01</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::DelayedCast</name>
    <filename>classtatami_1_1DelayedCast.html</filename>
    <templarg>typename Value_out_</templarg>
    <templarg>typename Index_out_</templarg>
    <templarg>typename Value_in_</templarg>
    <templarg>typename Index_in_</templarg>
    <base>Matrix&lt; Value_out_, Index_out_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>DelayedCast</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>affd97e917000abf9d782fdc9d6599c60</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_in_, Index_in_ &gt; &gt; p)</arglist>
    </member>
    <member kind="function">
      <type>Index_out_</type>
      <name>nrow</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>ac7cdb52624f338c4470990fd96a5ba25</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Index_out_</type>
      <name>ncol</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>adf4e97588cce683486827ab3325bb0ba</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>af0a4a7726a5b62c213c1d60470f026cb</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sparse_proportion</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>a70fbb58b425d6a2898dcbd22e6492795</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>prefer_rows</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>a1c131c401dfbd89d3908c522a8876554</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>prefer_rows_proportion</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>a1b56caede22f9ef0b5ea2165265f0f51</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>uses_oracle</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>a07167517245ed8f4d3c112e4f4f3edac</anchor>
      <arglist>(bool row) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_out_, Index_out_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>aa4ac0acd64fc490be115037e2192a409</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_out_, Index_out_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>a0afdc535aa0540ce1e0fce2beda40d1c</anchor>
      <arglist>(Index_out_ block_start, Index_out_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_out_, Index_out_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>a87aa25e476372f2b45080c050031e8d5</anchor>
      <arglist>(std::vector&lt; Index_out_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_out_, Index_out_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>a792d159e920c23cd07faf20a22cb6a26</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_out_, Index_out_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>a2791fdc5d90beacdf0b03e434896c307</anchor>
      <arglist>(Index_out_ block_start, Index_out_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_out_, Index_out_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>a482e7b7aadd3610e060e44c93dafc047</anchor>
      <arglist>(std::vector&lt; Index_out_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_out_, Index_out_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>a5fbaeb7a93a078983add7ff26c611e40</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_out_, Index_out_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>a1c77a0b18385f569df2fc7808e84fc92</anchor>
      <arglist>(Index_out_ block_start, Index_out_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_out_, Index_out_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>a476e3097dba559b475145350183593fa</anchor>
      <arglist>(std::vector&lt; Index_out_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_out_, Index_out_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>a52b9a5d62e3ce0b7618d763b82ba0416</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_out_, Index_out_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>a7ab6420842e1686361b23823a751fee8</anchor>
      <arglist>(Index_out_ block_start, Index_out_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_out_, Index_out_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedCast.html</anchorfile>
      <anchor>afad930739c2f81aa6e563201907c0b1b</anchor>
      <arglist>(std::vector&lt; Index_out_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedDivideScalarHelper</name>
    <filename>structtatami_1_1DelayedDivideScalarHelper.html</filename>
    <templarg>bool RIGHT</templarg>
    <templarg>typename T</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedDivideScalarHelper</name>
      <anchorfile>structtatami_1_1DelayedDivideScalarHelper.html</anchorfile>
      <anchor>a32dc69c234baba7c77f03a9ac41766fd</anchor>
      <arglist>(T s)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structtatami_1_1DelayedDivideScalarHelper.html</anchorfile>
      <anchor>a819a3adf111e3067178e92a6f7b2c8c5</anchor>
      <arglist>(size_t r, size_t c, T val) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>sparse_</name>
      <anchorfile>structtatami_1_1DelayedDivideScalarHelper.html</anchorfile>
      <anchor>aa52a33e8c62ec2fb4f6a2ad551fc2bd0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_row_</name>
      <anchorfile>structtatami_1_1DelayedDivideScalarHelper.html</anchorfile>
      <anchor>a004a96471b5b983ffdf72716f9ff9149</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_column_</name>
      <anchorfile>structtatami_1_1DelayedDivideScalarHelper.html</anchorfile>
      <anchor>a92771f2988221e29390cdf9bbdeb027a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedDivideVectorHelper</name>
    <filename>structtatami_1_1DelayedDivideVectorHelper.html</filename>
    <templarg>bool RIGHT</templarg>
    <templarg>int MARGIN</templarg>
    <templarg>typename T</templarg>
    <templarg>class V</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedDivideVectorHelper</name>
      <anchorfile>structtatami_1_1DelayedDivideVectorHelper.html</anchorfile>
      <anchor>adf8bad6ea3dd7523dee92f050ae33b13</anchor>
      <arglist>(V v)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structtatami_1_1DelayedDivideVectorHelper.html</anchorfile>
      <anchor>a63ca8047273ae80a88165d81be120146</anchor>
      <arglist>(size_t r, size_t c, T val) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>sparse_</name>
      <anchorfile>structtatami_1_1DelayedDivideVectorHelper.html</anchorfile>
      <anchor>afa05be70ca3b5ef7c8a5d9a2f479da39</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_row_</name>
      <anchorfile>structtatami_1_1DelayedDivideVectorHelper.html</anchorfile>
      <anchor>affbd7fdf4644051641b7af497a30e757</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_column_</name>
      <anchorfile>structtatami_1_1DelayedDivideVectorHelper.html</anchorfile>
      <anchor>af29325d80b8100a8a7f728e3d9050bb4</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedExpHelper</name>
    <filename>structtatami_1_1DelayedExpHelper.html</filename>
    <templarg>typename T</templarg>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structtatami_1_1DelayedExpHelper.html</anchorfile>
      <anchor>affeefc01a2a66d3704aef703e040d653</anchor>
      <arglist>(size_t r, size_t c, T val) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>sparse_</name>
      <anchorfile>structtatami_1_1DelayedExpHelper.html</anchorfile>
      <anchor>af7119e0283f7b49cfab0a1ebcaad1a63</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_row_</name>
      <anchorfile>structtatami_1_1DelayedExpHelper.html</anchorfile>
      <anchor>adebad4420563835a9b563939022f8caf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_column_</name>
      <anchorfile>structtatami_1_1DelayedExpHelper.html</anchorfile>
      <anchor>ad145ef0bef810261b962390f52fc0b74</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedLog1pHelper</name>
    <filename>structtatami_1_1DelayedLog1pHelper.html</filename>
    <templarg>typename T</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedLog1pHelper</name>
      <anchorfile>structtatami_1_1DelayedLog1pHelper.html</anchorfile>
      <anchor>a54be07a933d88fa1bfcd54d383a6b327</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DelayedLog1pHelper</name>
      <anchorfile>structtatami_1_1DelayedLog1pHelper.html</anchorfile>
      <anchor>ae1fbc1aaae2c404dd4a88b6a15034978</anchor>
      <arglist>(double base)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structtatami_1_1DelayedLog1pHelper.html</anchorfile>
      <anchor>a1b82fc162ab70e33f19542f93aa6d387</anchor>
      <arglist>(size_t r, size_t c, T val) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>sparse_</name>
      <anchorfile>structtatami_1_1DelayedLog1pHelper.html</anchorfile>
      <anchor>a8a0947ed6dbbec1f598367c2d5f305ce</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_row_</name>
      <anchorfile>structtatami_1_1DelayedLog1pHelper.html</anchorfile>
      <anchor>a61543ff7ed9e4985d57d8abeb8e8e44c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_column_</name>
      <anchorfile>structtatami_1_1DelayedLog1pHelper.html</anchorfile>
      <anchor>a8e5b43c2cc1b8409718e55798f89ea68</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedLogHelper</name>
    <filename>structtatami_1_1DelayedLogHelper.html</filename>
    <templarg>typename T</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedLogHelper</name>
      <anchorfile>structtatami_1_1DelayedLogHelper.html</anchorfile>
      <anchor>aad43deb006798a17b736adf7cbaa60ef</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DelayedLogHelper</name>
      <anchorfile>structtatami_1_1DelayedLogHelper.html</anchorfile>
      <anchor>a506aca4f5b63bcdaad3ddab4231fd0f7</anchor>
      <arglist>(double base)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structtatami_1_1DelayedLogHelper.html</anchorfile>
      <anchor>ab443629dd5741c6bba0deb07c8ff92c4</anchor>
      <arglist>(size_t r, size_t c, T val) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>sparse_</name>
      <anchorfile>structtatami_1_1DelayedLogHelper.html</anchorfile>
      <anchor>a7be77608c9405d255c03b3df1858979e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_row_</name>
      <anchorfile>structtatami_1_1DelayedLogHelper.html</anchorfile>
      <anchor>ad18c2630e67fa8faa2e66b8d4fed6d42</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_column_</name>
      <anchorfile>structtatami_1_1DelayedLogHelper.html</anchorfile>
      <anchor>a56a684bb72e50ee986b2c3ff71d0204d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedMultiplyScalarHelper</name>
    <filename>structtatami_1_1DelayedMultiplyScalarHelper.html</filename>
    <templarg>typename T</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedMultiplyScalarHelper</name>
      <anchorfile>structtatami_1_1DelayedMultiplyScalarHelper.html</anchorfile>
      <anchor>ae3e6b355cb53636374b5d229eb75c2e1</anchor>
      <arglist>(T s)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structtatami_1_1DelayedMultiplyScalarHelper.html</anchorfile>
      <anchor>aa4469882c45e9c6f37d478bdc9ae7984</anchor>
      <arglist>(size_t r, size_t c, T val) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>sparse_</name>
      <anchorfile>structtatami_1_1DelayedMultiplyScalarHelper.html</anchorfile>
      <anchor>aaaa889a6282ca9dcbdf4414da6798a27</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_row_</name>
      <anchorfile>structtatami_1_1DelayedMultiplyScalarHelper.html</anchorfile>
      <anchor>a4e86086797e8b2d0524bebe102ad8e1a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_column_</name>
      <anchorfile>structtatami_1_1DelayedMultiplyScalarHelper.html</anchorfile>
      <anchor>afa5fe94c939e25718a9515c2fc40c194</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedMultiplyVectorHelper</name>
    <filename>structtatami_1_1DelayedMultiplyVectorHelper.html</filename>
    <templarg>int MARGIN</templarg>
    <templarg>typename T</templarg>
    <templarg>class V</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedMultiplyVectorHelper</name>
      <anchorfile>structtatami_1_1DelayedMultiplyVectorHelper.html</anchorfile>
      <anchor>a264178e8902ffd7b728816ca0968e73d</anchor>
      <arglist>(V v)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structtatami_1_1DelayedMultiplyVectorHelper.html</anchorfile>
      <anchor>ac64849c5578ccfbb1262020421448f07</anchor>
      <arglist>(size_t r, size_t c, T val) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>sparse_</name>
      <anchorfile>structtatami_1_1DelayedMultiplyVectorHelper.html</anchorfile>
      <anchor>aecd7bcb7fcfb026308ab23aedf6bce4f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_row_</name>
      <anchorfile>structtatami_1_1DelayedMultiplyVectorHelper.html</anchorfile>
      <anchor>a6f767aebb25acf36159c92049dc42e1d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_column_</name>
      <anchorfile>structtatami_1_1DelayedMultiplyVectorHelper.html</anchorfile>
      <anchor>a900a5bee090779271fa72ac9a749b39c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedRoundHelper</name>
    <filename>structtatami_1_1DelayedRoundHelper.html</filename>
    <templarg>typename T</templarg>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structtatami_1_1DelayedRoundHelper.html</anchorfile>
      <anchor>aa6c9c7ab1d87872d51473361f25bb284</anchor>
      <arglist>(size_t r, size_t c, T val) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>sparse_</name>
      <anchorfile>structtatami_1_1DelayedRoundHelper.html</anchorfile>
      <anchor>a86a97924213721bcbdd7a60bf2e307b1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_row_</name>
      <anchorfile>structtatami_1_1DelayedRoundHelper.html</anchorfile>
      <anchor>ab0ef63f284ea59c04d87afb0e297ea23</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_column_</name>
      <anchorfile>structtatami_1_1DelayedRoundHelper.html</anchorfile>
      <anchor>a31e20fdc07a26688b906832f9c1d2325</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedSqrtHelper</name>
    <filename>structtatami_1_1DelayedSqrtHelper.html</filename>
    <templarg>typename T</templarg>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structtatami_1_1DelayedSqrtHelper.html</anchorfile>
      <anchor>ae0941c96f07180027d17918d6a8df941</anchor>
      <arglist>(size_t r, size_t c, T val) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>sparse_</name>
      <anchorfile>structtatami_1_1DelayedSqrtHelper.html</anchorfile>
      <anchor>a503a69a39bc995e83754e885d5728543</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_row_</name>
      <anchorfile>structtatami_1_1DelayedSqrtHelper.html</anchorfile>
      <anchor>a22cd75a6a73b81a3a864de25f5f62e24</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_column_</name>
      <anchorfile>structtatami_1_1DelayedSqrtHelper.html</anchorfile>
      <anchor>afe09f6ced951d8931418ee7ea6538178</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::DelayedSubset</name>
    <filename>classtatami_1_1DelayedSubset.html</filename>
    <templarg>int margin_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class IndexStorage_</templarg>
    <base>Matrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>DelayedSubset</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>a6b832655afd6ded8c5e5198dd576d604</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p, IndexStorage_ idx)</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>nrow</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>a363f65ce5c29b4a4b9c8e49037559f69</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>ncol</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>acf487dcd36c8612f8014bbe3b0362f4a</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>a2738d02a23fad9f5dca3125905533067</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sparse_proportion</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>acd7dfd1edbf3022c5140f5bc17330887</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>prefer_rows</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>ad7d937bb28815076108cf1d2cb7eb9ac</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>prefer_rows_proportion</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>a86417eb7b4ea08c49a03b91098b50568</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>uses_oracle</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>a3d6296d55f02ba1f39b5b9d66e50386c</anchor>
      <arglist>(bool row) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>a520f715166298b1a82e65f9a19f98686</anchor>
      <arglist>(const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>a37c70a9729008bf4c0295732edbee9ca</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>ad3d06f77fd7401cec1a497257655c280</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>aad418776eb05e3fe3c2f311cbcee2d5d</anchor>
      <arglist>(const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>a2bb1557b57f3048e9ad3046d637e3675</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>a72022dbd75040ebe21ffdf5483d3d932</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>a9604ab0813ca3a9827f7fdcb8276e9cb</anchor>
      <arglist>(const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>a2f37a9cf5db5736d85281274d5888998</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>adfa2f65489d5c90786416c3c714cf821</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>a82f97c7cd46e9e4a3ce9f20114553a81</anchor>
      <arglist>(const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>acf5f157d2a355448c02dba67e6bb17bb</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedSubset.html</anchorfile>
      <anchor>a0bb000a1ebd708a148785f8a904cc9f0</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;options) const</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::DelayedSubsetBlock</name>
    <filename>classtatami_1_1DelayedSubsetBlock.html</filename>
    <templarg>int margin_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <base>Matrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>DelayedSubsetBlock</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>abeaa6e9309f85b80dcb32f68c30e44c7</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p, Index_ f, Index_ l)</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>nrow</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>ad464fba67568e017244f1a5378262b1f</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>ncol</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>aeb6512ab3d6b70fa3ead75fd47995a60</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>afbbdf9fe3282eb0e157020b44812daf1</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sparse_proportion</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>a6d5db7efb0567a2447c82f1cefad7997</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>prefer_rows</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>a00ebfc270ed722fad12d9ca3349cbd33</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>prefer_rows_proportion</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>a9de5fba27b120ec0a07c008ab964524a</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>uses_oracle</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>ac989b4ed58351821d0b4e995db3ab658</anchor>
      <arglist>(bool row) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>a8cb0f716fed6ad1683a098accf4153b7</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>a0c0e2302475402e3240e07529a802c58</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>a03fa3bbc8570e89733b725ed9296c215</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>a489ddc00702fe3699f818a13d9707c59</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>a42d9195918c209e70759d5b8cb2444f9</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>aedc41486ad938840ec159620b5bc120a</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>a7b4bccf8b6e5d5473934499a1be8c4ba</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>ad99125d6c095e615f104f86d4e93633a</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>a495c8cb29e2a3fc049fe4468c0fbe8e4</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>ac1f92e628289749083b8459f05cbd84d</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>ad5be144b2bddc041681ce3b340844708</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetBlock.html</anchorfile>
      <anchor>a2bc07cc9c9319a8a4f4eaf16caf74bf3</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::DelayedSubsetSorted</name>
    <filename>classtatami_1_1DelayedSubsetSorted.html</filename>
    <templarg>int margin_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class IndexStorage_</templarg>
    <base>Matrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>DelayedSubsetSorted</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>acb54e427b4b971040fff8b0b327337bf</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p, IndexStorage_ idx, bool check=true)</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>nrow</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>a4244a48118361e58e3b8dcf26384df45</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>ncol</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>abffcd6776a07da8f14fa4ca2b3410e8d</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>ad3de98a3b6a0bd33081544b4757d014f</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sparse_proportion</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>a75b68a4d28ecc0073a81e6c407cb8930</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>prefer_rows</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>a685fbba599f884a4cf95a44e4fc4c7b9</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>prefer_rows_proportion</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>aa643b08f0d8443d1611665205fd0bfe1</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>uses_oracle</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>a810f2bab68978a1297e44fd6015dcb31</anchor>
      <arglist>(bool row) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>a66efaad2aa759c175c7d1e2b9d9d9919</anchor>
      <arglist>(const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>a0f5aae7a06ea022edba764e349be4ef1</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>a82f95e2c19acbaa7d7d5ce3045ec8ede</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>af61b1f12840591ac26292d79905bcd48</anchor>
      <arglist>(const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>ac296502d83bd8d623bac89f8afc708d9</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>a0d189857678a091bf263269e150217e0</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>a8f36a618c9495e2899a5dff11bd6c75b</anchor>
      <arglist>(const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>af314f2bb2c3d535e271f412d8c76f066</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>ac6497dc8677ec669a2bce95c166b02ef</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>a37f83309c4d561f45702596224778486</anchor>
      <arglist>(const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>a195bdde3ffd966b3a421a3ee263d4b1b</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetSorted.html</anchorfile>
      <anchor>a3a6143b184db99fc1ab30422fdaec5f9</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;options) const</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::DelayedSubsetSortedUnique</name>
    <filename>classtatami_1_1DelayedSubsetSortedUnique.html</filename>
    <templarg>int margin_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class IndexStorage_</templarg>
    <base>Matrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>DelayedSubsetSortedUnique</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>a4a56857aecad86ffbf3345bef61f99e6</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p, IndexStorage_ idx, bool check=true)</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>nrow</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>ae01ccb17c3afde476f89dfe6b6312b45</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>ncol</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>ad37c87845de22fd19218b32672f0a32f</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>a964965ffd2dbb46be8c2e43e43fe7021</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sparse_proportion</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>acc91dead096884f86ef77b76a20dd775</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>prefer_rows</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>abc50ddffc07bfe1075abcd7b4cc41e02</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>prefer_rows_proportion</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>ad55e79fa2702133b64af71e5b93a41e3</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>uses_oracle</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>afd0685e407df2d3c0c08419f06d08802</anchor>
      <arglist>(bool row) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>a95eaa479c3e83de0e1b53470253db651</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>afb987ea0309670ed53272bccfa99affd</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>a24c590cca14c2046c10f701ed0f402c0</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>af511de57e55eb7050552c6e92f4111ff</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>a69a4131fd196926170018947359e38ba</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>ab1d00c4b33bc6a8b70cf7b62744d69cf</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>a6964e114bbe219fdd321ca86b51b04c5</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>ad52e3e667f4824187dfd85287060cbdc</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>a41812757a117fd56cccd2268729b56c1</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>aa78dc1f08ae74964ecfdc6dbf3d1f4b8</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>a256fe6dc55e2b7317292f5ccac395a0e</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetSortedUnique.html</anchorfile>
      <anchor>a398c9a5518f32c1ecf3db3629e3c5a76</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::DelayedSubsetUnique</name>
    <filename>classtatami_1_1DelayedSubsetUnique.html</filename>
    <templarg>int margin_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class IndexStorage_</templarg>
    <base>Matrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>DelayedSubsetUnique</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>a6137eae868d3ca1c7116d3440278dbef</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p, IndexStorage_ idx, bool check=true)</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>nrow</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>a8cd3203978e03fc049acdb5c30609204</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>ncol</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>abc7eb3bb7c9f4024489deae9c1d3a78d</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>ab89df61da7d2dbda82198df41961f02d</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sparse_proportion</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>a54c985af221be5d60b406f61e741691c</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>prefer_rows</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>a8d030a0b68a377bc30962b0662acd856</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>prefer_rows_proportion</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>af303be4a8cf50222b93bd126a69920d7</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>uses_oracle</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>af9d1f0ab74ed37b79696f6fc993014f4</anchor>
      <arglist>(bool row) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>ab457721fe5480a252d4faf5f5fb6bf44</anchor>
      <arglist>(const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>a7b3d5379597fabc1bc2cb2be70ee5859</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>a2f3685a1092d6497bbdac10ab5578e05</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>a0ae48533d1627324e08ca464aaff25d7</anchor>
      <arglist>(const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>a976cc7b3e6808575f4b0eb3c3125a9f6</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>a0e1ee2c7754ccceb01ff02cce4340b6a</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>a3a2f7dc0c2816054c2ca48aca530be0a</anchor>
      <arglist>(const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>a680baef1a688e43eaf04e86f6c4d57b3</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>a190cd0ada220fa3664c2bfdcf8e9ea0a</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>a77d74e2c87b0ca638c5810cf8778369b</anchor>
      <arglist>(const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>a8823ff62bcb495f742a1943637196247</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;options) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedSubsetUnique.html</anchorfile>
      <anchor>a4eb8f770e84ca2cb67bedb287a0d27ad</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;options) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedSubtractScalarHelper</name>
    <filename>structtatami_1_1DelayedSubtractScalarHelper.html</filename>
    <templarg>bool RIGHT</templarg>
    <templarg>typename T</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedSubtractScalarHelper</name>
      <anchorfile>structtatami_1_1DelayedSubtractScalarHelper.html</anchorfile>
      <anchor>abf2792fab8a48f9a89ee50ad6c6d8da3</anchor>
      <arglist>(T s)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structtatami_1_1DelayedSubtractScalarHelper.html</anchorfile>
      <anchor>ad686b1974b4ff2023a140fef7b9b2620</anchor>
      <arglist>(size_t r, size_t c, T val) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>sparse_</name>
      <anchorfile>structtatami_1_1DelayedSubtractScalarHelper.html</anchorfile>
      <anchor>a22961a62ff0ca4a14e14fd58189bf4a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_row_</name>
      <anchorfile>structtatami_1_1DelayedSubtractScalarHelper.html</anchorfile>
      <anchor>af8863bfbf50f5d116194937d5340c424</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_column_</name>
      <anchorfile>structtatami_1_1DelayedSubtractScalarHelper.html</anchorfile>
      <anchor>a94a1f89fb9d8badc425504444bccca2c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedSubtractVectorHelper</name>
    <filename>structtatami_1_1DelayedSubtractVectorHelper.html</filename>
    <templarg>bool RIGHT</templarg>
    <templarg>int MARGIN</templarg>
    <templarg>typename T</templarg>
    <templarg>class V</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedSubtractVectorHelper</name>
      <anchorfile>structtatami_1_1DelayedSubtractVectorHelper.html</anchorfile>
      <anchor>ad9f4a928250416796f4ead2b451d0701</anchor>
      <arglist>(V v)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structtatami_1_1DelayedSubtractVectorHelper.html</anchorfile>
      <anchor>af7b4833df07f470bcbb7bf82c9c5efbe</anchor>
      <arglist>(size_t r, size_t c, T val) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>sparse_</name>
      <anchorfile>structtatami_1_1DelayedSubtractVectorHelper.html</anchorfile>
      <anchor>a611fa78d54d10f856bcfb6ba67db4c2c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_row_</name>
      <anchorfile>structtatami_1_1DelayedSubtractVectorHelper.html</anchorfile>
      <anchor>a72b6ea7ba4c03c81abb4f7d82510a2f3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const bool</type>
      <name>needs_column_</name>
      <anchorfile>structtatami_1_1DelayedSubtractVectorHelper.html</anchorfile>
      <anchor>ac14a39bd201052add2dab5830b03bcdf</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::DelayedTranspose</name>
    <filename>classtatami_1_1DelayedTranspose.html</filename>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <base>Matrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>DelayedTranspose</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>a1650e93e9ce61974a46be83ea67493df</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p)</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>nrow</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>ae861c999ade80dfaab4e00ea0a3c9f89</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>ncol</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>a495080ac6d15551d6816adf85ab2d481</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>a7e981425a1c2653be11ca4e1b14a2a64</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sparse_proportion</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>a778b7daddea3238c9aa8efe01b059fea</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>prefer_rows</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>a8d386b8a0d05d6bcdbc1dd80c7d8118e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>prefer_rows_proportion</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>af218bfcac74f403cc9d5bd6be84cc7a5</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>uses_oracle</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>a6e814addfd1b93d7999a5383d0e69104</anchor>
      <arglist>(bool row) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>a88e3074c6b7d05f49f6ae9a0ec660729</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>a5e8cedd000e22ff767c5acc9fb49c704</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>aac3f9816a6694fda0a957c77ba4f7539</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>aca3af5f469b31c8cdfd17110f46c1d68</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>a78584737c5b404e06898550a8234b86a</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>a6cddd31a66c2e566a1969fc1146099bb</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>aa8f0b02f3dc8f7c6207f8a90f703f81b</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>a94b76d50f8046fd5f3abc99f1dd84fef</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>af9bdca9700d767e78128e3dcbad9e9ff</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>ac9ee6ab5a2900718bfc929520324d7b3</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>ace4ddb430ed329ab3ddf12c45c6f18a2</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedTranspose.html</anchorfile>
      <anchor>aaa45e99399c24068984dc2159d1be726</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::DelayedUnaryIsometricOp</name>
    <filename>classtatami_1_1DelayedUnaryIsometricOp.html</filename>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class Operation_</templarg>
    <base>Matrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>DelayedUnaryIsometricOp</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>a308f8ac0b9e202bcc1865500f28ee0cf</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p, Operation_ op)</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>nrow</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>a27526af79dcb3b5a3b9460733355137e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>ncol</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>a57d5be64d07b2ca81dc232cbac702ea2</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>abf2131b00738a76e81afb260bc3e4ed4</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sparse_proportion</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>a8975153bf959ef935fb88e0558d0d3e6</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>prefer_rows</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>af5b5d90d807777c8ee966416b2fe44b8</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>prefer_rows_proportion</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>a228ccbff3b0bd609ada680f2ed419224</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>uses_oracle</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>a9cbdbefdafe60b00f400fac148339566</anchor>
      <arglist>(bool row) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>a8621b45c366fc253ace223469bcf8ff3</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>a1ec5c659f9ee6234c4c2cba95dcb9e0d</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>af6551d961b81422de6ef6945a383d6ba</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>ad7f567c47c9e952162dc0ea08ed993b5</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>af17b48717458ee16db84fb5c68d3aff0</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>ac977f1622ef3396eeadb53fd4211c45d</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>acd65b765fcf7d80bf3798f545ed4a6f9</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>abdf36bab2e9f3570338db28d7ae8ca2d</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>afec1641d42641341033df0a01f966579</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>a3a71b16bc5028eaaffc19c4eb6deb6a2</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>a0d17779f2558cb0325b3a13b3db7d63d</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedUnaryIsometricOp.html</anchorfile>
      <anchor>a0b8635c5cfdd8822e289248939b277a0</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::DenseExtractor</name>
    <filename>classtatami_1_1DenseExtractor.html</filename>
    <templarg>DimensionSelectionType selection_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <member kind="function" virtualness="pure">
      <type>virtual const Value_ *</type>
      <name>fetch</name>
      <anchorfile>classtatami_1_1DenseExtractor.html</anchorfile>
      <anchor>a809f1682a0ba66035066175e37bb0e83</anchor>
      <arglist>(Index_ i, Value_ *buffer)=0</arglist>
    </member>
    <member kind="function">
      <type>const Value_ *</type>
      <name>fetch_copy</name>
      <anchorfile>classtatami_1_1DenseExtractor.html</anchorfile>
      <anchor>a67bd8c1e45a792c99014e1ea625fa13f</anchor>
      <arglist>(Index_ i, Value_ *buffer)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Value_ &gt;</type>
      <name>fetch</name>
      <anchorfile>classtatami_1_1DenseExtractor.html</anchorfile>
      <anchor>a2a61bcd97df31e1f5061de8fd0c7a8e4</anchor>
      <arglist>(Index_ i)</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static constexpr bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1DenseExtractor.html</anchorfile>
      <anchor>a439bdce51f6ef1467f7973af5c76d5d1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static constexpr DimensionSelectionType</type>
      <name>selection</name>
      <anchorfile>classtatami_1_1DenseExtractor.html</anchorfile>
      <anchor>ae0fd341cc3f969a13ab88c6747d15a87</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::DenseMatrix</name>
    <filename>classtatami_1_1DenseMatrix.html</filename>
    <templarg>bool row_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class Storage_</templarg>
    <base>tatami::VirtualDenseMatrix</base>
    <member kind="function">
      <type></type>
      <name>DenseMatrix</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>a114b648506fa0802dfd07a416d52cc09</anchor>
      <arglist>(Index_ nr, Index_ nc, const Storage_ &amp;source)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DenseMatrix</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>a88972832ba613fa4ea6820fdbe5dcacd</anchor>
      <arglist>(Index_ nr, Index_ nc, Storage_ &amp;&amp;source)</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>nrow</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>a108c3fffbcdd9ca44cde3c5e2ec8ce9a</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>ncol</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>ae025690fa69ae5acc9c2c25931c30efb</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>prefer_rows</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>ab2dad56ac067668f2b845b7bcacf2f94</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>uses_oracle</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>a9b62638264c912ecb06cd8cf66bbe967</anchor>
      <arglist>(bool) const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>prefer_rows_proportion</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>a6e2b959fc2b6f7de296d244cbba5cfa8</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>afcb586f220c1df06c464def9c6c4f004</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>ae94aa68c8f754dbd45bad996a4d4af9d</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>ae9fe012376d25e7357d4590f200606fd</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>a8ace942278f3ba1638401d5fd864d28a</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>a08ca9bd8c2cd11fc553110486e7f7f22</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>a94146a6703b5315b87257acc5115b153</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::ExtractorBase</name>
    <filename>structtatami_1_1ExtractorBase.html</filename>
    <templarg>typename Index_</templarg>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>set_oracle</name>
      <anchorfile>structtatami_1_1ExtractorBase.html</anchorfile>
      <anchor>ab0e98cdedccd44333ab88cee7548221b</anchor>
      <arglist>(std::unique_ptr&lt; Oracle&lt; Index_ &gt; &gt; o)=0</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::FixedOracle</name>
    <filename>structtatami_1_1FixedOracle.html</filename>
    <templarg>typename Index_</templarg>
    <base>tatami::Oracle</base>
    <member kind="function">
      <type></type>
      <name>FixedOracle</name>
      <anchorfile>structtatami_1_1FixedOracle.html</anchorfile>
      <anchor>ac04dcf58c95882cf83bc5204b1944776</anchor>
      <arglist>(const Index_ *r, size_t n)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>predict</name>
      <anchorfile>structtatami_1_1FixedOracle.html</anchorfile>
      <anchor>a931ccde93f123c8c718f808e89cf8f36</anchor>
      <arglist>(Index_ *predicted, size_t number)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::FullExtractor</name>
    <filename>structtatami_1_1FullExtractor.html</filename>
    <templarg>typename Index_</templarg>
    <base>tatami::ExtractorBase</base>
    <member kind="variable">
      <type>Index_</type>
      <name>full_length</name>
      <anchorfile>structtatami_1_1FullExtractor.html</anchorfile>
      <anchor>a9dea647495eb0edbed7b77502ca528e3</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::MatrixMarket::HeaderDetails</name>
    <filename>structtatami_1_1MatrixMarket_1_1HeaderDetails.html</filename>
    <member kind="variable">
      <type>size_t</type>
      <name>nrow</name>
      <anchorfile>structtatami_1_1MatrixMarket_1_1HeaderDetails.html</anchorfile>
      <anchor>a5cdfdb10078f2f25cae2f2aa22279007</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>size_t</type>
      <name>ncol</name>
      <anchorfile>structtatami_1_1MatrixMarket_1_1HeaderDetails.html</anchorfile>
      <anchor>a35de5188ab362aaa52aa39d2cfd79188</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>size_t</type>
      <name>nlines</name>
      <anchorfile>structtatami_1_1MatrixMarket_1_1HeaderDetails.html</anchorfile>
      <anchor>a63af108a27977c4e626181e58f28cb95</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::IndexExtractor</name>
    <filename>structtatami_1_1IndexExtractor.html</filename>
    <templarg>typename Index_</templarg>
    <base>tatami::ExtractorBase</base>
    <member kind="function" virtualness="pure">
      <type>virtual const Index_ *</type>
      <name>index_start</name>
      <anchorfile>structtatami_1_1IndexExtractor.html</anchorfile>
      <anchor>a089566f519e13b933bc58c9cf48f14bf</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="variable">
      <type>Index_</type>
      <name>index_length</name>
      <anchorfile>structtatami_1_1IndexExtractor.html</anchorfile>
      <anchor>ace949af0d062f6d1a07c2771d072805a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::SomeNumericArray::Iterator</name>
    <filename>structtatami_1_1SomeNumericArray_1_1Iterator.html</filename>
    <member kind="typedef">
      <type>std::random_access_iterator_tag</type>
      <name>iterator_category</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>aaf0075ceb671a85fc42404298baeed07</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>std::ptrdiff_t</type>
      <name>difference_type</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>ae2439d56ee9c796afae3b07f70e1c3b1</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>T</type>
      <name>value_type</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>a452df712534d47fdb8543805affd038b</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>const T *</type>
      <name>pointer</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>adfb70cf5c5bbfdfdebe7f20244fc5f30</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>const T &amp;</type>
      <name>reference</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>ab83a68ab47e92b1265fd043bc4351be9</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Iterator</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>aaec76c06dc4b04937e4c0e16233a7fdb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Iterator</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>a47aa875df6ad992824acccefb811c4d4</anchor>
      <arglist>(const SomeNumericArray *p, size_t i)</arglist>
    </member>
    <member kind="function">
      <type>value_type</type>
      <name>operator*</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>adadc4fec22839e1d7f826594eb605fc1</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>value_type</type>
      <name>operator[]</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>a220df2018efe0af029200950bf8f2656</anchor>
      <arglist>(size_t i) const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>afd39145d4998f6c0fc3d2cdc1f0ecd64</anchor>
      <arglist>(const Iterator &amp;right) const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator!=</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>acab7e46b690f81ad367780bd6cb33ebd</anchor>
      <arglist>(const Iterator &amp;right) const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator&lt;</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>a3a0271b552cbd873c9598740d472e00e</anchor>
      <arglist>(const Iterator &amp;right) const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator&gt;=</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>a82077d04886a5a0088372b03707d2fca</anchor>
      <arglist>(const Iterator &amp;right) const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator&gt;</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>a4296f53f8ed91bb3324e9a681b2d3c35</anchor>
      <arglist>(const Iterator &amp;right) const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator&lt;=</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>aed25f7d5d5e8252746cac7c2ad2fccd6</anchor>
      <arglist>(const Iterator &amp;right) const</arglist>
    </member>
    <member kind="function">
      <type>Iterator &amp;</type>
      <name>operator+=</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>a631fc3c3234a2229f01944127ec58988</anchor>
      <arglist>(size_t n)</arglist>
    </member>
    <member kind="function">
      <type>Iterator &amp;</type>
      <name>operator++</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>ad7fcc6f55aa2b723c4f9b6b3505f27a9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Iterator</type>
      <name>operator++</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>a28163076266a7b4987d8ff0fbfbaf4c2</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>Iterator &amp;</type>
      <name>operator-=</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>a45ee8f9680f7184067f58cf8835ea0ce</anchor>
      <arglist>(size_t n)</arglist>
    </member>
    <member kind="function">
      <type>Iterator &amp;</type>
      <name>operator--</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>a259cf6cdfcfb7c89cf6a64262ae3af4c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Iterator</type>
      <name>operator--</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>a41dbc8abc6e9a4c5bc2493d549ee862c</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>Iterator</type>
      <name>operator+</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>af692b5785cec55d153341dfabe197bc3</anchor>
      <arglist>(size_t n) const</arglist>
    </member>
    <member kind="function">
      <type>Iterator</type>
      <name>operator-</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>acfb66b80475de1ad9aa1233ccf224691</anchor>
      <arglist>(size_t n) const</arglist>
    </member>
    <member kind="function">
      <type>std::ptrdiff_t</type>
      <name>operator-</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>a87e797f6c5fc8e7433ba17eb89c6d18b</anchor>
      <arglist>(const Iterator &amp;right) const</arglist>
    </member>
    <member kind="friend">
      <type>friend Iterator</type>
      <name>operator+</name>
      <anchorfile>structtatami_1_1SomeNumericArray_1_1Iterator.html</anchorfile>
      <anchor>af891cdee3c447bfccb98be31584ab4a8</anchor>
      <arglist>(size_t n, const Iterator &amp;it)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::LayeredMatrixData</name>
    <filename>structtatami_1_1LayeredMatrixData.html</filename>
    <templarg>typename T</templarg>
    <templarg>typename IDX</templarg>
    <member kind="variable">
      <type>std::shared_ptr&lt; Matrix&lt; T, IDX &gt; &gt;</type>
      <name>matrix</name>
      <anchorfile>structtatami_1_1LayeredMatrixData.html</anchorfile>
      <anchor>a117e250097c0219902fad9160e5d389c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; size_t &gt;</type>
      <name>permutation</name>
      <anchorfile>structtatami_1_1LayeredMatrixData.html</anchorfile>
      <anchor>ae9c9cfb331519d5c7e18f4e51fc2898c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::LruChunkCache</name>
    <filename>classtatami_1_1LruChunkCache.html</filename>
    <templarg>typename Id_</templarg>
    <templarg>class ChunkContents_</templarg>
    <member kind="function">
      <type></type>
      <name>LruChunkCache</name>
      <anchorfile>classtatami_1_1LruChunkCache.html</anchorfile>
      <anchor>a064612c107fa2ce3625dc512bbcb3a0f</anchor>
      <arglist>(size_t m)</arglist>
    </member>
    <member kind="function">
      <type>const ChunkContents_ &amp;</type>
      <name>find_chunk</name>
      <anchorfile>classtatami_1_1LruChunkCache.html</anchorfile>
      <anchor>a17636ac967438f3ca1579a33af523e03</anchor>
      <arglist>(Id_ id, Cfunction_ create, Pfunction_ populate)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::Matrix</name>
    <filename>classtatami_1_1Matrix.html</filename>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <member kind="typedef">
      <type>Value_</type>
      <name>value_type</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a5b783d3ab440b7696c635a4b16516ba8</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Index_</type>
      <name>index_type</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>adebf536caf3a3c1751eefc6e10f0b16e</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual Index_</type>
      <name>nrow</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a54fe7b4baf78069d35ff00357a1b6cc6</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual Index_</type>
      <name>ncol</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a1f6fec43d5bae3f831841646bc02377b</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a9de24a05a6aee62d4ba59799a73242e7</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual double</type>
      <name>sparse_proportion</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a8a1a29b602b903201fc36ad65f76c78c</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>prefer_rows</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a69382d2a4e66cda1a9eb48dc25ab113a</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual double</type>
      <name>prefer_rows_proportion</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a2ddcc730155ca894df1c58b0a963d5e7</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>uses_oracle</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a68da21245203cf9349648c5452cdb2ca</anchor>
      <arglist>(bool row) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>aa29b9428f53453369256d7febc9ed0fc</anchor>
      <arglist>(const Options &amp;opt) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a7d96127de1d22a5d976b9a2f46ae0ab4</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>aba9cc13a9df94e59d03772bde6ef565a</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a9f6c3453a5c0c28a764f300f45bf860f</anchor>
      <arglist>(const Options &amp;opt) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a3ddceca01de5d1985aafa6e8bfdcf604</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a28ddabd4babf3d162f0901578c520a9a</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>adee56bee86c500a4f4a3940e8ef18459</anchor>
      <arglist>(const Options &amp;opt) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a82e592bf93aa216045a93efb923107e6</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>aa452ec9e476c69cdd9c7079a4d382af9</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a6fcbdc787f4b4980f57c228dc89f6289</anchor>
      <arglist>(const Options &amp;opt) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a369b6b56d1f290d0b298ef840119468a</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a120c3a294a1ffd715c86b136902b8c51</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const =0</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a3913e6cfd25586f37c0833a21c8f0ca6</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a6456594fe017e5f61a2adce0ce152ee4</anchor>
      <arglist>(Index_ block_start, Index_ block_length) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>aa2b4e61e29787566fcd78209d4930c2e</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a94a543328188b8212c9267f9c5edc431</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a39764e844f1da2b4027b61eb92b28674</anchor>
      <arglist>(Index_ block_start, Index_ block_length) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a35af5355144252c421d3f4ee6c2aa51b</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a7274f8de1b6918f45291e10c054252a1</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>ae49828cfac71e469c569bd74a4ed006e</anchor>
      <arglist>(Index_ block_start, Index_ block_length) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a26f66f9ce349aabc63227aa4f1cad6ab</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a8bfc42ee94b0068d1542c579718b52b4</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a390f34d93c2b8d9dd561c6807aa8eac5</anchor>
      <arglist>(Index_ block_start, Index_ block_length) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1Matrix.html</anchorfile>
      <anchor>a80f0867067f1b1e5a235a99ecf58924b</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::Options</name>
    <filename>structtatami_1_1Options.html</filename>
    <member kind="variable">
      <type>bool</type>
      <name>sparse_extract_index</name>
      <anchorfile>structtatami_1_1Options.html</anchorfile>
      <anchor>a192a47606233c7f5c76bb1e7a62867a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>sparse_extract_value</name>
      <anchorfile>structtatami_1_1Options.html</anchorfile>
      <anchor>abe1ec7b8bfc82728742ca4a496e59548</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>sparse_ordered_index</name>
      <anchorfile>structtatami_1_1Options.html</anchorfile>
      <anchor>a4d2e6f2631a7e813bfbf3f215594a79b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>cache_for_reuse</name>
      <anchorfile>structtatami_1_1Options.html</anchorfile>
      <anchor>a8eda7dd2e46293d3f3fee30bdbddb9c8</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::Oracle</name>
    <filename>structtatami_1_1Oracle.html</filename>
    <templarg>typename Index_</templarg>
    <member kind="function" virtualness="pure">
      <type>virtual size_t</type>
      <name>predict</name>
      <anchorfile>structtatami_1_1Oracle.html</anchorfile>
      <anchor>a32ab54d21fc982580f18bea802d89f48</anchor>
      <arglist>(Index_ *predicted, size_t number)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::OracleChunkCache</name>
    <filename>classtatami_1_1OracleChunkCache.html</filename>
    <templarg>typename Id_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class ChunkContents_</templarg>
    <member kind="function">
      <type></type>
      <name>OracleChunkCache</name>
      <anchorfile>classtatami_1_1OracleChunkCache.html</anchorfile>
      <anchor>a0bd9bfcb80bcaa3c92ccde442c6c15a9</anchor>
      <arglist>(std::unique_ptr&lt; Oracle&lt; Index_ &gt; &gt; oracle, size_t per_iteration, size_t num_chunks)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const ChunkContents_ *, Index_ &gt;</type>
      <name>next_chunk</name>
      <anchorfile>classtatami_1_1OracleChunkCache.html</anchorfile>
      <anchor>a85322f0675c67217a450e061880e0065</anchor>
      <arglist>(Ifunction_ identify, Sfunction_ swap, Rfunction_ ready, Afunction_ allocate, Pfunction_ populate)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::OracleStream</name>
    <filename>structtatami_1_1OracleStream.html</filename>
    <templarg>typename Index_</templarg>
    <member kind="function">
      <type></type>
      <name>OracleStream</name>
      <anchorfile>structtatami_1_1OracleStream.html</anchorfile>
      <anchor>a293c82567a7c94b9dc86e39110b30688</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>OracleStream</name>
      <anchorfile>structtatami_1_1OracleStream.html</anchorfile>
      <anchor>abac31a84479939fec89f7037a5dd8f06</anchor>
      <arglist>(std::unique_ptr&lt; Oracle&lt; Index_ &gt; &gt; o)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>structtatami_1_1OracleStream.html</anchorfile>
      <anchor>aa870c1926cb411c6de10be96956459ea</anchor>
      <arglist>(std::unique_ptr&lt; Oracle&lt; Index_ &gt; &gt; o)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>next</name>
      <anchorfile>structtatami_1_1OracleStream.html</anchorfile>
      <anchor>a370450d2c7194730cd61a8cc8f8ec390</anchor>
      <arglist>(Index_ &amp;prediction)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>back</name>
      <anchorfile>structtatami_1_1OracleStream.html</anchorfile>
      <anchor>a3ce855b3a8ffa7887b20f489faf1a68e</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::SemiCompressedSparseMatrix</name>
    <filename>classtatami_1_1SemiCompressedSparseMatrix.html</filename>
    <templarg>bool row_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class IndexStorage_</templarg>
    <templarg>class PointerStorage_</templarg>
    <base>tatami::Matrix</base>
    <member kind="function">
      <type></type>
      <name>SemiCompressedSparseMatrix</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>ab0a0579fc173803dc856faa141ea7024</anchor>
      <arglist>(Index_ nr, Index_ nc, IndexStorage_ idx, PointerStorage_ ptr, bool check=true)</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>nrow</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>a0a168f9099fd7ffd63867af6365fd829</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>ncol</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>a6b2a17e27bc6a56f08096b6fbd147735</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>a61a520975cbbdc711dd9e32681b11d49</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sparse_proportion</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>a849c19bd6cb418ebc09647b9b67084b9</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>prefer_rows</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>a74c4e8ac60e55b1d62aa159b85806467</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>prefer_rows_proportion</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>a2f603458319450865dc1226f05319050</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>uses_oracle</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>a6d059ba00d1970499fe34579c6ca5120</anchor>
      <arglist>(bool) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>ae7ca8ec0934562765ac468936ac61dbe</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>ae2164a57ed55013cff3e21fba1b28af4</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>a2f7652e1874cbeda36b47de5c09c4847</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>a58fe6fd533b511ff1764a0c4969b90e8</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>ae51cb25e91b0f54d780db17f04422c06</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>aadf33c3b1c3198f272e994cdefbcf60b</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>a5fb0d0da1c3b59fb9794d94e85938b9f</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>a60c6510f93cf1345a74d2eb04796de33</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>a79c1bb50d009f097453942f9eb8ce1d7</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>ae6876383fce85a907a4324595f8b0370</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>ada2224d2142ca6b0837741c0cf2fa9be</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1SemiCompressedSparseMatrix.html</anchorfile>
      <anchor>a25a6d306f43ef68e2b7945be91330f1f</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::SomeNumericArray</name>
    <filename>structtatami_1_1SomeNumericArray.html</filename>
    <templarg>typename T</templarg>
    <class kind="struct">tatami::SomeNumericArray::Iterator</class>
    <member kind="enumeration">
      <type></type>
      <name>Type</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>a6e23ba528d0302095fd423e7777455bd</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SomeNumericArray</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>ad027a9c8e476b4f8ebe8d93ec0f4c0b6</anchor>
      <arglist>(void *x, size_t n, Type t)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SomeNumericArray</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>ad1aca8120bd6bcd704137da9604516d6</anchor>
      <arglist>(const int8_t *x, size_t n)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SomeNumericArray</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>acbcd9f0a5f4f09c48600ce99afedcdd6</anchor>
      <arglist>(const uint8_t *x, size_t n)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SomeNumericArray</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>a4083f387fa053c1d66768ca5101b4185</anchor>
      <arglist>(const int16_t *x, size_t n)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SomeNumericArray</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>a61c4c9def9cd3364a44a842735f6b858</anchor>
      <arglist>(const uint16_t *x, size_t n)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SomeNumericArray</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>a245702fd8b44250f36d97fc9ec983d3b</anchor>
      <arglist>(const int32_t *x, size_t n)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SomeNumericArray</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>a6a06d411ceb7f22bd1cb276674a5c7e2</anchor>
      <arglist>(const uint32_t *x, size_t n)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SomeNumericArray</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>a3aaf11fe3d96df29e2e1a8d9705db9a4</anchor>
      <arglist>(const int64_t *x, size_t n)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SomeNumericArray</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>a961ad134b8a10f90e704036ee8b2505b</anchor>
      <arglist>(const uint64_t *x, size_t n)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SomeNumericArray</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>ab6a6d31ce3313e940d36cc6489b44391</anchor>
      <arglist>(const float *x, size_t n)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SomeNumericArray</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>a5db2f8e0c9f1c5d026682d95be81d87d</anchor>
      <arglist>(const double *x, size_t n)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator[]</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>ac40db52de47708abb156ad33f24619a4</anchor>
      <arglist>(size_t i) const</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>size</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>a6b23cda3aa1555f3605bdc425985e161</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Iterator</type>
      <name>begin</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>a5b7999f71425029d4b4de32c2025fdd9</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Iterator</type>
      <name>end</name>
      <anchorfile>structtatami_1_1SomeNumericArray.html</anchorfile>
      <anchor>aaf5d7921f2356a9a3aca745723611d71</anchor>
      <arglist>() const</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::SparseExtractor</name>
    <filename>classtatami_1_1SparseExtractor.html</filename>
    <templarg>DimensionSelectionType selection_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <member kind="function" virtualness="pure">
      <type>virtual SparseRange&lt; Value_, Index_ &gt;</type>
      <name>fetch</name>
      <anchorfile>classtatami_1_1SparseExtractor.html</anchorfile>
      <anchor>a721c2b7e0bf9958e122570dea88bb41c</anchor>
      <arglist>(Index_ i, Value_ *vbuffer, Index_ *ibuffer)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual SparseRange&lt; Value_, Index_ &gt;</type>
      <name>fetch_copy</name>
      <anchorfile>classtatami_1_1SparseExtractor.html</anchorfile>
      <anchor>aa60fb36ca94bc0fad1e83c7836ad2ce0</anchor>
      <arglist>(Index_ i, Value_ *vbuffer, Index_ *ibuffer)</arglist>
    </member>
    <member kind="function">
      <type>SparseRangeCopy&lt; Value_, Index_ &gt;</type>
      <name>fetch</name>
      <anchorfile>classtatami_1_1SparseExtractor.html</anchorfile>
      <anchor>a3a719d7fc7c8f75448d0f56b0aa67ec1</anchor>
      <arglist>(Index_ i)</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static constexpr bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1SparseExtractor.html</anchorfile>
      <anchor>af47a70fd33b7524f043bc3a797a3ed7d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static constexpr DimensionSelectionType</type>
      <name>selection</name>
      <anchorfile>classtatami_1_1SparseExtractor.html</anchorfile>
      <anchor>a44e2cb904dccbb9fb36104dd761a080a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::SparseRange</name>
    <filename>structtatami_1_1SparseRange.html</filename>
    <templarg>typename Value</templarg>
    <templarg>typename Index</templarg>
    <member kind="function">
      <type></type>
      <name>SparseRange</name>
      <anchorfile>structtatami_1_1SparseRange.html</anchorfile>
      <anchor>af2711e7bdfa4cd93da0255b1d8c95e7a</anchor>
      <arglist>(Index n, const Value *v=NULL, const Index *i=NULL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SparseRange</name>
      <anchorfile>structtatami_1_1SparseRange.html</anchorfile>
      <anchor>a469d0435e9a8da9818dbd7cbc122faf9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>Index</type>
      <name>number</name>
      <anchorfile>structtatami_1_1SparseRange.html</anchorfile>
      <anchor>a46230be86a913041c8c04cbe0fee732c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const Value *</type>
      <name>value</name>
      <anchorfile>structtatami_1_1SparseRange.html</anchorfile>
      <anchor>a6bd0dcad11850847dd8dd336e1932732</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const Index *</type>
      <name>index</name>
      <anchorfile>structtatami_1_1SparseRange.html</anchorfile>
      <anchor>ae187f7d169bbdc0701fe7414dd051717</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::SparseRangeCopy</name>
    <filename>structtatami_1_1SparseRangeCopy.html</filename>
    <templarg>typename Value</templarg>
    <templarg>typename Index</templarg>
    <member kind="function">
      <type></type>
      <name>SparseRangeCopy</name>
      <anchorfile>structtatami_1_1SparseRangeCopy.html</anchorfile>
      <anchor>a3ddd45e45bb70ea81655a5a5c55de4e7</anchor>
      <arglist>(Index n)</arglist>
    </member>
    <member kind="variable">
      <type>Index</type>
      <name>number</name>
      <anchorfile>structtatami_1_1SparseRangeCopy.html</anchorfile>
      <anchor>adab7f22070406ccdf552d15c60d88f7b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Value &gt;</type>
      <name>value</name>
      <anchorfile>structtatami_1_1SparseRangeCopy.html</anchorfile>
      <anchor>a6f924c14d5830aedf2d175b850099b45</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Index &gt;</type>
      <name>index</name>
      <anchorfile>structtatami_1_1SparseRangeCopy.html</anchorfile>
      <anchor>a62432f2387296b4070ce917e3c2eaac8</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::VirtualDenseMatrix</name>
    <filename>classtatami_1_1VirtualDenseMatrix.html</filename>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <base>tatami::Matrix</base>
    <member kind="function">
      <type>bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1VirtualDenseMatrix.html</anchorfile>
      <anchor>a10ca176c496b318d7c84a73000d309ba</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sparse_proportion</name>
      <anchorfile>classtatami_1_1VirtualDenseMatrix.html</anchorfile>
      <anchor>ad4a10daeac4ec28e4d8b256c9772f414</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1VirtualDenseMatrix.html</anchorfile>
      <anchor>a6d39bc3ff038786f058d6180f1369252</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1VirtualDenseMatrix.html</anchorfile>
      <anchor>a956946c98a7d06e917c78b5dcdc836ab</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1VirtualDenseMatrix.html</anchorfile>
      <anchor>a4c9eb67e9aaf1bd533b7ad3b9074ce6f</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1VirtualDenseMatrix.html</anchorfile>
      <anchor>a3b44e951cc4c966f7f561fddbec0cb35</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1VirtualDenseMatrix.html</anchorfile>
      <anchor>a8de5cef8cbf819e538f2289f5f52e6c4</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1VirtualDenseMatrix.html</anchorfile>
      <anchor>a4c16857ea52bde64900581af4db2da97</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>tatami</name>
    <filename>namespacetatami.html</filename>
    <namespace>tatami::MatrixMarket</namespace>
    <class kind="class">tatami::ArrayView</class>
    <class kind="struct">tatami::BlockExtractor</class>
    <class kind="class">tatami::CompressedSparseMatrix</class>
    <class kind="struct">tatami::ConsecutiveOracle</class>
    <class kind="struct">tatami::DelayedAbsHelper</class>
    <class kind="struct">tatami::DelayedAddScalarHelper</class>
    <class kind="struct">tatami::DelayedAddVectorHelper</class>
    <class kind="class">tatami::DelayedBind</class>
    <class kind="class">tatami::DelayedCast</class>
    <class kind="struct">tatami::DelayedDivideScalarHelper</class>
    <class kind="struct">tatami::DelayedDivideVectorHelper</class>
    <class kind="struct">tatami::DelayedExpHelper</class>
    <class kind="struct">tatami::DelayedLog1pHelper</class>
    <class kind="struct">tatami::DelayedLogHelper</class>
    <class kind="struct">tatami::DelayedMultiplyScalarHelper</class>
    <class kind="struct">tatami::DelayedMultiplyVectorHelper</class>
    <class kind="struct">tatami::DelayedRoundHelper</class>
    <class kind="struct">tatami::DelayedSqrtHelper</class>
    <class kind="class">tatami::DelayedSubset</class>
    <class kind="class">tatami::DelayedSubsetBlock</class>
    <class kind="class">tatami::DelayedSubsetSorted</class>
    <class kind="class">tatami::DelayedSubsetSortedUnique</class>
    <class kind="class">tatami::DelayedSubsetUnique</class>
    <class kind="struct">tatami::DelayedSubtractScalarHelper</class>
    <class kind="struct">tatami::DelayedSubtractVectorHelper</class>
    <class kind="class">tatami::DelayedTranspose</class>
    <class kind="class">tatami::DelayedUnaryIsometricOp</class>
    <class kind="class">tatami::DenseExtractor</class>
    <class kind="class">tatami::DenseMatrix</class>
    <class kind="struct">tatami::ExtractorBase</class>
    <class kind="struct">tatami::FixedOracle</class>
    <class kind="struct">tatami::FullExtractor</class>
    <class kind="struct">tatami::IndexExtractor</class>
    <class kind="struct">tatami::LayeredMatrixData</class>
    <class kind="class">tatami::LruChunkCache</class>
    <class kind="class">tatami::Matrix</class>
    <class kind="struct">tatami::Options</class>
    <class kind="struct">tatami::Oracle</class>
    <class kind="class">tatami::OracleChunkCache</class>
    <class kind="struct">tatami::OracleStream</class>
    <class kind="class">tatami::SemiCompressedSparseMatrix</class>
    <class kind="struct">tatami::SomeNumericArray</class>
    <class kind="class">tatami::SparseExtractor</class>
    <class kind="struct">tatami::SparseRange</class>
    <class kind="struct">tatami::SparseRangeCopy</class>
    <class kind="class">tatami::VirtualDenseMatrix</class>
    <member kind="typedef">
      <type>DenseMatrix&lt; false, Value_, Index_, Storage_ &gt;</type>
      <name>DenseColumnMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac47a769e00660eb7e9b5fcd543bcf2d3</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>DenseMatrix&lt; true, Value_, Index_, Storage_ &gt;</type>
      <name>DenseRowMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a51122d20490b377cd3f4609cc044f314</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>typename std::conditional&lt; selection_==DimensionSelectionType::FULL, FullExtractor&lt; Index_ &gt;, typename std::conditional&lt; selection_==DimensionSelectionType::BLOCK, BlockExtractor&lt; Index_ &gt;, IndexExtractor&lt; Index_ &gt; &gt;::type &gt;::type</type>
      <name>ConditionalSelectionExtractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a85ed61a4f772a2f7be4a12f739554e6e</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>typename std::conditional&lt; sparse_, SparseExtractor&lt; selection_, Value_, Index_ &gt;, DenseExtractor&lt; selection_, Value_, Index_ &gt; &gt;::type</type>
      <name>Extractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ae9f8db5316521603085577d977a6955f</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Extractor&lt; DimensionSelectionType::FULL, false, Value_, Index_ &gt;</type>
      <name>FullDenseExtractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a6008dbced6de41e5619156b5335f5762</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Extractor&lt; DimensionSelectionType::BLOCK, false, Value_, Index_ &gt;</type>
      <name>BlockDenseExtractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ae75de1fc78b7d361ea8b59a5379ea4da</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Extractor&lt; DimensionSelectionType::INDEX, false, Value_, Index_ &gt;</type>
      <name>IndexDenseExtractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a47ce406c32c3914c2ecce187e21b6ced</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Extractor&lt; DimensionSelectionType::FULL, true, Value_, Index_ &gt;</type>
      <name>FullSparseExtractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a0fbb0624c8e1913a87e8fb5c975400e1</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Extractor&lt; DimensionSelectionType::BLOCK, true, Value_, Index_ &gt;</type>
      <name>BlockSparseExtractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac8d0024399a66ce61f6315f5f46ebb63</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Extractor&lt; DimensionSelectionType::INDEX, true, Value_, Index_ &gt;</type>
      <name>IndexSparseExtractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a4b67b4d1b6c00cd0bd449703432a5f7b</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Matrix&lt; double, int &gt;</type>
      <name>NumericMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a35c670894994f1d620abb55953f98441</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>CompressedSparseMatrix&lt; false, Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_ &gt;</type>
      <name>CompressedSparseColumnMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a18cee3a5d9734f0092b03d023cfe4b6a</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>CompressedSparseMatrix&lt; true, Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_ &gt;</type>
      <name>CompressedSparseRowMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a0111adeeb583aeb7e24e9e1e25be4aa0</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>SemiCompressedSparseMatrix&lt; false, Value_, Index_, IndexStorage_, PointerStorage_ &gt;</type>
      <name>SemiCompressedSparseColumnMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a6f03d0d880bc056e09c2cbb80eb2c2ec</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>SemiCompressedSparseMatrix&lt; true, Value_, Index_, IndexStorage_, PointerStorage_ &gt;</type>
      <name>SemiCompressedSparseRowMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a412cb6ee12f3ee81d404d6eb0e494e4d</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumeration">
      <type></type>
      <name>DimensionSelectionType</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a0a2ecaf58e2b69bb4a808e814aeb16a1</anchor>
      <arglist></arglist>
      <enumvalue file="namespacetatami.html" anchor="a0a2ecaf58e2b69bb4a808e814aeb16a1aba7de5bc6888294e5884b024a4c894f1">FULL</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="a0a2ecaf58e2b69bb4a808e814aeb16a1a4d34f53389ed7f28ca91fc31ea360a66">BLOCK</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="a0a2ecaf58e2b69bb4a808e814aeb16a1acb4ae3b37047fb4b2c0d16f8bf84f076">INDEX</enumvalue>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>extracted_length</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>af9d13ceaa112d2c091265510d741488d</anchor>
      <arglist>(const ConditionalSelectionExtractor&lt; selection_, Index_ &gt; &amp;ex)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>make_DelayedBind</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a5f8e0f69139575707aa9314174b415b3</anchor>
      <arglist>(std::vector&lt; std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt; &gt; ps)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_out_, Index_out_ &gt; &gt;</type>
      <name>make_DelayedCast</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aca62c3bf751cdd06a08e8e503b0b591a</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_in_, Index_in_ &gt; &gt; p)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>make_DelayedTranspose</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>afa35d8e9fe286967f327ec0eb6bd5005</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>make_DelayedSubsetBlock</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a8f4aafc0a1fbdc0c31bc122d24122a63</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p, Index_ f, Index_ l)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>make_DelayedSubset</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>abfc63177b00e6e3e2fa47754b8d87704</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p, IndexStorage_ idx)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>new_extractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a95b6638dceda82d82a7579dc88a45709</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *ptr, Args_ &amp;&amp;... args)</arglist>
    </member>
    <member kind="function">
      <type>DelayedAddVectorHelper&lt; MARGIN, T, V &gt;</type>
      <name>make_DelayedAddVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aa82419090e54a4674d971429e93e6af5</anchor>
      <arglist>(V v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedSubtractVectorHelper&lt; RIGHT, MARGIN, T, V &gt;</type>
      <name>make_DelayedSubtractVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>abb68a8d27e821f213965edab7f029c48</anchor>
      <arglist>(V v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedMultiplyVectorHelper&lt; MARGIN, T, V &gt;</type>
      <name>make_DelayedMultiplyVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a316f4cf0ae5a9e83c2dad6ccbbbcfb43</anchor>
      <arglist>(V v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedDivideVectorHelper&lt; RIGHT, MARGIN, T, V &gt;</type>
      <name>make_DelayedDivideVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ab104e5029db4d09afe0b47a125f457bf</anchor>
      <arglist>(V v)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>make_DelayedUnaryIsometricOp</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a07b072944c766d85dc4ce4dc734e0b75</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p, Operation_ op)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>column_medians</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a44202790861791b1ed9df6d480c69f6a</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>row_medians</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ab877f67694b3fa500f59469f044b1c01</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>column_maxs</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a8102a2423efab5e329543f4235ac2290</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>row_maxs</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a47890bfd538f65971f9e5dba5e1ed785</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>column_mins</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a60e44772f82183100c7c25c181479b56</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>row_mins</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a3bcdcb0499b36b5b844ebd88fc30eefe</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; std::vector&lt; Output_ &gt;, std::vector&lt; Output_ &gt; &gt;</type>
      <name>column_ranges</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ab5eff79824930c458e730d23f640c1ab</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; std::vector&lt; Output_ &gt;, std::vector&lt; Output_ &gt; &gt;</type>
      <name>row_ranges</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a2b2475a166ad744b31b1275e26e562be</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>column_sums</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a3aab6733d637b66abdf45f7300d2e1ba</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>row_sums</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a53a658059404691856bef57fb85d83d6</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>parallelize</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a29ce7a2219ea60d45de1aa3d4de66063</anchor>
      <arglist>(Function_ fun, size_t tasks, size_t threads)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>consecutive_extractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a36c6ecf33bcb87e1ed33c0a7d744dd82</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *mat, Index_ iter_start, Index_ iter_length, Args_ &amp;&amp;... args)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>column_variances</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac37af9a76d15f08a3634881696a27f65</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Output_ &gt;</type>
      <name>row_variances</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ab806279616e19f6150376a7c07d5b64b</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *p, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; std::shared_ptr&lt; Matrix &gt;, std::vector&lt; size_t &gt; &gt;</type>
      <name>bind_intersection</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a87c52cf8a974ad0025550b002d391d39</anchor>
      <arglist>(const std::vector&lt; std::shared_ptr&lt; Matrix &gt; &gt; &amp;inputs, const std::vector&lt; const Id * &gt; &amp;ids)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; size_t &gt;</type>
      <name>compress_sparse_triplets</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ab17e92414b0bff60f7b7a6431ac8a330</anchor>
      <arglist>(size_t nr, size_t nc, U &amp;values, V &amp;rows, W &amp;cols)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>convert_to_dense</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>afaccb18d22f35a4ed6b6a7091342da06</anchor>
      <arglist>(const Matrix_ *incoming, StoredValue_ *store, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index &gt; &gt;</type>
      <name>convert_to_dense</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>af80143ac537339fe8dafd892632e96de</anchor>
      <arglist>(const Matrix_ *incoming, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_dense</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a64e8745f8aaaa5ad93c6ce65a0a6591e</anchor>
      <arglist>(const Matrix_ *incoming, int order, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_sparse</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a220c47e470b3ade64ede7303e50ca92d</anchor>
      <arglist>(const InputMatrix_ *incoming, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_sparse</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a5847f3cd78ce89170ca43262d7f5f8e3</anchor>
      <arglist>(const InputMatrix_ *incoming, int order, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>process_consecutive_indices</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>af01c93a616eb99f6a17861a8b19f7ee0</anchor>
      <arglist>(const Index_ *indices, Index_ length, Function_ fun)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; const Matrix&lt; T, IDX &gt; &gt;</type>
      <name>wrap_shared_ptr</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a7518f5e8e09a6f6d7d3955b8ea286689</anchor>
      <arglist>(const Matrix&lt; T, IDX &gt; *ptr)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>tatami::MatrixMarket</name>
    <filename>namespacetatami_1_1MatrixMarket.html</filename>
    <class kind="struct">tatami::MatrixMarket::HeaderDetails</class>
    <member kind="function">
      <type>LayeredMatrixData&lt; T, IDX &gt;</type>
      <name>load_layered_sparse_matrix_from_file</name>
      <anchorfile>namespacetatami_1_1MatrixMarket.html</anchorfile>
      <anchor>a1957591502ac286ad6b1bcae78ee75f9</anchor>
      <arglist>(const char *filepath, int compression=0, size_t bufsize=65536)</arglist>
    </member>
    <member kind="function">
      <type>LayeredMatrixData&lt; T, IDX &gt;</type>
      <name>load_layered_sparse_matrix_from_buffer</name>
      <anchorfile>namespacetatami_1_1MatrixMarket.html</anchorfile>
      <anchor>af7655515516f20bf22da4311b7e8dcc5</anchor>
      <arglist>(const unsigned char *buffer, size_t n, int compression=0, size_t bufsize=65536)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; T, IDX &gt; &gt;</type>
      <name>load_sparse_matrix_from_file</name>
      <anchorfile>namespacetatami_1_1MatrixMarket.html</anchorfile>
      <anchor>acb8ab6a2dba7c9b9088245792d8dd4c5</anchor>
      <arglist>(const char *filepath, int compression=0, size_t bufsize=65536)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; T, IDX &gt; &gt;</type>
      <name>load_sparse_matrix_from_buffer</name>
      <anchorfile>namespacetatami_1_1MatrixMarket.html</anchorfile>
      <anchor>abca83cf0a5a51b93f78991719c8bd772</anchor>
      <arglist>(const unsigned char *buffer, size_t n, int compression=0, size_t bufsize=65536)</arglist>
    </member>
    <member kind="function">
      <type>HeaderDetails</type>
      <name>extract_header_from_file</name>
      <anchorfile>namespacetatami_1_1MatrixMarket.html</anchorfile>
      <anchor>afc195d253d64d0a67216cfd6d52e8639</anchor>
      <arglist>(const char *filepath, int compression=0, size_t bufsize=65536)</arglist>
    </member>
    <member kind="function">
      <type>HeaderDetails</type>
      <name>extract_header_from_buffer</name>
      <anchorfile>namespacetatami_1_1MatrixMarket.html</anchorfile>
      <anchor>ae0ced3181c20a3281a2ab4844b978313</anchor>
      <arglist>(const unsigned char *buffer, size_t n, int compression=0, size_t bufsize=65536)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>A C++ API for all sorts of matrices</title>
    <filename>index.html</filename>
    <docanchor file="index.html" title="A C++ API for all sorts of matrices">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
