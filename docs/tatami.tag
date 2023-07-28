<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.9.5">
  <compound kind="file">
    <name>Extractor.hpp</name>
    <path>/github/workspace/include/tatami/base/</path>
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
    <path>/github/workspace/include/tatami/base/</path>
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
    <path>/github/workspace/include/tatami/base/</path>
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
    <name>SparseRange.hpp</name>
    <path>/github/workspace/include/tatami/base/</path>
    <filename>SparseRange_8hpp.html</filename>
    <class kind="struct">tatami::SparseRange</class>
    <class kind="struct">tatami::SparseRangeCopy</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>convert_to_dense.hpp</name>
    <path>/github/workspace/include/tatami/dense/</path>
    <filename>dense_2convert__to__dense_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>void</type>
      <name>convert_to_dense</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a3d709db6be55e94d987d38c4c07c71c0</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, StoredValue_ *store, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index &gt; &gt;</type>
      <name>convert_to_dense</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a3a0e048d06c98ea3706466139c4c97dc</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_dense</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aa74ad14ba410177396121d272371dffd</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, int order, int threads=1)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>DenseMatrix.hpp</name>
    <path>/github/workspace/include/tatami/dense/</path>
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
    <path>/github/workspace/include/tatami/dense/</path>
    <filename>VirtualDenseMatrix_8hpp.html</filename>
    <class kind="class">tatami::VirtualDenseMatrix</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>arith_utils.hpp</name>
    <path>/github/workspace/include/tatami/isometric/</path>
    <filename>arith__utils_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="enumeration">
      <type></type>
      <name>DelayedArithOp</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aab44a37b3762de0c5b1ffbfceb25fa0f</anchor>
      <arglist></arglist>
      <enumvalue file="namespacetatami.html" anchor="aab44a37b3762de0c5b1ffbfceb25fa0fa9eeb52badb613229884838847294b90d">ADD</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="aab44a37b3762de0c5b1ffbfceb25fa0fa23ebcc4776b613af25dfbe7c8ce4813e">SUBTRACT</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="aab44a37b3762de0c5b1ffbfceb25fa0fa080aaf8d817ada96fca7096b7b55bd30">MULTIPLY</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="aab44a37b3762de0c5b1ffbfceb25fa0fa210c66d794cec40488f3f8f634d6c33b">DIVIDE</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="aab44a37b3762de0c5b1ffbfceb25fa0fac9c9c146c630ca5ef9197c73c032f4a6">POWER</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="aab44a37b3762de0c5b1ffbfceb25fa0fa928ab45d616dde447dbbbd0270db87ad">MODULO</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="aab44a37b3762de0c5b1ffbfceb25fa0fa051460a4a75d4d251a41a7c04bf49412">INTEGER_DIVIDE</enumvalue>
    </member>
  </compound>
  <compound kind="file">
    <name>DelayedBinaryIsometricOp.hpp</name>
    <path>/github/workspace/include/tatami/isometric/binary/</path>
    <filename>DelayedBinaryIsometricOp_8hpp.html</filename>
    <class kind="class">tatami::DelayedBinaryIsometricOp</class>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>make_DelayedBinaryIsometricOp</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>af39c0767d23c6bdc5f4267f062ef2e17</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; left, std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; right, Operation_ op)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>boolean_utils.hpp</name>
    <path>/github/workspace/include/tatami/isometric/</path>
    <filename>boolean__utils_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="enumeration">
      <type></type>
      <name>DelayedBooleanOp</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a2104862d4068933ea4cc805c92f82d07</anchor>
      <arglist></arglist>
      <enumvalue file="namespacetatami.html" anchor="a2104862d4068933ea4cc805c92f82d07a558ffc8f5770d8e4f95f51d822685532">AND</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="a2104862d4068933ea4cc805c92f82d07a1d00e7dce692e8dc3f6877f035e3a616">OR</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="a2104862d4068933ea4cc805c92f82d07a97675eb3f268048604dc5155511a2a4d">XOR</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="a2104862d4068933ea4cc805c92f82d07a969f331a87d8c958473c32b4d0e61a44">EQUAL</enumvalue>
    </member>
  </compound>
  <compound kind="file">
    <name>compare_utils.hpp</name>
    <path>/github/workspace/include/tatami/isometric/</path>
    <filename>compare__utils_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="enumeration">
      <type></type>
      <name>DelayedCompareOp</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac4fc175a57ace709941b5ca7ddb19708</anchor>
      <arglist></arglist>
      <enumvalue file="namespacetatami.html" anchor="ac4fc175a57ace709941b5ca7ddb19708a969f331a87d8c958473c32b4d0e61a44">EQUAL</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="ac4fc175a57ace709941b5ca7ddb19708a1625ef4fe09f68fa20d3ff6e02cd5c8e">GREATER_THAN</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="ac4fc175a57ace709941b5ca7ddb19708aa327176a0a845c117bdfadec134a95e9">LESS_THAN</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="ac4fc175a57ace709941b5ca7ddb19708aa6eac69202c3dc2978176801a84e4d1d">GREATER_THAN_OR_EQUAL</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="ac4fc175a57ace709941b5ca7ddb19708a8397780541b6289d2a0b991d1c28c432">LESS_THAN_OR_EQUAL</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="ac4fc175a57ace709941b5ca7ddb19708a4ea2d378cdec20f59330f113297bc1ce">NOT_EQUAL</enumvalue>
    </member>
  </compound>
  <compound kind="file">
    <name>arith_helpers.hpp</name>
    <path>/github/workspace/include/tatami/isometric/binary/</path>
    <filename>binary_2arith__helpers_8hpp.html</filename>
    <class kind="struct">tatami::DelayedBinaryArithHelper</class>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>DelayedBinaryArithHelper&lt; DelayedArithOp::ADD &gt;</type>
      <name>make_DelayedBinaryAddHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a1d138b1e7f6a26f814c025363cc3db80</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryArithHelper&lt; DelayedArithOp::SUBTRACT &gt;</type>
      <name>make_DelayedBinarySubtractHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac4c951d489cb0f2bb0e90f2cf4a25862</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryArithHelper&lt; DelayedArithOp::MULTIPLY &gt;</type>
      <name>make_DelayedBinaryMultiplyHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a2569c540083a24f92af8140358e1e9c4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryArithHelper&lt; DelayedArithOp::DIVIDE &gt;</type>
      <name>make_DelayedBinaryDivideHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aeb4cf766f850766b966f7121728a6af8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryArithHelper&lt; DelayedArithOp::POWER &gt;</type>
      <name>make_DelayedBinaryPowerHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ab6c73c04e7b08130ed7b46d93c4dfd11</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryArithHelper&lt; DelayedArithOp::MODULO &gt;</type>
      <name>make_DelayedBinaryModuloHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a04293eb1e8eefb7024cc192d75ac093e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryArithHelper&lt; DelayedArithOp::INTEGER_DIVIDE &gt;</type>
      <name>make_DelayedBinaryIntegerDivideHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aa8f54741424bef6a21225935a45e9d53</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>arith_helpers.hpp</name>
    <path>/github/workspace/include/tatami/isometric/unary/</path>
    <filename>unary_2arith__helpers_8hpp.html</filename>
    <class kind="struct">tatami::DelayedArithScalarHelper</class>
    <class kind="struct">tatami::DelayedArithVectorHelper</class>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>DelayedArithScalarHelper&lt; DelayedArithOp::ADD, true, Value_, Scalar_ &gt;</type>
      <name>make_DelayedAddScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a33b2b9798f2a43b62821131ba3a3f6bd</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithScalarHelper&lt; DelayedArithOp::SUBTRACT, right_, Value_, Scalar_ &gt;</type>
      <name>make_DelayedSubtractScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a9bdfe26c0f1426b611e91ee7569045de</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithScalarHelper&lt; DelayedArithOp::MULTIPLY, true, Value_, Scalar_ &gt;</type>
      <name>make_DelayedMultiplyScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aa653b51960498cc467ea69c18ea0c097</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithScalarHelper&lt; DelayedArithOp::DIVIDE, right_, Value_, Scalar_ &gt;</type>
      <name>make_DelayedDivideScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a60076096ee4296c0e187a226be5ca4b5</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithScalarHelper&lt; DelayedArithOp::POWER, right_, Value_, Scalar_ &gt;</type>
      <name>make_DelayedPowerScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a8c478145f87c37a9fe6d4e8490fbf05c</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithScalarHelper&lt; DelayedArithOp::MODULO, right_, Value_, Scalar_ &gt;</type>
      <name>make_DelayedModuloScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a5504f584d6db28c05d78f1917c40810d</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithScalarHelper&lt; DelayedArithOp::INTEGER_DIVIDE, right_, Value_, Scalar_ &gt;</type>
      <name>make_DelayedIntegerDivideScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a5ec92ba88a0f5e9d7bf52021ed28e859</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithVectorHelper&lt; DelayedArithOp::ADD, true, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedAddVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a628c631be4424b37c1c1e68e28c5e982</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithVectorHelper&lt; DelayedArithOp::SUBTRACT, right_, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedSubtractVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a800b92854282015d96308fb283aeb508</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithVectorHelper&lt; DelayedArithOp::MULTIPLY, true, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedMultiplyVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ab929acdeb634187cd08d1100bbbf1b29</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithVectorHelper&lt; DelayedArithOp::DIVIDE, right_, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedDivideVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a1a4629f17aa4f3a4be06f3479bb6f68f</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithVectorHelper&lt; DelayedArithOp::POWER, right_, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedPowerVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aafe489bd754c040491ebc753a4e5656e</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithVectorHelper&lt; DelayedArithOp::MODULO, right_, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedModuloVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a4220ce5709e46ebce77b240b572d97d9</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithVectorHelper&lt; DelayedArithOp::INTEGER_DIVIDE, right_, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedIntegerDivideVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a857fc6c11ec69895ac85991ac83c395d</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>boolean_helpers.hpp</name>
    <path>/github/workspace/include/tatami/isometric/binary/</path>
    <filename>binary_2boolean__helpers_8hpp.html</filename>
    <class kind="struct">tatami::DelayedBinaryBooleanHelper</class>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>DelayedBinaryBooleanHelper&lt; DelayedBooleanOp::EQUAL &gt;</type>
      <name>make_DelayedBinaryBooleanEqualHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a1063aea86897ef76e17b1772320b8f7d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryBooleanHelper&lt; DelayedBooleanOp::AND &gt;</type>
      <name>make_DelayedBinaryBooleanAndHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a31625dbb50d420fe1b0da406f19ef33d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryBooleanHelper&lt; DelayedBooleanOp::OR &gt;</type>
      <name>make_DelayedBinaryBooleanOrHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a65bad100b39d9372fbeb23d16dae6588</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryBooleanHelper&lt; DelayedBooleanOp::XOR &gt;</type>
      <name>make_DelayedBinaryBooleanXorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aeb091953ab0f55406935a52b7ecb7350</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>boolean_helpers.hpp</name>
    <path>/github/workspace/include/tatami/isometric/unary/</path>
    <filename>unary_2boolean__helpers_8hpp.html</filename>
    <class kind="struct">tatami::DelayedBooleanScalarHelper</class>
    <class kind="struct">tatami::DelayedBooleanNotHelper</class>
    <class kind="struct">tatami::DelayedBooleanVectorHelper</class>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>DelayedBooleanNotHelper&lt; Value_ &gt;</type>
      <name>make_DelayedBooleanNotHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac812e757b55c2145bc3d4cd113c712fd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanScalarHelper&lt; DelayedBooleanOp::AND, Value_ &gt;</type>
      <name>make_DelayedBooleanAndScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ae460e24124cef43d8a5d3b82e8f02325</anchor>
      <arglist>(bool s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanScalarHelper&lt; DelayedBooleanOp::OR &gt;</type>
      <name>make_DelayedBooleanOrScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac7dbc64e6b0a0c8c48b1423b4b6748f1</anchor>
      <arglist>(bool s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanScalarHelper&lt; DelayedBooleanOp::XOR &gt;</type>
      <name>make_DelayedBooleanXorScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a29f6f1b1d28c670f38fda651fe2aa03a</anchor>
      <arglist>(bool s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanScalarHelper&lt; DelayedBooleanOp::EQUAL &gt;</type>
      <name>make_DelayedBooleanEqualScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a6bc5d47e0c391c17c85bddf19cc5620a</anchor>
      <arglist>(bool s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanVectorHelper&lt; DelayedBooleanOp::AND, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedBooleanAndVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a1fd06669afab8565747e08bbce4b041f</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanVectorHelper&lt; DelayedBooleanOp::OR, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedBooleanOrVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a503f23fe3e5a2c5f5bc6c0db01a7f64d</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanVectorHelper&lt; DelayedBooleanOp::XOR, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedBooleanXorVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a4bcdace63565fc670b0393675e33bb11</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanVectorHelper&lt; DelayedBooleanOp::EQUAL, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedBooleanEqualVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a8b07dc8b7a71776bcc222d7c437e2ca3</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>compare_helpers.hpp</name>
    <path>/github/workspace/include/tatami/isometric/binary/</path>
    <filename>binary_2compare__helpers_8hpp.html</filename>
    <class kind="struct">tatami::DelayedBinaryCompareHelper</class>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>DelayedBinaryCompareHelper&lt; DelayedCompareOp::EQUAL &gt;</type>
      <name>make_DelayedBinaryEqualHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac3599db8f6d19fc9172da45a3b6244fc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryCompareHelper&lt; DelayedCompareOp::GREATER_THAN &gt;</type>
      <name>make_DelayedBinaryGreaterThanHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>abe4fb79c5a2cddab2f5c9cf0bd9fe3d0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryCompareHelper&lt; DelayedCompareOp::LESS_THAN &gt;</type>
      <name>make_DelayedBinaryLessThanHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ae6a24f5bf2c2c61382aef43f2c8d3362</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryCompareHelper&lt; DelayedCompareOp::GREATER_THAN_OR_EQUAL &gt;</type>
      <name>make_DelayedBinaryGreaterThanOrEqualHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ad31cd299427d73538dc00dcbffa9cd2b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryCompareHelper&lt; DelayedCompareOp::LESS_THAN_OR_EQUAL &gt;</type>
      <name>make_DelayedBinaryLessThanOrEqualHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a38c7cc26409e5c995b7264164b148f4c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryCompareHelper&lt; DelayedCompareOp::NOT_EQUAL &gt;</type>
      <name>make_DelayedBinaryNotEqualHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a3073c4e33dcf89416af292399ba55c29</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>compare_helpers.hpp</name>
    <path>/github/workspace/include/tatami/isometric/unary/</path>
    <filename>unary_2compare__helpers_8hpp.html</filename>
    <class kind="struct">tatami::DelayedCompareScalarHelper</class>
    <class kind="struct">tatami::DelayedCompareVectorHelper</class>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>DelayedCompareScalarHelper&lt; DelayedCompareOp::EQUAL, Value_, Scalar_ &gt;</type>
      <name>make_DelayedEqualScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a11bfea31e4524927274dc1d80c32fd96</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareScalarHelper&lt; DelayedCompareOp::GREATER_THAN, Value_, Scalar_ &gt;</type>
      <name>make_DelayedGreaterThanScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a1961ba2a59d705adf7b8627523a62b3a</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareScalarHelper&lt; DelayedCompareOp::LESS_THAN, Value_, Scalar_ &gt;</type>
      <name>make_DelayedLessThanScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a7d2d381ec9c0d87f5d561b4dd9d5a58e</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareScalarHelper&lt; DelayedCompareOp::GREATER_THAN_OR_EQUAL, Value_, Scalar_ &gt;</type>
      <name>make_DelayedGreaterThanOrEqualScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac3df12ccdf4692c0f660f1634a79b857</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareScalarHelper&lt; DelayedCompareOp::LESS_THAN_OR_EQUAL, Value_, Scalar_ &gt;</type>
      <name>make_DelayedLessThanOrEqualScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>af5e67b02534cbe1813be5e2d690af88b</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareScalarHelper&lt; DelayedCompareOp::NOT_EQUAL, Value_, Scalar_ &gt;</type>
      <name>make_DelayedNotEqualScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac5d5cef56d2119e73e3f670d49de6b48</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareVectorHelper&lt; DelayedCompareOp::EQUAL, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedEqualVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aa40c54d5c5cc1109d42a150f5cd7ed2f</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareVectorHelper&lt; DelayedCompareOp::GREATER_THAN, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedGreaterThanVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aaf1951e22415be4714495556f701a7fd</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareVectorHelper&lt; DelayedCompareOp::LESS_THAN, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedLessThanVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a6dca8b9072813f6e9e52aa0536817b82</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareVectorHelper&lt; DelayedCompareOp::GREATER_THAN_OR_EQUAL, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedGreaterThanOrEqualVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a70a3e0489fdaccc80c7cc0125e9cf906</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareVectorHelper&lt; DelayedCompareOp::LESS_THAN_OR_EQUAL, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedLessThanOrEqualVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a83999d2ac877c5b15f6fd5a0f33a894a</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareVectorHelper&lt; DelayedCompareOp::NOT_EQUAL, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedNotEqualVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a6c84e8b757cf5ba6479a6ea3e6e530e3</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>DelayedUnaryIsometricOp.hpp</name>
    <path>/github/workspace/include/tatami/isometric/unary/</path>
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
    <path>/github/workspace/include/tatami/isometric/unary/</path>
    <filename>math__helpers_8hpp.html</filename>
    <class kind="struct">tatami::DelayedAbsHelper</class>
    <class kind="struct">tatami::DelayedSignHelper</class>
    <class kind="struct">tatami::DelayedLogHelper</class>
    <class kind="struct">tatami::DelayedSqrtHelper</class>
    <class kind="struct">tatami::DelayedCeilingHelper</class>
    <class kind="struct">tatami::DelayedFloorHelper</class>
    <class kind="struct">tatami::DelayedTruncHelper</class>
    <class kind="struct">tatami::DelayedLog1pHelper</class>
    <class kind="struct">tatami::DelayedRoundHelper</class>
    <class kind="struct">tatami::DelayedExpHelper</class>
    <class kind="struct">tatami::DelayedExpm1Helper</class>
    <class kind="struct">tatami::DelayedAcosHelper</class>
    <class kind="struct">tatami::DelayedAcoshHelper</class>
    <class kind="struct">tatami::DelayedAsinHelper</class>
    <class kind="struct">tatami::DelayedAsinhHelper</class>
    <class kind="struct">tatami::DelayedAtanHelper</class>
    <class kind="struct">tatami::DelayedAtanhHelper</class>
    <class kind="struct">tatami::DelayedCosHelper</class>
    <class kind="struct">tatami::DelayedCoshHelper</class>
    <class kind="struct">tatami::DelayedSinHelper</class>
    <class kind="struct">tatami::DelayedSinhHelper</class>
    <class kind="struct">tatami::DelayedTanHelper</class>
    <class kind="struct">tatami::DelayedTanhHelper</class>
    <class kind="struct">tatami::DelayedGammaHelper</class>
    <class kind="struct">tatami::DelayedLgammaHelper</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>DelayedBind.hpp</name>
    <path>/github/workspace/include/tatami/other/</path>
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
    <path>/github/workspace/include/tatami/other/</path>
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
    <path>/github/workspace/include/tatami/other/</path>
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
    <path>/github/workspace/include/tatami/sparse/</path>
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
    <name>convert_to_compressed_sparse.hpp</name>
    <path>/github/workspace/include/tatami/sparse/</path>
    <filename>convert__to__compressed__sparse_8hpp.html</filename>
    <class kind="struct">tatami::CompressedSparseContents</class>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>CompressedSparseContents&lt; Value_, Index_ &gt;</type>
      <name>retrieve_compressed_sparse_contents</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ad1aee6313d02036d2f56703f372e1df6</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, bool two_pass, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_compressed_sparse</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a48dab0e55d3b852afcfd9e93835190df</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, bool two_pass=false, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_compressed_sparse</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>acf7ca799fac9cff0c16a77bf98e32ff3</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, int order, bool two_pass=false, int threads=1)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>convert_to_fragmented_sparse.hpp</name>
    <path>/github/workspace/include/tatami/sparse/</path>
    <filename>convert__to__fragmented__sparse_8hpp.html</filename>
    <class kind="struct">tatami::FragmentedSparseContents</class>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>FragmentedSparseContents&lt; Value_, Index_ &gt;</type>
      <name>retrieve_fragmented_sparse_contents</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a012194aba94c37703f15555247e02155</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_fragmented_sparse</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a953ad06ecc9a36c8d2c1f3febe84a787</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_fragmented_sparse</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ae46babb9d985cffcb54bcba89553a811</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, int order, int threads=1)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>FragmentedSparseMatrix.hpp</name>
    <path>/github/workspace/include/tatami/sparse/</path>
    <filename>FragmentedSparseMatrix_8hpp.html</filename>
    <class kind="class">tatami::FragmentedSparseMatrix</class>
    <namespace>tatami</namespace>
    <member kind="typedef">
      <type>FragmentedSparseMatrix&lt; false, Value_, Index_, ValueVectorStorage_, IndexVectorStorage_ &gt;</type>
      <name>FragmentedSparseColumnMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a6b9735bdc6e0a2856cd4747efbbf9f3c</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>FragmentedSparseMatrix&lt; true, Value_, Index_, ValueVectorStorage_, IndexVectorStorage_ &gt;</type>
      <name>FragmentedSparseRowMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a3d97fc9bd75052a79bb643d0fb5a0b98</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>SemiCompressedSparseMatrix.hpp</name>
    <path>/github/workspace/include/tatami/sparse/</path>
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
    <name>medians.hpp</name>
    <path>/github/workspace/include/tatami/stats/</path>
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
    <path>/github/workspace/include/tatami/stats/</path>
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
    <path>/github/workspace/include/tatami/stats/</path>
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
    <name>utils.hpp</name>
    <path>/github/workspace/include/tatami/stats/</path>
    <filename>stats_2utils_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="function">
      <type>void</type>
      <name>parallelize</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ad2eddeaad50bdad0bc1e5c8f2f8549e0</anchor>
      <arglist>(Function_ fun, Index_ tasks, size_t #if defined(_OPENMP)||defined(TATAMI_CUSTOM_PARALLEL) threads #endif)</arglist>
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
    <name>variances.hpp</name>
    <path>/github/workspace/include/tatami/stats/</path>
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
    <name>DelayedSubset.hpp</name>
    <path>/github/workspace/include/tatami/subset/</path>
    <filename>DelayedSubset_8hpp.html</filename>
    <class kind="class">tatami::DelayedSubset</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>DelayedSubsetBlock.hpp</name>
    <path>/github/workspace/include/tatami/subset/</path>
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
    <path>/github/workspace/include/tatami/subset/</path>
    <filename>DelayedSubsetSorted_8hpp.html</filename>
    <class kind="class">tatami::DelayedSubsetSorted</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>DelayedSubsetSortedUnique.hpp</name>
    <path>/github/workspace/include/tatami/subset/</path>
    <filename>DelayedSubsetSortedUnique_8hpp.html</filename>
    <class kind="class">tatami::DelayedSubsetSortedUnique</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>DelayedSubsetUnique.hpp</name>
    <path>/github/workspace/include/tatami/subset/</path>
    <filename>DelayedSubsetUnique_8hpp.html</filename>
    <class kind="class">tatami::DelayedSubsetUnique</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>make_DelayedSubset.hpp</name>
    <path>/github/workspace/include/tatami/subset/</path>
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
    <name>tatami.hpp</name>
    <path>/github/workspace/include/tatami/</path>
    <filename>tatami_8hpp.html</filename>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>ArrayView.hpp</name>
    <path>/github/workspace/include/tatami/utils/</path>
    <filename>ArrayView_8hpp.html</filename>
    <class kind="class">tatami::ArrayView</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>bind_intersection.hpp</name>
    <path>/github/workspace/include/tatami/utils/</path>
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
    <path>/github/workspace/include/tatami/utils/</path>
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
    <name>ElementType.hpp</name>
    <path>/github/workspace/include/tatami/utils/</path>
    <filename>ElementType_8hpp.html</filename>
    <namespace>tatami</namespace>
    <member kind="typedef">
      <type>typename std::remove_cv&lt; typename std::remove_reference&lt; decltype(std::declval&lt; Array_ &gt;()[0])&gt;::type &gt;::type</type>
      <name>ElementType</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a58a028d23a7be58854b2e60dfae1b04b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Oracles.hpp</name>
    <path>/github/workspace/include/tatami/utils/</path>
    <filename>Oracles_8hpp.html</filename>
    <class kind="struct">tatami::FixedOracle</class>
    <class kind="struct">tatami::ConsecutiveOracle</class>
    <class kind="struct">tatami::OracleStream</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>process_consecutive_indices.hpp</name>
    <path>/github/workspace/include/tatami/utils/</path>
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
    <path>/github/workspace/include/tatami/utils/</path>
    <filename>SomeNumericArray_8hpp.html</filename>
    <class kind="struct">tatami::SomeNumericArray</class>
    <class kind="struct">tatami::SomeNumericArray::Iterator</class>
    <namespace>tatami</namespace>
  </compound>
  <compound kind="file">
    <name>wrap_shared_ptr.hpp</name>
    <path>/github/workspace/include/tatami/utils/</path>
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
  <compound kind="struct">
    <name>tatami::CompressedSparseContents</name>
    <filename>structtatami_1_1CompressedSparseContents.html</filename>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <member kind="variable">
      <type>std::vector&lt; Value_ &gt;</type>
      <name>value</name>
      <anchorfile>structtatami_1_1CompressedSparseContents.html</anchorfile>
      <anchor>a42f014b8d8da15e6c49511db087bebf7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Index_ &gt;</type>
      <name>index</name>
      <anchorfile>structtatami_1_1CompressedSparseContents.html</anchorfile>
      <anchor>ad6aa16e860a82ebc6de6774e8b6905ea</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; size_t &gt;</type>
      <name>pointers</name>
      <anchorfile>structtatami_1_1CompressedSparseContents.html</anchorfile>
      <anchor>ad56fc74d757e9dcd8ffbdc48d5adc692</anchor>
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
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedAcosHelper</name>
    <filename>structtatami_1_1DelayedAcosHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedAcoshHelper</name>
    <filename>structtatami_1_1DelayedAcoshHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedArithScalarHelper</name>
    <filename>structtatami_1_1DelayedArithScalarHelper.html</filename>
    <templarg>DelayedArithOp op_</templarg>
    <templarg>bool right_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Scalar_</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedArithScalarHelper</name>
      <anchorfile>structtatami_1_1DelayedArithScalarHelper.html</anchorfile>
      <anchor>ad49bfd018e4e05434b52bc3eca1563b1</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedArithVectorHelper</name>
    <filename>structtatami_1_1DelayedArithVectorHelper.html</filename>
    <templarg>DelayedArithOp op_</templarg>
    <templarg>bool right_</templarg>
    <templarg>int margin_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Vector_</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedArithVectorHelper</name>
      <anchorfile>structtatami_1_1DelayedArithVectorHelper.html</anchorfile>
      <anchor>af5554f1c599eeabfa01f630da1e7004e</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedAsinHelper</name>
    <filename>structtatami_1_1DelayedAsinHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedAsinhHelper</name>
    <filename>structtatami_1_1DelayedAsinhHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedAtanHelper</name>
    <filename>structtatami_1_1DelayedAtanHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedAtanhHelper</name>
    <filename>structtatami_1_1DelayedAtanhHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedBinaryArithHelper</name>
    <filename>structtatami_1_1DelayedBinaryArithHelper.html</filename>
    <templarg>DelayedArithOp op_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedBinaryBooleanHelper</name>
    <filename>structtatami_1_1DelayedBinaryBooleanHelper.html</filename>
    <templarg>DelayedBooleanOp op_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedBinaryCompareHelper</name>
    <filename>structtatami_1_1DelayedBinaryCompareHelper.html</filename>
    <templarg>DelayedCompareOp op_</templarg>
  </compound>
  <compound kind="class">
    <name>tatami::DelayedBinaryIsometricOp</name>
    <filename>classtatami_1_1DelayedBinaryIsometricOp.html</filename>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class Operation_</templarg>
    <base>Matrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>DelayedBinaryIsometricOp</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>af35faa9b83b2534f4b9d4e6e118b5728</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; l, std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; r, Operation_ op)</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>nrow</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>ae948f15be02f19e00dd61056ea958a83</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>ncol</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>a28043844566505282c396dfa00b15881</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>a6c917133bab98a7937550bb359037d23</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sparse_proportion</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>a9a1adeab0a859ee84f7af45f213b5d2d</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>prefer_rows</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>a98ef4a94d64a2055688b9c736823a7d0</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>prefer_rows_proportion</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>a5963ab48765408c1ac25dd56d038ba41</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>uses_oracle</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>a7b09379690fd41d786df93427ffc816d</anchor>
      <arglist>(bool row) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>ab475b20177f2308e36ab16c2269f2306</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>ab08cb6e685230a45ffbcf62e72770f06</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>a76b120dcf5e8594b1f627fe8aecaa4f4</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>ad79a514a85a201920ca3da1a726bb6ba</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>a5d722d31354f593722cb86a5eda9a79c</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>ad9330217c788494fd60c00d1ea796642</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>a7111b21e62a770a91eac54a4c612b289</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>ad1c13faa0e69a8c0b4054928dd8217a9</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>a6d6563a8ae122d8771e8a9d74970b3af</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>adca14ab1b00e1fe2bb8a28d10d597e3e</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>a3b071d6a73682e496c0c662d7cffab48</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1DelayedBinaryIsometricOp.html</anchorfile>
      <anchor>a5875b9f788ebcdca296f7424342abeed</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
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
  <compound kind="struct">
    <name>tatami::DelayedBooleanNotHelper</name>
    <filename>structtatami_1_1DelayedBooleanNotHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedBooleanScalarHelper</name>
    <filename>structtatami_1_1DelayedBooleanScalarHelper.html</filename>
    <templarg>DelayedBooleanOp op_</templarg>
    <templarg>typename Value_</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedBooleanScalarHelper</name>
      <anchorfile>structtatami_1_1DelayedBooleanScalarHelper.html</anchorfile>
      <anchor>af70a6159d971454fef9eaa22c6ad6c2e</anchor>
      <arglist>(bool s)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedBooleanVectorHelper</name>
    <filename>structtatami_1_1DelayedBooleanVectorHelper.html</filename>
    <templarg>DelayedBooleanOp op_</templarg>
    <templarg>int margin_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Vector_</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedBooleanVectorHelper</name>
      <anchorfile>structtatami_1_1DelayedBooleanVectorHelper.html</anchorfile>
      <anchor>a8eae5e50b3c4fafa9c5f92a3c8b7c18b</anchor>
      <arglist>(Vector_ v)</arglist>
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
    <name>tatami::DelayedCeilingHelper</name>
    <filename>structtatami_1_1DelayedCeilingHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedCompareScalarHelper</name>
    <filename>structtatami_1_1DelayedCompareScalarHelper.html</filename>
    <templarg>DelayedCompareOp op_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Scalar_</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedCompareScalarHelper</name>
      <anchorfile>structtatami_1_1DelayedCompareScalarHelper.html</anchorfile>
      <anchor>a5b9e28a56fc093e389749a34c50e959c</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedCompareVectorHelper</name>
    <filename>structtatami_1_1DelayedCompareVectorHelper.html</filename>
    <templarg>DelayedCompareOp op_</templarg>
    <templarg>int margin_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Vector_</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedCompareVectorHelper</name>
      <anchorfile>structtatami_1_1DelayedCompareVectorHelper.html</anchorfile>
      <anchor>a907bed566e1c9cf7be5e917a491cc9f0</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedCosHelper</name>
    <filename>structtatami_1_1DelayedCosHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedCoshHelper</name>
    <filename>structtatami_1_1DelayedCoshHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedExpHelper</name>
    <filename>structtatami_1_1DelayedExpHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedExpm1Helper</name>
    <filename>structtatami_1_1DelayedExpm1Helper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedFloorHelper</name>
    <filename>structtatami_1_1DelayedFloorHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedGammaHelper</name>
    <filename>structtatami_1_1DelayedGammaHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedLgammaHelper</name>
    <filename>structtatami_1_1DelayedLgammaHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedLog1pHelper</name>
    <filename>structtatami_1_1DelayedLog1pHelper.html</filename>
    <templarg>typename Value_</templarg>
    <templarg>typename Base_</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedLog1pHelper</name>
      <anchorfile>structtatami_1_1DelayedLog1pHelper.html</anchorfile>
      <anchor>a72f0a36b37285d7db294261af4e97f9f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DelayedLog1pHelper</name>
      <anchorfile>structtatami_1_1DelayedLog1pHelper.html</anchorfile>
      <anchor>a2b0dbc27f431eda6d26a40872aa31a8a</anchor>
      <arglist>(Base_ base)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedLogHelper</name>
    <filename>structtatami_1_1DelayedLogHelper.html</filename>
    <templarg>typename Value_</templarg>
    <templarg>typename Base_</templarg>
    <member kind="function">
      <type></type>
      <name>DelayedLogHelper</name>
      <anchorfile>structtatami_1_1DelayedLogHelper.html</anchorfile>
      <anchor>a4f43ca158fd03c86d913c36deac643af</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DelayedLogHelper</name>
      <anchorfile>structtatami_1_1DelayedLogHelper.html</anchorfile>
      <anchor>a8f4ac06070ccc7c92b201b70d7a1f785</anchor>
      <arglist>(Base_ base)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedRoundHelper</name>
    <filename>structtatami_1_1DelayedRoundHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedSignHelper</name>
    <filename>structtatami_1_1DelayedSignHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedSinHelper</name>
    <filename>structtatami_1_1DelayedSinHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedSinhHelper</name>
    <filename>structtatami_1_1DelayedSinhHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedSqrtHelper</name>
    <filename>structtatami_1_1DelayedSqrtHelper.html</filename>
    <templarg>typename Value_</templarg>
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
      <anchor>a6f7dfa2450c9f1f8faf3e2da147fcdc5</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p, Index_ s, Index_ l)</arglist>
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
    <name>tatami::DelayedTanHelper</name>
    <filename>structtatami_1_1DelayedTanHelper.html</filename>
    <templarg>typename Value_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami::DelayedTanhHelper</name>
    <filename>structtatami_1_1DelayedTanhHelper.html</filename>
    <templarg>typename Value_</templarg>
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
  <compound kind="struct">
    <name>tatami::DelayedTruncHelper</name>
    <filename>structtatami_1_1DelayedTruncHelper.html</filename>
    <templarg>typename Value_</templarg>
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
      <anchor>a5e757df46d08e67de4cfc5eb8ad82617</anchor>
      <arglist>(const Options &amp;) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>a0d32ec294715928d00e31d9cf08727d4</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>a20787448e2d340dbaa8a432dcf42a400</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>a79a6302ec5b15d513216fbe36667d430</anchor>
      <arglist>(const Options &amp;) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>a817a447dba85286831ba9dfdd99b4831</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1DenseMatrix.html</anchorfile>
      <anchor>a8abfcfd6aaf5ba60071cb0445adf5fd3</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;) const</arglist>
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
    <name>tatami::FragmentedSparseContents</name>
    <filename>structtatami_1_1FragmentedSparseContents.html</filename>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Value_ &gt; &gt;</type>
      <name>value</name>
      <anchorfile>structtatami_1_1FragmentedSparseContents.html</anchorfile>
      <anchor>a24193eb145792effdc67d3e1425f1cd2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Index_ &gt; &gt;</type>
      <name>index</name>
      <anchorfile>structtatami_1_1FragmentedSparseContents.html</anchorfile>
      <anchor>a8c1a2bebb066b64451576dc356877091</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami::FragmentedSparseMatrix</name>
    <filename>classtatami_1_1FragmentedSparseMatrix.html</filename>
    <templarg>bool row_</templarg>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class ValueVectorStorage_</templarg>
    <templarg>class IndexVectorStorage_</templarg>
    <base>tatami::Matrix</base>
    <member kind="function">
      <type></type>
      <name>FragmentedSparseMatrix</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>a75a66fd646bab027ff762cd824e52fb0</anchor>
      <arglist>(Index_ nr, Index_ nc, ValueVectorStorage_ vals, IndexVectorStorage_ idx, bool check=true)</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>nrow</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>a713230b7d55a09d472e09e8363654183</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>ncol</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>a6b302efc9d37ff270780579781f259e5</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>sparse</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>a93dad37ca355923e02f5fa346dfd74c6</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sparse_proportion</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>ab1e3895c3ef0e0decf14fd05ec26b9bf</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>prefer_rows</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>ace970c2637bf0070cae8f85e65845fc9</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>prefer_rows_proportion</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>a55f6ca0df7adae11d5a194749803fb18</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>uses_oracle</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>ae712803cedf2672f33d46fb952c99d4a</anchor>
      <arglist>(bool) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>a4c0df3ff19af3218780c8e69fa48c794</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>aabb6feb08662d779c77bd42934762ba1</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_row</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>a4c88478fb6a11f15e5d02ca1673d1bac</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>a0fb0d3f53b5ae2c8bf3e4b60699ddf62</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>a07d3618a0d47e61599e37f2f55e77e4f</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexDenseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>dense_column</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>a9767cf1003a4418f72b2a3bdbbb4ec34</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>a69bcdaf32f0d8b62954693864eca0f28</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>afd3cb26fc21e2a7fd4ef8d8186c5ec5d</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_row</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>ac2b5165f7caf5d7972ee4f762b68a845</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FullSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>a36f443f694e812a8ee151d812c3df30e</anchor>
      <arglist>(const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; BlockSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>a17e339136971965ebbc9942f073eb309</anchor>
      <arglist>(Index_ block_start, Index_ block_length, const Options &amp;opt) const</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; IndexSparseExtractor&lt; Value_, Index_ &gt; &gt;</type>
      <name>sparse_column</name>
      <anchorfile>classtatami_1_1FragmentedSparseMatrix.html</anchorfile>
      <anchor>afb6e2f27ca5ca808926b0d7b8e52a34b</anchor>
      <arglist>(std::vector&lt; Index_ &gt; indices, const Options &amp;opt) const</arglist>
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
    <class kind="class">tatami::ArrayView</class>
    <class kind="struct">tatami::BlockExtractor</class>
    <class kind="struct">tatami::CompressedSparseContents</class>
    <class kind="class">tatami::CompressedSparseMatrix</class>
    <class kind="struct">tatami::ConsecutiveOracle</class>
    <class kind="struct">tatami::DelayedAbsHelper</class>
    <class kind="struct">tatami::DelayedAcosHelper</class>
    <class kind="struct">tatami::DelayedAcoshHelper</class>
    <class kind="struct">tatami::DelayedArithScalarHelper</class>
    <class kind="struct">tatami::DelayedArithVectorHelper</class>
    <class kind="struct">tatami::DelayedAsinHelper</class>
    <class kind="struct">tatami::DelayedAsinhHelper</class>
    <class kind="struct">tatami::DelayedAtanHelper</class>
    <class kind="struct">tatami::DelayedAtanhHelper</class>
    <class kind="struct">tatami::DelayedBinaryArithHelper</class>
    <class kind="struct">tatami::DelayedBinaryBooleanHelper</class>
    <class kind="struct">tatami::DelayedBinaryCompareHelper</class>
    <class kind="class">tatami::DelayedBinaryIsometricOp</class>
    <class kind="class">tatami::DelayedBind</class>
    <class kind="struct">tatami::DelayedBooleanNotHelper</class>
    <class kind="struct">tatami::DelayedBooleanScalarHelper</class>
    <class kind="struct">tatami::DelayedBooleanVectorHelper</class>
    <class kind="class">tatami::DelayedCast</class>
    <class kind="struct">tatami::DelayedCeilingHelper</class>
    <class kind="struct">tatami::DelayedCompareScalarHelper</class>
    <class kind="struct">tatami::DelayedCompareVectorHelper</class>
    <class kind="struct">tatami::DelayedCosHelper</class>
    <class kind="struct">tatami::DelayedCoshHelper</class>
    <class kind="struct">tatami::DelayedExpHelper</class>
    <class kind="struct">tatami::DelayedExpm1Helper</class>
    <class kind="struct">tatami::DelayedFloorHelper</class>
    <class kind="struct">tatami::DelayedGammaHelper</class>
    <class kind="struct">tatami::DelayedLgammaHelper</class>
    <class kind="struct">tatami::DelayedLog1pHelper</class>
    <class kind="struct">tatami::DelayedLogHelper</class>
    <class kind="struct">tatami::DelayedRoundHelper</class>
    <class kind="struct">tatami::DelayedSignHelper</class>
    <class kind="struct">tatami::DelayedSinHelper</class>
    <class kind="struct">tatami::DelayedSinhHelper</class>
    <class kind="struct">tatami::DelayedSqrtHelper</class>
    <class kind="class">tatami::DelayedSubset</class>
    <class kind="class">tatami::DelayedSubsetBlock</class>
    <class kind="class">tatami::DelayedSubsetSorted</class>
    <class kind="class">tatami::DelayedSubsetSortedUnique</class>
    <class kind="class">tatami::DelayedSubsetUnique</class>
    <class kind="struct">tatami::DelayedTanHelper</class>
    <class kind="struct">tatami::DelayedTanhHelper</class>
    <class kind="class">tatami::DelayedTranspose</class>
    <class kind="struct">tatami::DelayedTruncHelper</class>
    <class kind="class">tatami::DelayedUnaryIsometricOp</class>
    <class kind="class">tatami::DenseExtractor</class>
    <class kind="class">tatami::DenseMatrix</class>
    <class kind="struct">tatami::ExtractorBase</class>
    <class kind="struct">tatami::FixedOracle</class>
    <class kind="struct">tatami::FragmentedSparseContents</class>
    <class kind="class">tatami::FragmentedSparseMatrix</class>
    <class kind="struct">tatami::FullExtractor</class>
    <class kind="struct">tatami::IndexExtractor</class>
    <class kind="class">tatami::Matrix</class>
    <class kind="struct">tatami::Options</class>
    <class kind="struct">tatami::Oracle</class>
    <class kind="struct">tatami::OracleStream</class>
    <class kind="class">tatami::SemiCompressedSparseMatrix</class>
    <class kind="struct">tatami::SomeNumericArray</class>
    <class kind="class">tatami::SparseExtractor</class>
    <class kind="struct">tatami::SparseRange</class>
    <class kind="struct">tatami::SparseRangeCopy</class>
    <class kind="class">tatami::VirtualDenseMatrix</class>
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
      <type>FragmentedSparseMatrix&lt; false, Value_, Index_, ValueVectorStorage_, IndexVectorStorage_ &gt;</type>
      <name>FragmentedSparseColumnMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a6b9735bdc6e0a2856cd4747efbbf9f3c</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>FragmentedSparseMatrix&lt; true, Value_, Index_, ValueVectorStorage_, IndexVectorStorage_ &gt;</type>
      <name>FragmentedSparseRowMatrix</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a3d97fc9bd75052a79bb643d0fb5a0b98</anchor>
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
    <member kind="typedef">
      <type>typename std::remove_cv&lt; typename std::remove_reference&lt; decltype(std::declval&lt; Array_ &gt;()[0])&gt;::type &gt;::type</type>
      <name>ElementType</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a58a028d23a7be58854b2e60dfae1b04b</anchor>
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
    <member kind="enumeration">
      <type></type>
      <name>DelayedArithOp</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aab44a37b3762de0c5b1ffbfceb25fa0f</anchor>
      <arglist></arglist>
      <enumvalue file="namespacetatami.html" anchor="aab44a37b3762de0c5b1ffbfceb25fa0fa9eeb52badb613229884838847294b90d">ADD</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="aab44a37b3762de0c5b1ffbfceb25fa0fa23ebcc4776b613af25dfbe7c8ce4813e">SUBTRACT</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="aab44a37b3762de0c5b1ffbfceb25fa0fa080aaf8d817ada96fca7096b7b55bd30">MULTIPLY</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="aab44a37b3762de0c5b1ffbfceb25fa0fa210c66d794cec40488f3f8f634d6c33b">DIVIDE</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="aab44a37b3762de0c5b1ffbfceb25fa0fac9c9c146c630ca5ef9197c73c032f4a6">POWER</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="aab44a37b3762de0c5b1ffbfceb25fa0fa928ab45d616dde447dbbbd0270db87ad">MODULO</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="aab44a37b3762de0c5b1ffbfceb25fa0fa051460a4a75d4d251a41a7c04bf49412">INTEGER_DIVIDE</enumvalue>
    </member>
    <member kind="enumeration">
      <type></type>
      <name>DelayedBooleanOp</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a2104862d4068933ea4cc805c92f82d07</anchor>
      <arglist></arglist>
      <enumvalue file="namespacetatami.html" anchor="a2104862d4068933ea4cc805c92f82d07a558ffc8f5770d8e4f95f51d822685532">AND</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="a2104862d4068933ea4cc805c92f82d07a1d00e7dce692e8dc3f6877f035e3a616">OR</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="a2104862d4068933ea4cc805c92f82d07a97675eb3f268048604dc5155511a2a4d">XOR</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="a2104862d4068933ea4cc805c92f82d07a969f331a87d8c958473c32b4d0e61a44">EQUAL</enumvalue>
    </member>
    <member kind="enumeration">
      <type></type>
      <name>DelayedCompareOp</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac4fc175a57ace709941b5ca7ddb19708</anchor>
      <arglist></arglist>
      <enumvalue file="namespacetatami.html" anchor="ac4fc175a57ace709941b5ca7ddb19708a969f331a87d8c958473c32b4d0e61a44">EQUAL</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="ac4fc175a57ace709941b5ca7ddb19708a1625ef4fe09f68fa20d3ff6e02cd5c8e">GREATER_THAN</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="ac4fc175a57ace709941b5ca7ddb19708aa327176a0a845c117bdfadec134a95e9">LESS_THAN</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="ac4fc175a57ace709941b5ca7ddb19708aa6eac69202c3dc2978176801a84e4d1d">GREATER_THAN_OR_EQUAL</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="ac4fc175a57ace709941b5ca7ddb19708a8397780541b6289d2a0b991d1c28c432">LESS_THAN_OR_EQUAL</enumvalue>
      <enumvalue file="namespacetatami.html" anchor="ac4fc175a57ace709941b5ca7ddb19708a4ea2d378cdec20f59330f113297bc1ce">NOT_EQUAL</enumvalue>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>extracted_length</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>af9d13ceaa112d2c091265510d741488d</anchor>
      <arglist>(const ConditionalSelectionExtractor&lt; selection_, Index_ &gt; &amp;ex)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>new_extractor</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a95b6638dceda82d82a7579dc88a45709</anchor>
      <arglist>(const Matrix&lt; Value_, Index_ &gt; *ptr, Args_ &amp;&amp;... args)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>convert_to_dense</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a3d709db6be55e94d987d38c4c07c71c0</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, StoredValue_ *store, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index &gt; &gt;</type>
      <name>convert_to_dense</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a3a0e048d06c98ea3706466139c4c97dc</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_dense</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aa74ad14ba410177396121d272371dffd</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, int order, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryArithHelper&lt; DelayedArithOp::ADD &gt;</type>
      <name>make_DelayedBinaryAddHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a1d138b1e7f6a26f814c025363cc3db80</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryArithHelper&lt; DelayedArithOp::SUBTRACT &gt;</type>
      <name>make_DelayedBinarySubtractHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac4c951d489cb0f2bb0e90f2cf4a25862</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryArithHelper&lt; DelayedArithOp::MULTIPLY &gt;</type>
      <name>make_DelayedBinaryMultiplyHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a2569c540083a24f92af8140358e1e9c4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryArithHelper&lt; DelayedArithOp::DIVIDE &gt;</type>
      <name>make_DelayedBinaryDivideHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aeb4cf766f850766b966f7121728a6af8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryArithHelper&lt; DelayedArithOp::POWER &gt;</type>
      <name>make_DelayedBinaryPowerHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ab6c73c04e7b08130ed7b46d93c4dfd11</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryArithHelper&lt; DelayedArithOp::MODULO &gt;</type>
      <name>make_DelayedBinaryModuloHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a04293eb1e8eefb7024cc192d75ac093e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryArithHelper&lt; DelayedArithOp::INTEGER_DIVIDE &gt;</type>
      <name>make_DelayedBinaryIntegerDivideHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aa8f54741424bef6a21225935a45e9d53</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryBooleanHelper&lt; DelayedBooleanOp::EQUAL &gt;</type>
      <name>make_DelayedBinaryBooleanEqualHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a1063aea86897ef76e17b1772320b8f7d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryBooleanHelper&lt; DelayedBooleanOp::AND &gt;</type>
      <name>make_DelayedBinaryBooleanAndHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a31625dbb50d420fe1b0da406f19ef33d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryBooleanHelper&lt; DelayedBooleanOp::OR &gt;</type>
      <name>make_DelayedBinaryBooleanOrHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a65bad100b39d9372fbeb23d16dae6588</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryBooleanHelper&lt; DelayedBooleanOp::XOR &gt;</type>
      <name>make_DelayedBinaryBooleanXorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aeb091953ab0f55406935a52b7ecb7350</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryCompareHelper&lt; DelayedCompareOp::EQUAL &gt;</type>
      <name>make_DelayedBinaryEqualHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac3599db8f6d19fc9172da45a3b6244fc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryCompareHelper&lt; DelayedCompareOp::GREATER_THAN &gt;</type>
      <name>make_DelayedBinaryGreaterThanHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>abe4fb79c5a2cddab2f5c9cf0bd9fe3d0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryCompareHelper&lt; DelayedCompareOp::LESS_THAN &gt;</type>
      <name>make_DelayedBinaryLessThanHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ae6a24f5bf2c2c61382aef43f2c8d3362</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryCompareHelper&lt; DelayedCompareOp::GREATER_THAN_OR_EQUAL &gt;</type>
      <name>make_DelayedBinaryGreaterThanOrEqualHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ad31cd299427d73538dc00dcbffa9cd2b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryCompareHelper&lt; DelayedCompareOp::LESS_THAN_OR_EQUAL &gt;</type>
      <name>make_DelayedBinaryLessThanOrEqualHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a38c7cc26409e5c995b7264164b148f4c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBinaryCompareHelper&lt; DelayedCompareOp::NOT_EQUAL &gt;</type>
      <name>make_DelayedBinaryNotEqualHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a3073c4e33dcf89416af292399ba55c29</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>make_DelayedBinaryIsometricOp</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>af39c0767d23c6bdc5f4267f062ef2e17</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; left, std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; right, Operation_ op)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithScalarHelper&lt; DelayedArithOp::ADD, true, Value_, Scalar_ &gt;</type>
      <name>make_DelayedAddScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a33b2b9798f2a43b62821131ba3a3f6bd</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithScalarHelper&lt; DelayedArithOp::SUBTRACT, right_, Value_, Scalar_ &gt;</type>
      <name>make_DelayedSubtractScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a9bdfe26c0f1426b611e91ee7569045de</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithScalarHelper&lt; DelayedArithOp::MULTIPLY, true, Value_, Scalar_ &gt;</type>
      <name>make_DelayedMultiplyScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aa653b51960498cc467ea69c18ea0c097</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithScalarHelper&lt; DelayedArithOp::DIVIDE, right_, Value_, Scalar_ &gt;</type>
      <name>make_DelayedDivideScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a60076096ee4296c0e187a226be5ca4b5</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithScalarHelper&lt; DelayedArithOp::POWER, right_, Value_, Scalar_ &gt;</type>
      <name>make_DelayedPowerScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a8c478145f87c37a9fe6d4e8490fbf05c</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithScalarHelper&lt; DelayedArithOp::MODULO, right_, Value_, Scalar_ &gt;</type>
      <name>make_DelayedModuloScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a5504f584d6db28c05d78f1917c40810d</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithScalarHelper&lt; DelayedArithOp::INTEGER_DIVIDE, right_, Value_, Scalar_ &gt;</type>
      <name>make_DelayedIntegerDivideScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a5ec92ba88a0f5e9d7bf52021ed28e859</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithVectorHelper&lt; DelayedArithOp::ADD, true, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedAddVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a628c631be4424b37c1c1e68e28c5e982</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithVectorHelper&lt; DelayedArithOp::SUBTRACT, right_, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedSubtractVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a800b92854282015d96308fb283aeb508</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithVectorHelper&lt; DelayedArithOp::MULTIPLY, true, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedMultiplyVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ab929acdeb634187cd08d1100bbbf1b29</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithVectorHelper&lt; DelayedArithOp::DIVIDE, right_, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedDivideVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a1a4629f17aa4f3a4be06f3479bb6f68f</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithVectorHelper&lt; DelayedArithOp::POWER, right_, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedPowerVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aafe489bd754c040491ebc753a4e5656e</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithVectorHelper&lt; DelayedArithOp::MODULO, right_, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedModuloVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a4220ce5709e46ebce77b240b572d97d9</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedArithVectorHelper&lt; DelayedArithOp::INTEGER_DIVIDE, right_, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedIntegerDivideVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a857fc6c11ec69895ac85991ac83c395d</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanNotHelper&lt; Value_ &gt;</type>
      <name>make_DelayedBooleanNotHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac812e757b55c2145bc3d4cd113c712fd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanScalarHelper&lt; DelayedBooleanOp::AND, Value_ &gt;</type>
      <name>make_DelayedBooleanAndScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ae460e24124cef43d8a5d3b82e8f02325</anchor>
      <arglist>(bool s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanScalarHelper&lt; DelayedBooleanOp::OR &gt;</type>
      <name>make_DelayedBooleanOrScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac7dbc64e6b0a0c8c48b1423b4b6748f1</anchor>
      <arglist>(bool s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanScalarHelper&lt; DelayedBooleanOp::XOR &gt;</type>
      <name>make_DelayedBooleanXorScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a29f6f1b1d28c670f38fda651fe2aa03a</anchor>
      <arglist>(bool s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanScalarHelper&lt; DelayedBooleanOp::EQUAL &gt;</type>
      <name>make_DelayedBooleanEqualScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a6bc5d47e0c391c17c85bddf19cc5620a</anchor>
      <arglist>(bool s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanVectorHelper&lt; DelayedBooleanOp::AND, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedBooleanAndVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a1fd06669afab8565747e08bbce4b041f</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanVectorHelper&lt; DelayedBooleanOp::OR, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedBooleanOrVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a503f23fe3e5a2c5f5bc6c0db01a7f64d</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanVectorHelper&lt; DelayedBooleanOp::XOR, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedBooleanXorVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a4bcdace63565fc670b0393675e33bb11</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedBooleanVectorHelper&lt; DelayedBooleanOp::EQUAL, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedBooleanEqualVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a8b07dc8b7a71776bcc222d7c437e2ca3</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareScalarHelper&lt; DelayedCompareOp::EQUAL, Value_, Scalar_ &gt;</type>
      <name>make_DelayedEqualScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a11bfea31e4524927274dc1d80c32fd96</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareScalarHelper&lt; DelayedCompareOp::GREATER_THAN, Value_, Scalar_ &gt;</type>
      <name>make_DelayedGreaterThanScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a1961ba2a59d705adf7b8627523a62b3a</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareScalarHelper&lt; DelayedCompareOp::LESS_THAN, Value_, Scalar_ &gt;</type>
      <name>make_DelayedLessThanScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a7d2d381ec9c0d87f5d561b4dd9d5a58e</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareScalarHelper&lt; DelayedCompareOp::GREATER_THAN_OR_EQUAL, Value_, Scalar_ &gt;</type>
      <name>make_DelayedGreaterThanOrEqualScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac3df12ccdf4692c0f660f1634a79b857</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareScalarHelper&lt; DelayedCompareOp::LESS_THAN_OR_EQUAL, Value_, Scalar_ &gt;</type>
      <name>make_DelayedLessThanOrEqualScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>af5e67b02534cbe1813be5e2d690af88b</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareScalarHelper&lt; DelayedCompareOp::NOT_EQUAL, Value_, Scalar_ &gt;</type>
      <name>make_DelayedNotEqualScalarHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ac5d5cef56d2119e73e3f670d49de6b48</anchor>
      <arglist>(Scalar_ s)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareVectorHelper&lt; DelayedCompareOp::EQUAL, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedEqualVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aa40c54d5c5cc1109d42a150f5cd7ed2f</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareVectorHelper&lt; DelayedCompareOp::GREATER_THAN, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedGreaterThanVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>aaf1951e22415be4714495556f701a7fd</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareVectorHelper&lt; DelayedCompareOp::LESS_THAN, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedLessThanVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a6dca8b9072813f6e9e52aa0536817b82</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareVectorHelper&lt; DelayedCompareOp::GREATER_THAN_OR_EQUAL, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedGreaterThanOrEqualVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a70a3e0489fdaccc80c7cc0125e9cf906</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareVectorHelper&lt; DelayedCompareOp::LESS_THAN_OR_EQUAL, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedLessThanOrEqualVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a83999d2ac877c5b15f6fd5a0f33a894a</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>DelayedCompareVectorHelper&lt; DelayedCompareOp::NOT_EQUAL, margin_, Value_, Vector_ &gt;</type>
      <name>make_DelayedNotEqualVectorHelper</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a6c84e8b757cf5ba6479a6ea3e6e530e3</anchor>
      <arglist>(Vector_ v)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>make_DelayedUnaryIsometricOp</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a07b072944c766d85dc4ce4dc734e0b75</anchor>
      <arglist>(std::shared_ptr&lt; const Matrix&lt; Value_, Index_ &gt; &gt; p, Operation_ op)</arglist>
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
      <type>CompressedSparseContents&lt; Value_, Index_ &gt;</type>
      <name>retrieve_compressed_sparse_contents</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ad1aee6313d02036d2f56703f372e1df6</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, bool two_pass, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_compressed_sparse</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a48dab0e55d3b852afcfd9e93835190df</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, bool two_pass=false, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_compressed_sparse</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>acf7ca799fac9cff0c16a77bf98e32ff3</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, int order, bool two_pass=false, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>FragmentedSparseContents&lt; Value_, Index_ &gt;</type>
      <name>retrieve_fragmented_sparse_contents</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a012194aba94c37703f15555247e02155</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_fragmented_sparse</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>a953ad06ecc9a36c8d2c1f3febe84a787</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, int threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>convert_to_fragmented_sparse</name>
      <anchorfile>namespacetatami.html</anchorfile>
      <anchor>ae46babb9d985cffcb54bcba89553a811</anchor>
      <arglist>(const Matrix&lt; InputValue_, InputIndex_ &gt; *incoming, int order, int threads=1)</arglist>
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
      <anchor>ad2eddeaad50bdad0bc1e5c8f2f8549e0</anchor>
      <arglist>(Function_ fun, Index_ tasks, size_t #if defined(_OPENMP)||defined(TATAMI_CUSTOM_PARALLEL) threads #endif)</arglist>
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
  <compound kind="page">
    <name>index</name>
    <title>A C++ API for all sorts of matrices</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__github_workspace_README</docanchor>
  </compound>
</tagfile>
