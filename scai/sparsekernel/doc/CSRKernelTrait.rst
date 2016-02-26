CSRKernelTrait
==============

Conversion
----------

====================== ============================================================= ==== ==== ===
**Functionname**       **Description**                                               Host CUDA MIC
====================== ============================================================= ==== ==== ===
sortRowElements        sorts the elements of a row by increasing column indexes      *
sizes2offsets          computes offset array from sizes array                        *    *    *
offsets2sizes          computes sizes array from offset array                        *    *    *   
convertCSR2CSC         converts from CSR2CSC                                         *    *
compress               fill compresses CSR data in new data structures               *
====================== ============================================================= ==== ==== ===

Calculation
-----------

====================== ============================================================= ==== ==== ===
**Functionname**       **Description**                                               Host CUDA MIC
====================== ============================================================= ==== ==== ===
jacobi                 compute one jacobi iteration step                             *    *    *
jacobiHalo             compute one jacobi iteration step on halo values              *    *    *
jacobiHaloWithDiag     compute one jacobi iteration step on halo values              *    *    *
matrixAddSizes         computes row sizes for result of matrix addition              *    *    *
matrixMultiplySizes    computes row sizes for result of matrix multiplication        *    *    *
matrixMultiplyJA       computes column indexes for result of matrix multiplication   *
scaleRows              multiplies each row with an own value                         *    *    *
absMaxDiffVal          computes the maximal element-wise difference for two matrices *         *
normalGEMV             matrix-vector multiplication                                  *    *    *
normalGEVM             vector-matrix multiplication                                  *    *
sparseGEMV             matrix-vector multiplication with just non-zero rows          *    *    *
sparseGEVM             vector-matrix multiplication with just non-zero rows          *    *
gemm                   matrix-matrix multiplication (CSR * Dense)                    *
matrixAdd              matrix-matrix addition (all CSR)                              *    *    *
matrixMultiply         matrix-matrix multiplication  (all CSR)                       *    *    *
====================== ============================================================= ==== ==== ===

Properties
----------

====================== ============================================================= ==== ==== ===
**Functionname**       **Description**                                               Host CUDA MIC
====================== ============================================================= ==== ==== ===
validOffsets           checks for legal offset array                                 *         *
hasDiagonalProperty    checks if CSR data has diagonal property                      *    *    *
countNonZeros          count non-zero entries                                        *
====================== ============================================================= ==== ==== ===

