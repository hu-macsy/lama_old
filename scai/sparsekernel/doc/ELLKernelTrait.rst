ELLKernelTrait
==============

Conversion
----------

========================= ============================================================= ==== ==== ===
**Functionname**          **Description**                                               Host CUDA MIC
========================= ============================================================= ==== ==== ===
fillELLValues             fill up ja and values array                                   *    *
getCSRValues              ELL --> CSR: conversion ELL to CSR                            *    *    *
setCSRValues              CSR --> ELL: conversion CSR to ELL                            *    *    *
compressIA                compress the ia array by using values array and epsilon       *
compressValues            compress ja and values array by using epsilon                 *
getRow                    returns a row of the matrix                                   *    *    *
getValue                  get single element of matrix                                  *    *    *
countNonEmptyRowsBySizes  count non-empty rows by sizes array                           *    *    *
setNonEmptyRowsBySizes    set non-empty rows by sizes array                             *    *    *
========================= ============================================================= ==== ==== ===

Calculation
-----------

========================= ============================================================= ==== ==== ===
**Functionname**          **Description**                                               Host CUDA MIC
========================= ============================================================= ==== ==== ===
jacobi                    compute one jacobi iteration step                             *    *    *
jacobiHalo                compute one jacobi iteration step on halo values              *    *    *
normalGEMV                matrix-vector multiplication                                  *    *    *
sparseGEMV                matrix-vector multiplication with just non-zero rows          *    *    *
normalGEVM                vector-matrix multiplication                                  *    *
sparseGEVM                vector-matrix multiplication with just non-zero rows          *    *
absMaxVal                 compute the maximal absolute value                            *
scaleValue                scale with array                                              *    *    *
matrixMultiplySizes       computes row sizes for result of matrix multiplication        *
matrixAddSizes            computes row sizes for result of matrix addition              *
matrixAdd                 matrix-matrix addition (all ELL)                              *
matrixMultiply            matrix-matrix multiplication  (all ELL)                       *
========================= ============================================================= ==== ==== ===

Properties
----------

========================= ============================================================= ==== ==== ===
**Functionname**          **Description**                                               Host CUDA MIC
========================= ============================================================= ==== ==== ===
hasDiagonalProperty       Checks if the first n entries are the diagonal elements       *    *    *
check                     Checks integrity of ia array                                  *    *    *
========================= ============================================================= ==== ==== ===

