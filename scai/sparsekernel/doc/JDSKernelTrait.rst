JDSKernelTrait
==============

Conversion
----------

========================= ============================================================= ==== ==== ===
**Functionname**          **Description**                                               Host CUDA MIC
========================= ============================================================= ==== ==== ===
sortRows                  sorting of values in descending order                         *    *    *
setInversePerm            compute inverse permutation for a given permutation           *    *    *
ilg2dlg                   compute dlg array from ilg array                              *    *    *
getCSRValues              JDS --> CSR: conversion JDS to CSR                            *    *    *
setCSRValues              CSR --> JDS: conversion CSR to JDS                            *    *    *
getRow                    returns a row of the matrix                                   *    *    *
getValue                  get single element of matrix                                  *    *    *
========================= ============================================================= ==== ==== ===

Calculation
-----------

========================= ============================================================= ==== ==== ===
**Functionname**          **Description**                                               Host CUDA MIC
========================= ============================================================= ==== ==== ===
jacobi                    compute one jacobi iteration step                             *    *    *
jacobiHalo                compute one jacobi iteration step on halo values              *    *    *
normalGEMV                matrix-vector multiplication                                  *    *    *
normalGEVM                vector-matrix multiplication                                  *    *
scaleValue                scale with array                                              *    *    *
========================= ============================================================= ==== ==== ===

Properties
----------

========================= ============================================================= ==== ==== ===
**Functionname**          **Description**                                               Host CUDA MIC
========================= ============================================================= ==== ==== ===
checkDiagonalProperty     Checks if the first n entries are the diagonal elements       *    *    *
========================= ============================================================= ==== ==== ===
