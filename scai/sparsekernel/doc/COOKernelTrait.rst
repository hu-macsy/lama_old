COOKernelTrait
==============

Conversion
----------

====================== ============================================================= ==== ==== ===
**Functionname**       **Description**                                               Host CUDA MIC
====================== ============================================================= ==== ==== ===
getCSRSizes            COO --> CSR: get sparse row sizes                             *         *
offsets2ia             CSR --> COO: converting offset array (CSR) to ia array (COO)  *    *    *
getCSRValues           COO --> CSR: conversion COO to CSR                            *
setCSRData             CSR --> COO: conversion CSR to COO                            *    *    *
====================== ============================================================= ==== ==== ===

Calculation
-----------

====================== ============================================================= ==== ==== ===
**Functionname**       **Description**                                               Host CUDA MIC
====================== ============================================================= ==== ==== ===
scaleRows              multiplies each row with an own value                         *
normalGEMV             matrix-vector multiplication                                  *    *    *
normalGEVM             vector-matrix multiplication                                  *    *
jacobi                 compute one jacobi iteration step                             *
jacobiHalo             compute one jacobi iteration step on halo values
====================== ============================================================= ==== ==== ===

Properties
----------

====================== ============================================================= ==== ==== ===
**Functionname**       **Description**                                               Host CUDA MIC
====================== ============================================================= ==== ==== ===
hasDiagonalProperty    Checks if the first n entries are the diagonal elements       *
====================== ============================================================= ==== ==== ===

