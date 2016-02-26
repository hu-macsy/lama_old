DIAKernelTrait
==============

Conversion
----------

========================= ============================================================= ==== ==== ===
**Functionname**          **Description**                                               Host CUDA MIC
========================= ============================================================= ==== ==== ===
getCSRSizes               DIA --> CSR: get sparse row sizes                             *
getCSRValues              DIA --> CSR: conversion DIA to CSR                            *
========================= ============================================================= ==== ==== ===

Calculation
-----------

========================= ============================================================= ==== ==== ===
**Functionname**          **Description**                                               Host CUDA MIC
========================= ============================================================= ==== ==== ===
normalGEMV                matrix-vector multiplication                                  *    *    *
normalGEVM                vector-matrix multiplication                                  *    *
jacobi                    compute one jacobi iteration step                             *         *
absMaxVal                 compute the maximal absolute value                            *
========================= ============================================================= ==== ==== ===

