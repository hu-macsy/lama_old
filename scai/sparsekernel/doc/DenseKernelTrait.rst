DenseKernelTrait
================

Conversion
----------

========================= ============================================================= ==== ==== ===
**Functionname**          **Description**                                               Host CUDA MIC
========================= ============================================================= ==== ==== ===
nonZeroValues             count non-zero values                                         *
getCSRSizes               Dense --> CSR: get sparse row sizes                           *
getCSRValues              Dense --> CSR: conversion Dense to CSR                        *
setCSRValues              CSR --> Dense: conversion CSR to Dense                        *
copyDenseValues           copy values of dense matrix                                   *
getRow                    returns a row of the matrix                                   *
getDiagonal               returns diagonal of a matrix                                  *
setDiagonal               sets the diagonal of a matrix                                 *
setValue                  sets a single value in a matrix                               *
setDiagonalValue          sets diagonal to a single value                               *
========================= ============================================================= ==== ==== ===

Calculation
-----------

========================= ============================================================= ==== ==== ===
**Functionname**          **Description**                                               Host CUDA MIC
========================= ============================================================= ==== ==== ===
scaleValue                scale all elements with a single value                        *
scaleRows                 multiplies each row with an own value                         *
========================= ============================================================= ==== ==== ===

