.. _sparsekernel_DIA:

Diagonal Storage Format (DIA)
=============================

The DIA format keeps the matrix in order, but extends every diagonal of
the matrix, that contains a non-zero-element. The other diagonals are ignored completely. The extension is done by
adding zeroes "outside" the matrix until all diagonals have the same specific length (numElementsPerDiagonal). This
specific length is either the number of rows (numRows) or the number of columns (numColumns) and depends on which of
these two integer holds the larger value. The number of diagonals is saved as well (numDiagonals) and the total number
of values (including the added zeroes, excluding the zeroes that were "deleted" in order to ignore the unnecessary
diagonals) is calculated like this: *numValues* = *numDiagonals* * *numElementsPerDiagonal*. All the elements
are stored in diagonal major order in an array (values) and another array shows the offset of the main diagonal
(offset). Negative values in the offset array represent diagonals "below" the main diagonal (its original position),
positive values represent diagonals "above" or "right" from the main diagonal.

Example
-------

.. math::

  A = \left(\begin{matrix} 6 & 0  & 0 & 4 \\
    7 & 0 & 0 & 0 \\
    0 & 0 & -9 & 4 \\
    2 & 5 & 0 & 3 \\
    2 & 0 & 0 & 1 \\
    0 & 0 & 0 & 0 \\
    0 & 1 & 0 & 2 \end{matrix}\right) 

The DIA format without diagonal element shifting looks like this:

.. math::
  values_{A} = \left(\begin{matrix} 
    0 & 0 & 0 & 0 & 0 & 6 & 0 &  - & 4 & - & - & - & -  & - & - \\
    - & 0 & 0 & 0 & 0 & 7 & 0 &  0 & - & 0 & - & - & -  & - & - \\
    - & - & 0 & 0 & 0 & 0 & 0 & -9 & 4 & - & 0 & - & -  & - & - \\
    - & - & - & 0 & 0 & 2 & 5 &  0 & 3 & 0 & - & 0 & -  & - & - \\
    - & - & - & - & 0 & 2 & 0 &  0 & 1 & 0 & 0 & - & 0  & - & - \\
    - & - & - & - & - & 0 & 0 &  0 & 0 & 0 & 0 & 0 & -  & 0 & - \\
    - & - & - & - & - & - & 1 &  0 & 2 & 0 & 0 & 0 & 0  & - & 0  \end{matrix}\right) 

Here are the corresponding arrays for the representation:

.. math::
    
    \begin{align}
    numRows &= 7 \\
    numColums &= 4 \\
    numDiagonals &= 8 \\
    values &= \left( \begin{matrix}
                     0 & 0 & 0 & 0 & 0 & 0 & 1 \\
                     0 & 0 & 0 & 0 & 2 & 0 & 0 \\
                     0 & 0 & 0 & 2 & 0 & 0 & 2 \\
                     0 & 0 & 0 & 5 & 1 & 0 & 0 \\
                     0 & 7 & 0 & 0 & 1 & 0 & 0 \\
                     6 & 0 & 9 & 3 & 0 & 0 & 0 \\
                     0 & 0 & 4 & 0 & 0 & 0 & 0 \\
                     4 & 0 & 0 & 0 & 0 & 0 & 0 \\
                     \end{matrix}\right) \\
    offset &= \left( \begin{matrix}
                      -5 & -4 & -3 & -2 & -1 & 0 & 1 & 3 
                     \end{matrix}\right) \\
    \end{align}

Remarks
-------

 * Using the DIA format is only helpful if there is only a limited number of diagonals
   that are really used.
 * The ``values`` array contains a lot of zero entries that stand for out-of-range 
   elements of the original matrix.

DIAKernelTrait
--------------

Conversion
^^^^^^^^^^

========================= ============================================================= ==== ====
**Functionname**          **Description**                                               Host CUDA
========================= ============================================================= ==== ====
getCSRSizes               DIA --> CSR: get sparse row sizes                             *
getCSRValues              DIA --> CSR: conversion DIA to CSR                            *
========================= ============================================================= ==== ====

Calculation
^^^^^^^^^^^

========================= ============================================================= ==== ====
**Functionname**          **Description**                                               Host CUDA
========================= ============================================================= ==== ====
normalGEMV                matrix-vector multiplication                                  *    *
jacobi                    compute one jacobi iteration step                             *
absMaxVal                 compute the maximal absolute value                            *
========================= ============================================================= ==== ====

