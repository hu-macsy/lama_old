.. _sparsekernel_JDS:

Jagged Diagonal Storage Format (JDS)
====================================

The JDS format avoids the need of filling each row with dummy data.
Therefore it sorts the rows by length, so the longest row stands on top of the matrix and the shortest at the bottom. 
Like the ELL matrix the elements in the *values*-array are entered in column major order. The JDS comes with the integer
*numValues*, *numRows*, *numColumns* and the number of jagged diagonals (which is equal to the number of
columns in the "jagged" Matrix): *numDiagonals*. It contains 5 arrays: One array for the length of each column
(dlg) and one for the length of each row (ilg), the permutation array which shows, where the lines were supposed to
be before the assorting (perm) and the arrays for the elements (values) and their original column indices (ja) as in
the CSR and ELL format.

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

Starting from this compressed representation

.. math::

  ia_{A} = \left(\begin{matrix} 2 \\
    1  \\
    2 \\
    3 \\
    2 \\
    0 \\
    2 \end{matrix}\right) 
  values_{A} = \left(\begin{matrix} 6 & 4  \\
    7 \\
    -9 & 4 \\
    2 & 5 & 3 \\
    2 & 1 \\
     \\
    1 & 2 \end{matrix}\right) 
  ja_{A} = \left(\begin{matrix} 0 & 3  \\
    0  \\
    2 & 3  \\
    0 & 1 & 3 \\
    0 & 3 \\
     \\
    1 & 3 \end{matrix}\right) 

the arrays are sorted corresponding to its sizes.

.. math::

  perm_{A} = \left(\begin{matrix} 3 \\
    0  \\
    1 \\
    2 \\
    4 \\
    6 \\
    5 \end{matrix}\right) 
  ilg_{A} = \left(\begin{matrix} 3 \\
    2  \\
    2 \\
    2 \\
    1 \\
    1 \\
    0 \end{matrix}\right) 
  values_{A} = \left(\begin{matrix} 
    2 & 5 & 3 \\
    6 & 4 \\
    -9 & 4  \\
    2 & 1 \\
    1 & 2   \\
    7 \\
        \end{matrix}\right) 
  ja_{A} = \left(\begin{matrix}
    0 & 1 & 3 \\
    0 & 3 \\
    2 & 3 \\
    0 & 3 \\
    1 & 3  \\
    0 \\
     \end{matrix}\right) 

Here are the JDS arrays for the representation:

.. math::
    
    \begin{align}
    numRows &= 7 \\
    numColums &= 4 \\
    numValues &= 12 \\
    numDiagonals &= 3 \\
    ilg &= \left[\begin{matrix} 3 & 2 & 2 & 2 & 2 & 1 & 0 \end{matrix}\right] \\
    perm &= \left[\begin{matrix} 3 & 0 & 2 & 4 & 6 & 1 & 5 \end{matrix}\right] \\
    ja     &= [ \begin{matrix} 0 & 0 &  2 & 0 & 1 & 0 & 1 & 3 & 3 & 3 & 3 & 3 \end{matrix} ] \\
    values &= [ \begin{matrix} 2 & 6 & -9 & 2 & 1 & 7 & 5 & 4 & 4 & 1 & 2 & 3 \end{matrix} ] \\
    dlg &= [ \begin{matrix} 6 & 5 & 1 \end{matrix} ] \\
    \end{align}

JDSKernelTrait
--------------

Conversion
^^^^^^^^^^

========================= ============================================================= ==== ====
**Functionname**          **Description**                                               Host CUDA
========================= ============================================================= ==== ====
sortRows                  sorting of values in descending order                         *    *
ilg2dlg                   compute dlg array from ilg array                              *    *
getCSRValues              JDS --> CSR: conversion JDS to CSR                            *    *
setCSRValues              CSR --> JDS: conversion CSR to JDS                            *    *
getRow                    returns a row of the matrix                                   *    *
getValue                  get single element of matrix                                  *    *
========================= ============================================================= ==== ====

Calculation
^^^^^^^^^^^

========================= ============================================================= ==== ====
**Functionname**          **Description**                                               Host CUDA
========================= ============================================================= ==== ====
jacobi                    compute one jacobi iteration step                             *    *
jacobiHalo                compute one jacobi iteration step on halo values              *    *
normalGEMV                matrix-vector multiplication                                  *    *
scaleValue                scale with array                                              *    *
========================= ============================================================= ==== ====

Properties
^^^^^^^^^^

========================= ============================================================= ==== ====
**Functionname**          **Description**                                               Host CUDA
========================= ============================================================= ==== ====
checkDiagonalProperty     Checks if the first n entries are the diagonal elements       *    *
========================= ============================================================= ==== ====
