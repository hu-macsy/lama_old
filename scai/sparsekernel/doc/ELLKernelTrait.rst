.. _sparsekernel_ELL:

ELLPACK-R Storage Format (ELL)
==============================

The ELL format is much like the CSR format, though the compressed matrix is filled with zeroes to obtain a "shortened"
version of the original. The ELL matrix saves the number of rows it (and the original matrix) have (numRows), the
number of columns the ELL format has which is equal to the length of its longest rows (numValuesPerRow), the original
number of columns (numColumns). The total number of values, including the zeroes that are used as a filler, can be
calculated with *numRows* * *numValuesPerRow*. Additionally the ELL format saves three arrays as well, one for
all the values in the ELL matrix, which are stored in column major order (values), one for the number of non-zero
values (plus the main diagonal zeroes, if needed) (ia) and one for their associated columns (ja). In case of the
"filler"-zeroes, ja points at the last element's column of this row.

The ELL format is used for matrices with about equivalent numbers of non-zero-values in each row.

In contrary to the CSR format the ELL format uses the array ``ia`` with the sizes of each row. 
The offsets become redundant as the position of an entry can be directly determined.

.. code-block:: c++

   #pragma omp parallel for
   for ( IndexType i = 0; i < numRows; ++i )
   {
       // code that operates on row i of the storage

       for ( IndexType k = 0; k < ia[i]; ++k )
       {
           ... csrJA[k * numRows + i] for column pos, csrValues[k * numRows + i] for value
       }
   }

Example
-------

Matrix:

-------

.. math::

  A = \left(\begin{matrix} 6 & 0  & 0 & 4 \\
    7 & 0 & 0 & 0 \\
    0 & 0 & -9 & 4 \\
    2 & 5 & 0 & 3 \\
    2 & 0 & 0 & 1 \\
    0 & 0 & 0 & 0 \\
    0 & 1 & 0 & 2 \end{matrix}\right) 

The ELL representation of the sparse matrix is:

.. math::

  values_{A} = \left(\begin{matrix} 6 & 4 & - \\
    7 & - & - \\
    -9 & 4 & - \\
    2 & 5 & 3 \\
    2 & 1 & - \\
    - & - & - \\
    1 & 2 & - \end{matrix}\right) 
  ja_{A} = \left(\begin{matrix} 0 & 3 & - \\
    0 & - & - \\
    2 & 3  & - \\
    0 & 1 & 3 \\
    0 & 3 & - \\
    - & - & - \\
    1 & 3 & - \end{matrix}\right) 
  ia_{A} = \left(\begin{matrix} 2 \\
    1  \\
    2 \\
    3 \\
    2 \\
    0 \\
    2 \end{matrix}\right) 

Here are the corresponding arrays for the representation:

.. math::
    
    \begin{align}
    numRows &= 7 \\
    numColums &= 4 \\
    numValuesPerRow &= 3 \\
    ia &= [\begin{matrix} 2 & 1 & 2 & 3 & 2 & 0 & 2 \end{matrix} ] \\
    ja &= \left[\begin{matrix} 
                 0 & 0 & 2 & 0 & 0 & - & 1 & 
                 3 & - & 3 & 1 & 3 & - & 3 &
                 - & - & - & 3 & - & - & -   \end{matrix} \right] \\
    values &= \left[\begin{matrix} 
                 6 & 7 & -9 & 2 & 2 & - & 1 & 
                 4 & - & 4 & 5 & 1 & - & 2 &
                 - & - & - & 3 & - & - & -   \end{matrix} \right]
    \end{align}

Remarks
-------

 * LAMA uses always zero-based indexing within the array ``ja`` for the column indexes.
 * The arrays ``ja`` and ``values`` are always filled with zero values that allows for 
   faster matrix-vector multiplication on devices that work in an SIMD mode.
 * An explicit offset array is never needed as the offset can be computed by a closed formula.
 * The arrays ``ja`` and ``values`` are stored column-wise. 

    
ELLKernelTrait
--------------

Conversion
^^^^^^^^^^

========================= ============================================================= ==== ====
**Functionname**          **Description**                                               Host CUDA
========================= ============================================================= ==== ====
fillELLValues             fill up ja and values array                                   *    *
getCSRValues              ELL --> CSR: conversion ELL to CSR                            *    *
setCSRValues              CSR --> ELL: conversion CSR to ELL                            *    *
compressIA                compress the ia array by using values array and epsilon       *
compressValues            compress ja and values array by using epsilon                 *
getRow                    returns a row of the matrix                                   *    *
getValue                  get single element of matrix                                  *    *
countNonEmptyRowsBySizes  count non-empty rows by sizes array                           *    *
setNonEmptyRowsBySizes    set non-empty rows by sizes array                             *    *
========================= ============================================================= ==== ====

Calculation
^^^^^^^^^^^

========================= ============================================================= ==== ====
**Functionname**          **Description**                                               Host CUDA
========================= ============================================================= ==== ====
jacobi                    compute one jacobi iteration step                             *    *
jacobiHalo                compute one jacobi iteration step on halo values              *    *
normalGEMV                matrix-vector multiplication                                  *    *
sparseGEMV                matrix-vector multiplication with just non-zero rows          *    *
normalGEVM                vector-matrix multiplication                                  *    *
sparseGEVM                vector-matrix multiplication with just non-zero rows          *    *
absMaxVal                 compute the maximal absolute value                            *
scaleValue                scale with array                                              *    *
matrixMultiplySizes       computes row sizes for result of matrix multiplication        *
matrixAddSizes            computes row sizes for result of matrix addition              *
matrixAdd                 matrix-matrix addition (all ELL)                              *
matrixMultiply            matrix-matrix multiplication  (all ELL)                       *
========================= ============================================================= ==== ====

Properties
^^^^^^^^^^

========================= ============================================================= ==== ====
**Functionname**          **Description**                                               Host CUDA
========================= ============================================================= ==== ====
hasDiagonalProperty       Checks if the first n entries are the diagonal elements       *    *
check                     Checks integrity of ia array                                  *    *
========================= ============================================================= ==== ====

