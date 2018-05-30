.. _sparsekernel_CSR:

Compressed Sparse Row Format (CSR)
==================================

The compressed sparse row (CSR) format represents the nonzero values of a matrix by three arrays:

- the CSR ``values`` array contains the non-zero values stored row-wise, its size is exactly
  ``numValues`` 
- the CSR ja array contains the corresponding column indices for the values
- the CSR ia array contains the running sum (offsets) of the  row sizes; it has numRows + 1 entries, the first one is always 0 and the last one
  is numValues

In contrary to an array ``ia`` that contains only the sizes of each row, the offset array is very helpful for 
more efficient operations. Its main advantage is that it allows the row parallelization of many operations, 
i.e. each thread can operate independently on a bundle of rows.

.. code-block:: c++

   #pragma omp parallel for
   for ( IndexType i = 0; i < numRows; ++i )
   {
       // code that operates on row i of the storage

       for ( IndexType k = ia[i]; k < ia[i+1]; ++k )
       {
           ... csrJA[k] for column pos, csrValues[k] for value
       }
   }

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

The CSR representation of the sparse matrix is:

.. math::

  values_{A} = \left(\begin{matrix} 6 & 4 \\
    7  \\
    -9 & 4 \\
    2 & 5 & 3 \\
    2 & 1 \\
    \\
    1 & 2 \end{matrix}\right) 
  ja_{A} = \left(\begin{matrix} 0 & 3 \\
    0  \\
    2 & 3 \\
    0 & 1 & 3 \\
    0 & 3 \\
    \\
    1 & 3 \end{matrix}\right) 
  sizes_{A} = \left(\begin{matrix} 2 \\
    1  \\
    2 \\
    3 \\
    2 \\
    0 \\
    2 \end{matrix}\right) 
  offsets_{A} = \left(\begin{matrix} 0 \\
    2 \\
    3  \\
    5 \\
    8 \\
    10 \\
    10 \\
    12 \end{matrix}\right) 

Here are the corresponding arrays for the representation:

.. math::
    
    \begin{align}
    numRows &= 7 \\
    numColums &= 4 \\
    numValues &= 12 \\
    ia &= [\begin{matrix} 0 & 2 & 3 & 5 & 8 & 10 & 10 & 12 \end{matrix}] \\
    ja &= [\begin{matrix} 0 & 3 & 0 & 2 & 3 & 0 & 1 & 3 & 0 & 3 & 1 & 3 \end{matrix}] \\
    values &= [\begin{matrix} 6 & 4 & 7 & -9 & 4 & 2 & 5 & 3 & 2 & 1 & 1 & 2 \end{matrix}]
    \end{align}

Remarks
-------

 * LAMA uses always zero-based indexing.
 * The entries for one row do not have to be sorted by column indices. But 
   there are some operations (matrix add, matrix-multiplication, element-wise binary
   operations) that will be faster if they are sorted. 
 * It is assumed that there will be never two entries for the same matrix position (i,j).

CSRKernelTrait
--------------

The CSRKernelTrait contains various function which operates on the introduced CSR-format. 
These functions are grouped into conversion, caculation und properties. For the CSR-format
we offer the most functionality compared to the other formats. The following tables show 
which function has been implemented on which back-end.

Conversion
^^^^^^^^^^

====================== ============================================================= ==== ====
**Functionname**       **Description**                                               Host CUDA
====================== ============================================================= ==== ====
sortRowElements        sorts the elements of a row by increasing column indexes      *
sizes2offsets          computes offset array from sizes array                        *    *
offsets2sizes          computes sizes array from offset array                        *    *
convertCSR2CSC         converts from CSR2CSC                                         *    *
compress               fill compresses CSR data in new data structures               *
====================== ============================================================= ==== ====

Calculation
^^^^^^^^^^^

====================== ============================================================= ==== ====
**Functionname**       **Description**                                               Host CUDA
====================== ============================================================= ==== ====
jacobi                 compute one jacobi iteration step                             *    *
jacobiHalo             compute one jacobi iteration step on halo values              *    *
matrixAddSizes         computes row sizes for result of matrix addition              *    *
matrixMultiplySizes    computes row sizes for result of matrix multiplication        *    *
matrixMultiplyJA       computes column indexes for result of matrix multiplication   *
scaleRows              multiplies each row with an own value                         *    *
absMaxDiffVal          computes the maximal element-wise difference for two matrices *
normalGEMV             matrix-vector multiplication                                  *    *
normalGEVM             vector-matrix multiplication                                  *    *
sparseGEMV             matrix-vector multiplication with just non-zero rows          *    *
sparseGEVM             vector-matrix multiplication with just non-zero rows          *    *
gemm                   matrix-matrix multiplication (CSR * Dense)                    *
matrixAdd              matrix-matrix addition (all CSR)                              *    *
matrixMultiply         matrix-matrix multiplication  (all CSR)                       *    *
====================== ============================================================= ==== ====

Properties
^^^^^^^^^^

====================== ============================================================= ==== ====
**Functionname**       **Description**                                               Host CUDA
====================== ============================================================= ==== ====
validOffsets           checks for legal offset array                                 *
====================== ============================================================= ==== ====

