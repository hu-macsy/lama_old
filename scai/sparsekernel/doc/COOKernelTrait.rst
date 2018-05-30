.. _sparsekernel_COO:

Coordinate Storage Format (COO)
===============================

The coordinate format represents the nonzero values of a matrix by tuples (row, column, value). Instead of
tuples one can also use three arrays that have all the same size.

- the COO ``values`` array contains the non-zero values, its size is exactly ``numValues``
- the COO ``ja`` array contains the corresponding column indices for the values
- the COO ``ia`` array contains the corresponding row indices for the values

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

The COO representation of the sparse matrix is:

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
  ia_{A} = \left(\begin{matrix} 0 & 0 \\
    1  \\
    2 & 2 \\
    3 & 3 & 3 \\
    4 & 4 \\
    \\
    6 & 6 \end{matrix}\right) 

Here are the corresponding arrays for the representation:

.. math::
    
    \begin{align}
    numRows &= 7 \\
    numColums &= 4 \\
    numValues &= 12 \\
    ia &= \left[ \begin{matrix} 0 & 0 & 1 & 2 & 2 & 3 & 3 & 3 & 4 & 4 & 6 & 6 \end{matrix} \right] \\
    ja &= \left[ \begin{matrix}  0 & 3 & 0 & 2 & 3 & 0 & 1 & 3 & 0 & 3 & 1 & 3 \end{matrix}\right] \\
    values &= \left[ \begin{matrix} 6 & 4 & 7 & -9 & 4 & 2 & 5 & 3 & 2 & 1 & 1 & 2 \end{matrix}\right] 
    \end{align}

Remarks
-------

 * LAMA uses always zero-based indexing within the arrays ``ia`` and ``ja`` for the row and the column indexes.
 * The entries will always be sorted, first by row-index and then by column-index.
   Most operations for COO data can be implemented much faster if the data is sorted.
 * For operations on COO data is assumed that there will be never two entries for the same matrix position (i,j).
   There is a method that eliminates double entries in COO data, where entries might either be replaced
   or combined by a binary operation (usually summing it up).
 * Conversion beween COO and CSR is very simple if entries are sorted. Only the 
   COO ``ia`` array has to be converted to an CSR ``ia`` offset array and vice versa.

Within LAMA the COO format is also used for assembling matrix entries.

COOKernelTrait
--------------

The COOKernelTrait contains various function which operates on the introduced COO-format. 
These functions are grouped into conversion, caculation und properties. The following tables show 
which function has been implemented on which back-end.

Conversion
^^^^^^^^^^

====================== ============================================================= ==== ====
**Functionname**       **Description**                                               Host CUDA
====================== ============================================================= ==== ====
offsets2ia             CSR --> COO: converting offset array (CSR) to ia array (COO)  *    *
sort                   sorts the entries, first row-wise then by columns             *
unique                 eliminates double entries in sorted COO data                  *    *
====================== ============================================================= ==== ====

Calculation
^^^^^^^^^^^

====================== ============================================================= ==== ====
**Functionname**       **Description**                                               Host CUDA
====================== ============================================================= ==== ====
scaleRows              multiplies each row with an own value                         *
normalGEMV             matrix-vector multiplication                                  *    *
jacobi                 compute one jacobi iteration step                             *
jacobiHalo             compute one jacobi iteration step on halo values
====================== ============================================================= ==== ====

Properties
^^^^^^^^^^

====================== ============================================================= ==== ====
**Functionname**       **Description**                                               Host CUDA
====================== ============================================================= ==== ====
isSorted               Checks if the entries are sorted by (rows, columns)           *
====================== ============================================================= ==== ====

