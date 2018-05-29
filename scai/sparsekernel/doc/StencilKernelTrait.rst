.. _sparsekernel_Stencil:

Stencil Storage Format (Stencil)
================================

The stencil storage format stands for a matrix that represents a linear mapping on
an n-dimenisonal grid (n = 1, 2, 3, or 4).
Each point of the grid is updated in the same way by the values of certain neighbored points.

Example
-------

.. math::

  A = \left(\begin{matrix} 
    2  & -1 &  0 &  0 &  0 &  0 &  0 \\
   -1  &  2 & -1 &  0 &  0 &  0 &  0 \\
    0  & -1 &  2 & -1 &  0 &  0 &  0 \\
    0  &  0 & -1 &  2 & -1 &  0 &  0 \\
    0  &  0 &  0 & -1 &  2 & -1 &  0 \\
    0  &  0 &  0 &  0 & -1 &  2 & -1 \\
    0  &  0 &  0 &  0 &  0 & -1 &  2 \\
    \end{matrix}\right) 

Here are the corresponding arrays for the representation:

.. math::
    
    \begin{align}
    grid &= [ 7 ]
    stencil &= [ -1:-1, 0:2, 1:1 ]
    \end{align}

Remarks
-------

 * A stencil storage can never be updated. 
 * Be careful about non-absorbing boundary conditions.
 * Conversion into a stencil storage is NOT supported.
 * Conversion of a stencil storage in other formats (especially CSR) is supported
 * Conversion into the DIA format is useful as the number of diagonals corresponds to the number
   of stencil points and the values array can be updated to deal with special bounary conditions.

StencilKernelTrait
------------------

The StencilKernelTrait contains various functions for operations on a stencil storage.
The most important ones are the conversion into the CSR format and the matrix-vector 
multiplication.

Conversion
^^^^^^^^^^

====================== ============================================================= ==== ====
**Functionname**       **Description**                                               Host CUDA
====================== ============================================================= ==== ====
convert2CSR            converts to CSR                                               *    *
====================== ============================================================= ==== ====

Calculation
^^^^^^^^^^^

====================== ============================================================= ==== ====
**Functionname**       **Description**                                               Host CUDA
====================== ============================================================= ==== ====
normalGEMV             matrix-vector multiplication                                  *    *
====================== ============================================================= ==== ====

