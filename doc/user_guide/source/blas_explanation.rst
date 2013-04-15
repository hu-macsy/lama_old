.. _blas_explanation:

Brief Introduction into BLAS
============================
 
`Wikipedia Original`_.

.. _Wikipedia Original: http://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms

**Basic Linea Algebra Subprograms** (BLAS) is a de facto application programming interface standard
for publishing libraries to perform basic linear algebra operations such as vector and matrix multiplication.
They were first published in 1979, and are used to build larger packages such as LAPACK. Heavily used in
high-performance computing, highly optimized implementations of the BLAS interface have been developed by
hardware vendors such as Intel and AMD, as well as by other authors, e.g. Goto BLAS and ATLAS (a portable
self-optimizing BLAS). The LINPACK benchmark relies heavily on DGEMM, a BLAS subroutine, for its performance.

Functionality
-------------

Level 1
^^^^^^^^
Level 1 contains vector operations of the form

.. math::
    y := a x + y 

as well as scalar dot products and vector norms, among other things.

Level 2
^^^^^^^^
Level 2 contains vector-matrix operations of the form

.. math::
    y := a A * x + b y


Level 3
^^^^^^^^
Level 3 contains matrix-matrix operations of the form

.. math::
   C := a A * B + b C

as well as solving 

.. math::
   B := a T^{-1} * B 
    
for triangular matrices T, among other things.
This level contains the widely used General Matrix Multiply operation.