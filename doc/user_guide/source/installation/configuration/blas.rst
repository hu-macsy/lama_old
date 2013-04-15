.. _blas:

BLAS and LAPACK
---------------

`Wikipedia BLAS`_.

.. _Wikipedia BLAS: http://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms

Basic Linear Algebra Subroutine (BLAS) is a de facto application programming interface standard for 
publishing libraries to perform basic linear algebra operations such as vector and matrix multiplication.

`Wikipedia LAPACK`_.

.. _Wikipedia LAPACK: http://en.wikipedia.org/wiki/LAPACK

LAPACK (Linear Algebra PACKage) is a software library for numerical linear algebra. It provides routines 
for solving systems of linear equations and linear least squares, eigenvalue problems, and singular value decomposition.
Usually implementations of LAPACK rely on efficient implementations of BLAS routines.

LAMA uses BLAS and LAPACK routines, mainly for operations on dense matrices, but also for operations on vectors.
Libraries implementing this standard are available on nearly all platforms and are mandatory for building LAMA.

Highly optimized implementations of the BLAS interface have been developed by hardware vendors 
such as Intel (MKL) and AMD, as well as by other authors, e.g. Goto BLAS and ATLAS (a portable self-optimizing BLAS).

We recommend using Intel MKL libraries for BLAS, LAPACK and (optional) ScaLAPACK routines. 
If the environment variable MKL_ROOT is set to the corresponding installation directory,
the configuration with cmake should find all needed libraries.
On some machines it was necessary to guide the configuration with the hint whether a 64-bit or a 32-bit version of
the MKL is needed; this is possible by the CMake variable ``MKL_Is64``:

.. code-block:: bash

    export MKL_ROOT=<path/to/mkl/installation>

or

.. code-block:: bash

    cmake -D MKL_ROOT=<path/to/mkl/installation> -D MKL_Is64=TRUE ...

If the MKL is not available, the LAMA configuration tries to find other 
BLAS or LAPACK libraries. If it fails, you will see a typical message like this::

    Could NOT find LAMA_BLAS (missing: LAMA_BLAS_LIBRARIES)

Unfortunately, the CMake modules have some problems in finding the BLAS and LAPACK libraries. 
Especially the module FindLAPACK.cmake has the problem that in some versions the CMAKE_FIND_LIBRARY_SUFFIXES 
variable was misspelled as CMAKE_FIND_LIBRRAY_SUFFIXES. Here is a solution of how to specify the corresponding libraries directly::

  cmake -D BLAS_blas_LIBRARY=/usr/lib64/libblas.so.3  \
        -D LAPACK_lapack_LIBRARY=/usr/lib64/liblapack.so.3

ScaLAPACK is optional and is used for building the inverse of a distributed dense matrix.
If it is not available, the inverse is calculated on a replicated matrix. 
