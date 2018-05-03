.. _BLAS:

BLAS and LAPACK
^^^^^^^^^^^^^^^

Basic Linear Algebra Subroutine (BLAS) is a de facto application programming interface (API) standard for publishing libraries to perform basic linear algebra operations such as matrix-vector-multiplication.

LAPACK (Linear Algebra PACKage) is a software library for numerical linear algebra. It provides routines for solving systems of linear equations and linear least squares, eigenvalue problems, and singular value decomposition. Usually implementations of LAPACK rely on efficient implementations of BLAS routines.

LAMA uses BLAS and LAPACK routines for operations on dense matrices and dense vectors.
Libraries implementing this standard are available on nearly all platforms and we recommend to use one of them for building LAMA.
Alternatively we have an internal implementation if these libraries are not available on your system.

CMake will look automatically for available libraries in our recommended order:

* |MKL|
* Generic BLAS (CMake supported: |NetLibBLAS| with |NetLibLAPACK|)
* internal LAMA implementation

You may also make the choice for CMake by defining SCAI_BLAS_LIBRARY to MKL, BLAS or INTERNALBLAS for the respective installation. NOTE: If you choose MKL or BLAS and it is not found by CMake the configuration is incomplete.

There are also many other BLAS implementations, e.g. |GotoBLAS| (development has been abandoned), |OpenBLAS| (fork of Goto BLAS) and |ATLAS| (a portable self-optimizing BLAS). Due to standardized API they also should work, but sometimes its a little bit difficult that CMake will find it. Look for the instructions below for :ref:`generic_blas`.

Intel MKL
"""""""""

If your Intel MKL is not installed in the system path, you can give CMake hints to find it with **MKL_ROOT** via environment variable or by passing it with -D to the cmake call.

On some machines it is necessary to guide the configuration with the hint whether a 64-bit or a 32-bit version of the MKL is needed; this is possible by the CMake variable ``MKL_Is64``:

.. code-block:: bash

    export MKL_ROOT=<path/to/mkl/installation>

or

.. code-block:: bash

    cmake -DMKL_ROOT=<path/to/mkl/installation> -DMKL_Is64=TRUE ...

.. _generic_blas:

Generic BLAS
""""""""""""

If CMake does not find a generic BLAS library on its own, you can give hints by passing **BLAS_blas_LIBRARY** and **LAPACK_lapack_LIBRARY**. This may be the case if you have another library than the Netlib implementation, but also with older CMake installation (unfortunately the module FindLAPACK.cmake has the problem that in some versions the CMAKE_FIND_LIBRARY_SUFFIXES variable was misspelled as CMAKE_FIND_LIBRRAY_SUFFIXES). In all cases try to specify the libraries directly:

.. code-block:: bash

  cmake -DBLAS_blas_LIBRARY=/usr/lib64/libblas.so.3 -DLAPACK_lapack_LIBRARY=/usr/lib64/liblapack.so.3


.. |MKL| raw:: html

		 <a href="https://software.intel.com/en-us/intel-mkl" target="_blank">Intel MKL</a>
		 
.. |OpenBLAS| raw:: html

		 <a href="http://www.openblas.net/" target="_blank">OpenBLAS</a>

.. |GOTOBLAS| raw:: html

		 <a href="http://c2.com/cgi/wiki?GotoBlas" target="_blank">Goto BLAS</a>


.. |ATLAS| raw:: html

		 <a href="http://math-atlas.sourceforge.net/" target="_blank">ATLAS</a>

.. |NetLibBLAS| raw:: html

		 <a href="http://www.netlib.org/blas/index.html" target="_blank">NetLib BLAS</a>

.. |NetLibLAPACK| raw:: html

		 <a href="http://www.netlib.org/lapack/index.html" target="_blank">NetLib LAPACK</a>

.. _blaskernel-tested-versions:

Tested BLAS Libraries
"""""""""""""""""""""

The following blas libraries has been tested. The blaskernel library contains a own implementation of the blas function
which is called OpenMP-BLAS. 

==============    ================   ===========================================================================   ======================================
BLAS-library      Backend            Supported ValueTypes                                                          Status 
==============    ================   ===========================================================================   ======================================
Fortran-BLAS      Host               float, double, ComplexFloat, ComplexDouble                                    successfully tested
MKL               Host               float, double, ComplexFloat, ComplexDouble                                    successfully tested
OpenMP-BLAS       Host               float, double, long double, ComplexFloat, ComplexDouble, ComplexLongDouble    successfully tested
cuBLAS            CUDA               float, double, ComplexFloat, ComplexDouble                                    successfully tested
==============    ================   ===========================================================================   ======================================

Tested LAPACK Libraries
"""""""""""""""""""""""

The following lapack libraries has been tested. The blaskernel library contains a own implementation of the lapack function
which is called OpenMP-LAPACK. 

===============    ================   ===========================================================================   ======================================
LAPACK-library      Backend            Supported ValueTypes                                                          Status 
===============    ================   ===========================================================================   ======================================
Fortran-LAPACK      Host               float, double, ComplexFloat, ComplexDouble                                    successfully tested
MKL                 Host               float, double, ComplexFloat, ComplexDouble                                    successfully tested
OpenMP-LAPACK       Host               float, double, long double, ComplexFloat, ComplexDouble, ComplexLongDouble    successfully tested
===============    ================   ===========================================================================   ======================================
