.. _main-page_blaskernel:

###############
SCAI BLASKernel 
###############

***********
Description
***********

The BLASKernel library provides wrapper for various BLAS and LAPACK routines in such a way that these routines can 
be called as kernel functions on different devices with the same signature.

The BLASKernel library supports both approaches: either the BLAS or LAPACK routine is implemented by the library
itself (internal version) or the BLAS or LAPACK routine of an efficient vendor library is wrapped (external version).
The latter approach is usually the more efficient one as vendors provide optimized implementations on their
hardware. This is not only the case for usual host CPUs but also for CUDA devices, where the cuBLAS library is available.

The internal version of the BLAS and LAPACK routines is a fallback implementation on Host CPUs if no external BLAS or
LAPACK library is available. In contrary to most vendor libraries the routines can be instantiated for any arithmetic value
type, i.e. also for long or long complex value types. All routines of the internal version will be registered in 
the kernel registry, but will be replaced with wrapper routines to external vendor libraries if available.

The current version of the BLASKernel library supports only a subset of the BLAS and LAPACK routines that are available.
This subset contains all routines that are currently used in other libraries like the LAMA and solver library. 
But it is easy to extend the BLASKernel library for further routines.

********
Contents
********

.. toctree::
   :titlesonly:
   :maxdepth: 1
   
   Kernels
   Supported
   Wrapping
   
*********
Relations
*********

The provided functions for the different back-ends register theirselfs in the KernelRegistry through the BLASKernelTrait. 

*******
Example
*******

How to use registered functions from the KernelRegistry:

.. toctree::
   :titlesonly:
   :maxdepth: 1

   Usage

How to extend the supported blas functions:

.. toctree::
   :titlesonly:
   :maxdepth: 1

   Extension

*********************
Environment-Variables
*********************

- SCAI_USE_BLAS=[0|1|2|3] (default=0)

  - usage of external blas library on host
   
  - parameters:

    - *0* --> don't use external blas library
      
    - *1* --> use level 1-3 from external library
      
    - *2* --> use level 2-3 from external library ( matrix-vector and matrix-matrix multiplication)
      
    - *3* --> use level 3 from external library (matrix-matrix multiplication)
      
- SCAI_USE_LAPACK=[0|1] (default=0)

  - usage of external lapack library on host

************
Dependencies
************

Internal dependencies:

- :ref:`SCAI Common - Basic Concepts <scaicommon:main-page_common>`
- :ref:`SCAI Logging - Logging Macros <scailogging:main-page_logging>`
- :ref:`SCAI Tracing - Tracing Macros <scaitracing:main-page_tracing>`
- :ref:`SCAI Tasking - Asynchronous Tasks <scaitasking:main-page_tasking>`
- :ref:`SCAI Hmemo - Heterogeneous Memory Architecture <scaihmemo:main-page_hmemo>` (currently just used in testing)
- :ref:`SCAI Kregistry - Generic Kernel Registry <scaikregistry:main-page_kregistry>`

External dependencies:

- :ref:`BLAS`
- :ref:`CUBLAS as part of the CUDA SDK toolkit <scaicommon:CUDA>`

Tested versions of external dependencies

.. toctree::
   :hidden:
   
   blas

************
Related Work
************

- Software

  - |MKL|
  
    - BLAS/LAPACK implementation optimized for Intel x86 CPUs
  
  - |CUBLAS|
  
    - BLAS implementation for Nvidia GPUs
  
  - |clBLAS|
  
    - BLAS implementation for OpenCL devices
  
  - |OpenBLAS|
  
    - free and parallel BLAS implementation for various CPU architectures
  
  - |NetLibBLAS|
  
    - free but not parallel BLAS implementation for various CPU architectures
  
  - |NetLibLAPACK|
  
    - free but not parallel LAPACK implementation for various CPU architectures
  
- Papers / Books

  - Basic Linear Algebra Subprograms for Fortran Usage by Lawson et. al.
  
  - An Extended Set of FORTRAN Basic Linear Algebra Subprograms by Dongarra et. al.
  
  - A Set of Level 3 Basic Linear Algebra Subprograms by Dongarra et. al.
  

.. Links
  
.. |MKL| raw:: html

    <a href="https://software.intel.com/en-us/intel-mkl" target="_blank">MKL</a>
    
.. |CUBLAS| raw:: html

    <a href="http://docs.nvidia.com/cuda/cublas/" target="_blank">CUBLAS</a>   
   
.. |clBLAS| raw:: html

    <a href="https://github.com/clMathLibraries/clBLAS" target="_blank">clBLAS</a>

.. |OpenBLAS| raw:: html

    <a href="http://www.openblas.net/" target="_blank">OpenBLAS</a>

.. |NetLibBLAS| raw:: html

    <a href="http://www.netlib.org/blas/index.html" target="_blank">NetLib BLAS</a>

.. |NetLibLAPACK| raw:: html

    <a href="http://www.netlib.org/lapack/index.html" target="_blank">NetLib LAPACK</a>


.. |GOTOBLAS| raw:: html

    <a href="http://c2.com/cgi/wiki?GotoBlas" target="_blank">Goto BLAS</a>


.. |ATLAS| raw:: html

    <a href="http://math-atlas.sourceforge.net/" target="_blank">ATLAS</a>
