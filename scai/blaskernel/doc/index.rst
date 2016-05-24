.. _main-page_blaskernel:

###############
SCAI BLASKernel 
###############

***********
Description
***********

The BLASKernel library is a wrapper for various blas and lapack implementations. Lots of vendors provide optimized implementations
for their hardware. For example the MKL is optimized for execution on Intel processors. This library can be used with different
libraries. If no library is found a integrated implementation is used. Currently not all blas and lapack calls are wrapped, but the
extension is very easy. 
   
********
Contents
********

.. toctree::
   :titlesonly:
   :maxdepth: 1
   
   Wrapping
   
*************
Relations
*************

The provided functions for the different back-ends register theirselfs in the KernelRegistry through the
BLASKernelTrait. 

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

- SCAI_USE_BLAS=[0|1|2|3]

  - usage of external blas library on host
   
  - parameters:

    - *0* --> don't use external blas library
      
    - *1* --> use level 1-3 from external library
      
    - *2* --> use level 2-3 from external library ( matrix-vector and matrix-matrix multiplication)
      
    - *3* --> use level 3 from external library (matrix-matrix multiplication)
      
- SCAI_USE_LAPACK=[0|1]

  - usage of external lapack library on host

************
Dependencies
************

- Internal:

  - :ref:`SCAI Common - Basic Concepts <scaicommon:main-page_common>`

  - :ref:`SCAI Logging - Logging Macros <scailogging:main-page_logging>`

  - :ref:`SCAI Tracing - Tracing Macros <scaitracing:main-page_tracing>`

  - :ref:`SCAI Tasking - Asynchronous Tasks <scaitasking:main-page_tasking>`

  - :ref:`SCAI Hmemo - Heterogeneous Memory Architecture <scaihmemo:main-page_hmemo>` (currently just used in testing)

  - :ref:`SCAI Kregistry - Generic Kernel Registry <scaikregistry:main-page_kregistry>`

- External:

  - CUBLAS
   
  - Intel MKL

Tested versions of external dependencies

.. toctree::
   :titlesonly:
   :maxdepth: 1
   
   Tested_Versions

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

