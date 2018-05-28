.. _main-page_sparsekernel:

#################
SCAI SparseKernel 
#################

***********
Description
***********

The SparseKernel library provides kernel routines for different sparse matrix formats.
Currently the COO, CSR, DIA, ELL and JDS formats are supported. The best supported format is 
the CSR format. For each other format there are routines to convert to and from CSR, 
this is also supported for the dense format.
Besides operations for conversion, operations 
for matrix-vector and matrix-matrix computations and for property checkings are included.
Every format supports the matrix-vector multiplication. 
The matrix-matrix multiplication is currently only supported for the CSR, ELL, and dense format.

Like all other kernel routines in the LAMA framework, the sparse matrix kernels are registered 
in the kernel registry and are always implemented on the Host device but also on all other backends wherever possible. 
For every sparse format and for the dense format a corresponding KernelTrait structure groups the
supported routines of the format and specifies name and signature of each routine. Within the library,
these KernelTraits are used for the registration of the implementations. When using the SparseKernel library,
these KernelTraits are used to access the different kernel implementations on the supported backends.


********
Contents
********

.. toctree::
   :titlesonly:
   :maxdepth: 1
   
   COOKernelTrait
   CSRKernelTrait
   DenseKernelTrait
   DIAKernelTrait
   ELLKernelTrait
   JDSKernelTrait
   StencilKernelTrait

*******
Example
*******

Calculate matrix-vector multiplication with a matrix stored in CSR format.

.. math::

  A = \left(\begin{matrix} 
    6  &  1 &  0 &  4  \\
    -2 &  8 &  3 &  0  \\
     0 &  1 &  9 &  4  \\
     -1 &  0 &  1 &  3  \\
    \end{matrix}\right) 

.. code-block:: c++

    ContextPtr ctx = Context::getContextPtr();  // default context
	
    HArray<IndexType> csrIA(  {  0,          3,          6,          9,         12 } );
    HArray<IndexType> csrJA(  {  0,  1,  3,  1,  0,  2,  2,  1,  3,  3,  0,  2 } );
    HArray<double> csrValues( {  6,  1,  4,  8, -2,  3,  9,  1,  4,  3, -1,  1 } );

    HArray<double> x( 4, 4.0 );
    HArray<double> b( 4, 2.0 );

    // compute b = b + A * x 

    CSRUtils::gemv( b, 1, x, 1, b, csrIA, csrJA, csrValues, ctx );

The CSR arrays might also be computed from the dense array.

.. code-block:: c++

    ContextPtr ctx = Context::getContextPtr();  // default context
    
    HArray<double> values( {  6,  1,  0, 4,  -2,  8,  3, 0 , 0,  1,  9,  4, -1, 0, 1, 3 }, ctx );

    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<double> csrValues;

    DenseUtils::compressSparseRows( csrIA, csrJA, csrValues, denseValues, ctx );

    HArray<double> x( 4, 4.0 );
    HArray<double> b( 4, 2.0 );

    // compute b = b + A * x 

    CSRUtils::gemv( b, 1, x, 1, b, csrIA, csrJA, csrValues, ctx );

*********************
Environment-Variables
*********************

- ``SCAI_CUDA_USE_CUSPARSE [bool, default = 1]``

  - if enabled, wrapped routines of the cuSPARSE library are taken for CSR kernel routines
    on the CUDA device wherever possible
  - if disabled, own CUDA kernel implementations are taken for CSR kernel routines

- ``SCAI_USE_MKL [:bool, default = 1]``

  - if enabled, wrapped routines of the MKL library are taken for CSR kernel routines
    on the Host wherever possible
  - if disabled, own OpenMP implementations are taken for CSR kernel routines

The environment variable ``SCAI_USE_MKL`` has no effect for the selection of BLAS or LAPACK kernels.

  
************
Dependencies
************

Internal dependencies:

- :ref:`SCAI Common - Basic Concepts <scaicommon:main-page_common>`

- :ref:`SCAI Logging - Logging Macros <scailogging:main-page_logging>`:
  The logging features of the LAMA framework is used in all subprojects and also here. The general
  strategy of logging is here to have at least one INFO logging message for each kernel routine.

- :ref:`SCAI Tracing - Tracing Macros <scaitracing:main-page_tracing>`:
  Each of the kernels is traced as a region.

- :ref:`SCAI Tasking - Asynchronous Tasks <scaitasking:main-page_tasking>`:
  Most of the computational kernels  of the SparseKernel library can run asynchronously.

- :ref:`SCAI Hmemo - Heterogeneous Memory Architecture <scaihmemo:main-page_hmemo>` (currently just used in testing)

- :ref:`SCAI Kregistry - Generic Kernel Registry <scaikregistry:main-page_kregistry>`:
  All kernels  of the SparseKernel library are registered in the kernel registry.

- :ref:`SCAI UtilsKernel - Utils Kernel <scaiutilskernel:main-page_utilskernel>`:
  Some kernels of the SparseKernel library use for their implementation kernels of the UtilsKernel library
  (e.g. running sums of arrays). 

External dependencies:

- The :ref:`Intel MKL <scaiblaskernel:BLAS>` does not only provide implementations for BLAS and LAPACK
  routines but also routines for the CSR storage format. If the MKL is not available, own OpenMP implementations
  are used. The MKL can be used on the Host device.

- :ref:`cuSPARSE library of the CUDA toolkit <scaicommon:CUDA>` is used for the implementation of CSR kernel
  routines on NVidia GPU devices.

************
Related Work
************

There are different libraries which provides similar functionality. Some of them can be used
inside of SparseKernel.

- Software

  - |CUSPARSE|

    - COO, CSR, CSC, HYB (mixed ELL/COO), BSR, BSRX
  
    - provides just functions for Nvidia GPUs
    
  - |CUSP|
  
    - COO, CSR, DIA, ELL, HYB, Permutation
    
    - useable on x86 CPUs and Nvidia GPUs
    
  - |clSPARSE|
  
    - COO, CSR 
  
    - useable on devices that can be utilized using OpenCL
  
  - |MKL|
 
    - BSR, COO, CSC, CSR, DIA, SKY
  
    - provides just functions on x86 CPUs where most of the functions are OpenMP parallelized.

- Papers / Books

  - Yousef Saad |SAAD| 
   
    - describes various sparse storage formats ( COO, CSC, CSR, MSR, DIA, ELL, JDS ) 
   
   
.. |CUSPARSE| raw:: html

        <a href="http://docs.nvidia.com/cuda/cusparse" target="_blank">CUSPARSE</a>   
   
.. |MKL| raw:: html

		 <a href="https://software.intel.com/en-us/articles/intel-mkl-sparse-blas-overview" target="_blank">MKL</a>
		 
.. |CUSP| raw:: html

		 <a href="https://github.com/cusplibrary/cusplibrary" target="_blank">CUSP</a>
		 		 
.. |clSPARSE| raw:: html

		 <a href="https://github.com/clMathLibraries/clSPARSE" target="_blank">clSPARSE</a>
		 
.. |SAAD| raw:: html

		 <a href="http://www-users.cs.umn.edu/~saad/books.html" target="_blank">Iterative Methods for Sparse Linear Systems</a>
		 
