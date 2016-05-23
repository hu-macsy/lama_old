.. _main-page_sparsekernel:

#################
SCAI SparseKernel 
#################

**********
Motivation
**********

The SparseKernel library is used to group functionality for different sparse matrix formats.
Currently the COO, CSR, DIA, ELL and JDS formats are supported. Additionally some routines 
for converting dense to CSR are included. The functions are provided for different backends. 
Every format supports the matrix-vector multiplication. The matrix-matrix multiplication
is currently just supported by CSR and ELL. 

*********
Contents:
*********

.. toctree::
   :titlesonly:
   :maxdepth: 1
   
   COOKernelTrait
   CSRKernelTrait
   DenseKernelTrait
   DIAKernelTrait
   ELLKernelTrait
   JDSKernelTrait

*************
Relationships
*************

The best supported format is the CSR format. Every of the other formats contains operations to
copy to and from CSR. To avoid code duplications some operations rely on code provided by the
UtilsKernel-library. For every format a corresponding KernelTrait exists. Through this the
functions are registered in the KernelRegistry. 

*******
Example
*******

Calculate matrix-vector multiplication with given ELL input-set on host.

.. code-block:: c++

	// Get Function from Registry
	KernelTraitContextFunction<ELLKernelTrait::normalGEMV<ValueType> > gemv;
	
	gemv[context::Host]( b, one, x, zero, b, m, max_nnz, sizes, ja, values );
	

*********************
Environment-Variables
*********************

- SCAI_CUDA_USE_CUSPARSE

  - used to turn on/off the usage of cusparse for the csr format on cuda
  
  - default value: off 

- SCAI_USE_MKL [host and mic]

  - used to turn on/off the usage of mkl for the csr format on host and mic
  
  - default value: off

************
Dependencies
************

- Internal:

  - common

  - logging

  - tracing

  - tasking

  - hmemo (currently just used in testing)

  - kregistry

  - utilskernel

- External:

  - Intel MKL

  - CUSPARSE

************
Related Work
************

There are different libraries which provides similar functionality. Some of them can be used
inside of SparseKernel.

- Software

  - `CUSPARSE`_

    - COO, CSR, CSC, HYB (mixed ELL/COO), BSR, BSRX
  
    - provides just functions for Nvidia GPUs
    
  - `CUSP`_
  
    - COO, CSR, DIA, ELL, HYB, Permutation
    
    - useable on x86 CPUs and Nvidia GPUs
    
  - `clSPARSE`_
  
    - COO, CSR 
  
    - useable on devices that can be utilized using OpenCL
  
  - `MKL`_
 
    - BSR, COO, CSC, CSR, DIA, SKY
  
    - provides just functions on x86 CPUs 

- Books

  - Yousef Saad `Iterative Methods for Sparse Linear Systems`_ 
   
    - describes various sparse storage formats ( COO, CSC, CSR, MSR, DIA, ELL, JDS ) 
   
   
.. _CUSPARSE : http://docs.nvidia.com/cuda/cusparse
.. _MKL : https://software.intel.com/en-us/articles/intel-mkl-sparse-blas-overview
.. _CUSP : https://github.com/cusplibrary/cusplibrary
.. _clSPARSE : https://github.com/clMathLibraries/clSPARSE
.. _Iterative Methods for Sparse Linear Systems : http://www-users.cs.umn.edu/~saad/books.html