.. _main-page_sparsekernel:

#################
SCAI SparseKernel 
#################

***********
Description
***********

The SparseKernel library is used to group functionality for different sparse matrix formats.
Currently the COO, CSR, DIA, ELL and JDS formats are supported. The best supported format is 
the CSR format. Every of the other formats contains operations to copy to and from CSR. 
They are additional provided for dense storage. Besides operations for conversion, operations 
for calculations and property checking are included. The functions are 
provided for different backends. Every format supports the matrix-vector multiplication. 
The matrix-matrix multiplication is currently just supported by CSR and ELL. 

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

*********
Relations
*********

To avoid code duplications some operations rely on code provided by the
UtilsKernel-library. For every format a corresponding KernelTrait exists. Through this the
functions are registered in the KernelRegistry. 

*******
Example
*******

Calculate matrix-vector multiplication with given ELL input-set on host.

.. code-block:: c++

	// Get function from Registry
	kregistry::KernelTraitContextFunction<ELLKernelTrait::normalGEMV<ValueType> > gemv;
	
	// Call function on host
	gemv[context::Host]( b, one, x, zero, b, m, max_nnz, sizes, ja, values );
	

*********************
Environment-Variables
*********************

- SCAI_CUDA_USE_CUSPARSE=[0|1]

  - used to turn on/off the usage of cusparse for the csr format on cuda
  
  - default value: off 

- SCAI_USE_MKL=[0|1]

  - used to turn on/off the usage of mkl for the csr format on host and mic
  
  - default value: off

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

  - :ref:`SCAI UtilsKernel - Utils Kernel <scaiutilskernel:main-page_utilskernel>`

- External:

  - Intel MKL

  - CUSPARSE

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
  
    - provides just functions on x86 CPUs 

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
		 
		 