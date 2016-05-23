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

