.. _main-page_sparsekernel:

#################
SCAI SparseKernel 
#################

*************
Specification
*************

* Provides wrapper for different SparseBLAS-libraries
* Contains own implementations for different sparse matrix formats
* Supports different backends (currently Host, CUDA, MIC)
* Provides a set of KernelTraits, grouped by sparse matrix formats
* Uses the Kernel Registry to register kernels
* Internal dependencies: common, logging, tracing, tasking, hmemo, kregistry, utilskernel
* External dependencies: [cuda --> cusparse], [mic --> mkl] 

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

*******
Classes
*******

===================   =======================================================
Name                  Description
===================   =======================================================
COOKernelTrait        routines for COO matrix storage
CSRKernelTrait        routines for CSR matrix storage
DenseKernelTrait      routines for Dense matrix storage
DIAKernelTrait        routines for DIA matrix storage
ELLKernelTrait        routines for ELL matrix storage
JDSKernelTrait        routines for JDS matrix storage
===================   =======================================================