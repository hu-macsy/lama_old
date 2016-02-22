.. _main-page_blaskernel:

SCAI BLASKernel 
===============

Contents:

.. toctree::
   :titlesonly:
   :maxdepth: 1
   
   naming_scheme_wrapping
   extension
   tested_version
   
   
************
Introduction
************

* Provides wrapper for different BLAS-libraries
* Provides wrapper for different LAPACK-libraries
* Supports different backends (currently Host, CUDA, MIC)
* Provides the BLASKernelTrait
* Uses the Kernel Registry to register kernels
* Internal dependencies: common, logging, tracing, tasking, hmemo, kregistry
* External dependencies: [cuda --> cublas], [mic --> mkl], [blas], [lapack] 
* Tested BLAS- and LAPACK-Libraries are listed here: :ref:`blaskernel-tested-versions`

The library doesn't cover all blas functionallity, yet. Currently just the ones are 
accessible which are used inside LAMA.  

