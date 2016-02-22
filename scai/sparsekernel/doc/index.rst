.. _main-page_sparsekernel:

SCAI SparseKernel 
===============

Contents:

.. toctree::
   :titlesonly:
   :maxdepth: 1
   
   
   
************
Introduction
************

* Provides wrapper for different SparseBLAS-libraries
* Contains own implementations for different sparse matrix formats
* Supports different backends (currently Host, CUDA, MIC)
* Provides a set of KernelTraits, grouped by sparse matrix formats
* Uses the Kernel Registry to register kernels
* Internal dependencies: common, logging, tracing, tasking, hmemo, kregistry, utilskernel
* External dependencies: [cuda --> cusparse], [mic --> mkl] 

