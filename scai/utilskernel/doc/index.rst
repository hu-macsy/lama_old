.. _main-page_utilskernel:

################
SCAI UtilsKernel 
################

********
Contents
********

.. toctree::
   :titlesonly:
   :maxdepth: 1
   
   LArray
   LAMAKernel
   
*************
Specification
*************

* Contains utilities for the different backends
* Supports different backends (currently Host, CUDA, MIC)
* Provides the UtilKernelTrait
* Uses the Kernel Registry to register kernels
* Internal dependencies: common, logging, tracing, tasking, hmemo, kregistry
