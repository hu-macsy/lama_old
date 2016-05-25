.. _main-page_utilskernel:

################
SCAI UtilsKernel 
################

***********
Description
***********

The UtilsKernel library provides various functions for array operations. Additionally it contains an extension to the 
HArray provided by the HMemo library called LArray which allows to execute some of these functions directly on an array.

********
Contents
********

.. toctree::
   :titlesonly:
   :maxdepth: 1
   
   LArray
   LAMAKernel
   ReductionOp

*********
Relations
*********

The LArray class is a extension of the HArray class which is located in the HMemo library. 
Besides this the UtilsKernel library extends the KernelTraitContextFunction which comes from KRegistry to work with ContextPtr
insted of the ContextType. The provided array functions are registered in the KernelRegistry. 

*******
Example
*******

The following example shows how you can operate on the LArray:

.. code-block:: c++

  LArray<float> x( n );
  
  /*
   * ... fill x with data
   */ 
   
   float l1 = x.l1Norm(); // not possible on _HArray / HArray
   float l2 = x.l2Norm();

This examples shows how to get a ContextPtr from the Context-factory and use it to
get a function from the KernelRegistry:

.. code-block:: c++

  // Get ContextPtr
  hmemo::ContextPtr loc = hmemo::Context::getContextPtr( common::context::Host );
  
  // Get function from Registry
  LAMAKernel<UtilKernelTrait::setVal<double> > setVal;
  
  // Call function
  setVal[loc]( b, size, double( 2 ), reduction::COPY );  

Another example of the LAMAKernel is shown :ref:`here <utilskernel_lamakernel_example>`.

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

  - :ref:`SCAI BLASKernel - BLAS Kernel <scaiblaskernel:main-page_blaskernel>`
  
- External:

  - Thrust

************
Related Work
************
   
