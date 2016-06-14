.. _main-page_utilskernel:

################
SCAI UtilsKernel 
################

***********
Description
***********

The UtilsKernel library provides various kernel implementations for array operations. 
These implementations are registered in the KernelRegistry. 

Additionally it contains a class LArray that is derived from the 
class HArray provided by the HMemo library. This derived class provides the above array operations
as operators or methods of this class and makes so the use of heterogeneous arrays more convenient.

Besides this, the UtilsKernel library extends the class KernelTraitContextFunction from the KRegistry library.
This extension uses directly a context (ContextPtr) insted of the context type. 

********
Contents
********

.. toctree::
   :titlesonly:
   :maxdepth: 1
   
   LArray
   LAMAKernel
   ReductionOp

*******
Example
*******

LArray
------

The following example shows how you can operate on local arrays:

.. code-block:: c++

  const int N = 10;

  LArray<float> x( N );
  LArray<float> y( N );
  
  for ( IndexType i = 0; i < N; ++i )
  {
      x[i] = i;
  }

  y = x;
  y *= y;
  y *= 2.0f;
   
  float l1 = x.l1Norm(); // not possible on _HArray / HArray
  float l2 = x.l2Norm();

  float sum = x.sum();   // build the sum of all elements
  float min = y.min();   // minimal value of y
  float max = x.max();   // maximal value of x

Using the [] operator for reading and writing single values of a local array is very convenient
but not recommended for a large number of values as it includes the overhead of getting the access
to an array.

.. code-block:: c++

   // more efficient                            // for ( .. ) x[i] = i  implies

   WriteAccess<float> writeX( x );              for ( IndexType i = 0; i < N; ++i )
   for ( IndexType i = 0; i < N; ++ i )         {
   {                                                WriteAccess<float> writeX( x );
       writeX[i] = i;                               writeX[i] = i;
   }                                            }


LAMAKernel
----------

This examples shows how to get a ContextPtr from the Context-factory and use it to
get a function from the KernelRegistry:

.. code-block:: c++

  // Get ContextPtr
  hmemo::ContextPtr loc = hmemo::Context::getContextPtr( common::context::Host );
  
  // Get function from Registry
  static LAMAKernel<UtilKernelTrait::setVal<double> > setVal;

  // update the context to make sure that setVal is available
  setval.getSupportedContext( loc );

  // Call function
  setVal[loc]( b, size, double( 2 ), reduction::COPY );  

The above scheme is used in other libraries, e.g. in the LAMA core, where kernel routines are used.
Another example of the LAMAKernel is shown :ref:`here <utilskernel_lamakernel_example>`.

************
Dependencies
************

Internal dependencies:

- :ref:`SCAI Common - Basic Concepts <scaicommon:main-page_common>`
- :ref:`SCAI Logging - Logging Macros <scailogging:main-page_logging>`
- :ref:`SCAI Tracing - Tracing Macros <scaitracing:main-page_tracing>`
- :ref:`SCAI Tasking - Asynchronous Tasks <scaitasking:main-page_tasking>`
- :ref:`SCAI Hmemo - Heterogeneous Memory Architecture <scaihmemo:main-page_hmemo>` (currently just used in testing)
- :ref:`SCAI Kregistry - Generic Kernel Registry <scaikregistry:main-page_kregistry>`
- :ref:`SCAI BLASKernel - BLAS Kernel <scaiblaskernel:main-page_blaskernel>`
  
External dependencies:

- :ref:`Thrust as part of the CUDA SDK toolkit <scaicommon:CUDA>`

************
Related Work
************
   
