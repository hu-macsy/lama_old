.. _main-page_utilskernel:

################
SCAI UtilsKernel 
################

***********
Description
***********

The UtilsKernel library provides various kernel implementations for operations on heterogeneous arrays.
The implementations are registered in the KernelRegistry. 

********
Contents
********

.. toctree::
   :titlesonly:
   :maxdepth: 1
   
   HArrayUtils
   Functions
   LAMAKernel

*******
Example
*******

HArray
------

The following example shows how you can operate on heterogeneous arrays:

.. code-block:: c++

  const int N = 10;

  HArray<float> x( N );
  HArray<float> y( N );
  
  for ( IndexType i = 0; i < N; ++i )
  {
      x[i] = i;
  }

  HArrayUtils::assign( y, x );                           // y = x;
  HArrayUtils::setArray( y, y, common::BinaryOp::MULT ); // y *= y;;
  HArrayUtils::setScalar( y, 2, common::BinaryOp::MULT ); // y *= 2;;
   
  float l1 = HArrayUtils::l1Norm( x ); 
  float l2 = HArrayUtils::l2Norm( x ); 

  float sum = HArrayUtils::sum( x );   // build the sum of all elements
  float min = HArrayUtils::min( y );   // minimal value of y
  float max = HArrayUtils::max( x );   // maximal value of x

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
   
