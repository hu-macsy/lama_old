.. _main-page_kregistry:

##############
SCAI kregistry
##############

*************
Specification
*************

 * Handles kernel registration
 * Supports different backends
 * Provides the ContextFunction, KernelContextFunction
 * Internal dependencies: common, logging
 * External dependencies: none (C++ template library)
 * Different contextes are provided as enum from ‚common‘, the Context class (of hmemo) is not needed
 * Uses own exception class
 * Tests and examples are available.

********
Contents
********

.. toctree::
   :titlesonly:
   :maxdepth: 2
   
   KernelRegistry
   ContextFunction
   KernelFunction
   KernelContextFunction
   KernelTrait
   Example

**********
Motivation
**********

* Writing ‚kernel‘ functions on LAMA arrays that run on different devices (context)
* Kernel functions have usually basic functionality (no side effects) and might be called at different places 
   (e.g. initialization of an array, scaling of an array)
* Decision about context (location of execution) is separated from kernel implementations
* Read- Write- Accesses on LAMA arrays on the chosen location take care that data will be available
* One implementation on CPU host should always be available
* Further implementations can be added dynamically to improve performance (faster device) or to reduce memory transfer

.. code-block:: c++

  // method computes: a = b + c

  template <typename T>
  void add ( HArray<T>& a, const HArray<T>& b, const HArray<T>& c )
  {
       int n = <get size of a, checkt for sufficient sizes of b and c>
       ContextPtr ctx =  choose best context for execution
       WriteAccess<T> wa (a, ctx )
       ReadAccess<T> rb( b, ctx )
       ReadAccess<T> rc (c, ctx )
  
       void ( kernel_add* ) ( T* a, const T* b, const T *c, int n ) = implementation for ctx
  
       // get routine gets pointer to the data for the selected context

       kernel_add( wa.get(), rb.get(), rc.get(), n )
  }

The library ''kregistry'' manages different kernel implementations dependent on the 
context where the routine is executed.

 * Use of switch statement makes code not extendable
 * Virtual functions for an interface class provided for each context would be great: but not possible for template functions
 * Class KernelRegistry: provides function pointers for different implementations 
 * Registry allows for adding new implementations dynamically
 * Wanted: fallback if routine is not implemented on specified context
