.. _main-page_kregistry:

##############
SCAI KRegistry
##############

***********
Description
***********

The library KRegistry manages different kernel implementations dependent on the 
context where the routine is executed. The implementations have to register theirselfs
at the KernelRegistry. This can be done by using the Registrator, which is shipped within. 
To support new hardware just the kernels have to be added. Algorithm which are 
formulated using the KernelRegistry don't have to be altered for this. The KernelRegistry 
supports KernelTraits which are a set of kernels grouped by affinity.

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
   Registrator

*********
Relations
*********

The KernelRegistry is the counterpart to the HArray. The HArray is abstraction of different memory locations. While the 
KernelRegistry is the abstraction of the execution location. 

*******
Example
*******

The following example shows how a kernel ``add`` can be registered in the KernelRegistry using the ``ExampleKernelTrait``. 
The registration is done through static instantiation of a global boolean. The function will be registered for the valuetypes
int, float and double.  

.. code-block:: c++

  // method computes: a = b + c
  template<typename T>
  void add ( T* a, const T* b, const T* c, const IndexType n )
  {
  	for( IndexType i = 0; i < n; ++i )
  	{
  		a[i] = b[i] + c[i];
  	}
  }
  
  // KernelTrait
  struct ExampleKernelTrait
  {
    template<typename T>
  	struct add
  	{
  		typedef ( *FuncType )( T* a, const T* b, const T* c, const IndexType n );
  		
  		static const char* getId() { return "Example.add"; }
  	};
  };
  
  // Function to register ExampleKernelTrait in the KernelRegistry
  bool register_example()
  {
  	common::context::ContextType ctx = common::context::Host;
  	
  	KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD;
  
  	KernelRegistry::set<ExampleKernelTrait::add<int> >( add, ctx, flag ); 
  	KernelRegistry::set<ExampleKernelTrait::add<float> >( add, ctx, flag ); 
  	KernelRegistry::set<ExampleKernelTrait::add<double> >( add, ctx, flag ); 
  	
  	return true;
  }
  
  // static instantiation of the registration
  bool registered_example = register_example();
  

*********************
Environment-Variables
*********************

No additional environmental variables are defined in the project. 

************
Dependencies
************

- Internal:

  - :ref:`SCAI Common - Basic Concepts <scaicommon:main-page_common>`

  - :ref:`SCAI Logging - Logging Macros <scailogging:main-page_logging>`

************
Related Work
************

 