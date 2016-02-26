KernelRegistry
==============

Functionality
-------------

 * Kernel registry is large table of function pointers
 * Keys are: 
   - Name of the routine
   - signature of the routine (very often template routines with parameters for the different value types)
   - context where the function is executed
 * Methods working on LAMA arrays use kernel registry to find the correct implementations
 * Implementation and registration of kernel functions are completely independent.
 * KernelRegistry not implemented as (common::) Factory
 * one derived class for each kernel routine is too much overhead
 * Grouping of routines for context not really supported

Data structures for kernel registry:

 * Kernel registry is implemented by using std::map

   * extensible, i.e. new kernel routines can be added without changing data structures
   * rather fast

 * Function pointers for the different contextes are put together in one array

   * Some memory overhead as many pointers might be NULL
   * Static array implies maximal number of supported contexts
   * But: all implementations of one kernel routine on all contexts are available with one access in routines using the corresponding function
   * Use of static variable makes implementation very efficient as routine is only looked up once
   
Implementation
--------------

.. code-block:: c++

   class KernelRegistry
   {
       template<typename FunctionType>
       static void set( FunctionType fn, const char* name, ContextType ctx, 
                        bool replace = false );
   
       template<typename FunctionType>
       static void get( FunctionType& fn, const char* name, ContextType ctx );
   
       template<typename KernelTrait>
       static void set( typename KernelTrait::FuncType fn, ContextType ctx, 
                        bool replace = false );
   
       template<typename KernelTrait>
       static void get( ContextFunction<KernelTrait::FunctionType>& contextFunction );
   };

Some remarks:

 * Set routine using a kernel trait doesn not use the string argument as it is provided by the trait
 * kernel trait provides also the function type (signature)
 * the class is implemented as a singleton, i.e. there is only one incarnation
 * i.e. all kernel routines are registered in one map, as multiple maps make handling more complicated
 * For unique keys it is useful to prefix names of kernel routines, e.g. ELL.getCSRValues and JDS.getCSRValues.
   