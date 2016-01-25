.. _main-page:

SCAI kregistry
==============

Motivation
----------

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

The library ''kregistry'' managages different kernel implementations dependent on the 
context where the routine is executed.

 * Use of switch statement makes code not extendable
 * Virtual functions for an interface class provided for each context would be great: but not possible for template functions
 * Class KernelRegistry: provides function pointers for different implementations 
 * Registry allows for adding new implementations dynamically
 * Wanted: fallback if routine is not implemented on specified context

Kernel Functions
----------------

A kernel function is a function that is executed on a certain device and operates on data that has been allocated on it.

.. code-block:: c++

   void kernel_add( double* a, const double* b, const double* c, int n )

   template<typename T>
   void kernel_add( T* a, const T* b, const T* c, int n );

The routine might be executed either on the CPU host or on a CUDA device or on any other context.
Each pointer argument points to data that has been allocated on the same device where the function 
is called.

A kernel function has the following properties:

 * it takes full advantage of available parallelism on the device where they are implemented (OpenMP on CPU, CUDA on GPU device, OpenMP on Intel MIC)
 * it can assume that all data is available
 * is not synonymous with CUDA kernel routines as here the kernel routine is responsible for launching the corresponding kernel
 * it can be implemented by using devie specific libraries like BLAS, MKL on the CPU, Thrust library or cuSPARSE on the GPU
 * it should never use communication in kernel routines (possible but makes encapsulation and porting difficult, e.g. SAMG setup, ScaL, ...)
 * it might be executed asynchronously

Usually, a kernel function should be available for all supported context types.

Functionality of Kernel Registry
--------------------------------

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

Class ContextFunction
---------------------

.. code-block:: c++

  template<typename FunctionType>
  class ContextFunction : public _ContextFunction
  {
  public:
      ContextFunction(); /** Only default Constructor at this time */
      ContextFunction( const ContextFunction& other );
      FunctionType get( ContextType ctx ) const; 
      void set( ContextType ctx, FunctionType fn );
      using _ContextFunction::validContext;
      ContextType validContext( ContextType preferedCtx );
      ContextType validContext( const _ContextFunction& other, ContextType preferedCtx );
  };

The class ''ContextFunction'' has the following properties:

 * an object of the class ''ContextFunction'' is just an array of typed function pointers
 * the array of typed function pointers is indexed by the context type (enum value)
 * base class _ContextFunction is used as an array of untyped function pointers (for internal purpose only)
 * validContext returns preferred context if function pointer has been set, 
   otherwise the first context type for which it has been registered (usually Host).

Kernel Trait
------------

A kernel trait is a struct that has an entry for the name of the kernel routine and a type 
definition for its signature.

.. code-block:: c++

   struct function_name
   {
       typedef res_type ( *FuncType ) ( type1, ..., typen );
       static const char* getId() { return "function_name"; }
   };

Template arguments can be used to define traits for different value types.

.. code-block:: c++

   template <typename T1, typename T2, typename T3>
   struct template_function_name
   {
       typedef res_type ( *FuncType ) ( T1*, const T2*, T3 );
       static const char* getId() { return „template_function_name"; }
   };

A kernel trait is a struct that is introduced for the following reasons:

 * it avoids misspelling of function namesr; different strings used for registration and access
   will not be detected at compile time, but different struct identifiers.
 * The same is true for the correct signature, i.e. functions are never registered and
   used with different signatures.
 * They are used for Doxygen documentation of the function behavior

Multiple kernel traits might be grouped for kernel routines used in a certain module. The 
name of the group should appear as prefix in the id of the kernel routines.

.. code-block:: c++

  struct UtilKernelTrait
  {
      template <typename ValueType>
      struct maxval
      {   /** @brief Find maximal value of n contiguously stored values.
           *  @param[in] array is an array of values
           *  @param[in] n is the size of array
           *  @return maximum of all values in array                      
          */
          typedef ValueType ( *FuncType ) ( const ValueType array[], const IndexType n );

          static const char* getId() { return "Utils.maxval"; }
      };
  
      // Other traits of Utils routines
  
      ...

  }; // struct UtilKernelTrait


Class KernelRegistry
--------------------

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

Class KernelContextFunction
---------------------------

The class KernelContextFunction is derived from the class ContextFunction and provides the following
additional functionality:

 * The constructor initializes the array of function pointers with the corresponding entries of
   the kernel registry, i.e. all registered functions are available.
 * The operator [] is introduced to make its use more convenient.

.. code-block:: c++

   template<typename FunctionType> 
   class KernelContextFunction : public ContextFunction<FunctionType>
   {
   public:
        KernelContextFunction( const char* name );
        FunctionType operator[] ( ContextType ctx )
   };

The class KernelTraitContextFunction uses for its constructor the corresponding kernel trait.

.. code-block:: c++

   template<typename KernelTrait> 
   class KernelTraitContextFunction : public KernelContextFunction<typename KernelTrait::FuncType>
   {
   public:
       typedef typename KernelTrait::FuncType ContextFunctionType;
       KernelTraitContextFunction() : 
       KernelContextFunction<ContextFunctionType>( KernelTrait::getId() )
       using KernelContextFunction<ContextFunctionType>::operator[];
   };

Class LamaKernel
----------------

The class LAMAKernel is a further extension with the only difference that it uses a Context object for
the methods and for indexing and not the context type. This makes it use more convenient.
It is not introduced in the library kregistry.

.. code-block:: c++

  template<typename KernelTrait>
  class LAMAKernel : public kregistry::KernelTraitContextFunction<KernelTrait>
  {
  public:
      typedef typename KernelTrait::FuncType ContextFunctionType;
      LAMAKernel();
      ContextFunctionType operator[] ( hmemo::ContextPtr context );
      hmemo::ContextPtr getValidContext( hmemo::ContextPtr defaultContext );
      hmemo::ContextPtr getValidContext( kregistry::_ContextFunction other, 
                                         hmemo::ContextPtr defaultContext )
  }; 

Example
-------

The following example shows a typical use of this class LAMAKernel:

.. code-block:: c++

  template<typename ValueType>
  Scalar DenseVector<ValueType>::max() const
  {
    static LAMAKernel<UtilKernelTrait::maxval<ValueType> > maxval;
    ContextPtr loc = maxval.getValidContext( mLocalValues.getValidContext() );
    ReadAccess<ValueType> localValues( mLocalValues, loc );
    ValueType localMax = maxval[loc]( localValues.get(), localValues.size() );
    return getDistribution().getCommunicator().max( localMax );
  }

Some remarks:

 * Preferred context for execution is where valid values of the array with the local part of the vector is available (mLocalValues)
 * Context will be set to host if no kernel function has been registered for the prefered context
 * After the decision about the context, the corresponding read and write accesses for the heterogeneoues 
 * maxval[ contexType ] throws Exception if function pointer is NULL, i.e. no function has been registered for the given context type.

Here is an OpenMP implementation of the kernel routine maxval:

.. code-block:: c++

  template<typename ValueType>
  ValueType OpenMPUtils::maxval( const ValueType array[], const IndexType n )
  {
      ValueType val = Trait<ValueType>.minval;
      #pragma omp parallel
      {
          ValueType threadVal = 0;
          #pragma omp for schedule( SCAI_OMP_SCHEDULE )
          for( IndexType i = 0; i < n; ++i )
          {
              ValueType elem = abs( array[i] );
              if( elem > threadVal ) threadVal = elem;
          }
          #pragma omp critical
          {
              if( threadVal > val ) val = threadVal;
          }
      }
      return val;
  }

This method is registered as follows:

.. code-block:: c++

  using common::context::Host;
  KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD;
  KernelRegistry::set<UtilsInterface::maxval<int> >( maxval, Host, flag );
  KernelRegistry::set<UtilsInterface::maxval<float> >( maxval, Host, flag );
  KernelRegistry::set<UtilsInterface::maxval<double> >( maxval, Host, flag );
  ...

The registration of the kernel routine ''maxval'' also implies instantiation for the corresponding
value type. LAMA uses always traits for registration and the class LAMAKernel ith trait as template parameter for access to kernel routines.

Library kregistry
-----------------

 * Internal dependencies: common, logging
 * External dependencies: none (C++ template library)
 * Different contextes are provided as enum from ‚common‘, the Context class (of hmemo) is not needed
 * Uses own exception class
 * Tests and examples are available.

LAMA itself uses later different groups for kernel traits, but these are not part of this library here.

===================   =======================================================
Name                  Description
===================   =======================================================
UtilKernelTrait       general routines
BLASKernelTrait       routines with counterparts in BLAS,LAPACK
CSRKernelTrait        routines for CSR matrix storage
ELLKernelTrait        routines for ELL matrix storage
DenseKernelTrait      routines for Dense matrix storage
JDSKernelTrait        routines for JDS matrix storage
COOKernelTrait        routines for COO matrix storage
DIAKernelTrait        routines for DIA matrix storage
===================   =======================================================

