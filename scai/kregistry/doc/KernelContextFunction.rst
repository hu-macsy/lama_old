Class KernelContextFunction
===========================

The class ``KernelContextFunction`` is derived from the class ``ContextFunction`` and provides the following
additional functionality:

- The constructor initializes the array of function pointers with the corresponding entries of
   the kernel registry, i.e. all registered functions are available.
   
- The operator [] is introduced to make its use more convenient.


.. code-block:: c++

   template<typename FunctionType> 
   class KernelContextFunction : public ContextFunction<FunctionType>
   {
   public:
        KernelContextFunction( const char* name );
        FunctionType operator[] ( ContextType ctx )
   };

The class ``KernelTraitContextFunction`` uses for its constructor the corresponding ``KernelTrait``.

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