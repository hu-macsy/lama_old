Class LamaKernel
================

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