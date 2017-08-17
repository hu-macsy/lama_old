Class LAMAKernel
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
  
.. _utilskernel_lamakernel_example:  
  
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

- Preferred context for execution is where valid values of the array with the local part of the vector is available (mLocalValues)

- ``Context`` will be set to host if no kernel function has been registered for the prefered context

- After the decision about the context, the corresponding read and write accesses for the heterogeneoues 

- ``maxval[ contexType ]`` throws ``Exception`` if function pointer is NULL, i.e. no function has been registered for the given context type.

Here is an OpenMP implementation of the kernel routine maxval:

.. code-block:: c++

  template<typename ValueType>
  ValueType OpenMPUtils::maxval( const ValueType array[], const IndexType n )
  {
      ValueType val = Trait<ValueType>.minval;
      #pragma omp parallel
      {
          ValueType threadVal = 0;
          #pragma omp for 
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

The registration of the kernel routine ``maxval`` also implies instantiation for the corresponding value type. 
