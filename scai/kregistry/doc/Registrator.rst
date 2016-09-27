Registrator
===========

.. image:: _images/Registrator.png
    :align: center
    :width: 700px

The ``Registrator`` is used for the registration of functions that have one or two template
arguments. It uses meta-programming techniques to register these functions for a certain number
of types.

The template class ``mepr::RegistratorV`` is used to call the method with one template argument and 
``mepor::RegistratorVO`` for the method with two template arguments. This meta-programming
 is based on template programming that where only classes or structs can be used as arguments, but not functions. This is the 
reason why the ``registerKernels`` functions are encapsulated in a struct.

.. code-block:: c++

	kregistry::mepr::RegistratorV<ArrayKernels, SCAI_ARRAY_TYPES_HOST_LIST>::registerKernels( flag );
	kregistry::mepr::RegistratorVO<BinOpKernels, SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::registerKernels( flag );

The class ``RegistratorV`` calls the method ``registerKernels`` of the template class ``ArrayKernels`` with one template argument.

.. code-block:: c++

    template<typename ValueType>
    struct ArrayKernels
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag )
        {
            using kregistry::KernelRegistry;
            const common::context::ContextType ctx = common::context::Host;
            KernelRegistry::set<UtilKernelTrait::vectorScale<ValueType> >( vectorScale, ctx, flag );
            KernelRegistry::set<UtilKernelTrait::reduce<ValueType> >( reduce, ctx, flag );
            ....
        }
    };

The class ``RegistratorVO`` calls the method ``registerKernels`` of the template class ``BinOpKernels`` that has two template arguments.
It is called for for each possible pair of types of is two type list arguments.

.. code-block:: c++

    template<typename ValueType, OtherValueType>
    struct BinOpKernels
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag )
        {
            using kregistry::KernelRegistry;
            const common::context::ContextType ctx = common::context::Host;
            KernelRegistry::set<UtilKernelTrait::setScale<ValueType, OtherValueType> >( setScale, ctx, flag );
            KernelRegistry::set<UtilKernelTrait::setGather<ValueType, OtherValueType> >( setGather, ctx, flag );
            ....
        }
    };

The registration of kernel functions without template argument does not require any kind of meta-prorgramming.
For consistency with the other notation, the corresponding ``registerKernels`` function is also wrapped with
a corresponding structure.

.. code-block:: c++

    struct IndexKernels
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag )
        {
            using kregistry::KernelRegistry;
            const common::context::ContextType ctx = common::context::Host;
            KernelRegistry::set<UtilKernelTrait::validIndexes>( validIndexes, ctx, flag );
            KernelRegistry::set<UtilKernelTrait::countBuckets<IndexType> >( countBuckets, ctx, flag );
            KernelRegistry::set<UtilKernelTrait::sortInBuckets<IndexType> >( sortInBuckets, ctx, flag );
            ....
        }
    };

The following code shows the typical outline of a KernelClass that provides different kernel routines 
and registers them. A kernel class is always a singleton class. The single object registers all kernel routines 
during the static initialization when it is constructed.

.. code-block:: c++

    class KernelClass
    {
        // Here we have all methods implementing kernels

        ...

        struct Kernels
        {
            // register all kernels without template argument
            static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
        }

        template<typename ValueType>
        struct KernelsV
        {
            // register all kernels with one template argument
            static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
        }

        template<typename ValueType, typename OtherValueType>
        struct KernelsVO
        {
            // register all kernels with two template arguments
            static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
        }

        KernelClass::KernelClass()
        {
            Kernels::registerKernels( flag );
	        kregistry::mepr::RegistratorV<KernelsV, SCAI_ARRAY_TYPES_HOST_LIST>::registerKernels( flag );
	        kregistry::mepr::RegistratorVO<KernelsVO, SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::registerKernels( flag );
        }
    }

Here are some further remarks:

* Kernel functions are never bound to an object, so they are always static and the same is 
  true for the methods ``registerKernels``.

* In some situations a kernel routine with one template argument is only instantiated for
  one type, e.g. for the ``IndexType``. In this case no meta-programming is required and
  this routine might be register in the ``registerKernels`` method without template
  argument.

* The above approach allows also to register kernel routines for different type lists. E.g.
  one kernel routine might be registered for all array types, e.g. another kernel routine
  only for numeric types and not for the IndexType. 

.. code-block:: c++

   kregistry::mepr::RegistratorV<ArithmeticKernels, SCAI_ARITHMETIC_TYPES_HOST_LIST>::registerKernels( flag );
   kregistry::mepr::RegistratorV<ArrayKernels, SCAI_ARRAY_TYPES_HOST_LIST>::registerKernels( flag );

* The destructor can be used in the same way as the constructor to unregister the kernel methods.
  This is the reason why the mehtods ``registerKernels`` have one ``KernelRegistryFlag`` argument.

