.. _Factory:

Factory 
=======

LAMA makes like many other object-oriented software heavy use of inheritance. Typically for each
base class there are many derived classes available.

Wanted: applications that can deal with all kind of derived classes even with dynamically added ones (by library module)

.. code-block:: c++

   ContextPtr cudaContext = Context::create( "CUDA" );
   shared_ptr<_Harray> myArray ( _HArray::create( ScalarType::FLOAT ) );

A Factory is a class that allows the dynamic creation of derived objects for a base class without
knowing at compile time which derived classes might be available.

.. code-block:: c++

   class Base : public Factory<InputType, OutputType>


.. literalinclude:: ../examples/DemoFactory.cpp
   :language: c++
   :lines: 52-53,62

.. literalinclude:: ../examples/DemoFactory.cpp
   :language: c++
   :lines: 181-209

This is the ways to defined derived classes that register in the factory for the base class.

.. literalinclude:: ../examples/DemoFactory.cpp
   :language: c++
   :lines: 106-122,130

A guard for the registration should be instantiated explicitly in the code.

.. code-block:: c++

    Base::Register<Derived1>::RegisterGuard Base::Register<Derived1>::register

For template classes this mechanism is as follows:

.. literalinclude:: ../examples/DemoFactory.cpp
   :language: c++
   :lines: 139-156,166

Guard for registration should be initiated explicitly.

.. literalinclude:: ../examples/DemoFactory.cpp
   :language: c++
   :lines: 170-173

Factory1 is similiar to Factory but allows one additional argument for the creation of objects.

.. code-block:: c++

  class Base : public Factory1<InputType, OtherType, OutputType>

