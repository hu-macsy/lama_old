.. _TypeTrait:

TypeTrait
=========

When writing template functions that might be instantiated for different value
types, it is sometimes the case that this code contains some specific part that
really depends on the value type.

The TypeTrait provided by the common library
is a struct that encapsulates all type-specific routines and values
so that the code can be written without specific knowledge for which types
the code will be instantiated later. It is very similar to the ``numeric-limits`` class
that provides in C++ a standardized way to query various properties of arithmetic types.

.. literalinclude:: ../TypeTraits.hpp 
   :language: c++
   :lines: 52-103

The member variable stype of the TypeTrait is set with the enum value of 
ScalarType that stands for the corresponding type. It can be used to implement
a query operation for the value type.

.. code-block:: c++

  template <typename ValueType>
  class Data : public _Data
  {
      ...
      virtual ScalarType getValueType()
      {
          return TypeTrait<ValueType>::scalarType;
      }
  }

The id routine is helpful to print the name of the value in logging messages. It might
be more convenient than using the name method of typeinfo.

.. code-block:: c++

  template <typename ValueType>
  void f( ... ) 
  {
      // std::cout << "f<" typeid( ValueType ).name() << ">" << std::endl;

      std::cout << "f<" << TypeTrait<ValueType>::id() << ">" << std::endl;
      ....
  }

For each supported arithmetic type a template-specific implementation is provided.

