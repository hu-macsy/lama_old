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

Math
----

Math is a wrapper for functions from cmath. Due to the compatibility of cmath to C it is not the case that mathematical operations are overloaded
for all arithmetic types. E.g. there are different routines like sqrtf, sqrt and sqrtl
for the square root function. The Math wrapper enables you
just to write Math::sqrt. This is very helpful to write templated code. Besides sqrt Math
provides wrapping functions for abs and conj. Currently int, long, long, float, double, long double,
Complex<float>, Complex<double> and Complex<long double> are supported.  

The following template function returns the L2-norm for an array that contains
values of a certain type. If you call this code with ValueType=float it will lead to compiler warning, 
because the sqrt-function expects an argument of type double. 

.. code-block:: c++

  template<typename ValueType>
  ValueType l2norm( const ValueType vals[], const int n )
  {
      ValueType sum = 0;
      for ( int i = 0; i < n; ++i )
      {
          sum += vals[i] * vals[i];
      }
      return sqrt( sum );
  }


For this and some other situations, a struct Math is used that provides 
type-specific functions or values that can be used. In the example above, the 
code becomes:

.. code-block:: c++

  template<typename ValueType>
  ValueType l2norm( const ValueType vals[], const int n )
  {
      ValueType sum = 0;
      for ( int i = 0; i < n; ++i )
      {
          sum += vals[i] * vals[i];
      }
      return Math::sqrt( sum );
  }
