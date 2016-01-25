TypeTrait
---------

When writing template functions that might be instantiated for different value
types, it is sometimes the case that this code contains some specific part that
really depends on the value type.

The following template function returns the L2-norm for an array that contains
values of a certain type.

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

Usually, operator and function overloading guarantees that the code can be instantiated
for different ValueTypes like float, double, complex. 
Unfortunately, in C++ it is not the case that mathematical operations are overloaded
for all arithmetic types. E.g. there are different routines like sqrtf, sqrt and sqrtl
for the square root function.

For this and some other situations, a struct TypeTrait is used that provides 
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
      return TypeTrait<ValueType>::sqrt( sum );
  }

The TypeTrait provided by the common library
is a struct that encapsulates all type-specific routines and values
so that the code can be written without specific knowledge for which types
the code will be instantiated later. It is very similar to the ``numeric-limits`` class
that provides in C++ a standardized way to query various properties of arithmetic types.

.. code-block:: c++

  template <typename ValueType>
  struct TypeTrait
  {
      // square root function
      static inline ValueType sqrt( const ValueType x );
  
      // absolute value function
      static inline ValueType abs( const ValueType y );

      // conversion ValueType -> enum ScalarType
      static scalar::ScalarType stype = scalar::UNKNOWN; 
  
      // convenient name of type for output
      static inline const char* id();
  
      // Get value-specific epsilon for comparison.
      static inline ValueType getEps();
  
      // Get maximal value of a ValueType, used for reductions on arrays. */
      static inline ValueType getMax();
  };

The member variable scalar of the TypeTrait is set with the enum value of 
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
The following code shows the implementation for the ValueType double.

.. code-block:: c++

  template<>
  class TypeTrait<double>
  {
  public:
      static inline double sqrt( double x ) 
      {
          return ::sqrt( x );  // from cmath
      }
      static inline double abs( double x ) 
      {
          return ::fabs( x );  // from cmath
      }
      static inline double getEps()
      {
          return std::numeric_limits<double>::epsilon();
      }
      static inline double getMax() 
      {
          return std::numeric_limits<double>::max();
      }
      static scalar::ScalarType stype = scalar::DOUBLE;
      static inline const char* id()
      {
          return "double";
      }
  }

