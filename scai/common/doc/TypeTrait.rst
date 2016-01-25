TypeTraits
----------

In LAMA, nearly all routines with mathematical operations use a template
argument for the arithmetic type of the involved variables. The TypeTraits
is a struct that encapsulates all type-specific routines and values
so that the code can be written without specific knowledge for which types
the code will be instantiated later.

.. code-block:: c++

  template <typename ValueType>
  struct TypeTraits
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

For each supported arithmetic type a template-specific implementation must be provided.

.. code-block:: c++

  template<>
  class TypeTraits<double>
  {
  public:
      static inline double sqrt( double x ) 
      {
          return ::sqrt( x );
      }
      static inline double abs( double x ) 
      {
          return ::fabs( x );
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
          return scalar2str( stype );  // use same string as ScalarType
      }
  }

Unfortunately, in C++ it is not the case that mathematical operations are overloaded
for all arithmetic types. E.g. there are different routines like sqrtf, sqrt and sqrtl
for the square root function.
The following examples show typical calls of the routines sqrt or abs. 

.. code-block:: c++

  template<typename ValueType>
  ValueType l2norm( const ValueType vals[], const int n )
  {
      ValueType sum = 0;
      ...
      return TypeTraits<ValueType>sqrt( sum );
  }

  template<typename ValueType>
  ValueType l1norm( const ValueType vals[], const int n )
  {
      ValueType maxval = 0;
      for ( int i = 0; i < n; ++i )
      {
          ValueType v = TypeTraits<ValueType>::abs( vals[i] );
          if ( v > maxval ) maxval = v;
      }
      return maxval;
  }

