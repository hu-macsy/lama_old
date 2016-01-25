ContextType
-----------

 * Enumeration for different devices / locations
 * Used as key for ContextFactory to get a certain device

.. code-block:: c++

  namespace context
  {
      enum ContextType
      {
          Host,          //!< context for cpu + main memory
          CUDA,          //!< CUDA GPU device
          OpenCL,        //!< OpenCL GPU device, currently not supported
          MIC,           //!< Intel MIC
          UserContext,   //!< can be used for a new derived Context class
          MaxContext     //!< used for dimension of ContextType arrays
      }; 
  
      std::ostream& operator<<( std::ostream& stream, const ContextType& type );

  } // namespace context


ScalarType
----------

 * ScalarType is an enumeration type that has a value at least for each supported arithmetic type
 * Some additional values are used for internal purposes
 * This enum type can be used for operations where type ids itself are not allowed
 * Conversion of type ids to enum values by TypeTraits

.. code-block:: c++

  namespace scalar
  {
      typedef enum
      {
          INDEX_TYPE,          //!<  synonymous for IndexType
          FLOAT,               //!<  synonymous for float
          DOUBLE,              //!<  synonymous for double
          LONG_DOUBLE,         //!<  synonymous for long double
          COMPLEX,             //!<  synonymous for complex
          DOUBLE_COMPLEX,      //!<  synonymous for double complex
          LONG_DOUBLE_COMPLEX, //!<  synonymous for long double complex
          PATTERN,             //!<  dummy type of size 0
          INTERNAL,            //!<  for internal use in I/O
          UNKNOWN
      } ScalarType;

      std::ostream& operator<<( std::ostream& stream, const ScalarType& type );
  }

Example of Use:

.. code-block:: c++

  class _HArray
  {
      ...
      virtual ScalarType getValueType() = 0;
  }

  template <typename ValueType>
  class HArray
  {
      ...
      virtual ScalarType getValueType()
      {
          return TypeTraits<ValueType>::scalarType;
      }
  }

  void sub( _Harray& array )
  {
      ScalarType type = array.getValueType();

      std::cout << "array has type " << type << std::endl;
      switch ( type ) :
          case scalar::FLOAT :
             ...
          case scalar::DOUBLE ::
             ...
      ...
  }
