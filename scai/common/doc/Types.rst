ContextType
-----------

A context specifies a device like CPU or GPU accelerator where memory can be
allocated and where operations on it might be executed.
The enumeration type for context might be used for registration and 
searching of objects belonging to a certain context.

.. literalinclude:: ../ContextType.hpp 
   :language: c++
   :lines: 54-55,63-72,82,85

..	.. code-block:: c++
	
	  namespace context
	  {
	      typdef enum 
	      {
	          Host,          //!< context for cpu + main memory
	          CUDA,          //!< CUDA GPU device
	          OpenCL,        //!< OpenCL GPU device, currently not supported
	          MIC,           //!< Intel MIC
	          UserContext,   //!< can be used for a new derived Context class
	          MaxContext     //!< used for dimension of ContextType arrays
	      } ContextType; 
	  
	      std::ostream& operator<<( std::ostream& stream, const ContextType& type );
	
	  } // namespace context


ScalarType
----------

ScalarType is an enumeration type that has a value at least for each supported arithmetic type.
Some additional values are used for internal purposes. This enum type can be used for operations 
where type ids itself are not allowed or inefficient.

.. literalinclude:: ../ScalarType.hpp 
   :language: c++
   :lines: 54-55,70-82,92-94

..	.. code-block:: c++
	
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

The following example shows how this enumeration type is used in LAMA or might
be used within other applications. Its main purpose is to provide some
kind of logic to guarantee the safety of the ``static_cast`` operation instead of
using the expensive ``dynamic_cast``.

The example consists of common abstract base class ``_Data`` and 
a derived template class ``Data`` that might be instantiated for different arithmetic types.

.. code-block:: c++

  class _Data
  {
      ...
      virtual ScalarType getValueType() = 0;
  }

  template <typename ValueType>
  class Data : public _Data
  {
      ...
      virtual ScalarType getValueType()
      {
          return TypeTraits<ValueType>::scalarType;
      }
  }

Often, for an arbitrary object of the base class it is important to find
out from which derived class it is. As the virtual method ``getValueType``
returns a corresponding enum value, this value identifies exactly to which 
template class the object belongs.

.. code-block:: c++

  void sub( _Data& array )
  {
      ScalarType type = array.getValueType();

      std::cout << "array has type " << type << std::endl;
      switch ( type ) :
          case scalar::FLOAT :
             Data<float>& f_array = static_cast< Data<float>& >( array );
             ...
          case scalar::DOUBLE ::
             Data<double>& d_array = static_cast< Data<double>& >( array );
             ...
      ...
  }

The alternative solution is using a ``dynamic_cast`` operator and to check for each derived
class whether the corresponding cast was successful or not. This is rather expensive and
using the enumeration type provides the logic to guarantee the safety of the
static_cast operation.

ReductionOp
-----------

The following enumeration type specifies the different kind of binary operators
that can be used in reduction operators.

.. literalinclude:: ../ReductionOp.hpp 
   :language: c++
   :lines: 49-50,62-71,76,105

..	.. code-block:: c++
	
	  namespace reduction
	  {
	      typedef enum
	      {
	          COPY,     // for assign   x = y
	          ADD,      // for operator x += y
	          MULT,     // for operator x *= y
	          MIN,      // for operator x = min( x, y )
	          MAX,      // for operator x = max( x, y )
	          ABS_MAX   // for operator x = max( x, abs(y) )
	  
	      } ReductionOp;
	
	      std::ostream& operator<<( std::ostream& stream, const ReductionOp& op );
	  }

The enum class is used in ordeer to have one common function with an addtional op argument instead
of individual functions for each kind of operator.

