.. _ScalarType:

ScalarType
==========

ScalarType is an enumeration type that has a value at least for each supported arithmetic type.
Some additional values are used for internal purposes. This enum type can be used for operations 
where type ids itself are not allowed or inefficient.

=====================    =============================
Value                    Meaning
=====================    =============================
INDEX_TYPE               synonymous for IndexType
FLOAT                    synonymous for float
DOUBLE                   synonymous for double
LONG_DOUBLE              synonymous for long double
COMPLEX                  synonymous for complex
DOUBLE_COMPLEX           synonymous for double complex
LONG_DOUBLE_COMPLEX      synonymous for long double complex
PATTERN                  dummy type of size 0
INTERNAL                 for internal use in I/O
UNKNOWN                  internal use for the number of enum values
=====================    =============================

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
          return TypeTraits<ValueType>::stype;
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
          case ScalarType::FLOAT :
             Data<float>& f_array = static_cast< Data<float>& >( array );
             ...
          case ScalarType::DOUBLE ::
             Data<double>& d_array = static_cast< Data<double>& >( array );
             ...
      ...
  }

The alternative solution is using a ``dynamic_cast`` operator and to check for each derived
class whether the corresponding cast was successful or not. This is rather expensive and
using the enumeration type provides the logic to guarantee the safety of the
static_cast operation.


