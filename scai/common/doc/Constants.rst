.. _Constants:

Constants
=========

There are sometimes situations where a specific value of an argument
allow  more efficient implementations.
In the following example a lhs vector is replaced with a scaled
rhs vector. If the scale factor is 0, the multiplication can be replaced with just setting
lhs vector to 0, if it is 1, no scaling is required at all.

.. code-block:: c++

  template <ValueType>
  void mult( ValueType lhs[], const ValueType rhs[], const ValueType scale, int n )
  {
      if ( scale == ValueType( 0 ) )
      {
          // lhs = 0
      }
      else if ( scale == ValueType( 1 ) )
      {
          // lhs = rhs 
      }
      else
      {
          // lhs = scale * rhs
      }
  }
 
The enum class ``Constants`` provides an enum type for specific values:

.. code-block:: c++
	
	    enum class Constants
	    {
	        ONE,
	        ZERO
	    };

Furthermore, it has an overloaded comparison operator for equality:

.. code-block:: c++

   template<typename ValueType>
   bool operator==( const ValueType& x, const enum ConstantType& c );

The ``operator==`` does not test for exact equality but for equality in a certain 
range where the range is also specific for each value type as here the eps value
of the TypeTrait is used.

.. code-block:: c++

  template <ValueType>
  void mult( ValueType lhs[], const ValueType rhs[], const ValueType scale, int n )
  {
      if ( scale == Constants::ZERO )
      {
          // lhs = 0
      }
      else if ( scale == Constants::ONE )
      {
          // lhs = rhs 
      }
      else
      {
          // lhs = scale * rhs
      }
  }

Therefore the class Constants provides a standardized way to compare a variable
of an arbitrary value type with a constant value.
