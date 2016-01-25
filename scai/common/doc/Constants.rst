Constants
---------

.. code-block:: c++

    enum ConstantType
    {
        ONE,
        ZERO
    };


.. code-block:: c++

   template<typename ValueType>
   bool operator==( const ValueType& x, const enum ConstantType& c )
   {
       return    TypeTraits<ValueType>::abs( x - getConstant<ValueType>( c ) )
               < TypeTraits<ValueType>::eps ;
   }


.. code-block:: c++

   ValueType v = ...

   if ( v == Constants::ZERO )
   {
       ...
   }
