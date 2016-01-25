Exception
---------

Within LAMA one base class is provided for all LAMA relevant exceptions.

.. code-block:: c++

  class Exception
  {
      ...
  }


.. code-block:: c++

   COMMON_THROWEXCEPTION( msg )


.. code-block:: c++

   SCAI_THROWEXCEPTION( exceptionClass, msg )

Assertions
----------

.. code-block:: c++

    void sub( int val )
    {
        if ( val < 0 ) 
        {
            COMMON_THROWEXCEPTION( "sub: val must not be negative, val = " << val )
        }
        COMMON_ASSERT ( val % 2 == 0, "val = " << val << " must be even" )
        COMMON_ASSERT_LT( val, 10, "val = " << val << " must be less than 10" )
        COMMON_ASSERT_EQUAL( val, 4, "None" )
    }
    
    try
    {
        sub( vals[i] );
        std::cout << "Call of sub( " << vals[i] << ") terminated correctly" << std::endl;
    }
    catch ( const std::exception& exception )
    {
        // Importation: exception is a reference, so we get the routine of common::Exception
     
        std::cout << "Got exception: " << exception.what() << std::endl;
    }

There is one further macro provided that verifies assertions only for the debug mode.

.. code-block:: c++

    COMMON_ASSERT_DEBUG( msg )
