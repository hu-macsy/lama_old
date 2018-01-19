.. _Assertion:

Assertion
=========

It is very common to throw an exception if a certain condition is given.

.. code-block:: c++

  #include <scai/common/macros/throw.hpp>

  if ( val < 0 ) 
  {
      COMMON_THROWEXCEPTION( "sub: val must not be negative, val = " << val )
  }

Instead of this construct the ``ASSERT`` macro can be used in the following way
(be aware that the condition has negated):

.. code-block:: c++

  #include <scai/common/macros/assert.hpp>

  SCAI_ASSERT_GE( val, 0, "sub: val must not be negative, val = " << val )

Beside the more convenient usage, one other advantage of the assert macro is
that it can be compiled out of the code by using the flag:

.. code-block:: c++

  SCAI_ASSERT_ERROR( val >= 0, "sub: val must not be negative, val = " << val )
  SCAI_ASSERT_DEBUG( val >= 0, "sub: val must not be negative, val = " << val )

The following flags can be used to enable or disable ASSERT macros:

 * -DSCAI_ASSERT_LEVEL_OFF, all SCAI_ASSERT_ERROR and SCAI_ASSERT_DEBUG are ignored
 * -DSCAI_ASSERT_LEVEL_DEBUG, all macros are handled like SCAI_ASSERT
 * -DSCAI_ASSERT_LEVEL_ERROR, macros SCAI_ASSERT_ERROR are enabled, SCAI_ASSSERT_DEBUG are disabled

.. code-block:: c++

    SCAI_ASSERT( n1 == n2, "n1 and n2 have different values" )
    SCAI_ASSERT_EQ( n1, n2, "size mismatch, both vectors must have same length" )

The SCAI_ASSERT-macro relies on getting a boolean value which will be evaluated. You can
enter every possible term which results in a boolean value. Besides this some additional
macros are provided for comparison like equal, unequal, less equal, less than, greater 
equal or greater than. These macros and their writing with the normal SCAI_ASSERT macro
are listed below: 

.. code-block:: c++

    SCAI_ASSERT( n1 == n2, .. )     SCAI_ASSERT_EQ( n1, n2, .. )
    SCAI_ASSERT( n1 != n2, .. )     SCAI_ASSERT_NE( n1, n2, .. )
    SCAI_ASSERT( n1 <= n2, .. )     SCAI_ASSERT_LE( n1, n2, .. )
    SCAI_ASSERT( n1 < n2, .. )      SCAI_ASSERT_LT( n1, n2, .. )
    SCAI_ASSERT( n1 >= n2, .. )     SCAI_ASSERT_GE( n1, n2, .. )
    SCAI_ASSERT( n1 > n2, .. )      SCAI_ASSERT_GT( n1, n2, .. )

Every of the SCAI_ASSERT-macros is provided in the three assertion levels. 

.. code-block:: c++

    SCAI_ASSERT_<op>( val1, val2, ... )
    SCAI_ASSERT_<op>_DEBUG( val1, val2, ... )
    SCAI_ASSERT_<op>_ERROR( val1, val2, ... )

