Assertions
==========

Boost.Test offers a lot of assertions to check a statement. Here is a list of mostly needed assertions:

- BOOST_CHECK( expression ) : Checks if expression is true  
- BOOST_CHECK_CLOSE( val1, val2, tolerance) : Checks if two values are equal with a percentual tolerance
- BOOST_CHECK_EQUAL( val1, val2 ) : Checks if two values are equal
- BOOST_CHECK_THROW( expression, ExceptionType) : Checks if a expression throws an exception
- BOOST_CHECK_MESSAGE( expression, msg ) : Checks if expression is true; if not: print out message 

You can find a summary of available assertions in Boost.test here_

.. _here: http://www.boost.org/doc/libs/1_45_0/libs/test/doc/html/utf/testing-tools/reference.html

In addition to assertions of Boost there are some selfmade assertions, which fill some missing functions of Boost. These Assertion are defined in TestMacros.hpp.

Selfmade assertions of LAMA:

- SCAI_CHECK_SCALAR_CLOSE( x, y, ValueType, tolerance ) : Checks if two Scalars are equal with a percentual tolerance. The Scalar are transformed to ValueType.
- SCAI_CHECK_SCALAR_SMALL( x, ValueType, eps ) : Checks if a Scalar is equal to 0 with a absolute tolerance
- SCAI_CHECK_SCALAR_SMALL_EPS( x, ValueType ) : Checks if a Scalar is equal to 0 with a default absolute tolerance given by LAMA. This default tolerance is 1E-5.


These two macros are used by nearly all tests for objects, which inherit by class printable.

- LAMA_WRITEAT_TEST( printable ) : Checks if an object of type printable prints something (the content is not compared). 
- LAMA_WRITEAT_TEST_PTR( printable ) : Checks if an object of type printable prints something (the content is not compared).


Testcase macros: 

- LAMA_AUTO_TEST_CASE_TDUMMY( functionname, classname ) : This macro creates a functioncall of a templated testmethod with a dummy template
- LAMA_AUTO_TEST_CASE_T( functionname, classname ) : This macro creates a call of a single templated testmethod for each supported ValueType
- LAMA_AUTO_TEST_CASE_TT( functionname, classname ) : This macro creates a call of a double templated testmethod for each supported ValueType. This is helpful to test methods with casts of different valuetypes.

All theses three macros implement the CONTEXTLOOP() to call testmethods on different contexts automatically. The argument classname is used to call the correct testmethod of the given namespace of a testclass.
Those macro are especially used in Utiltest-classes, like ELLUtilsTest, JDSUtilsTest, etc. 
