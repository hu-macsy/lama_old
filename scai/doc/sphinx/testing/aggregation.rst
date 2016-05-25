:orphan:

Aggregation Unit Testing
========================

**WARNING: deprecated**

Building and Executing Unit Tests
---------------------------------

Assuming you are in your lama build directory you can build the unit tests with

.. code-block:: bash
	
	make

Or you can execute

.. code-block:: bash

	./lama_tests

to run all tests.

LAMA uses the testframwork Boost.Test from the Boost-Library. So all runtime parameter given by this
testframework are available.

Runtime parameter given by Boost.Test
-------------------------------------

You can run a single test by using :

.. code-block:: bash

	./lama_tests --run_test=Testsuitename/Testcasename

In LAMA all coherent tests were grouped in testsuites. All tests in a testsuites have unique names.
RegEx-Expressions are useful to run subgroups of tests with a common expression.

.. code-block:: bash

	./lama_tests --run_test=*Vector*/*Multiplication*

In this example all testcases with the expression "Multiplication" in all testsuites with the expression
"Vector" are executed.

There are some conventions with the names of testsuites:
All tests, which are made for running calculations parallel, begin with "P\_", so you can execute all
parallel tests by calling:

.. code-block:: bash

	./lama_tests --run_test=P_*/*


All tests, which are made for running calculations on CUDA, begin with "CUDA\_", so you can execute
these tests by calling:

.. code-block:: bash

	./lama_tests --run_test=CUDA_*/*

Receiving more information about the testrun, you can change the loglevel of Boost.Test. For example
like this:

.. code-block:: bash

	./lama_tests --log_level=test_suite

There are different possible levels: all, test_suite, message, warning and many more ...

You can have a look in a summary with all runtime parameters `here`__.

__ http://www.boost.org/doc/libs/1_45_0/libs/test/doc/html/utf/user-guide/runtime-config/reference.html

Creating new tests
------------------

Boost.Test offers different Makros to create testsuites and testcases. In LAMA each class has a related
testclass.

Basic structure of a test case
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The structure of a !TestClass.cpp looks like this:

.. code-block:: c++

  /*(1)*/ #include <boost/test/unit_test.hpp>

  /*(2)*/ typedef boost::mpl::list<double,float>  test_types;

          using namespace boost;
          using namespace lama;

          /* --------------------------------------------------------------------- */

  /*(3)*/ struct TestClassConfig
          {
              TestClassConfig()
              {
                  //SetUp
              }
  
              ~TestClassConfig() 
              {
                  //TearDown
              }
  
              //Declarations of objects
          };

  /*(4)*/ BOOST_FIXTURE_TEST_SUITE( TestClass, TestClassConfig );         //TestSuiteName is the filename "TestClass"
                                                                          //TestFixtureName is the filename+"Config"
  
          /* --------------------------------------------------------------------- */
  
  /*(5)*/ BOOST_AUTO_TEST_CASE( test1 )
          {
              //Test Case 1
          }
      
          BOOST_AUTO_TEST_CASE( test2 )
          {
              //Test Case 2
          }
  
  /*(6)*/ BOOST_AUTO_TEST_CASE_TEMPLATE( test3, T, test_types )
          {
              //Test Case 3 (template)
          }
  
          /* --------------------------------------------------------------------- */
  
  /*(7)*/ BOOST_AUTO_TEST_SUITE_END();

(1) :   Each testclass has to inlude this headerfile from Boost.Test.
(2) :   Defining a list of test-types. Each templated test runs as often as the number of given test_types.
(3) :   The possibility of !SetUp() and !TearDown-Methods can be realized by a struct.
(4) :   This struct can be committed in the makro BOOST_FIXTURE_TEST_SUITE(suitename, structname). If it is not neccessary to create those common used objects, an alternative is to use the makro BOOST_AUTO_TEST_SUITE( suitename ).
        There are some conventions about name-giving:
                (a) :   Each !TestSuiteName should be equal to the !FileName.
                (b) :   Each !TestFixtureName is created by adding "Config" to the !TestSuiteName and/or !FileName.
                        It seems that Boost.Test has some problems with equal !FixtureNames of different !TestClasses.
(5) :   Testcases are created by using the makro BOOST_AUTO_TEST_CASE( testcasename ).
(6) :   To parameterize testcases, you can use the makro BOOST_AUTO_TEST_CASE_TEMPLATE( casename, T, test_types). This case will run as often as the count of objects in the collection test_types. 
(7) :   The makro BOOST_AUTO_TEST_SUITE_END() will close the testsuite

Assertions
^^^^^^^^^^

You can find a summary of available assertions in Boost.test `here`__ .
In addition to assertions of Boost there are some selfmade assertions, which fill some missing functions
of Boost. These Assertion are defined in TestHelper.h 

__ http://www.boost.org/doc/libs/1_45_0/libs/test/doc/html/utf/testing-tools/reference.html

The assertion LAMA_BOOST_CHECK_CLOSE is made for comparing two Scalars. The epsilon is given in percentage units. This assertion transforms the scalars into value of type float or double and calls BOOST_CHECK_CLOSE from Boost.Test.

.. code-block:: c++

	LAMA_BOOST_CHECK_CLOSE( Scalar x, Scalar y, eps )

The assertion LAMA_BOOST_CHECK is made for comparing two Scalars. The epsilon is given as a floatingpoint
number. This brings an advantage by testing very small values. It takes the absolute value of the
difference of x and y and calls BOOST_CHECK from Boost.Test.

.. code-block:: c++

	LAMA_BOOST_CHECK( x, y, eps ) 

The assertion LAMA_BOOST_CHECK_EPS is made for comparing two Scalar. This is the same function as
LAMA_BOOST_CHECK, but in this case the epsilon is taken from TestHelper::eps<ValueType>().

.. code-block:: c++

	LAMA_BOOST_CHECK_EPS( x, y )

Some helpfull Classes
---------------------

TestSparseMatrices
^^^^^^^^^^^^^^^^^^

In the class TestSparseMatrices there are precalculated sparse matrices, which can be useful for some
tests. If other precalculated Matrices should be added this file is the right place for it.

EquationHelper
^^^^^^^^^^^^^^

In the class EquationHelper there are predefined solutionsystems. A system consists of an object of type
Matrix, a solutionvector and a rhs-vector.

Using CTest
-----------

CTest is part of Cmake and helps to integrate various testing executable with various runtime
configurations to the existing buildsystem.  

In an existing CMakeLists.txt file there are two Makros to use:

To activate a testing facility of CMake/CTest you have to use : 

.. code-block:: bash

	ENABLE_TESTING()


After this Makro you are able to add different test executables by using: 

.. code-block:: bash

	ADD_TEST(<NAME> <COMMAND>)

e.g.:

.. code-block:: bash

	ADD_TEST( TestRunName ./testrun )

Sources:

- `Cmake`_

- `LinuxMagazin`_

- `Boost`_

.. _Cmake: http://www.cmake.org/Wiki/CMake_Testing_With_CTest 
.. _LinuxMagazin: http://www.linux-magazin.de/Heft-Abo/Ausgaben/2007/02/Mal-ausspannen 
.. _Boost: https://svn.boost.org/trac/boost/wiki/CMakeTesting 

Archive
-------

Here are some ideas of testing for LAMA :doc:`brainstorming`
