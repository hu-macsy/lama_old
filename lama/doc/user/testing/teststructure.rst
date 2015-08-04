Teststructure
=============

Boost.Test offers different Makros to create testsuites and testcases.

Standard test class
-------------------

The structure of a classic TestClass.cpp in LAMA looks like this:

::

	#include <boost/test/unit_test.hpp>

	typedef boost::mpl::list<double,float>  test_types;

        using namespace boost;
        using namespace lama;

    /* --------------------------------------------------------------------- */

	struct TestClassConfig
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

	BOOST_FIXTURE_TEST_SUITE( TestClass, TestClassConfig );

    /* --------------------------------------------------------------------- */

	BOOST_AUTO_TEST_CASE( test1 )
	{ 
    	//Test Case 1 
	}
    
    BOOST_AUTO_TEST_CASE( test2 )
    { 
    	//Test Case 2 
    }

	BOOST_AUTO_TEST_CASE_TEMPLATE( test3, T, test_types )
    { 
    	//Test Case 3 (template-testcase) 
    }

    /* --------------------------------------------------------------------- */

	BOOST_AUTO_TEST_SUITE_END();

Each testclass has to inlude the headerfile unit_test.hpp from Boost.Test. Testobjects, that are used by many testcases in a testsuite can be created and deleted in a struct. 
We call is a TestClassConfig. The constructor and destructor creates and deletes these common used objects. This struct is the second argument of the makro "BOOST_FIXTURE_TEST_SUITE".
In Boost.Test there are two different testcases: BOOST_AUTO_TEST_CASE(name) and BOOST_AUTO_TEST_CASE_TEMPLATE( name, T, types). Examples for both types of testcases are shown above. 
If it is not neccessary to create those common used objects, the alternative for creating a testhierarchy, is to use the makro BOOST_AUTO_TEST_SUITE( suitename ).
Both kinds of testsuite have to be closed with the makro BOOST_AUTO_TEST_SUITE_END();

Inherited testclasses
---------------------

A lot of classes in LAMA inherits from abstract classes, e.g. all concrete storages: CSRStorage, COOStorage, ... inherits from class MatrixStorage. You can find these inheritances again in the testhierarchies.
For each inheritance we have created a common used test class, which tests all functions, that are given by the abstract class. CSRStorageTest.cpp, COOStorageTest.cpp, ... use the testclass MatrixStorageTest.cpp to test all common methods.

This code from a concrete testclass, like CSRStorageTest.cpp shows, how to invoke those inheritances:

::

	BOOST_AUTO_TEST_CASE_TEMPLATE( commonTestCases, T, test_types )
	{
	    typedef T ValueType;
	
	    CSRStorage<ValueType> csrStorage;
	    MatrixStorageTest<ValueType> matrixstorageTest( csrStorage );
	
	    if ( base_test_case )
	    {
	        MATRIXSTORAGE_COMMONTESTCASES( storageTest );
	    }
	    else
	    {
	    	matrixstorageTest.runTests();
	    }
	}

In those cases it is neccessary to create two objects: one concrete testobject ( here CSRStorage ) and a object 
of your common used test class ( here MatrixStorageTest ), using the concrete testobject as an argument.
Calling the method matrixstorageTest.runTests() will invoke all methods from the base testclass with this testobject.

The corresponding testclass (MatrixStorageTest.hpp) looks like this:

::
	
	static std::string storagetestclasses[] = { "CSRStorageTest", "COOStorageTest",
	                                            //other concrete testclasses 
	                                          };
	
	static std::string storagetestmethods[] = { "purgeTest", "setIdentityTest", "setCSRDataTest",
	                                            //other testmethods 
	                                          };
	
	template<typename T>
	class MatrixStorageTest
	{
	public:
	
	    typedef T ValueType;
	
		MatrixStorageTest( lama::MatrixStorage<T>& storage ) : mMatrixStorage( storage ) {};
	
	    void purgeTest();
		//all other definitions of testmethods here
	
	    void runTests();
	
		lama::MatrixStorage<T>& mMatrixStorage;

	};

	#define MATRIXSTORAGE_COMMONTESTCASES( testinstance )                   	\
	{   COMMONTESTCASEINVOKER( testinstance, purgeTest );                   	\
	    COMMONTESTCASEINVOKER( testinstance, /*all_other_testmethods here*/ ); 	\																				
 	}
 	 	

The makros MATRIXSTORAGE_COMMONTESTCASES, COMMONTESTCASEINVOKER and the two std::strings storagetestclasses 
and storagetestmethods are neccessary to invoke single testmethods from the common used test class, using the runtime parameter --run_test from Boost.Test. All those testmethods (e.g. purgeTest) are not registered automatically in the testhierarchy by Boost.Test.

The file MatrixStorageTest.cpp looks like this:

::

	//include headerfiles & declare namespaces
	
	LAMA_LOG_DEF_TEMPLATE_LOGGER(template<typename T>, MatrixStorageTest<T>::logger, "Test.MatrixStorageTest" );
	
	
	LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, StorageType, purgeTest )
		//test here
	LAMA_COMMON_TEST_CASE_END();
	
	
	LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, StorageType, writeAtTest )
		//test here
	LAMA_COMMON_TEST_CASE_TEMPLATE_END();
	
	/* ------------------------------------------------------------------------- */
	
	LAMA_COMMON_TEST_CASE_RUNNER_TEMPLATE( MatrixStorageTest )
	{
	    purgeTest();
		//calling all other testmethods
	}

The makros LAMA_COMMON_TEST_CASE_TEMPLATE, or LAMA_COMMON_TEST_CASE for a non-templated class, are useful to get some extra output,
if you invoke the test run with loglevel=test_suite. LAMA_COMMON_TEST_CASE_RUNNER_TEMPLATE or LAMA_COMMON_TEST_CASE_RUNNER encapsulates the invokes of all testmethods.

Examples of common base test class are:

- NormTest
- SparseMatrixTest
- MatrixStorageTest
- DistributionTest
- CommunicatorTest

Regular Expressions in common test classes
------------------------------------------

Boost.Test includes the facility to use regular expressions to invoke a subgroup of tests. Because we have implemented 
these common base classes, it is not possible to invoke them with regular expressions and the logic of Boost.Test.
The test executables are adapted, to invoke even these common base classes by expressions. The following example demonstrates 
a call of a testmethod, which is part of a common test class.

.. code-block:: bash

	./lama_test --run_test=MaxNormT*/Zero*

This command will invoke ZeroVectorTest of class MaxNormTest.
