Using different contexts
========================

LAMA can use different hardware architectures. This is implemented by using the context classes. Especially for testing it is profitable to have a construct to run a testcase on every available context.
This construct is called a Contextloop, which can be found in TestMacros.hpp. Using a CONTEXTLOOP() { ... } in combination with the makro GETCONTEXT( contextptr ) provides a facility to iterate 
over all available contexts and test the same function on different hardware architectures.

The usage of this loop is shown here ( part of CGTest.cpp ):

:: 

	BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithoutPreconditioning, T, test_types )
	{
	    typedef T ValueType;
	
	    CONTEXTLOOP()
	    {
	        GETCONTEXT( contextptr );
	        testSolveWithoutPreconditionmethod< CSRSparseMatrix<ValueType> >( contextptr );
			// ...
	    }
	}

The makro GETCONTEXT(contextptr) creates the required contextpointer, which will be used for the setContext-method of the testobject (here CGSolver).
The method testSolveWithoutPreconditionmethod contains the implementation of this test and will be invoked by the loop on all contexts.
