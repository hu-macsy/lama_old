Using Boost.Test
================

Boost.Test offers different runtime parameters, which simplify the usage of this testframework.

- run_test: Use this runtime parameter to run a single testcase:
 
Invoke of a single testcase:

.. code-block:: bash

	./lama_test --run_test=Testsuitename/Testcasename
	
Invoke a subgroup of tests by using regular expressions :

.. code-block:: bash

	./lama_test --run_test=*Vector*/*Multiplication*  
	
In this example all testcases with the expression "Multiplication" in all testsuites with the expression "Vector" are
executed.

- log_level: Receiving more information about the testrun, you can change the loglevel of Boost.Test.
 
E.g.: 

.. code-block:: bash 

	./lama_test --log_level=test_suite

There are different possible levels: all, test_suite, message, warning and many more ...

You can take a look at the summary with all runtime parameters and environment variables given by Boost.Test here_.

.. _here: http://www.boost.org/doc/libs/1_45_0/libs/test/doc/html/utf/user-guide/runtime-config/reference.html

- context: Additionally we have implement this runtime parameter. LAMA can execute calculations on different hardware 
architectures. For testing it is helpful to run a test just on a specific architecture.

Valid contexts are:
 
 - Host
 - CUDA
 - OpenCL ( not yet implemented )

The following example shows, how to run the testcase "clearTest" in the testsuite CSRStorageTest just on CUDA:

.. code-block:: bash

	./lama_tests --run_test=CSRStorageTest/clearTest --context=CUDA
	
For our runtime parameter "context" there will be a environment variable aswell. The following example shows how to use it:

.. code-block:: bash

	export LAMA_TEST_CONTEXT=CUDA 
