General information about testing in LAMA
-----------------------------------------

All testhierarchies in LAMA were made by using the testframework Boost.Test, which is part of Boost_.

At the moment there are three different test hierarchies: 

 - test/lama_test:
	this executable includes in general all seriell tests. Furthermore there are many testcases, which
	contains a contextloop, to execute a test on different hardware architectures.

 - test/distributed/lama_dist_test:
	All tests, which are made for a parallel execution or are part of the parallel environment, like
	different distributions are grouped in this executable.

 - test/cuda/lama_cuda_test:
	Testcases, that verify specific implementations to execute something on CUDA, are summarized here.
	To reduce code duplications you can find testscases, whose implementations are equal to testcases for
	Host in the executable for seriell tests. Method setContext will invoke there many times to run tests,
	like all Solvertests, Matrixtests or Storagetests even on CUDA and other architectures. The aim is to
	reduce the content of the CUDA-testdirectory as much as possible.

.. _Boost: www.boost.org 

For each class in scai/lama/ there is a corresponding testclass with the same name (e.g. MaxNorm.cpp -> MaxNormTest.cpp).
 
In every testclass there is exactly one testsuite, which groups testcases, that are made for testing functions from
this sourceclass. All names of testcases in a testsuite should be unique.

There are some conventions about name-giving as well:

 - Each TestSuiteName should be equal to the FileName. This makes it easy and intuitive to call testcases by using the runtime parameter --run_test. 
 - Each TestFixtureName is created by adding "Config" to the TestSuiteName. It seems that Boost.Test has some problems with equal Fixturenames of different Testclasses. With this convention we make sure that there are just unique names of testclasses and testfixtures.
