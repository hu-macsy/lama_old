:orphan:

Brainstorming
=============

The following list contains ideas of a brainstorming for "testing":

- Given situation:

  - CPPUNIT is used by everyone; everyone knows how to handle with it (perfectly or not).
  
  - CPPUNIT is the only testing-framework here, which is used for unit testing.
  
  - nobody has used or uses a different framework and nobody has the perfect knowledge about testing.
  
  - actually everyone is satisfied with CPPUNIT in general.
  
  - so the question is: do we need a second/different framework ?!

- Better enclosure of running tests: There should be a possibility to ...

  - ... run all tests (CPPUNIT already includes this feature).
  
  - ... run groups of test classes.
  
  - ... run just one test class (CPPUNIT already includes this feature).
  
  - ... run just one test method of a class.
  
  - ... see the progress of tests (a test with an endless loop is very difficult to fix)
  
  - ... have a verbose mode to see the test that is currently running

- Possibility of preinitialised fixtures:

  - CPPUNIT includes this feature in its fixture-class with setUp() and tearDown()-Methods...
  
  - but maybe its better and useful to have a global possibility (e.g. like Equation_Helper-class or maybe
  	some makros) to have initialised matrix-systems with prepaired testdata.

- Size of tests:

  - Do we need such a huge amount of test data ? or isn't a small amount of data (e.g. 3x3-matrix instead of
  	10x10-matrix or bigger) enough for testing functionalities ?
  
  - Do we really need a test or some tests for EACH function we develop? The test-classes/methods and the
  	build-time will grow extremly...

- General testing:

  - maybe its useful to run regularly all tests (in the night) for reporting defects.
  
    - Continuous Integration (Nightly Builds + Running all Test is definitely needed in the long term)
    
  - its a response for each developer and the correctness of the functions.

- Miscellaneous:

  - The most maybe all testclasses use "sblas". I have heard that it's out of date ...
  
  - Switching to boost Test would remove the dependency to cppunit without really introducing a new dependency
  	(because we are using boost anyhow)
  
  - It should be checked regularly that all codes compiles and runs with diferent compilers (gcc, intel, cl,
  	pgi, ...) on different platforms (linux, win32, macos?)
  
  - Is it possible to integrate the execution of parallel Tests (mpi) seamlessly, especially with automated
  	tests.
  
  - Is it possible to have an integration of valgrind checks?
  
  - It should be easy to add and (de)activate tests in development time

- Code coverage / Tests

  - There are not many tests that test for failures, i.e. wrong arguments give well-explained exceptions and
  	do not end in system crashes
  
  - Code coverage would identify branches in the code that have not been tested at all, also failures
  
  - Gnu and Intel compiler provide facilities for code coverage (measuring int)

The desired situation would be that all the code compiles on a set of compilers and a set of platforms warning
free and all tests run without any failure and valgrind warning.
