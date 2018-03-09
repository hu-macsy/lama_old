Build
-----

After a successful configuration a Makefile is generated in the build directory. 
The make utility takes this makefile as input file for building LAMA.

To build LAMA just invoke make in your build directory. Parallel builds are supported.

Libraries
^^^^^^^^^

.. code-block:: bash

   make [-j <number-of-build-processes>]

By this command, a library for each SCAI module is built.

For building unit tests, running checks, compiling example programs and building
documentation, separate make commands must be issued.

Tests
^^^^^

Building and running tests are optional features. It is recommended to
execute this step to verfiy a correct compilation of the LAMA libraries.

All unit tests (and also some other test scripts) are built and run by this command:

.. code-block:: bash

   make [-j <number_of_processes>] check

Each subproject comes with a test. Mostly it is one unit test executable, but there
are some exceptions. The test executables run via **CTest**, the testing tool
distributed as a part of CMake. The output of the above command should be as follows:

.. code-block:: none

   Test project <build-director>
         Start  1: commonTest
    1/14 Test  #1: commonTest .......................   Passed    6.07 sec
         Start  2: logging_test
    2/14 Test  #2: logging_test .....................   Passed    0.35 sec
         Start  3: tracing_test
    3/14 Test  #3: tracing_test .....................   Passed    2.05 sec
         Start  4: taskingTest
    4/14 Test  #4: taskingTest ......................   Passed    2.64 sec
         Start  5: kregistryTest
    5/14 Test  #5: kregistryTest ....................   Passed    0.00 sec
         Start  6: hmemoTest
    6/14 Test  #6: hmemoTest ........................   Passed    0.01 sec
         Start  7: blaskernelTest
    7/14 Test  #7: blaskernelTest ...................   Passed    0.08 sec
         Start  8: utilskernelTest
    8/14 Test  #8: utilskernelTest ..................   Passed    0.30 sec
         Start  9: sparsekernelTest
    9/14 Test  #9: sparsekernelTest .................   Passed    0.16 sec
         Start 10: dmemoTest
   10/14 Test #10: dmemoTest ........................   Passed    0.12 sec
         Start 11: lamaTest
   11/14 Test #11: lamaTest .........................   Passed    1.14 sec
         Start 12: lamaStorageTest
   12/14 Test #12: lamaStorageTest ..................   Passed    0.95 sec
         Start 13: lamaMatrixTest
   13/14 Test #13: lamaMatrixTest ...................   Passed    4.68 sec
         Start 14: solverTest
   14/14 Test #14: solverTest .......................   Passed   15.75 sec
   
   100% tests passed, 0 tests failed out of 14

   Total Test time (real) =  34.31 sec

If the tests are built the first time using multiple build processes is highly recommended.

Building the test executables and running the test script can be done in two separate steps:

.. code-block:: bash

   make [-j <number_of_processes>] tests
   make test

Please keep in mind that the target ``test`` will not build the tests but only run them.

We provide also a more detailed test script that will run example programs on different devices
and uses multiple processes that communicate or interact with each other.

Test programs are only exploited during the build phase, they will not be installed later.
Please keep also in mind, that the unit tests are not built if the Boost unit test framework
has not been found or explicitly disabled.

Examples
^^^^^^^^

All examples are built by the following command:

.. code-block:: bash

   make [-j <number_of_processes>] examples

Note: this command only compiles the example programs but will not run them. This might be done in the different
subdirectories **<module>/examples**. Here you also find hints which arguments can be used to run the example programs.

Documentation
^^^^^^^^^^^^^

LAMA provides two kind of documentation, one is the user documentation (built by Sphinx) and one
the system or API documentation (built by doxygen).

.. code-block:: bash

   make doc

Sphinx User Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^

If Sphinx is found the user documentation can be built as follows:

.. code-block:: bash

   make doc_libama
   <browser> doc/user/libama/html/index.html

Doxygen API Documentation
"""""""""""""""""""""""""

If Doxygen is found the API documentation can be built as follwos:

.. code-block:: bash

   make doxygendoc
   <browser> doc/system/html/index.html

In contrary to the user documentation, the API documentation is always built for
all SCAI modules.

Build All
"""""""""

If you want to build all libraries, tests, examples and documentation at once, 
you can do it and you might benefit of the full potential of parallel build:

.. code-block:: bash

   make [-j <number_of_build_processes] check examples doc


