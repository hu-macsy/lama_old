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

.. code-block:: bash

   make [-j <number_of_processes>] check

Note: Building the test executables and running the test script can be done in two separate steps:

.. code-block:: bash

   make [-j <number_of_processes>] tests
   make test

Please keep in mind that the target ``test`` will not build the unit tests but only run them.

We provide also a more detailed test script that will run example programs on different devices
and uses multiple processes that communicate or interact with each other.

Test programs are only exploited during the build phase, they will not be installed later.

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

   make doc_lama_all
   <browser> doc/user/lama_all/html/index.html

Doxygen API Documentation
"""""""""""""""""""""""""

If Doxygen is found the API documentation can be build as follwos:

.. code-block:: bash

   make doxygendoc
   <browser> doc/system/html/index.html

In contrary to the user documentation, the API documentation is always built for
all SCAI modules.
