Build
-----

After a successful configuration a Makefile is generated in the build directory. The make utility takes this makefile as input file for building LAMA.

Building Libraries
------------------

To build LAMA just invoke make in your build directory. Parallel builds are supported.

.. code-block:: bash

   make [-j <number-of-build-processes>]

In contrary to previous versions of LAMA, installation is completely separated. So you must call for installation:

.. code-block:: bash

   make install

The installation will be placed in the directory that has been specified during the configuration by the variable **CMAKE_INSTALL_PREFIX**. This installation directory will be refered later as **SCAI_ROOT**.

Header
^^^^^^

The directory **${SCAI_ROOT}/include** will contain a lot of include files organized in seperate directories for the subprojects. These include files are needed when using the LAMA library.

Libraries
^^^^^^^^^

The directory **${SCAI_ROOT}/lib** should contain the following libraries:

- libscai_common.so
- libscai_logging.so
- libscai_tracing.so
- libscai_tasking.so
- libscai_hmemo.so
- libscai_kregistry.so
- libscai_blaskernel.so
- libscai_utilskernel.so
- libscai_sparsekernel.so
- libscai_dmemo.so
- libscai_lama.so
- libscai_solver.so

Testing
-------

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
--------

All examples are built by the following command:

.. code-block:: bash

   make [-j <number_of_processes>] examples

Note: this command only compiles the example programs but will not run them. This might be done in the different
subdirectories **<module>/examples**. Here you also find hints which arguments can be used to run the example programs.


In the build directory the examples will be build directly by CMake. But you can adapt, rebuild and execute them there, if you want to with:

.. code-block:: bash

   cd <build_directory>/lama/examples/tutorial
   make
   ./simple.exe


Installing Libraries and Executables
------------------------------------

In contrary to previous versions of LAMA, installation is completely separated. So you must call for installation:

.. code-block:: bash

   make install

The installation will be placed in the directory that has been specified during the configuration by the variable **CMAKE_INSTALL_PREFIX**. This installation directory will be refered later as **SCAI_ROOT**.

Header
^^^^^^

The directory **${SCAI_ROOT}/include** will contain a lot of include files organized in seperate directories for the subprojects. These include files are needed when using the LAMA library.

Libraries
^^^^^^^^^

The directory **${SCAI_ROOT}/lib** should contain the following libraries:

- libscai_common.so
- libscai_logging.so
- libscai_tracing.so
- libscai_tasking.so
- libscai_hmemo.so
- libscai_kregistry.so
- libscai_blaskernel.so
- libscai_utilskernel.so
- libscai_sparsekernel.so
- libscai_dmemo.so
- libscai_lama.so
- libscai_solver.so

Examples
^^^^^^^^

All example programs will be installed in the folder **${SCAI_ROOT}/share/examples** where they are grouped by subprojects.

Note: Here you will also find ``Makefile``s  that might be used to compile the example programs without explicit use of cmake.

ToDo: we will provide very soon also CMake configuration files to compile LAMA programs using a ready LAMA installation.

Documentation
^^^^^^^^^^^^^

If Sphinx is found the user documentation can be built as follows:

.. code-block:: bash

   make doc_lama_all
   <browser> doc/user/lama_all/html/index.html

If also Doxygen is found the API documentation can be build by additionally calling

.. code-block:: bash

   make doxygendoc

User Documentation
""""""""""""""""""

This user documentation contains subdirectories for each subproject with the naming **scai-<libname>-<lib-version>**. They should give an overview on how to use the respective library.

The directory **scai-documentation** contains the overall documentations for the scai_libs (installation, tutorials, faqs).
  
API Documentation
"""""""""""""""""

You can open the main API documentation at **${SCAI_ROOT}/share/doc/system/index.html**. It will give you insights of all LAMA classes and their relations.

Examples
^^^^^^^^
The same examples in **${SCAI_ROOT}/share/examples** and can be built by a corresponding **Makefile**.
All these Makefiles will include a file **make.inc** that has been generated in such a way that it fits your installation. 
If the build directory is no more available you can copy the examples in your own directory from the installation directory.

.. code-block:: bash

   mkdir myExamples
   cp -r ${installation_directory}/share/examples/* .
   make

Take a look at the README files, where available, to find some more information about how to run the example programs.

You can use one of the Makefiles together with a make.inc file to compiler you own LAMA applications.
