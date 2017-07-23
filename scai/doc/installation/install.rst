Install
-------

In contrary to previous versions of LAMA, installation is completely separated. So you must call for installation:

.. code-block:: bash

   make install

The installation will be placed in the directory that has been specified during the configuration by the variable **CMAKE_INSTALL_PREFIX**. 
This installation directory will be refered later as **SCAI_ROOT**.

- Header files, libraries and example programs will always be installed.
- Test programs are only used during the build phase and will never be installed.
- Documentation is only installed if it has been built before.

Header
^^^^^^

The directory **${SCAI_ROOT}/include** will contain a lot of include files organized in seperate directories for the subprojects. 
These include files are needed when using the LAMA library.

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

Note: Here you will also find ``Makefile`` files that might be used to compile the example programs without explicit use of cmake.

ToDo: we will provide very soon also CMake configuration files to compile LAMA programs using a ready LAMA installation.

The same examples in **${SCAI_ROOT}/share/examples** and can be built by a corresponding **Makefile**.

All these Makefiles will include a file **make.inc** that has been generated in such a way that it fits your installation. 
You can copy the examples in your own directory from the installation directory.

.. code-block:: bash

   mkdir myExamples
   cp -r ${installation_directory}/share/examples/* .
   make

Take a look at the README files, where available, to find some more information about how to run the example programs.
You can use one of the Makefiles together with a make.inc file to compiler you own LAMA applications.

Documentation
^^^^^^^^^^^^^

This user documentation contains subdirectories for each subproject with the naming **scai-<libname>-<lib-version>**. 
They should give an overview on how to use the respective library.

The directory **scai-documentation** contains the overall documentations for the scai_libs (installation, tutorials, faqs).
  
API Documentation
"""""""""""""""""

You can open the main API documentation at **${SCAI_ROOT}/share/doc/system/index.html**. It will give you insights of all LAMA classes and their relations.

