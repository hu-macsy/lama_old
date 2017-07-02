Build
-----

After a successful configuration a Makefile is generated in the build directory and the whole directory tree. The make utility takes this first makefile as input file for building LAMA.

Building and Installing Libraries and Executables
-------------------------------------------------

To build LAMA just invoke make in your build directory. Parallel builds are supported.

.. code-block:: bash

   make [-j <number-of-build-processes>]

Due to the use of ExternalProjects for the LAMA subprojects everything is installed as well with the build process, so **NO** additional

.. code-block:: bash

   make install

call is needed. The installation will be placed in the directory that has been specified during the configuration by the variable **CMAKE_INSTALL_PREFIX**. This installation directory will be refered later as **SCAI_ROOT**.

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


Documentation
^^^^^^^^^^^^^

If Sphinx is found the directory **${SCAI_ROOT}/share/doc** contains automatically the subdirectory **user/html/scai-x.x.x** containing this user documentation. If also Doxygen is found the API documentation can be build by additionally calling

.. code-block:: bash

   make doxygendoc

in your build directory and you find the second subdirectory **system** with the API documentation.

If you do not want CMake to do build the user documentation you have to pass -DBUILD_DOC=OFF to CMake.

User Documentation
""""""""""""""""""

This user documentation contains subdirectories for each subproject with the naming **scai-<libname>-<lib-version>**. They should give an overview on how to use the respective library.

The directory **scai-documentation** contains the overall documentations for the scai_libs (installation, tutorials, faqs).
  
API Documentation
"""""""""""""""""

You can open the main API documentation at **${SCAI_ROOT}/share/doc/system/index.html**. It will give you insights of all LAMA classes and their relations.

Examples
^^^^^^^^

Every subproject is shipped with its own examples. You will find them in your build directory in the subdirectory **<project>/examples**. Additionally they will be installed in the folder **${SCAI_ROOT}/share/examples** where they are grouped by subprojects.  

In the build directory the examples will be build directly by CMake. But you can adapt, rebuild and execute them there, if you want to with:

.. code-block:: bash

   cd <build_directory>/lama/examples/tutorial
   make
   ./simple.exe

The same examples in **${SCAI_ROOT}/share/examples** and can be built by a corresponding **Makefile**.
All these Makefiles will include a file **make.inc** that has been generated in such a way that it fits your installation. 
If the build directory is no more available you can copy the examples in your own directory from the installation directory.

.. code-block:: bash

   mkdir myExamples
   cp -r ${installation_directory}/share/examples/* .
   make

Take a look at the README files, where available, to find some more information about how to run the example programs.

You can use one of the Makefiles together with a make.inc file to compiler you own LAMA applications.
