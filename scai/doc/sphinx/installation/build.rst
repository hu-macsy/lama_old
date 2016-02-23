Build
=====

After a successful configuration a Makefile is generated in the build directory.
The make utility takes this input file for building LAMA.

Building and Installing Libraries and Executables
-------------------------------------------------

To build LAMA just invoke make in your build directory. Parallel builds are
supported.

.. code-block:: bash 

   make [-j <number-of-build-processes>]

Due to the use of ExternalProjects for the LAMA subprojects everything is installed as well with the build process, so no additional

.. code-block:: bash 

   make install

call is needed. The installation will be placed in the directory that has been specified during the configuration by the variable
CMAKE_INSTALL_PREFIX. This installation directory will be refered later as ``SCAI_ROOT``.

The directory ``{SCAI_ROOT}/include`` will contain a lot of include files organized in seperate directories for the subprojects.
These include files are needed when using the LAMA library.

The directory ``{SCAI_ROOT}/lib`` should contain the following libraries:

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

User Documentation
------------------

If Sphinx is found by CMake the user documentation will be automatically built. 
If you don't want CMake to do this you have to pass -DBUILD_DOC=OFF to CMake.
You find the documentation in your build directory and find it in doc/sphinx/html/index.html or 
in the installation dir in share/doc/user/html/scai-<verion>/index.html.
  
API Documentation
-----------------

To build LAMAs API documentation call

.. code-block:: bash 

   make doxygendoc

in your build directory and find it in doc/doxygen/html/index.html or
in the installation dir in share/doc/system/index.html.
