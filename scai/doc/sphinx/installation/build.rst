Build
=====

After successful configuration a Makefile is generated in the build directory.
The make utility takes this input file for building LAMA.

The Makefile includes all source code dependencies, so if you change any source
file of LAMA, a further call of make will recompile everything (and only this)
that depends on the changed source files. 

Building libraries and executables
----------------------------------

To build LAMA just invoke make in your build directory. Parallel builds are
supported.

.. code-block:: bash 

   make [-j <number-of-build-processes>]

API Documentation
-----------------

To build LAMAs doxygen API documentation call

.. code-block:: bash 

   make doxygendoc

in your build directory and find it in doc/doxygen/html/index.html or
in the installation dir in share/doc/system/index.html.

User Documentation
------------------

If Sphinx is found by CMake the documentation will be automatically built. 
If you don't want CMake to do this you have to pass -DBUILD_DOC=OFF to cmake.
  
To build LAMAs user documentation call

.. code-block:: bash 

   make doc

in your build directory and find it in doc/sphinx/html/index.html or 
in the installation dir in share/doc/user/html/scai-<verion>/index.html.

Install
=======

To install LAMA to the configured installation path just call

.. code-block:: bash 

   make install

in your build directory.
   
The installation will be placed in the directory that has been specified during the configuration by the variable
CMAKE_INSTALL_PREFIX. This installation directory will be refered later as ``SCAI_ROOT``.

The directory ``{SCAI_ROOT}/include`` will contain a lot of include files
organized in different directories. These include files are needed when 
using the LAMA library.

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