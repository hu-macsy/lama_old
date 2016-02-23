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
