Install
=======

To install LAMA to the configured installation path just call

.. code-block:: bash 

   make install

in your build directory.
   
The installation will be placed in the directory that has been specified during the configuration by the variable
CMAKE_INSTALL_PREFIX. This installation directory will be refered later as ``LAMA_ROOT``.

The diretory ``{LAMA_ROOT}/include`` will contain a lot of include files
organized in different directories. These include files are needed when 
using the LAMA library.

The diretory ``{LAMA_ROOT}/lib`` should contain the following libraries:

- liblama.so  
- liblog4lama.so
- liblamacuda.so   (available if CUDA is available)
- libmpistubs.so   (available if MPI is available)
- many other libraries for experimental features
