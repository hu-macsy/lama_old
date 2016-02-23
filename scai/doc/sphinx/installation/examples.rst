Examples
========

Every subproject is shipped with its own examples. You will find them in the subdirectory <project>/examples
in your build directory. Additionally they will be installed in ``SCAI_ROOT`` in the folder share/examples
where they are grouped by subprojects.  

In the build dir the examples will be directly build by CMake. But you can adapt, rebuild and execute them here, if you want to.

.. code-block:: bash 

   cd <build_directory>/lama/examples/tutorial
   make
   ./simple.exe

The same examples can be found in the ``SCAI_ROOT/share/examples`` and can be built by a corresponding ``Makefile``.
All these Makefiles will include a file ``makefile.inc`` that has been generated in such a way that it fits your installation. 
If the build directory is no more available you can copy the examples in your own directory from the installation directory.

.. code-block:: bash 

   mkdir myExamples
   cp -r ${installation_directory}/share/examples/* .
   make

Take a look at the README files where you can find some more information about how to run the example programs.
