Examples
========

Every subproject is shipped with it's own examples. You will find them in the subdirectory <project>/examples
in your build directory. Additionally they will be installed in the installation path in the folder share/examples
where they are grouped by subprojects.  

.. You will find some example programs in the subdirectory ``examples`` of the build directory.

In the build dir the examples will be directly build by cmake. The ones in the installation path can
be built by a corresponding ``Makefile``.

.. The example programs will be built by cmake; they are directly built by a corresponding ``Makefile``.
 
All these Makefiles will include a file ``makefile.inc``
that has been generated in such a way that it fits your installation. 

.. The example programs can only be compiled after successful installation of LAMA.

.. code-block:: bash 

   cd <build_directory>/lama/examples
   ls
   cd tutorial
   make
   ./simple.exe

If the build directory is no more available you can copy the examples in your own
directory from the installation directory; they are in the subdirectory ``share/examples``.

.. code-block:: bash 

   mkdir myExamples
   cp -r ${installation_directory}/share/examples/* .
   make

Take a look at the README files where you can find some more information about how to
run the example programs.
