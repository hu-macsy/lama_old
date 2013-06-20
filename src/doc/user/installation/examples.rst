Examples
========

You will find some example programs in the subdirectory ``examples`` of the build directory.

The example programs will not be built by cmake; they are directly built by a 
corresponding ``Makefile``. All these Makefiles will include a file ``makefile.inc``
that has been generated in such a way that it fits your installation. The example
programs can only be compiled after successful installation of LAMA.

.. code-block:: bash 

   cd <build_directory>/examples
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
