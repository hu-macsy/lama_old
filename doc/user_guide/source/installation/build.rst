Build
-----

After successful configuration a Makefile is generated in the build directory.
The make utility takes this input file for building LAMA.

The Makefile includes all source code dependencies, so if you change any source
file of LAMA, a further call of make will recompile everything (and only this)
that depends on the changed source files. 

Building libraries and executables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To build LAMA just invoke make in your build directory. Parallel builds are
supported.

.. code-block:: bash 

   make -j <number-of-build-processes>

API Documentation
^^^^^^^^^^^^^^^^^
To build LAMAs doxygen API documentation call

.. code-block:: bash 

   make doc

in your build directory.

