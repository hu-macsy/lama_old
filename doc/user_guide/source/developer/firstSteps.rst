First Steps
===========

This description should give you a good start in programming with and in LAMA. It will show you the requirements for a 
smooth implementation process, how to work with the versioning tool git and how to properly build LAMA.

Requirements
^^^^^^^^^^^^

- Java SDK and Eclipse IDE C/C++ have to be both either 32bit or 64bit (esp. Linux)! 
- Please, also refer to the Eclipse documentation in order to avoid problems with non-suitable java versions
- The environment variable **JAVA_HOME** has to exist and point to the right location

.. code-block:: bash

	JAVA_HOME=<install_dir_java>/bin

- CMake (Version > 2.8) and cppunit have to be installed and be present in the standard executable path
- Doxygen has to be present on the system, if API documentation is supposed to be generated

SVN
---

Installation
^^^^^^^^^^^^

In order to check out this project's source code, it is advisable to rely on `Subclipse`_ for Ecplise or
`Tortoisesvn`_ for Windows. Please refer to the documentation on the respective sites for help with the
installation process.

.. _Subclipse: http://subclipse.tigris.org/
.. _Tortoisesvn: http://tortoisesvn.tigris.org/

Working with SVN
^^^^^^^^^^^^^^^^

Creating a Makefile with CMake
-------------------------------

The subfolder **build** is going to contain the complete build of this project.
In order to be able to use the normal **make** a makefile has to be generated.
This project is using CMake for this purpose. In order to create the Makefile with all dependencies for
a given system CMake has to be called in the following way:

.. code-block:: bash

	$ cmake [-D CMAKE_INCLUDE_PATH=<install_dir_cmake>/include \]
        [-D CMAKE_LIBRARY_PATH=<install_dir_cmake>/lib \]
        [-D DOXYGEN_EXECUTABLE=<doxygen_bin> \]
        <src_dir>

in a UNIX-like environment where *<src_dir>* points to the folder containing the source.

The option

.. code-block:: bash

	-D DOXYGEN_EXECUTABLE

defines the path to the doxygen executable, if not present in the standard path.

In order to execute CMake with settings for the Intel Compiler the following has to be specified (assuming a
bash shell):

.. code-block:: bash

	$ export CC=icc
	$ export CXX=icpc

before calling *cmake*

Overview: Makefile targets
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Compiling**

Compile everyting

.. code-block:: bash

	make

Compile the C-part of the library

.. code-block:: bash

	make lama

Compile the C++ and the C-part of the library

.. code-block:: bash

	$ make lama++

**Unit Tests**

Compile unit tests

.. code-block:: bash

	$ make lama_unit_test

Run the unit tests

.. code-block:: bash

	make test

Run tests with valgrind

.. code-block:: bash

	make test_valgrind

Run a specific test


.. code-block:: bash

	make <testname>

for a complete list of testnames refer to :ref:`testList`.
Run that specific test with valgrind

.. code-block:: bash

	make <testname>_valgrind

Create the API documentation with Doxygen

.. code-block:: bash

	make doc

**Cleaning up**

Remove executables and dependencies which are no longer needed:

.. code-block:: bash

	make clean

Remove everything which has been created after checkout

.. code-block:: bash

	make distclean
