Extension
=========

This tutorial describes how to add new project module in LAMA.

Let us assume that you want to add a new project module with the name ``mymodule`` 
(according to the LAMA coding conventions module names should only have lowercase letters). 

These are the necessary steps:

 * Create the directory ``scai/mymodule``
 * Add a configuration file ``CMakeLists.txt`` and of course the source files
   in this new directory
 * Optionally you can add the sub-directory examples (examples using the functionality of the
   new module, will be installed)
 * Optionally you can add the sub-directory test (Boost unit tests) 
 * Optionally you can add the sub-directory doc (Sphinx user documentation) 
 * Add the module name in the main configuration file ``scai/CMakeLists.txt``
   in the list variable ``SCAI_ALL_MODULES``. It should be in the list
   after all modules that are used (internal dependencies).
 * That's all.


Configuration File of a Module
------------------------------

The best way is to start by creating some directories and by copying some 
configuration files from an exsiting project, e.g. ``tasking``.

.. code-block:: bash

   cd ${LAMA_DIR}/scai

   mkdir mymodule
   mkdir mymodule/doc

   # copy an existing CMake configuration that you have to edit later

   cp tasking/CMakeLists.txt mymodule/CMakeLists.txt

   # copy an existing Sphinx configuration that you have to edit later
   
   cp tasking/doc/index.rst mymodule/doc/index.rst


The following changes of the cmake configuration file ``CMakeLists.txt``
are required:

 * Set the variable ``INTERNAL_DEPS``
 * Set the variable ``EXTERNAL_DEPS``
 * Set class and include files in the ``scai_project`` command

Most of all the other commands should just work fine.

User Documentation
------------------

The following changes of the sphinx configuration file ``index.rst`` are
required:

 * Set the name
 * Fill in the minimal chapters

If you want to add addtional pages you have to add the corresponding
files in the doc directory.

Boost Unit Test
---------------

For setting up a unit test of a new module, the recommended action is to
use one main source file and a source file for each class. The main file
is easily adapted from any main file of an existing module, the unit tests
of a class should have one test for each method.

The CMake configuration file for the test can be set up as follows:

.. code-block:: none

   set ( CXX_SOURCES

            mymoduleTest

            MyClass1Test  
            MyClass2Test  
   )

The following macro defines the unit test. It links it with all
required libraries. The option ``UNIT_TEST`` links the Boost
unit test framework, while the option ``RUN`` enables the
test to be run by CTest.

.. code-block:: none

    scai_test( EXECUTABLE mymoduleTest 
               FILES      ${CXX_SOURCES}
               UNIT_TEST
               RUN )

Optionally you can add a script to get your code be tested by a CI (Continuous Integration) system.

.. code-block:: none

    scai_test_scripts( SCRIPTS       xmltest.sh
                       CODE_COVERAGE ${USE_CODE_COVERAGE} )


Examples
--------

Just what you like to show how to use new module classes.

