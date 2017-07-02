Extension
=========

This tutorial describes how to add new project module in LAMA.

Let us assume that you want to add a new project module with the name ``mymodule`` 
(module names should only have lowercase letters). 
The best way is to start by creating some directories and by copying some 
configuration files from an exsiting project, e.g. ``tasking``.

.. code-block:: c++

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

The following changes of the sphinx configuration file ``index.rst`` are
required:

 * Set the name
 * Fill in the minimal chapters

If you want to add addtional pages you have to add the corresponding
files in the doc directory.


