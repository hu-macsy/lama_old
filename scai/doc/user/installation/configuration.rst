Configuration
=============

The configuration of LAMA tries to find required and optional software needed by LAMA and to identify correct
compilation and link flags. At the end of the configuration, a ``Makefile`` or ``Visual Studio project`` is generated
that can be used to build libraries and executables of LAMA.

This section will guide you through the Cmake configuration. If all packages are installed in the system path there is
less to do. If not or you have more than one option of a package you need to give Cmake some hints where to find the
packages or which one to choose. For a specific package please refer to the belonging subpage.

.. toctree::
   :titlesonly:

   configuration/c++
   configuration/boost
   configuration/blas
   CUDA <configuration/cuda>
   MPI <configuration/mpi>
   Tracing <configuration/tracing>

The following steps are for generating Makefiles under Linux. For the creation of a Visual Studio project under Windows
see :doc:`here <windowsTipps>`.

For its configuration LAMA uses the CMake utility. This utility is called in the build directory. Although in-source
builds are supported we recommend an out-of-source build. So create a build directory and change to it.

.. code-block:: bash 

   mkdir build
   cd build

To configure LAMA, call cmake in your build directory and specify the path to the LAMA source (**<project-root>/src**).
There are a couple of advanced options for defining special build-types or to switch on/off particular features.
The general cmake-call looks like this:

.. code-block:: bash

   cmake [options] <project-root>/src

For defining the options you need to write ''-D <VARIABLE_NAME>=value''.

If you have not installed CMake in the system path you need to specify to CMake include and library path by:

.. code-block:: bash

   cmake [-D CMAKE_INCLUDE_PATH=<install/dir/cmake>/include \]
         [-D CMAKE_LIBRARY_PATH=<install/dir/cmake>/lib \]
         <project-root>/src
  
All relevant LAMA options and their significations are listed in the table below.

.. csv-table:: 

   :header: "CMake variable", "Signification"

   "CMAKE_BUILD_TYPE", "Debug (default), or Release (for an optimized build)"
   "CMAKE_INSTALL_PREFIX", "Installation path of LAMA library"
   "LAMA_ADDITIONAL_LINK_LIBRARIES", "Paths to libraries that should be linked to lama, separated by ;"
   "LAMA_ADDITIONAL_LINK_FLAGS", "Additional linker Flags"
   "SCAI_ASSERT_LEVEL", "ERROR recommended, DEBUG adds much more assertions into the code"
   "LAMA_CMAKE_VERBOSE", "set to 1 or TRUE for more detailed infos about configuration"
   "LAMA_BLAS_LIBRARY", "Selected the used BLAS library: auto, MKL, ACML, GOTOBLAS, BLAS"
   "BUILD_TEST", "ON or OFF, explicitly enable/disable building the tests"
   "SCAI_LOG", "INFO recommended, alternatively DEBUG or even more detailed TRACE"
   "USE_CODE_COVERAGE", "ON or OFF, explicitly enable/disable use of Code Coverage"
   "USE_CUDA", "ON or OFF, explicitly enable/disable use of CUDA"
   "USE_MPI", "ON or OFF, explicitly enable/disable use of MPI"
   "USE_SCALAPACK", "ON or OFF, explicitly enable/disable use of SCALAPACK"
   
For an optimized installation of LAMA call

.. code-block:: bash

   cmake -D CMAKE_INSTALL_PREFIX=<path/to/install/dir> -D CMAKE_BUILD_TYPE=Release <project-root>/src

If all packages are installed in the system path, Cmake should find them automatically if they are available. If they are
located at different paths or you have additional packages installed on your system and want to use the alternative, you
have to give Cmake hints where to find the packages. For passing these hints you need to define Cmake-Variables by
-D <VARIABLE_NAME>=<path/to/package/dir>. A summary of possible Build-Variables for the packages are given in the
following table. For advanced help refer to the separate section.

.. csv-table:: CMake hint variables
   :header: "Package", "CMake variable", "Signification"
   
   "ACML (for BLAS and LAPack)", "ACML_ROOT or ACML_INCLUDE_DIR and ACML_LIBRARY_PATH", "path to ACML root dir, or paths to MKL include and library dir"
   "Boost", "BOOST_ROOT", "path to boost root dir"
   "CUDA", "CUDA_TOOLKIT_ROOT", "path to CUDA toolkit"
   "DOXYGEN", "DOXYGEN_EXECUTABLE", "path to doxygen binary dir"
   "FortranBLAS (for BLAS and LAPack)", "", ""
   "MKL (for BLAS, LAPack and ScaLAPack)", "MKL_ROOT or MKL_INCLUDE_DIR and MKL_LIBRARY_PATH", "path to MKL root dir, or paths to MKL include and library dir"
   "MPI", "MPI_ROOT", "path to MPI root dir"
   "Vampir Trace", "VAMPIRTRACE_ROOT" or "VT_ROOT", "path to vampir root dir"
   "experimental", "------------------" "------------------"
   "GASNet", "GASNET_LIBRARY_PATH", "path to GASNet library"
   "Metis", "METIS_ROOT or METIS_INCLUDE_DIR and METIS_LIBRARY_PATH", "path to Metis root dir, or paths to Metis include and library dir"
   "OpenCL", "OPENCL_ROOT or OPENCL_INCLUDE_DIR and OPENCL_LIBRARY_PATH", "path to OpenCL root dir, or paths to OpenCL include and library dir"
   "OpenShmem", "OSHMEM_INCLUDE_PATH and OSHMEM_LIBRARY_PATH", "paths to Open shared memory include and library dir"
   "ParMetis", "PARMETIS_ROOT or PARMETIS_INCLUDE_DIR and PARMETIS_LIBRARY_PATH", "path to ParMetis root dir, or paths to ParMetis include and library di"

For changing the C compiler there is a different way for telling cmake. Please see :doc:`here <configuration/c++>`. 

Your cmake-call may look like this then:

.. code-block:: bash

   cmake -D MKL_ROOT=/opt/intel/mkl -D VAMPIRTRACE_ROOT=/opt/vampir_trace -D CMAKE_INSTALL_PREFIX=<path/to/install/dir> -D CMAKE_BUILD_TYPE=Release <project-root>/src

The mkl library from '/opt/intel/mkl' and vampir library from '/opt/vampir_trace' will be used (instead of the
libraries in the system path). All other libraries will be taken from the system path. To check the configuration you may
call in the build directory:

.. code-block:: bash

   ccmake .

It will show you the advanced LAMA configuration definitions first, after toggeling (t) you can view all settings.


More information about useful CMake variables can be found in the `CMake Wiki`_.

.. _CMake Wiki: http://www.cmake.org/Wiki/CMake_Useful_Variables

A description of the find modules provided by CMake can be found in the `CMake documentation`_. We are using FindBoost,
FindMPI, FindCUDA.

.. _CMake documentation : http://www.cmake.org/cmake/help/v2.8.8/cmake.html

In short:

.. code-block:: bash 

   #Create build directory
   mkdir build
   cd build
   #Call cmake
   cmake [-D MKL_ROOT=/opt/intel/mkl] [-D ... ] -D CMAKE_INSTALL_PREFIX=<path/to/install/dir> [-D CMAKE_BUILD_TYPE=Release] <project-root>/src
