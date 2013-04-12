.. _configuration:

Configuration
-------------

.. toctree::
   :titlesonly:

   c++
   boost
   blas
   cuda
   mpi
   tracing

The configuration of LAMA tries to find required and optional software 
needed by LAMA and to identify correct compilation ang link flags.
At the end of the configuration, a ``Makefile`` is generated
that can be used to build libraries and executables of LAMA.

LAMA uses for its configuration the CMake utility. This utility
is called in the build directory.
Although in-source builds are supported we recommend an out-of-source build.
So create a build directory and change to it.

.. code-block:: bash 

   mkdir build
   cd build
   cmake [options] <project-root>/src

To configure LAMA, call cmake in your build directory and specify your
installation path with **CMAKE_INSTALL_PREFIX** and the path to the LAMA source
(**<project-root>/src**). If you want an optimized build set
**CMAKE_BUILD_TYPE=Release**. If **CMAKE_BUILD_TYPE** is not set the build
system is configured to do a debug build.

.. code-block:: bash

   cmake -DCMAKE_INSTALL_PREFIX=<path/to/install/dir> [-DCMAKE_BUILD_TYPE=Release] <project-root>/src

CMake selects the compiler to use automatically if you want to override this
behavior you can set the compiler with

.. code-block:: bash 

   CC=icc CXX=icpc cmake [...]

More information about useful CMake variables can be found in the `CMake Wiki`_.

.. _CMake Wiki: http://www.cmake.org/Wiki/CMake_Useful_Variables

CMake should find required and optional dependency automatically if they are
available. If CMake does not find a available dependency you can pass hints to
the CMake call:

.. code-block:: bash

   cmake -DVAMPIRTRACE_ROOT=/opt/vampir_trace -DMKL_ROOT=/opt/intel/mkl -DCMAKE_BUILD_TYPE=Release <project-root>/src

A description of the CMake provided find modules can be found in the `CMake
documentation`_.

.. _CMake documentation : http://www.cmake.org/cmake/help/cmake-2-8-docs.html

We are using FindBoost, FindMPI, FindCUDA.

.. csv-table:: 
   :header: "CMake variable", "Signification"
   :widths: 450, 560

   "LAMA_DEBUG_CMAKE", "set to 1 or TRUE for more detailed infos about configuration"
   "BOOST_ROOT", "Installation directory of Boost"
   "CMAKE_BUILD_TYPE", "Release, or Debug"
   "CUDA_TOOLKIT_ROOT", "optional, directory of the CUDA installation"
   "LAMA_ASSERT_LEVEL", "ERROR recommended, DEBUG adds much more assertions into the code"
   "LAMA_LOG_LEVEL", "INFO recomemended, alternatively DEBUG or more detailled TRACE"
   "LAMA_USE_CUDA", "ON or OFF, explicitly enable/disable use of CUDA"
   "LAMA_USE_MPI", "ON or OFF, explicitly enable/disable use of MPI"

For the CMake Modules we are providing with LAMA the following variables can be
used to pass hints:

.. csv-table:: 
   :header: "Requirement", "CMake variable"
   :widths: 450, 560

   "MKL (for BLAS, LAPack and ScaLAPack)", "MKL_ROOT or MKL_INCLUDE_DIR and MKL_LIBRARY_PATH"
   "ACML( for BLAS and LAPack)", "ACML_ROOT or ACML_INCLUDE_DIR and ACML_LIBRARY_PATH"
   "MPI", "MPI_ROOT"
   "Metis", "METIS_ROOT or METIS_INCLUDE_DIR and METIS_LIBRARY_PATH"
   "ParMetis", "PARMETIS_ROOT or PARMETIS_INCLUDE_DIR and PARMETIS_LIBRARY_PATH"
   "OpenCL", "OPENCL_ROOT or OPENCL_INCLUDE_DIR and OPENCL_LIBRARY_PATH"
   "Vampir Trace", "VAMPIRTRACE_ROOT" or "VT_ROOT"
   "GASNet", "GASNET_LIBRARY_PATH"
   "OpenShmem", "OSHMEM_INCLUDE_PATH and OSHMEM_LIBRARY_PATH"

If you have checked out from svn, you might have problems with different svn-version. You can skip this check with

.. code-block:: bash

   -DSKIP_SVN=true

In short:

.. code-block:: bash 

   #Create build directory
   mkdir build
   cd build
   #Call cmake
   cmake -DCMAKE_INSTALL_PREFIX=<path/to/install/dir> [ -DCMAKE_BUILD_TYPE=Release -DMKL_ROOT=/opt/intel/mkl -D... ] <project-root>/src

