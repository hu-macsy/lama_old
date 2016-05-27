Configuration
-------------

The configuration of LAMA tries to find required and optional software needed by LAMA and to identify correct
compilation and link flags. At the end of the configuration, a ``Makefile`` or ``Visual Studio project`` is generated
that can be used to build libraries and executables of LAMA. For its configuration LAMA uses the CMake utility.

This section will guide you through the CMake configuration. The following steps are for generating Makefiles under Linux.
For the creation of a Visual Studio project under Windows see :doc:`here <windowsTipps>`.

This CMake utility should be called in a seperate build directory. Although in-source builds are supported we recommend an out-of-source build. So after extracting your tar ball and changing into the new directory create a build directory and change to it.

.. code-block:: bash

   tar -xzvf libama-x.x.x.tar.gz
   cd libama-x.x.x
   mkdir build
   cd build

To configure LAMA, call cmake in your build directory and specify the path to the LAMA source directory (**../scai**).
If you do not install as root you have to set the CMAKE_INSTALL_PREFIX, because it is /usr/local by default.
Due to the use of ExternalProjects for LAMAs sublibraries installation will automatically take place after compilation.
There are a couple of advanced :ref:`options <options>` for defining special build-types or to switch on/off particular features, but the general cmake-call looks like this:

.. code-block:: bash

   cmake -DCMAKE_INSTALL_PREFIX=<path/to/install/dir> [options] ../scai

.. _options:

Options
^^^^^^^

For defining the options you need to write ''-D<VARIABLE_NAME>=value''. A whitespace between -D and the VARIABLE_NAME is optional, but a call without a whitespace is our recommendation with CMake.

CMake
"""""

If you have not installed CMake in the system path you need to specify to CMake include and library path by:

.. code-block:: bash

   cmake -DCMAKE_INSTALL_PREFIX=<path/to/install/dir>
         [-DCMAKE_INCLUDE_PATH=<install/dir/cmake>/include \]
         [-DCMAKE_LIBRARY_PATH=<install/dir/cmake>/lib \]
         ../scai

Build_Type
""""""""""

For an optimized installation of LAMA call define CMAKE_BUILD_TYPE=Release. Otherwise you build a debug version of LAMA.

.. code-block:: bash

   cmake -DCMAKE_INSTALL_PREFIX=<path/to/install/dir> -DCMAKE_BUILD_TYPE=Release ../scai

Build Specific Options
""""""""""""""""""""""

All relevant LAMA options and their significations are listed in the table below.

.. csv-table:: 
   :header: "CMake variable", "Signification"

   "CMAKE_BUILD_TYPE", "Debug (default), or Release (for an optimized build)"
   "CMAKE_INSTALL_PREFIX", "Installation path of LAMA library"
   "BUILD_DOC", "ON or OFF, explicitly enable/disable building the user doc"
   "BUILD_TEST", "ON or OFF, explicitly enable/disable building the tests"
   "SCAI_ASSERT_LEVEL", "ERROR recommended, DEBUG adds much more assertions into the code"
   "SCAI_BLAS_LIBRARY", "Selected the used BLAS library: auto, MKL, BLAS, INTERNALBLAS"
   "SCAI_CMAKE_VERBOSE", "set to 1 or TRUE for more detailed infos about configuration"
   "SCAI_LIBRARY_TYPE", "SHARED(default) or STATIC"
   "SCAI_LOGGING_LEVEL", "INFO recommended, alternatively DEBUG or even more detailed TRACE"
   "SCAI_TRACING", "ON or OFF, enable/disable region tracing"
   "USE_CODE_COVERAGE", "ON or OFF, explicitly enable/disable use of Code Coverage"
   "USE_CUDA", "ON or OFF, explicitly enable/disable use of CUDA"
   "USE_GPI", "ON or OFF, explicitly enable/disable use of GPI"
   "USE_GRAPHPARTITIONING", "ON or OFF, explicitly enable/disable use of graph partitioning with Metis"
   "USE_MIC", "ON or OFF, explicitly enable/disable use of MIC"
   "USE_MPI", "ON or OFF, explicitly enable/disable use of MPI"
   "USE_OPENMP", "ON or OFF, explicitly enable/disable use of OpenMP"
..   "USE_SCALAPACK", "ON or OFF, explicitly enable/disable use of SCALAPACK"
..   "SCAI_ADDITIONAL_LINK_LIBRARIES", "Paths to libraries that should be linked to lama, separated by ;"
..   "SCAI_ADDITIONAL_LINK_FLAGS", "Additional linker Flags"

Packages
""""""""

If all packages are installed in the system path, CMake should find them automatically - if available - and there is less to do. If they are installed somewhere else or you have more than one option of a package, you need to give CMake some hints where to find the packages or which one to choose. For a specific package please refer to the belonging subpage.
Especially for changing the C++ compiler please see the subpage, because there is a different way for telling CMake.

.. toctree::
   :titlesonly:

   configuration/cmake
   configuration/c++
   configuration/boost
   configuration/blas
   CUDA <configuration/cuda>
   MIC <configuration/mic>

External:

   * :ref:`MPI <scaidmemo:MPI>`
   * :ref:`GPI-2 <scaidmemo:GPI>`
   * :ref:`(Par)Metis <scaidmemo:Metis>`

A summary of possible hint-variables for the packages are given in the following table. For advanced help refer to the separate section.

.. csv-table:: CMake hint variables
   :header: "Package", "CMake variable", "Signification"
   
   "Boost", "BOOST_ROOT", "path to the boost root dir"
   "MKL", "MKL_ROOT or MKL_INCLUDE_DIR and MKL_LIBRARY_PATH", "path to the mkl root dir or include and library dir"
   "BLAS", "BLAS_blas_LIBRARY and LAPACK_lapack_LIBRARY", "path to the blas and lapack library"
   "CUDA", "CUDA_TOOLKIT_ROOT", "path to CUDA toolkit"
   "MIC", "-", "-"
   "DOXYGEN", "DOXYGEN_EXECUTABLE", "path to doxygen binary dir"
   "MPI", "MPI_ROOT", "path to MPI root dir"
   "GPI-2", "GPI2_ROOT", "path to GPI root dir"
   "Metis", "METIS_ROOT or METIS_INCLUDE_DIR and METIS_LIBRARY_PATH", "path to Metis root dir or paths to Metis include and library dir"
   "ParMetis", "PARMETIS_ROOT or PARMETIS_INCLUDE_DIR and PARMETIS_LIBRARY_PATH", "path to ParMetis root dir or paths to ParMetis include and library dir"
..   "experimental", ,
..   "OpenCL", "OPENCL_ROOT or OPENCL_INCLUDE_DIR and OPENCL_LIBRARY_PATH", "path to OpenCL root dir or paths to OpenCL include and library dir"

Your cmake-call may look like this then:

.. code-block:: bash

   cmake -DMKL_ROOT=/opt/intel/mkl -DCMAKE_INSTALL_PREFIX=<path/to/install/dir> -DCMAKE_BUILD_TYPE=Release ../src

The mkl library from '/opt/intel/mkl' will be used (instead of the libraries in the system path). All other libraries will be taken from the system path. 

Check
"""""

To check the configuration you may call in the build directory:

.. code-block:: bash

   ccmake .

It will show you the advanced LAMA configuration definitions first, after toggeling (t) you can view all settings.

Additional Information
^^^^^^^^^^^^^^^^^^^^^^

More information about useful CMake variables can be found in the `CMake Wiki`_.

.. _CMake Wiki: http://www.cmake.org/Wiki/CMake_Useful_Variables
