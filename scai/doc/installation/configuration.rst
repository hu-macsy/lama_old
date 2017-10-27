.. _configuration:

Configuration
-------------

The configuration of LAMA tries to find required and optional software needed by LAMA and to identify correct compilation and link flags. At the end of the configuration, a **Makefile** (default) is generated that can be used to build libraries and executables of LAMA. For its configuration LAMA uses the CMake utility. This section will guide you through the CMake configuration.

.. or **Visual Studio project**

.. The following steps are for generating Makefiles under Linux. For the creation of a Visual Studio project under Windows see :doc:`here <windowsTipps>`.

The CMake utility should be called in a seperate build directory. Although in-source builds are supported we recommend an out-of-source build - as CMake does as well. So after extracting your tar ball and changing into the new directory create a build directory and change to it.

.. code-block:: bash

   tar -xzvf libama-x.x.x.tar.gz
   cd libama-x.x.x
   mkdir build
   cd build

To configure LAMA, call ``cmake`` in your build directory and specify the path to the LAMA source directory (should be ``../scai``). If you do not install as root user you also have to set the CMAKE_INSTALL_PREFIX, because it is /usr/local by default. Due to the use of ExternalProjects for LAMAs sublibraries, installation will automatically take place after compilation.

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

For an optimized installation of LAMA define ``CMAKE_BUILD_TYPE=Release``. In our official releases (downloaded tar balls) this is the default, so you do not need to set this explicitly. If you have a working copy of git repository ``Debug`` is the default.

.. code-block:: bash

   cmake -DCMAKE_INSTALL_PREFIX=<path/to/install/dir> -DCMAKE_BUILD_TYPE=Release ../scai

Build Specific Options
""""""""""""""""""""""

All relevant LAMA options and their significations are listed in the table below. 

=========================      ======================================================================================
CMake variable                 Signification
=========================      ======================================================================================
``CMAKE_BUILD_TYPE``           Debug (build with debug flag) or Release (build with optimization flags)
``CMAKE_INSTALL_PREFIX``       Installation path of LAMA library
``SCAI_ASSERT_LEVEL``          ERROR (recommended) or DEBUG (adds much more assertions into the code)
``SCAI_BLAS_LIBRARY``          Select a BLAS library: auto, MKL, BLAS, INTERNALBLAS
``SCAI_CMAKE_VERBOSE``         Set to TRUE for more detailed informations about the configuration
``SCAI_HOST_TYPES_LIST``       Set data types for instatitions
                               (default: float;double;long double;ComplexFloat;ComplexDouble;ComplexLongDouble (all))
``SCAI_INDEX_TYPE``            Set data type used as IndexType, valid for all libraries (default: int)
``SCAI_LIBRARY_TYPE``          SHARED (default) or STATIC
``SCAI_LOG_LEVEL``             INFO recommended, alternatively DEBUG or even more detailed TRACE
``SCAI_TRACE``                 ON or OFF, enable/disable region tracing
``USE_CODE_COVERAGE``          ON or OFF, explicitly enable/disable use of Code Coverage
``USE_CUDA``                   ON or OFF, explicitly enable/disable use of CUDA
``USE_METIS``                  ON or OFF, explicitly enable/disable use redistributions with Metis
``USE_JAVA``                   ON or OFF, explicitly enable/disable use of Java for tracing GUI
``USE_MPI``                    ON or OFF, explicitly enable/disable use of MPI
``USE_OPENMP``                 ON or OFF, explicitly enable/disable use of OpenMP
=========================      ======================================================================================

.. USE_SCALAPACK            "ON or OFF, explicitly enable/disable use of SCALAPACK"
.. SCAI_ADDITIONAL_LINK_LIBRARIES "Paths to libraries that should be linked to lama, separated by ;"
.. SCAI_ADDITIONAL_LINK_FLAGS "Additional linker Flags"

Packages
""""""""

If all packages are installed in the system path, CMake should find them automatically - if available - and there is less to do. If they are installed somewhere else or you have more than one option of a package, you need to give CMake some hints where to find the packages or which one to choose. For a specific package please refer to the belonging subpage. Especially for changing the C++ compiler please see the subpage, because there is a different way for telling CMake.

Required dependencies:

* :ref:`Compiler`
* :ref:`Boost <scaicommon:Boost>` for the Common project and all unit tests 
* :ref:`scaicommon:PThread`

Optional dependencies:

* :ref:`BLAS`
* :ref:`CUDA <scaicommon:CUDA>`
* :ref:`MPI <scaidmemo:MPI>`
* :ref:`(Par)Metis <scaidmemo:Metis>`

A summary of possible hint-variables for the packages are given in the following table. For advanced help refer to the separate section.

========== ================================================================  ===========================================================================
Package    CMake variable                                                    Signification
========== ================================================================  ===========================================================================
Boost      BOOST_ROOT                                                        path to the boost root dir
MKL        MKL_ROOT or MKL_INCLUDE_DIR and MKL_LIBRARY_PATH                  path to the mkl root dir or include and library dir
BLAS       BLAS_blas_LIBRARY and LAPACK_lapack_LIBRARY                       path to the blas and lapack library
CUDA       CUDA_TOOLKIT_ROOT                                                 path to CUDA toolkit
DOXYGEN    DOXYGEN_EXECUTABLE                                                path to doxygen binary dir
MPI        MPI_ROOT                                                          path to MPI root dir
JAVA       JAVA_HOME                                                         path to JAVA installation directory, must be set as environment variable
Metis      METIS_ROOT or METIS_INCLUDE_DIR and METIS_LIBRARY_PATH            path to Metis root dir or paths to Metis include and library dir
ParMetis   PARMETIS_ROOT or PARMETIS_INCLUDE_DIR and PARMETIS_LIBRARY_PATH   path to ParMetis root dir or paths to ParMetis include and library dir
========== ================================================================  ===========================================================================

Your cmake-call may look like this then:

.. code-block:: bash

   export JAVA_HOME=/home/lama/development/jdk/jdk1.6.0_35
   cmake -DMKL_ROOT=/opt/intel/mkl -DCMAKE_INSTALL_PREFIX=<path/to/install/dir> -DCMAKE_BUILD_TYPE=Release ../src

The mkl library from '/opt/intel/mkl' will be used (instead of the libraries in the system path). All other libraries will be taken from the system path. 

Summary
"""""""

If no fatal cmake errors occured during the configuration in the end we print a full summary of the configuration. It should look slightly the same as this and give you an overview of all used libraries and all set options:

.. code-block:: none

   -- ==============================
   -- Summary of SCAI Configuration:
   -- ==============================
   -- 
   --  External Software/Packages:
   -- 
   --      CXX ENABLED
   --          FOUND      : Version GNU 4.8.5
   --          Executable : /usr/bin/g++-4.8
   --  
   --      C++11 support ENABLED
   --          FOUND      : Version via C++ compiler
   --  
   --      Boost Unit Test ENABLED
   --          FOUND      : Version 1.58.0
   --          Include(1) : /usr/include
   --          Libs(1)    : /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.so
   --  
   --      Sphinx ENABLED
   --          FOUND      : Version 1.3.6
   --  
   --      OpenMP ENABLED
   --          FOUND      : Version 3.1
   --          CXX Flags  : -fopenmp
   --  
   --      CUDA ENABLED
   --          FOUND      : Version 8.0
   --          CXX Flags  : optimize for compute capability 52
   --          Executable : /usr/local/cuda/bin/nvcc
   --          Include(1) : /usr/local/cuda/include
   --          Libs(5)    : /usr/lib/x86_64-linux-gnu/libcuda.so, ...
   --  
   --      Thread ENABLED
   --          FOUND      : Version 2.23
   --          Libs(1)    : -pthread
   --  
   --      Java ENABLED
   --          FOUND      : Version 1.8.0_91
   --          Executable : /usr/bin/javac
   --  
   --      SCAI BLAS ENABLED
   --          FOUND      : Version BLAS 3 Lapack 
   --          Libs(2)    : /usr/lib/libblas.so, /usr/lib/liblapack.so
   --  
   --      MPI ENABLED
   --          FOUND      : Version 1.10.2
   --          Include(4) : .., /usr/lib/openmpi/include, /usr/lib/openmpi/include/openmpi
   --          Libs(2)    : /usr/lib/openmpi/lib/libmpi_cxx.so, /usr/lib/openmpi/lib/libmpi.so
   --  
   --      Metis DISABLED
   --          NOT FOUND
   --  
   --      ZLIB ENABLED
   --          FOUND      : Version zlib 1.2.8
   --          Include(1) : /usr/include
   --          Libs(1)    : /usr/lib/x86_64-linux-gnu/libz.so
   --  
   --      PNG (image library) ENABLED
   --          FOUND      : Version 1.2.54
   --          Include(2) : /usr/include, /usr/include
   --          Libs(2)    : /usr/lib/x86_64-linux-gnu/libpng.so, /usr/lib/x86_64-linux-gnu/libz.so
   --  
   --      Doyxgen ENABLED
   --          FOUND      : Version 1.8.11
   --          Executable : /usr/bin/doxygen
   --  
   --  SCAI lama_all Version 2.2.0
   -- 
   --  Set SCAI modules    : lama
   --  Used SCAI modules   : common, logging, tracing, .., dmemo, lama
   --  Unused SCAI modules : solver
   -- 
   --  Build Type          : Debug
   --  Library Type        : SHARED
   --  Numeric Types       : float, double, long double, ComplexFloat, ComplexDouble, ComplexLongDouble
   --  IndexType           : int
   --  ASSERT Level        : DEBUG ( -DSCAI_ASSERT_LEVEL_DEBUG )
   --  LOG Level           : DEBUG ( -DSCAI_LOG_LEVEL_DEBUG )
   --  TRACING             : ON ( -DSCAI_TRACE_ON )
   --  CODE COVERAGE       : OFF
   -- 
   -- Configuring done
   -- Generating done
   -- Build files have been written to: [your_build_dir]

External software packages that are not found might reduce the functionality. Building tests or
documentation will be disabled during the build process if the corresponding functionality
is not available. Other missing packages reduces the runtime functionilty of LAMA.

Change
""""""

You may also change some of the configuration by calling ``ccmake`` in the build directory:

.. code-block:: bash

   ccmake .

It will show you the advanced LAMA configuration definitions first, after toggeling (t) you can view all settings.

Alternatively you can do your changes working with the cmake curses gui or by a second cmake-call. 

NOTE: For changing the C++ compiler start with an empty (!) build directory, otherwise CMake would skip essential tests.

Additional Information
^^^^^^^^^^^^^^^^^^^^^^

More information about useful CMake variables can be found in the |CMake_Wiki|.

.. |CMake_Wiki| raw:: html

  <a href="http://www.cmake.org/Wiki/CMake_Useful_Variables" target="_blank">CMake Wiki</a>

