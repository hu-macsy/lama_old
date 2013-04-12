.. _windowsTipps:

Windows Requirement Installation Tipps
======================================

.. [#f1] Make sure that cmake.exe is in your (system) path, subsequent cmake calls depend on it. You can also have the
         cygwin cmake in the path. It will cooperate with the cmake-gui. 

Boost
-----

Download, compile and install Boost. During LAMA's CMake configure specify BOOST_ROOT (-DBOOST_ROOT) pointing to the
installation path, e.g.

.. code-block:: bash 

   BOOST_ROOT=C:\boost\boost_1_49_0

How to build boost for Windows in general can be find `here`__ with the following modifications

__ http://stackoverflow.com/questions/2629421/how-to-use-boost-in-visual-studio-2010

Choose the target install instead of stage, so that the headers are copied to include (use --prefix to specify the path).
    
LAMA only depends on the boost libraries thread, program-options, and test 

Hence, after issuing **.\bootstrap** the subsequent command for bjam is for a 64 Bit architecture:

.. code-block:: bash 

   .\b2.exe --abbreviate-paths --toolset=msvc-10.0 --with=thread --with=test architecture=x86 address-model=64

The option **--abbreviate-paths** needs to be used, if your paths in the shell get very long. Building from a location
with a short path like C:\boost\boost_1_49_0 also helps preventing path name errors during the build MPI

Download and install HPC Pack 2008 R2 MS-MPI Redistributable Package. CMake should find MS-MPI automatically. If you
have installed the HPC Pack to a none default path specify MPI_ROOT during CMake configure, e.g.

.. code-block:: bash 
	
   MPI_ROOT=C:/Program Files/Microsoft HPC Pack 2008 R2

If you want to use a different MPI follow the documentation of the CMake Module FindMPI.


BLAS/LAPack
-----------

Alternative 1 Intel MKL
^^^^^^^^^^^^^^^^^^^^^^^

Download and install Intel MKL. During CMake configure specify MKL_ROOT pointing to the installation path, e.g.

.. code-block:: bash 

   MKL_ROOT=C:/Program Files (x86)/Intel/ComposerXE-2011/mkl

Sometimes the Architecture detection does not work correctly in our FindMKL CMake Module if your MKL installation is not
found specify MKL_Is64=true during CMake configure.

Alternative 2 ACML
^^^^^^^^^^^^^^^^^^

Download and install ACML . During CMake configure specify ACML_ROOT pointing to the installation path, e.g.

.. code-block:: bash 

   ACML_ROOT=C:/acml/acml-4.4.0/ifort64_mp

CUDA
----

Download and install the CUDA Toolkit from NVidia CUDA Downloads. If the Toolkit is installed at its default location
CMake should find it.

OpenCL
------

In the case of NVIDIA cards, OpenCL is installed together with CUDA. Here, you may want to disable CUDA by commenting
out the OpenCL line in the main CMakeLists.txt file:

.. code-block:: bash

   # Find optional packages
   ...
   set ( OPTIONAL_PACKAGES_TO_FIND
         OpenMP
         CppUnit
         Doxygen
         Threads
         CUDA
         #CUDA needs to be found before OpenCL because CUDA_INCLUDE_DIRS is used as a hint to find the OpenCL Headers
         OpenCL # comment this line out, if OpenCL is not wanted
         CBE
         DirectX
         OSHMEM REQUIRES GASNET
         GPI
         #add optional packages here 
       )
     ...

Otherwise, download and install the OpenCL SDK for your platform. You may need to set **OPENCL_LIBRARIES** and
**OPENCL_INCLUDE_DIRS**, e.g.

.. code-block:: bash 

   OPENCL_LIBRARIES="C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v3.2/lib/x64/OpenCL.lib"
   OPENCL_INCLUDE_DIRS="C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v3.2/include"
 