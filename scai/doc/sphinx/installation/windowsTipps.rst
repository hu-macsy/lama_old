:orphan:

Windows Installation
--------------------

CMake also provides the creation of Visual Studio projects, so you can build LAMA in Visual Studio after configuraition.
For this, make sure that cmake.exe is in your (system) path, subsequent cmake calls depend on it. You can also have the
cygwin cmake in the path. It will cooperate with the cmake-gui.

In general the configuration of LAMA under windows works as it does under linux by defining specific CMake variables.
The sections below will give you special tipps for the installation and usage of required and optional packages under
windows and will advise you of possible installion paths, because windows does not come with most of the prerequisites.
The variable can be readily defined in cmake-gui.   

Boost
^^^^^

Download, compile and install Boost. During LAMA's CMake configure specify BOOST_ROOT (-DBOOST_ROOT) pointing to the
installation path, e.g.

.. code-block:: bash 

   BOOST_ROOT=C:\boost\boost_1_49_0

How to build boost for Windows in general can be find `here`__ with the following modifications:
Choose the target install instead of stage, so that the headers are copied to include (use --prefix to specify the path).

__ http://stackoverflow.com/questions/2629421/how-to-use-boost-in-visual-studio-2010
 
LAMA only depends on the boost libraries thread, program-options (optional), and test (optional). So it
is not required to install the full library, which may cause some time. 

Hence, after issuing **.\bootstrap** the subsequent command for bjam is for a 64 Bit architecture:

.. code-block:: bash 

   .\b2.exe --abbreviate-paths --toolset=msvc-10.0 --with=thread --with=test architecture=x86 address-model=64

The option **--abbreviate-paths** needs to be used, if your paths in the shell get very long. Building from a location
with a short path like C:\boost\boost_1_49_0 also helps preventing path name errors during the build MPI.

Download and install HPC Pack 2008 R2 MS-MPI Redistributable Package. CMake should find MS-MPI automatically. If you
have installed the HPC Pack to a none default path specify MPI_ROOT during CMake configure, e.g.

.. code-block:: bash 
	
   MPI_ROOT=C:/Program Files/Microsoft HPC Pack 2008 R2

If you want to use a different MPI follow the documentation of the CMake Module FindMPI.


BLAS/LAPack
^^^^^^^^^^^

Alternative 1 Intel MKL
"""""""""""""""""""""""

Download and install Intel MKL. During CMake configure specify MKL_ROOT pointing to the installation path, e.g.

.. code-block:: bash 

   MKL_ROOT=C:/Program Files (x86)/Intel/ComposerXE-2011/mkl

Sometimes the architecture detection does not work correctly in our FindMKL CMake Module. If your MKL installation is not
found specify MKL_Is64=true during CMake configure.

Alternative 2 ACML
""""""""""""""""""

Download and install ACML . During CMake configure specify ACML_ROOT pointing to the installation path, e.g.

.. code-block:: bash 

   ACML_ROOT=C:/acml/acml-4.4.0/ifort64_mp

CUDA
^^^^

Download and install the CUDA Toolkit from NVidia CUDA Downloads. If the Toolkit is installed at its default location
CMake should find it.

OpenCL
^^^^^^

(only in branch feature/opencl)

In the case of NVIDIA cards, OpenCL is installed together with CUDA. Here, you may want to disable OpenCL by switching
it off by USE_OPENCL=OFF.

Otherwise, download and install the OpenCL SDK for your platform. You may need to set **OPENCL_LIBRARIES** and
**OPENCL_INCLUDE_DIRS**, e.g.

.. code-block:: bash 

   OPENCL_LIBRARIES="C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v3.2/lib/x64/OpenCL.lib"
   OPENCL_INCLUDE_DIRS="C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v3.2/include"
 