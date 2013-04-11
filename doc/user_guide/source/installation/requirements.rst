Software Prerequisites
======================

Required software
^^^^^^^^^^^^^^^^^
First of all you need to install the following software for building LAMA:

- :ref:`CMake Version 2.8 or later <configuration>`

- :ref:`C/C++ Compiler <c++>`

- :ref:`Boost (version 1.34 or later) <section_boost>`

- :ref:`BLAS <blas>`

- :ref:`LAPack <blas>`

On linux install the corresponding software via your package manager or simply use apt-get.
For Windows follow the helping description :ref:`here <windowsTipps>`.

.. [#f1] Make sure that cmake.exe is in your (system) path, subsequent cmake calls depend on it. You can also have the
         cygwin cmake in the path. It will cooperate with the cmake-gui. 

Running Configurations
^^^^^^^^^^^^^^^^^^^^^^
We have tested the installation of LAMA with various versions of the required software.
Have a look at the following table whether you configuration is supported.

.. csv-table:: 
   :header: "compiler", "boost", "BLAS", "cuda", "communication", "build"
   :widths: 200, 150, 150, 150, 150, 150 

   "gcc 4.4", "1.46", "MKL ??", "4.0", "OpenMPI 1.4.3", "ok"
   "gcc 4.6.1 / 4.6.2", "1.49.0", "blas 3", "4.2 / 5.0", "mpich2 1.2.1p1 (Parastation 5.0.25-2)", "ok"
   "gcc 4.6.1 / 4.6.2", "1.49.0", "blas 3", "4.2", "mvapich2", "ok"
   "", "", "", "", "", ""
   "", "", "", "", "", ""

Recommended software
^^^^^^^^^^^^^^^^^^^^
For using the support of parallel/hybrid machines in LAMA, you need the following:

- OpenMP 2.0 capable :ref:`C++ Compiler <c++>`

  - If no OpenMP Capable C++ Compiler is available shared memory parallelization is not used.

- :ref:`CUDA (version 4.0 or later) <section_cuda>`

  - Needed to utilize CUDA capable GPUs from NVidia.

- :ref:`MPI <section_mpi>`

  - Needed for Distributed Memory Parallelism.

Optional software
^^^^^^^^^^^^^^^^^
For optional features you may also install the following:

- Doxygen

  - Needed to build the API documentation.

- ScaLAPack

  - Needed for Distributed Memory Parallel LU-Decompostion.

- VampirTrace

  - Needed for performance analysis with `Vampir`_.
  
  .. _Vampir: http://http://vampir.eu/
  
Experimental software 
^^^^^^^^^^^^^^^^^^^^^
For experimental features (only available in the svn /trunk not in the release!)
you may also install the following:

- OpenCL
 
  - Needed to utilize OpenCL Accelerators.
   
- OpenShmem
 
  - Needed for Distributed Memory Parallelism with the PGAS Layer OpenShmem.
   
  - OpenShmem needs GASNet
   
- Metis/ParMetis
 
  - Needed to generate optimized Distributions based on SparseMatrix structure.
