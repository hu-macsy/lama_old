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
