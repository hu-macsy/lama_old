Software Prerequisites
======================

For the installation of LAMA you will need some general software packages. In the following sections we list the
:ref:`required <required>`, :ref:`recommended <recommended>`,  :ref:`optional <optional>` and
:ref:`experimental <experimental>` software for building LAMA.
You will get precise help to tell LAMA the installation path of the individual packages in seperate chapters.
At the end we will summarize our experience with different :ref:`configurations <supported>` regarding various versions of packages.

On linux install the corresponding software via your package manager or simply use apt-get.

On Windows you need to download the package from the precise page and install them on you machine.

.. 	_required:

Required software
-----------------

First of all you need to install the following software for building LAMA:

- CMake Version 2.8 or later

  - Needed for creating platform independent makefiles or projects.  

- :ref:`C/C++ Compiler <c++>` with OpenMP 2.0 capable

  - Supported Compilers are: GNU, Intel, PGI. 

- :ref:`Boost (version 1.34 or later) <section_boost>`

  - For LAMA we need the headers only packages *thread* for the internal thread management. *unit_test_framework* und *regex* are used for the tests (optional).

- :ref:`BLAS and LAPack <blas>`

   - Needed for the implementation of the dense blas functions. LAMA only prepares own implementations for sparse blas routines, for dense matrices we rely on given implementations. Supported libraries are: Intel MKL, GOTO, ATLAS, FortranBLAS. 

.. _recommended:

Recommended software
--------------------

For using the support of parallel/hybrid machines in LAMA, you need the following:

- :ref:`CUDA (version 4.0 or later) <section_cuda>`

  - Needed to utilize CUDA capable GPUs from NVidia.

- :ref:`MPI <section_mpi>`

  - Needed for distributed memory parallelism.

.. _optional:

Optional software
-----------------

For optional features you may also install the following:

- Doxygen

  - Needed to build the API documentation.

- ScaLAPack

  - Needed for Distributed Memory Parallel LU-Decompostion.

- VampirTrace

  - Needed for performance analysis with `Vampir`_.
  
  .. _Vampir: http://http://vampir.eu/

.. _experimental:

Experimental software
---------------------

For experimental features (only available in the corresponding git branches not in the release!) you may also install
the following:

- OpenCL (feature/opencl)
 
  - Needed to utilize OpenCL Accelerators, eg. AMD GPU's.
   
- OpenShmem (feature/pgas)
 
  - Needed for distributed memory parallelism with the PGAS Layer OpenShmem.
   
  - OpenShmem needs GASNet
   
- Metis/ParMetis (feature/loadbalancing)
 
  - Needed to generate optimized Distributions based on SparseMatrix structure.

.. _supported:

Configurations
--------------
We have tested the installation of LAMA with various versions of the required software.
Have a look at the following table whether you configuration is explicitly supported.
All tests are configured as release build.

.. csv-table:: 
   :header: "compiler", "boost", "BLAS", "cuda", "communication", "build"
   :widths: 150, 100, 100, 100, 200, 100                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               

   "gcc 4.4", "1.46", "MKL composerxe-2011.2.137", "4.2", "OpenMPI 1.4.3", "ok"
   "gcc 4.4", "1.46", "ACML", "4.2", "OpenMPI 1.4.3", "ok"
   "gcc 4.6.0", "1.41", "", "5.0", "OpenMPI 1.4.3", "ok"
   "gcc 4.6.1 / 4.6.2", "1.49.0", "blas 3", "4.2 / 5.0", "mpich2 1.2.1p1 (Parastation 5.0.25-2)", "ok"
   "gcc 4.6.1 / 4.6.2", "1.49.0", "blas 3", "4.2", "mvapich2", "ok"
   "-----", "-----", "-----", "-----", "-----", "-----", 
   "icc 12.1.3/13.0.0", "1.46", "MKL composerxe-2011.2.137", " not supported", "OpenMPI 1.4.3", "ok"
   "-----", "-----", "-----", "-----", "-----", "-----"
   "pgcpp ", "1.46", "ACML", " not supported", "OpenMPI 1.4.3", "ok"
   "", "", "", "", "", ""
 
If you have problems with the installation of supported configuration, do not hesitate to `contact`_ us.
If you have tested not listed configurations we are pleased to get new input for the list.

.. _`contact`: mailto:lama@scai.fraunhofer.de