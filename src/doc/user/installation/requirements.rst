Software Prerequisites
======================

For the installation of LAMA you will need some general software packages. In the following sections we list the
:ref:`required <required>`, :ref:`recommended <recommended>`,  :ref:`optional <optional>` and
:ref:`experimental <experimental>` software for building LAMA.
You will get precise help to tell LAMA the installation path of the individual packages in seperate chapters.
At the end we will summarize our experience with different :ref:`configurations <supported>` regarding various versions of packages.

On linux you can install the corresponding software via your package manager or simply use apt-get.

On Windows you need to download the package from the precise page and install them on your machine.

.. 	_required:

Required software
-----------------

First of all you need to install the following software for building LAMA:

- CMake Version 2.8 or later

  - Needed for creating platform independent makefiles or projects.  

- :doc:`C/C++ Compiler <configuration/c++>` with OpenMP 2.0 capable

  - Supported Compilers are: GNU, Intel, PGI. 

- :doc:`Boost (version 1.34 or later) <configuration/boost>`

  - For LAMA we need the headers only packages *thread* for the internal thread management. *unit_test_framework* and *regex* are used for the tests (optional).

- :doc:`BLAS and LAPack <configuration/blas>`

  - Needed for the implementation of the dense blas functions. LAMA only prepares own implementations for sparse blas routines, for dense matrices we rely on given implementations. Supported libraries are: Intel MKL, GOTO, ACML, FortranBLAS. 

.. _recommended:

Recommended software
--------------------

For using the support of parallel/hybrid machines in LAMA, you need the following:

- :doc:`CUDA (version 4.0 or later) <configuration/cuda>`

  - Needed to utilize CUDA capable GPUs from NVidia.

- :doc:`MPI <configuration/mpi>`

  - Needed for distributed memory parallelism.

.. _optional:

Optional software
-----------------

For optional features you may also install the following:

- Doxygen

  - Needed to build the API documentation.

- ScaLAPack

  - Needed for Distributed Memory Parallel LU-Decomposition.

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
Have a look at the following table whether your configuration is explicitly supported.
All tests are configured as release build.

.. csv-table:: tested configuration settings
   :header: "compiler", "boost", "BLAS", "cuda", "communication", "build"
   :widths: 150, 100, 100, 100, 200, 100 

   "gcc", "-----", "-----", "-----", "-----", "-----"
   "gcc 4.4", "1.46", "MKL composerxe-2011.2.137", "4.2", "OpenMPI 1.4.3", "ok"
   "gcc 4.4", "1.46", "ACML", "4.2", "OpenMPI 1.4.3", "ok"
   "gcc 4.6.0", "1.41", "blas 3", "5.0", "OpenMPI 1.4.3", "ok"
   "gcc 4.6.1 / 4.6.2", "1.49.0", "blas 3", "4.2 / 5.0", "mpich2 1.2.1p1 (Parastation 5.0.25-2)", "ok"
   "gcc 4.6.1 / 4.6.2", "1.49.0", "blas 3", "4.2", "mvapich2", "ok"
   "gcc 4.8.2" "1.55.0", "MKL composer-xe-2013.1.117" "-" "OpenMPI 1.7.3" "(1)"
   "icc", "-----", "-----", "-----", "-----", "-----"
   "icc 12.1.0 / 13.0.0", "1.46", "MKL composerxe-2011.2.137 / MKL 10.3.1 / MKL 11.0.0", "5.0", "OpenMPI 1.4.3 / OpenMPI 1.6.1 / IntelMPI 4.0.3.008 / IntelMPI 4.1.0.024", "ok"
   "pgi", "-----", "-----", "-----", "-----", "-----"
   "pgcpp ", "1.46", "ACML", " not supported", "OpenMPI 1.4.3", "ok"

(1) With new gcc-Versions we had problems with boost: for compiling with CUDA support you need a `nvcc-boost-patch`_.
For the combination of gcc 4.8.2 and boost 1.55.0 we needed some more patches (see. `boost-patches`_)
for known issues and needed to add the following to boost/tuple/detail/tuple_basic.hpp for suppression:

.. code-block:: c++
   :emphasize-lines: 1,2,3,4,5,14,15,16,17

    +#if BOOST_GCC >= 40700
    +#pragma GCC diagnostic push
    +#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
    +#endif
    +

    namespace boost
    {

    ...

    } // namespace boost
 
    +#if BOOST_GCC >= 40700
    +#pragma GCC diagnostic pop
    +#endif
    +

.. _`nvvc-boost_patch`: https://svn.boost.org/trac/boost/ticket/9392
.. _`boost-patches`: http://gcc.gnu.org/ml/gcc/2013-07/msg00237.html


If you have problems with the installation of supported configuration, do not hesitate to `contact`_ us.
If you have tested not listed configurations we are pleased to get new input for the list.

.. _`contact`: mailto:lama@scai.fraunhofer.de
