Software Prerequisites
----------------------

For the installation of LAMA you will need some general software packages. In the following sections we list the
:ref:`required <required>`, :ref:`recommended <recommended>` and :ref:`optional <optional>` software for building LAMA.

For installing dependencies please refer to the precise project pages. On linux it is the prefered way to install via the package manager of your system. But also custom build should be no problem. On Windows you need to download the package from the precise page and install it on your machine. You get precise help to explicitly tell LAMA the installation path of the individual packages in seperate chapters if they are not found implicitly.

At the end of the page, our experience with different :ref:`configurations <supported>` regarding various versions of packages is summarized.

..	_required:

Required Software
^^^^^^^^^^^^^^^^^

All you need for the LAMA installation are the following software packages:

- :doc:`CMake <configuration/cmake>` Version 2.8 or later

  - Needed for creating platform independent makefiles or projects.

- :doc:`C++ Compiler <configuration/c++>`

  - Supported Compilers are: GNU, Intel, PGI.

  - optionally: with OpenMP 2.0 capable

  - recommended but not mandatory: capable of C++11 (e.g. since gcc 4.7 or icc 11.1)

- :ref:`Boost <scaicommon:Boost>` (version 1.34 or later)

  - required if C++11 is not supported by the compiler
  - required if unit tests for the projects should be generated

- :ref:`PThread <scaicommon:PThread>` (version 1.34 or later)

  - required if C++11 features are not supported by the C++ compiler

.. _recommended:

Recommended Software
^^^^^^^^^^^^^^^^^^^^

The following software packages are not mandatory, but will give you a performance boost through optimized routines from BLAS and accelerator backends (CUDA, MIC) as well as distributed memory support

- :doc:`BLAS and LAPack <configuration/blas>`

  - Needed for the implementation of the dense blas functions. Supported libraries are: Intel MKL, BLAS. 

- :doc:`CUDA (version 4.0 or later) <configuration/cuda>`

  - Needed to utilize CUDA capable GPUs from Nvidia.
  
- :doc:`Intel MIC <configuration/mic>`

  - Needed to utilize Intel Xeon Phi

- :doc:`MPI <configuration/mpi>`

  - Needed for distributed memory parallelism through the Message Passing Interface.

- :doc:`GPI-2 <configuration/gpi>`

  - Needed for distributed memory parallelism through the Global Address Space Programming Interface.

.. _optional:

Optional Software
^^^^^^^^^^^^^^^^^

For optional features you may also install the following:

- Metis/ParMetis
 
  - Needed to generate optimized Distributions based on SparseMatrix structure.

- Sphinx

  - Needed to build the user documentation (this).

- Doxygen

  - Needed to build the API documentation.
 
.. - ScaLAPack

..   - Needed for Distributed Memory Parallel LU-Decomposition.

.. - VampirTrace

..   - Needed for performance analysis with `Vampir`_.
  
.. .. _Vampir: http://http://vampir.eu/

.. _supported:

Configurations
^^^^^^^^^^^^^^

We have tested the installation of LAMA with various versions of the required software.
Have a look at the following table whether your configuration is explicitly supported.
All tests are configured as release build.

.. csv-table:: tested configuration settings
   :header: "compiler", "boost", "BLAS", "cuda", "communication", "build"
   :widths: 150, 100, 100, 100, 200, 100 

   "**gcc**", "-", "-", "-", "-", "-"
   "gcc 4.4", "1.46", "MKL composerxe-2011.2.137", "4.2", "OpenMPI 1.4.3", "ok"
   "gcc 4.4", "1.46", "ACML", "4.2", "OpenMPI 1.4.3", "ok"
   "gcc 4.6.0", "1.41", "blas 3", "5.0", "OpenMPI 1.4.3", "ok"
   "gcc 4.6.1 / 4.6.2", "1.49.0", "blas 3", "4.2 / 5.0", "mpich2 1.2.1p1 (Parastation 5.0.25-2)", "ok"
   "gcc 4.6.1 / 4.6.2", "1.49.0", "blas 3", "4.2", "mvapich2", "ok"
   "gcc 4.8.2", "1.55.0", "MKL composer-xe-2013.1.117", "-", "OpenMPI 1.7.3", ":sup:`1`"
   "gcc 4.8.5", "1.53.0", "blas 3.4.2", "-", "OpenMPI 1.4.3", "-"
   "gcc 4.8.5", "1.53.0", "OpenBLAS :sup:`2`", "7.5", "-", "ok"
   "gcc 4.9.1", "1.58.0", "MKL composerxe-2015.3.187", "7.0", "OpenMPI 1.10.1", "ok"
   "**icc**", "-", "-", "-", "-", "-"
   "icc 12.1.0 / 13.0.0", "1.46", "MKL composerxe-2011.2.137 / MKL 10.3.1 / MKL 11.0.0", "5.0", "OpenMPI 1.4.3 / OpenMPI 1.6.1 / IntelMPI 4.0.3.008 / IntelMPI 4.1.0.024", "ok"
   "**pgi**", "-", "-", "-", "-", "-"
   "pgcpp ", "1.46", "ACML", " not supported", "OpenMPI 1.4.3", "ok"

:sup:`1` With new gcc-Versions we had problems with boost: for compiling with CUDA support you need a `nvcc-boost-patch`_.
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


:sup:`2` OpenBLAS (develop, commit: 6e7be06e072cddc5d34617f28f60a32484b9e910, date: 02/22/2016)

.. _`nvcc-boost-patch`: https://svn.boost.org/trac/boost/ticket/9392
.. _`boost-patches`: http://gcc.gnu.org/ml/gcc/2013-07/msg00237.html


If you have problems with the installation of supported configuration, do not hesitate to `contact`_ us.
If you have tested not listed configurations we are pleased to get new input for the list.

.. _`contact`: mailto:lama@scai.fraunhofer.de
