.. _requirements:

Software Prerequisites
----------------------

For the installation of LAMA you will need some software packages that are used by LAMA but not part of it. 
In the following sections we list the
:ref:`required <required>`, :ref:`recommended <recommended>` and :ref:`optional <optional>` software packages for building LAMA.

Some software packages used by LAMA might already be available on your system, otherwise they have to be installed.
On Linux it is the prefered way to install the additional software packages via the package manager of your system. 
But also custom-built installations should cause no problems as far as the installation path is known during the
configuration.
You get precise help to explicitly tell LAMA during the configuration the installation path of the individual packages;
this help you will find in the corresponding section for each software package that describe its use in LAMA.

At the end of the page, our experience with different :ref:`configurations <supported>` regarding various versions of packages is summarized.

..	_required:

Required Software
^^^^^^^^^^^^^^^^^

All you need for the LAMA installation are the following software packages:

- :ref:`cmake` Version 2.8.8 or later

  - Needed for creating platform independent makefiles or projects.

- :ref:`Compiler`

  - Supported Compilers are: GNU, Intel, Clang.

  - Optional: implementation of the |OpenMP_API|, version 2.0 or higher

  - Recommended but not mandatory: support of the ISO-C++ standard C++11 (e.g. since gcc 4.7, icc 11.1, clang 2.9)

- :ref:`Boost <scaicommon:Boost>` (version 1.34 or later)

  - Required if the C++11 standard is not supported by the C++ compiler

  - Required if unit tests of the projects should be built

- :ref:`PThread <scaicommon:PThread>` (version 1.34 or later)
  
  - Required for internal thread handling


.. |OpenMP_API| raw:: html

  <a href="http://www.openmp.org/wp/openmp-compilers/" target="_blank"> OpenMP API </a>

.. toctree::
   :hidden:

   configuration/cmake
   configuration/compiler

.. _recommended:

Recommended Software
^^^^^^^^^^^^^^^^^^^^

The following software packages are not mandatory, but will give you a performance boost through optimized routines from BLAS and accelerator backends (CUDA, MIC) as well as distributed memory support:

- :ref:`scaiblaskernel:BLAS`

  - Needed for the implementation of the dense blas functions. Supported libraries are: Intel MKL, BLAS. 

- :ref:`CUDA <scaicommon:CUDA>` (version 4.0 or later)

  - Needed to utilize CUDA capable GPUs from Nvidia (support for compute capability since cc 1.3).
  
- :doc:`Intel MIC <configuration/mic>`

  - Needed to utilize Intel Xeon Phi

- :ref:`MPI <scaidmemo:MPI>`

  - Needed for distributed memory parallelism through the Message Passing Interface.

- :ref:`GPI-2 <scaidmemo:GPI>`

  - Needed for distributed memory parallelism through the Global Address Space Programming Interface.

.. toctree::
   :hidden:

   MIC <configuration/mic>

.. _optional:

Optional Software
^^^^^^^^^^^^^^^^^

For optional features you may also install the following:

- :ref:`Metis/ParMetis <scaidmemo:Metis>` 
 
  - Needed to generate optimized Distributions based on SparseMatrix structure.

- Sphinx

   |Sphinx| is used for this user documentation.

.. |Sphinx| raw:: html

  <a href="http://sphinx-doc.org/" target="_blank">Sphinx</a>

- Doxygen

    |Doxygen| is used for our system documentation. You can also find it |apidoc|.
  
.. |Doxygen| raw:: html

  <a href="http://www.doxygen.org" target="_blank">Doxygen</a>

.. |apidoc| raw:: html

  <a href="https://test.libama.org/doxygen/index.html" target="_blank">here</a>

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

    #if BOOST_GCC >= 40700
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-local-typedefs"
    #endif

    namespace boost
    {

    ...

    } // namespace boost
 
    #if BOOST_GCC >= 40700
    #pragma GCC diagnostic pop
    #endif

:sup:`2` OpenBLAS (develop, commit: 6e7be06e072cddc5d34617f28f60a32484b9e910, date: 02/22/2016)

.. _`nvcc-boost-patch`: https://svn.boost.org/trac/boost/ticket/9392
.. _`boost-patches`: http://gcc.gnu.org/ml/gcc/2013-07/msg00237.html


If you have problems with the installation of supported configuration, do not hesitate to `contact`_ us.
If you have tested not listed configurations we are pleased to get new input for the list.

.. _`contact`: mailto:lama@scai.fraunhofer.de
