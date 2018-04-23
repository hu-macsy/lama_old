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

  - Now mandatory: support of the ISO-C++ standard C++11 (e.g. since gcc 4.7, icc 11.1, clang 2.9)

- :ref:`Boost <scaicommon:Boost>` (version 1.59 or later)

  - Required if unit tests of the projects should be built

.. |OpenMP_API| raw:: html

  <a href="http://www.openmp.org/wp/openmp-compilers/" target="_blank"> OpenMP API </a>

.. toctree::
   :hidden:

   configuration/cmake
   configuration/compiler

.. _recommended:

Recommended Software
^^^^^^^^^^^^^^^^^^^^

The following software packages are not mandatory, but will give you a performance boost through optimized routines from BLAS and accelerator backends (currently LAMA supports CUDA) as well as distributed memory support:

- :ref:`scaiblaskernel:BLAS`

  - Needed for the implementation of the dense blas functions. Supported libraries are: Intel MKL, BLAS. 

- :ref:`CUDA <scaicommon:CUDA>` (version 7.0 or later)

  - Needed to utilize CUDA capable GPUs from Nvidia (support for compute capability cc 2.0 or higher).
  - C++11 support required for the CUDA compiler (since nvcc 7.0)

- :ref:`MPI <scaidmemo:MPI>`

  - Needed for distributed memory parallelism through the Message Passing Interface.

.. toctree::
   :hidden:

.. _optional:

Optional Software
^^^^^^^^^^^^^^^^^

For optional features you may also install the following:

- :ref:`Metis<scaidmemo:Metis>` 
 
  - Needed to generate optimized distributions based on SparseMatrix structure.

- Sphinx

   |Sphinx| is used for this user documentation.

.. |Sphinx| raw:: html

  <a href="http://sphinx-doc.org/" target="_blank">Sphinx</a>

- Doxygen

    |Doxygen| is used for our system documentation. You can also find it |apidoc|.
  
.. |Doxygen| raw:: html

  <a href="http://www.doxygen.org" target="_blank">Doxygen</a>

.. |apidoc| raw:: html

  <a href="https://libama.org/doxygen/index.html" target="_blank">here</a>

.. _supported:

- Java

   The tracing module comes with a GUI written in Java that allows the visualization of
   trace files that can be generated when running LAMA applications. For the compilation of
   this GUI a Java compiler is needed and a **jar** tool to build the corresponding package.


Configurations
^^^^^^^^^^^^^^

We have tested the installation of LAMA with various versions of the required and optional software.

Nevertheless we encountered different problems, and here are some general advices:

- CXX Compiler (support of C++11 standard)

  We have very good experience with all versions of the GNU C++ compiler suite and the Intel C++ compiler
  tools. We now require using compiler versions that support the C++11 features, i.e. GCC 4.8.1 or higher,
  or Intel C++ 15.0 or higher. In this case, Boost is no more mandatory that also avoids a lot of problems.

- Boost (Version 1.61 or higher)

  Most Linux distributions come with a Boost installation that fits well the default GNU C++ compiler. 
  Generally speaking, the Boost library should have been built with the same compiler that is exploited for the
  compilation of C++ code using the boost functionality. 
  It is very likely that you run in serious problems
  when using the Boost library compiled with a later compiler version (e.g. g++ 5.4 ) than the compiler used for the C++ code
  (e.g. g++ 5.3 ) using the boost functionality. In these cases, we strongly recommend to build a Boost installation with 
  the corresponding C++ compiler.

.. code-block:: bash

   tar xvfz boost_1_61_0.tar.gz 
   cd boost_1_61_0
   bootstrap.sh --prefix=$BOOST_ROOT --with-libraries=test
   ./b2

- CUDA (Version 7.0 or higher)

  Unfortunately the CUDA toolkit versions cannot keep up with the latest GNU C++ compiler or Intel compiler version.
  In this case, you have to install/use a previous compiler version. E.g. CUDA Version 8.0 does not support gcc versions
  later than 5.3. Futhermore, C++11 support is only given with CUDA Version 7.0 or higher.
  And when using an older gcc or icc version, you might run in problems when using an installed Boost version.


.. _`contact`: mailto:lama@scai.fraunhofer.de
