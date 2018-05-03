.. _Compiler:

C++ Compiler
^^^^^^^^^^^^

LAMA is written in C++, therefore a C++ compiler is mandatory.

On Linux systems, the use of the |GNU| C++ compiler or the |Intel| C++ compiler is recommended. We also support |LLVM| (especially to support MAC).

.. |GNU| raw:: html

  <a href="https://gcc.gnu.org/" target="_blank">GNU</a>

.. |Intel| raw:: html

  <a href="https://software.intel.com/en-us/c-compilers" target="_blank">Intel</a>

.. |LLVM| raw:: html

  <a href="http://llvm.org/" target="_blank">Clang/LLVM</a>

.. On Window systems, the use of Visual Studio is recommended.

LAMA exploits OpenMP parallelization for the use of multiple cores of one CPU and/or for multiple CPUs sharing a global address space. Therefore a C++ compiler supporting OpenMP is strongly recommended, but not mandatory.

Furthermore, a C++ compiler supporting the C++11 standard is required. Most compilers will use the flag
``-std=c++11`` that enables the corresponding C++11 features.

Selection of Compiler for the Configuration
-------------------------------------------

During the configuration with cmake a valid C++ compiler will be identified on your machine. You can force the use of a
certain C and C++ compiler in the following way

.. code-block:: bash

  export CC=/usr/bin/gcc-4.8
  export CXX=/usr/bin/g++-4.8
  cmake ...

or

.. code-block:: bash

  CC=/opt/intel/bin/icc CXX=/opt/intel/bin/icpc cmake ...

You also can use the name of the compiler executable directly if it is in your path:

.. code-block:: bash

  CC=icc;CXX=icpc cmake ...
  CC=gcc-4.8;CXX=g++-4.8 cmake ...

Alternatively, you can pass the compiler executable path with -D as any other option:

.. code-block:: bash

   cmake -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=/opt/intel/bin/icc ...

Dependencies with Other Packages
--------------------------------

The C++ compiler interacts with nearly all other packages needed for LAMA, either the mandatory or the optional ones.

* The CUDA compiler driver ``nvcc`` uses a general purpose C++ host compiler that should be the same as used for LAMA. Usually each version of the CUDA toolkit that contains the ``nvcc`` compiler supports only GNU C++ compilers up to   a certain version, or even worse, supports only some versions of the Intel C++ compiler. Furthermore, the C++11 features used in LAMA are only supported by CUDA version 7.0 or higher. Therefore with enabled CUDA support, the number of supported C++ compilers is limited.

.. code-block:: bash

       // Error messages when using the NVIDIA (R) Cuda compiler driver release 6.5

       error -- unsupported GNU version! gcc 4.9 and up are not supported!
       error -- unsupported ICC configuration! Only ICC 14.0 on Linux x86_64 is supported!

* Other software packages like Boost, MPI, BLAS, or MKL should be available as developer versions that can be used   with the chosen C++ compiler. Even if this is usually the case, it happens very often that the header files of these   packages give compiler warnings that cause troubles for the compilation of LAMA. In such a case, compiler warnings   should be disabled (run ccmake and delete -Werror in ADDITIONAL_WARNING_FLAGS).

C Compiler
^^^^^^^^^^

CMake also needs to find a C compiler to look for some additional packages (e.g. PThreads, BLAS). It is not used in the further compilation steps, so it need not to be of the same family as the C++ compiler.
