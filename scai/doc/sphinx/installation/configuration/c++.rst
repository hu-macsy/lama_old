C++ Compiler
^^^^^^^^^^^^

LAMA is written in C++, therefore a C++ compiler is mandatory.

On Linux systems, the use of the GNU C++ compiler (version 4.4.x or higher) or the Intel C++ compiler is recommended.
On Window systems, the use of Visual Studio is recommended.

LAMA exploits OpenMP parallelization for the use of multiple cores of one CPU and/or for multiple CPUs sharing a global
address space. Therefore a C++ compiler supporting OpenMP is strongly recommended.

Furthermore, a C++ compiler supporting the C++11 features is helpful. In this case, LAMA can be built without any
installation of the Boost package.

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

You can use directly the name of the compiler executable if it is in your path:

.. code-block:: bash

  CC=icc;CXX=icpc cmake ...
  CC=gcc-4.8;CXX=g++-4.8 cmake ...

Dependencies with Other Packages
--------------------------------

The C++ compiler interacts with nearly all other packages needed for LAMA, either the mandatory or the optional ones.

* The CUDA compiler driver ``nvcc`` uses a general purpose C++ host compiler that should be the same as used for LAMA.
  Usually each version of the CUDA toolkit that contains the ``nvcc`` compiler supports only GNU C++ compilers up to 
  a certain version, or even worse, supports only some versions of the Intel C++ compiler. 
  Furthermore, the C++11 features 
  used in LAMA are only supported by CUDA version 7.0 or higher. Therefore with enabled
  CUDA support, the number of supported C++ compilers is limited.

.. code-block:: bash

       // Error messages when using the NVIDIA (R) Cuda compiler driver release 6.5

       error -- unsupported GNU version! gcc 4.9 and up are not supported!
       error -- unsupported ICC configuration! Only ICC 14.0 on Linux x86_64 is supported!

* If Intel MIC support is enabled, an Intel compiler supporting the offload programming model is required. A GNU C++
  compiler cannot be used.

* Other software packages like Boost, MPI, BLAS, or MKL should be available as developer versions that can be used
  with the chosen C++ compiler. Even if this is usually the case, it happens very often that the header files of these 
  packages give compiler warnings that cause troubles for the compilation of LAMA. In such a case, compiler warnings 
  should be disabled.


