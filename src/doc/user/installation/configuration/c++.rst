C/C++ Compiler
--------------

LAMA is written in C++, therefore a C++ compiler is mandatory.

On Linux systems, the use of the GNU C++ compiler (version 4.4.x or higher) or the Intel C++ compiler is recommended.
On Window systems, the use of Visual Studio is recommended.

LAMA exploits OpenMP parallelization for the use of multiple cores of one CPU and/or for multiple CPUs sharing a global
address space. Therefore a C++ compiler supporting OpenMP is strongly recommended.

During the configuration with cmake a valid C++ compiler will be identified on your machine. You can force the use of a
certain C and C++ compiler in the following way::

  export CC=<path/to/C-compiler>
  export CXX<path/to/C++ compiler>
  cmake ..

or::

  CC=/opt/intel/bin/icc CXX=/opt/intel/bin/icpc cmake ...

You can use directy the name of the compiler executable if it is in your path::

  CC=icc;CXX=icpc cmake ...

Using g++ version 4.3.4 caused compilation errors as some boost Spirit grammars could not be identified
(maybe there is a better workaround than using a newer compiler version)::

   comparisons like ‘X<=Y<=Z’ do not have their mathematical meaning

