Directory Structure
-------------------

LAMA itself is the largest SCAI project coming with the distribution, but there are many others:

- common for some basic stuff
- logging is a logging library developed by SCAI
- tracing is a tracing library developed by SCAI
- tasking is a library supporting thread pools and asynchronous execution
- hmemo is a library for hybrid architecutes where an array might use memory on different devices

The SCAI LAMA distribution contains the following directories

- scai is the directory containing the sources of the different projects
- tools contains some tools for source code formatting, etc.
- CMake contains SCAI-specific extensions of cmake
- res contains some test cases

::

   ./CMake       # SCAI cmake module files
   ./res         # input test files
   ./scai        # different SCAI project
   ./tools       # tools that might be used for code development

Each project is like lama in a separate directory of the scai directory. 

::

   ./scai 
   ./scai/common
   ./scai/hmemo
   ./scai/lama
   ./scai/logging
   ./scai/tasking
   ./scai/tracing

Each of the project directories scai/project contains the hpp and cpp source files and its cmake
configuration file CMakeLists.txt. Generally there are one subdirectory with example programs
and one with test programs. The example programs demonstrate the typical use of the library
provided by the project while the test programs are intended to be used for automated testing 
of the functionality.

::

    ./scai/project
    ./scai/project/CMakeLists.txt
    ./scai/project/class.hpp
    ./scai/project/class.cpp
    ...
    ./scai/project/examples
    ./scai/project/examples/CMakeLists.txt
    ./scai/project/examples/example.cpp
    ./scai/project/test
    ./scai/project/test/CMakeLists.txt
    ./scai/project/test/Test1.cpp
    ./scai/project/test/Test2.cpp

LAMA is a mathematical library for heterogeneous architectures. This implies that code for
different architectures is provided, e.g. CUDA code for NVIDIA Tesla graphic cards, OpenCL code
for OpenCL devices and so on. As this code is intended to be optional, it is always contained in a separate
subdirectory, e.g. cuda, and is only included if the corresponding device should be supported.

::

    ./scai/project
    ...
    ./scai/project/cuda
    ./scai/project/cuda/CMakeLists.txt
    ./scai/project/cuda/CUDAclass.hpp
    ./scai/project/cuda/CUDAclass.cpp
    ...
    ./scai/project/opencl
    ./scai/project/opencl/CMakeLists.txt
    ./scai/project/opencl/OpenCLclass.hpp
    ./scai/project/opencl/OpenCLclass.cpp
    ...

Each of the projects has its own namespace ``scai::<project>``. Further subdirectories with code for
specific devices do not introduce an additional namespace.

