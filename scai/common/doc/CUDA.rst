.. _CUDA:

Using CUDA in LAMA
==================

CUDA is needed to utilize CUDA capable GPUs from NVIDIA. The CUDA toolkit can be downloaded |CUDA_Download|.

.. |CUDA_Download| raw:: html

  <a href="https://developer.nvidia.com/cuda-downloads" target="_blank">here</a>

CUDA is optional and LAMA can be built without it. But you will not be able to run any code on 
the GPU.

Use of CUDA in LAMA
-------------------

* Some CUDA classes in the common project
* As a context for heterogeneus memory
* For the implemenation of kernels on the CUDA context

For the implementation of LAMA kernel routines on a GPU, any of the following possibilities can be exploited:

- Hand-written CUDA kernels,

- Using Thrust, a C++ template library for CUDA based on the Standard Template Library (STL). 
  LAMA uses it to implement some operations on matrices and vectors with minimal programming effort
  through a high-level interface that is fully interoperable with CUDA C.
  Since Thrust is a template library of header files, no further installation is necessary for using Thrust.

- Using routines of the NVIDIA CUDA Basic Linear Algebra Subroutines (|cuBLAS|) library as a 
  a GPU-accelerated version of the complete standard BLAS library like the MKL BLAS for the CPU.

- Using routines of the NVIDIA CUDA Sparse Matrix library (|cuSPARSE|) as a collection of basic linear algebra subroutines 
  used for sparse matrices with the same interfaces as the MKL for the CPU.

.. |cuSPARSE| raw:: html

  <a href="https://developer.nvidia.com/cusparse" target="_blank">cuSPARSE</a>

.. |cuBLAS| raw:: html

  <a href="https://developer.nvidia.com/cusparse" target="_blank">cuBLAS</a>

Due to the interoperability of all these possibilities with CUDA, the CUDA approach provides
higher functionality and higher performance for GPU devices than an approach based on OpenCL does.

CUDA Version
------------

The CUDA version must be at least 4.0 or higher. 

It is recommended to use always the latest CUDA version:

* Some CUDA versions do not work correctly with the Intel C++ compiler

* The C++11 features can only be used with CUDA version 7.0 or higher. So even if the native
  C++ compiler supports it you have to disable it for the whole configuration.

.. code-block:: bash

  -DCXX_SUPPORTS_C11=False

* Latest NVIDIA GPUs can still be used with older CUDA versions (upward compability)
  but might not give the best possible performance.

CUDA Driver and CUDA Runtime
----------------------------

CUDA is composed of two APIs:

* A low-level API called the CUDA driver API,
* A higher-level API called the CUDA runtime API that is implemented on top of the CUDA driver API.

The CUDA runtime is part the CUDA toolkit and used for the compilation of CUDA code.

Note: Due to the multi-threaded execution model supported by LAMA, the CUDA driver API is used
as only this API provides sufficient control over one CUDA context used by multiple threads.

Compute Capability
------------------

The compute capability of the GPU device must be at least 1.3 or higher. 

CMake Configuration
-------------------

The configuration of LAMA usually finds an available CUDA installation on your system.
If not, you can give it a hint where to find it

.. code-block:: bash

   cmake -D CUDA_TOOLKIT_ROOT=<path/to/cuda/installation>

If CUDA is available on your system but you do not want to use it, you can switch off its use as follows

.. code-block:: bash

   cmake -D USE_CUDA=OFF

Furthermore, you can change relevant CMake variables for CUDA by using the ccmake utility.

Via ccmake you can verify that the Boost variables needed for LAMA have correct values. 
Here are the important ones:

.. code-block:: bash

   CUDA_COMPUTE_CAPABILITY      30
   CUDA_CUDART_LIBRARY          /home/lama/cuda/lib64/libcudart.so
   CUDA_CUDA_LIBRARY            /usr/lib64/libcuda.so   
   CUDA_NVCC_EXECUTABLE         /home/lama/cuda/bin/nvcc
   CUDA_TOOLKIT_ROOT_DIR        /home/lama/cuda
   CUDA_VERSION                 6.5
   CUDA_cublas_LIBRARY          /home/lama/cuda/lib64/libcublas.so 
   CUDA_cusparse_LIBRARY        /home/lama/cuda/lib64/libcusparse.so

The CMake configuration for LAMA tries to identify the compute capability of the GPU
available as device 0 on the machine where LAMA is built. It uses this information 
to optimize the CUDA code for this device. The value of the variable
``CUDA_COMPUTE_CAPABILITY`` must be edited in the following situations:

* If you want to use LAMA on several devices, please use the lowest compute capability for
  configuration. The code will run also on devices with higher compute capabilities.
* It can be NOTFOUND if no device has been found during configuration. Please set it to the capability
  similiar as above.
* If the CUDA toolkit version is too old to support the compute capability (e.g. compute capability 5.3
  is supported only since CUDA version 7.0), use the latest one
  that is supported by that version.

Note: The current LAMA release 2.0 does not support the optimization for different compute capabilities.

Known Problems
^^^^^^^^^^^^^^

On some Linux systems, compilation fails when using Thrust as it might be confused about system files

.. code-block:: bash

   error: kernel launches from templates are not allowed in system files

The following workaround helps to get rid of this compilation problem:

.. code-block:: bash

   unset CPLUS_INCLUDE_PATH

