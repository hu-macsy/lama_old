.. _CUDAError:

CUDAError
=========

The include file ``CUDAError.hpp`` provides some macros and utilities for calling
functions of the CUDA driver or runtime API, but also for calling cuBLAS or cuSPARSE
routines. These macros will not only print meaningful error messages in case
of failurues but also the file and line number where the call has been made.

.. code-block:: c++

    SCAI_CUDA_DRV_CALL( call, msg ) 
    SCAI_CUDA_RT_CALL( call, msg ) 
    SCAI_CUBLAS_CALL( call, msg ) 
    SCAI_CUSPARSE_CALL( call, msg ) 

The following macro is very useful to verify if a CUDA context has been set
(see also discussions for class CUDAAccess).

.. code-block:: c++

    SCAI_CHECK_CUDA_ACCESS

CUDA Driver API
---------------

Many routines can be called either via the CUDA driver API or via the CUDA runtime API.
Whereever possible, the use of the CUDA driver API is the preferred solution.

CUDA Runtime API
----------------

The runtime API is used mainly for launching CUDA kernels.

cuBLAS API
----------

The NVIDIA CUDA Basic Linear Algebra Subroutines (cuBLAS) library is a GPU-accelerated version of the 
complete standard BLAS library that delivers 6x to 17x faster performance than the latest MKL BLAS.

When using the cuBLAS library, a handle must be created that is returned by an initialization routine
(``cublasInit``). This handle must be passed in each call of a cuBLAS function.

.. code-block:: c++

    CUDACtx context( nr );

    CUDAAccess access( context );

    ...
    float dot;
    cublasHandle_t handle = context.getcuBLASHandle();
    SCAI_CUBLAS_CALL( cublasSdot( getcuBLASHandle(), ..., &dot ), "cublasSDot for float" );


cuSPARSE API
------------

The cuSPARSE library contains a set of basic linear algebra subroutines used for handling sparse matrices. 
It is implemented on top of the NVIDIA CUDA runtime 
and is designed to be called from C and C++. 

Like for the cuBLAS library, a handle is created and returned by an initialization routine (``cusparseCreate``).
This  handle must be passed to every cuSPARSE function.

.. code-block:: c++

    void sub( ... )
    {
        cusparseHandle_t handle = CUDAAccess::getCurrentCUDACtx().getcuSparseHandle();
        SCAI_CUSPARSE_CALL( cusparseScsr2csc( handle, ... ), "csr2csc<float>" )
    }

    CUDACtx context( nr );
    {
        CUDAAccess access( context );
        sub( ... );
    }

For an object of the class ``CUDACtx`` the method ``getcuSparseHandle`` returns the handle to be used
for the corresponding context. The example demonstrates how this handle can be accessed globally in
subroutines.
