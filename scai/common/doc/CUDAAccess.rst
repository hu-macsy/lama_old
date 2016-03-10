CUDAAccess
==========

When a CUDA context has been created with a CUDA device it has been disabled and must
be enabled explicitly via a push routine before further CUDA calls can be executed.

.. code-block:: c++

    CUDADevice device( nr );
    ...
    SCAI_CUDA_DRV_CALL( cuCtxPushCurrent( device.getCUcontext(), "" )
    // now cuMemAlloc can be called to reserve memory
    SCAI_CUDA_DRV_CALL( cuMemAlloc( .. ), "" )
    ...
    SCAI_CUDA_DRV_CALL( cuCtxPopCurrent( ... ), "" )

Most CUDA applications do not use these push and pop context routines but they are mandatory
if another thread will use the same device within one application.

The class CUDAAccess has been introduced to make sure that a context is always disabled by
the destructor of the corresponding access object.

.. code-block:: c++

    using namespace scai::common;

    CUDADevice device( nr );

    void* pointer;
    size_t size = ...

    // first code block that might also be executed by any thread

    {
        CUDAAccess access( device );
        SCAI_CUDA_DRV_CALL( cuMemAlloc( &pointer, size ), "cuMemAlloc( size = " << size << " ) failed." )
        ...
        // context is disabled at the end of this scope
    }

    // code block that might also be executed by another thread

    {
        CUDAAccess access( device );
        SCAI_CUDA_DRV_CALL( cuMemFree( pointer ), "cuMemFree( " << pointer << " ) failed" )
    }

A complete CUDA example is this one:

.. literalinclude:: ../examples/cuda/CUDAExample.cu

