.. _CUDAAccess:

CUDAAccess
==========

Before calling most of the CUDA driver or CUDA runtime routines, a CUDA context must
be enabled.  Especially the launching of CUDA kernels requires an active CUDA context.
As this context is not an explicit argument it is set globally but thread-specific,
i.e. different threads can have different context enabled. 

Applications that use only one CUDA context with the same thread (most likely the
main thread) use always the same context and do not have to disable or enable it
But when switching between different contexts or different threads, the context has to 
be disabled (pop routine) and another one has to be enabled (push routine).

When a CUDA context has been created via the constructor of a CUDACtx it has been disabled and must
be enabled explicitly via a push routine of the CUDA driver API before further CUDA calls can be executed.
It should be disabled after these other CUDA calls.

.. code-block:: c++

    CUDACtx context( nr );
    ...
    SCAI_CUDA_DRV_CALL( cuCtxPushCurrent( context.getCUcontext(), "" )
    // now cuMemAlloc can be called to reserve memory
    SCAI_CUDA_DRV_CALL( cuMemAlloc( .. ), "" )
    ...
    SCAI_CUDA_DRV_CALL( cuCtxPopCurrent( ... ), "" )

The class CUDAAccess has been introduced to make sure that a context is always enabled
by the constructor and disabled by
the destructor of the corresponding access object.

.. code-block:: c++

    using namespace scai::common;

    CUDACtx context( nr );

    void* pointer;
    size_t size = ...

    // first code block that might also be executed by any thread

    {
        CUDAAccess access( context );
        SCAI_CUDA_DRV_CALL( cuMemAlloc( &pointer, size ), "cuMemAlloc( size = " << size << " ) failed." )
        ...
        // context is disabled at the end of this scope
    }

    // code block that might also be executed by another thread

    {
        CUDAAccess access( context );
        SCAI_CUDA_DRV_CALL( cuMemFree( pointer ), "cuMemFree( " << pointer << " ) failed" )
    }

When a context is accessed via CUDAAccess, it is possible to query for the currently accessed context.
By this way, the latest accessed context is always available globally.

.. code-block:: c++

    void sub()
    {
        // query the currently used CUDA context

        const CUDACtx& context = CUDAAccess::getCurrentCUDACtx();
        std::cout << "actual context is " << context.getDeviceNr() << std::endl;
    }

    ....

    {
        CUDAAccess access( context );
        sub();                        // context is not be passed as argument
    }

A complete CUDA example is this one:

.. literalinclude:: ../examples/cuda/CUDAExample.cu

