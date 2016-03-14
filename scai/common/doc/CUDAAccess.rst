CUDAAccess
==========

Before calling most of the CUDA driver or CUDA runtime routines, a CUDA context must
be set globally. This context is not an explicit argument. Especially the launching
of CUDA kernels requires an active CUDA context.

Applications that use only one CUDA device with the same thread (most likely the
main thread) use always the same context and do not have to disable or enable it
But when switching between different devices or different threads, the context has to 
be disabled (pop routine) and another one has to be enabled (push routine).

When a CUDA context has been created via the constructor of a CUDADevice it has been disabled and must
be enabled explicitly via a push routine before further CUDA calls can be executed.
It should be disabled after these other CUDA calls.

.. code-block:: c++

    CUDADevice device( nr );
    ...
    SCAI_CUDA_DRV_CALL( cuCtxPushCurrent( device.getCUcontext(), "" )
    // now cuMemAlloc can be called to reserve memory
    SCAI_CUDA_DRV_CALL( cuMemAlloc( .. ), "" )
    ...
    SCAI_CUDA_DRV_CALL( cuCtxPopCurrent( ... ), "" )

The class CUDAAccess has been introduced to make sure that a context is always enabled
by the constructor and disabled by
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

When a device is accessed via CUDAAccess, it is possible to query for the currently accessed device.
By this way, the latest accessed device is always available globally.

.. code-block:: c++

    void sub()
    {
        // query the currently accessed device

        const CUDADevice& device = CUDAAccess::getCurrentCUDADevice();
        std::cout << "actual device is " << device.getDeviceNr() << std::endl;
    }

    ....

    {
        CUDAAccess access( device );
        sub();                        // device is not be passed as argument
    }

A complete CUDA example is this one:

.. literalinclude:: ../examples/cuda/CUDAExample.cu

