CUDAStreamSyncToken
===================

 * Class is a derived class of SyncToken 
 * Encapsulates the reservation of a CUDA stream with the constructor and the
   wait on the stream, relase with the destructor.
 * The scheduling of computations or memory transfers on the CUDA stream is
   done after the constructor call with other method calls.

.. code-block:: c++

    using namespace scai;
    using namespace tasking;

    common::CUDADevice device( deviceNr );

    {
        CUDAStreamSyncToken token( device, CUDAStreamSyncToken::ComputeStream );

        common::CUDAAccess access( device );

        cudaStream_t stream = token.getCUDAStream();
     
        cuda_kernel<<< dimGrid, dimBlock, 0, stream >>>( ... );

        // run other stuff asynchronously with the CUDA kernel

        ...

         // wait is done implicitly at end of scope
     }

Some remarks:

 * The constructor of a CUDAStreamSyncToken pushes and pops the CUcontext, so a new
   CUDAAccess is required afterwards before starting the kernel.
 * Be careful when starting multiple kernels in one CUDA stream. The wait method of the
   CUDAStreamSyncToken waits for the completion of all kernels. If synchronization with
   less kernels is needed, an additional event must be added and set for the token.

.. code-block:: c++

        cudaStream_t stream = token.getCUDAStream();
     
        cuda_kernel<<< dimGrid, dimBlock, 0, stream >>>( ... );

        token.setEvent();

        cuda_kernel1<<< dimGrid, dimBlock, 0, stream >>>( ... );

        token.wait();  // only waits for completion of first kernel
