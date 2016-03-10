CUDAStreamPool
===============

 * Class to define a pool of CUDA streams
 * Each CUDA device has its own pool of streams
 * An object of type CUstream can be used to run CUDA kernels or memory transfers
   asynchronously

.. code-block:: c++

    common::CUDADevice device( deviceNr );

    tasking::CUDAStreamPool& pool = tasking::CUDAStreamPool::getPool( device );

    CUstream str1 = pool.reserveStream( CUDAStreamSyncToken::ComputeStream );
    CUstream str2 = pool.reserveStream( CUDAStreamSyncToken::TransferStream );
    CUstream str3 = pool.reserveStream( CUDAStreamSyncToken::TransferStream );
    ...
    ...
    pool.releaseStream( str1 );
    pool.releaseStream( str2 );
    pool.releaseStream( str3 );

    // free of pool can be done explicitly as done here

    tasking::CUDAStreamPool::freePool( device );

Some remarks:

 * The kind of stream (ComputeStream or TransferStream) should be specified as 
   compute kernels on CUDA and memory transfers can overlap with computations on Host CPU
   at the same time.
 * The reserved streams must be released after its use 
 * The pool of streams will be freed implicitly with the destructor of a CUDA device
   (otherwise a pool of streams might be reused incorrectly with a new constructed
   CUDA device).
