.. _MPICommunicator:

MPICommunicator
===============

MPICommunicator is a derived Communicator class that implements
the communication routines via MPI.

.. code-block:: c++

   CommunicatorPtr mpiComm = Communicator::getCommunicatorPtr( Communicator::MPI );

Note: It is not necessary to call ``MPI_Init`` in your application. It will be
called implicitly with the first call to get the MPI communicator from the factory.

.. code-block:: c++

   const IndexType N = 1000;

   // Constructor of distribution uses default communicator; in case of MPI it calls MPI_Init
   DistributonPtr dist( new BlockDistribution>( N ) );
   ...

``MPI_Finalize`` will be called automatically at program exit. 
Sometimes it might be necessary to call the finalize method of the communicator that will
call ``MPI_Finalize`` in case of a MPI communicator. By this way it is possible to run explicitly
some code after MPI has terminated.

.. code-block:: c++

    auto comm = Communicator::getCommunicatorPtr();
    ...
    comm->finalize();  // use of MPI is no more possible after this call
    ...

CUDA-Aware MPI
==============

Running an MPI application on a node with multiple GPUs can be done as follows:

.. code-block:: c++
    mpirun -np 2 dmemo/examples/BenchComm --SCAI_CONTEXT=CUDA --SCAI_DEVICE=0,1

If the MPI implementation (e.g. here mvapich) supports 

.. code-block:: c++
    mpirun -np 2 -env MV2_USE_CUDA 1 -env SCAI_MPI_CUDA dmemo/examples/BenchComm --SCAI_CONTEXT=CUDA --SCAI_DEVICE=0,1

.. code-block:: bash
    mpirun -np 2 -env MV2_USE_CUDA 1 -env SCAI_MPI_CUDA dmemo/examples/BenchComm --SCAI_CONTEXT=CUDA --SCAI_DEVICE=0,1

