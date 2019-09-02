.. _MPICommunicator:

MPICommunicator
===============

MPICommunicator is a derived Communicator class that implements
the communication routines via MPI.

.. code-block:: c++

   CommunicatorPtr mpiComm = Communicator::getCommunicatorPtr( CommunicatorType::MPI );

Note: It is not necessary to call ``MPI_Init`` in your application. It will be
called implicitly with the first call to get the MPI communicator from the factory.
``MPI_Free`` will be called implicitly as soon as all shared pointers to this object
are destructed.

.. code-block:: c++

   const IndexType N = 1000;

   // Constructor of distribution uses default communicator; in case of MPI it calls MPI_Init
   DistributonPtr dist( new BlockDistribution>( N ) );
   ...
   // calls MPI_Free
   dist = DistributionPtr();   
