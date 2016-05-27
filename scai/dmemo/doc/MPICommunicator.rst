.. _MPICommunicator:

MPICommunicator
===============

MPICommunicator is a derived Communicator class that implements
the communication routines via MPI.


Communicator Factory
^^^^^^^^^^^^^^^^^^^^

.. code-block:: c++

   CommunicatorPtr mpiComm = Communicator::get( Communicator::MPI );
   CommunicatorPtr gpiComm = Communicator::get( Communicator::GPI );

