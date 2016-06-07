.. _NoCommunicator:

NoCommunicator
==============

NoCommunicator is a derived Communicator class that provides
the communication routines on a single machine.


Communicator Factory
^^^^^^^^^^^^^^^^^^^^

.. code-block:: c++

   CommunicatorPtr noComm = Communicator::getCommunicatorPtr( Communicator::NO );

