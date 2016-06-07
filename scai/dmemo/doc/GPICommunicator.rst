.. _GPICommunicator:

GPICommunicator
===============

GPICommunicator is a derived Communicator class that implements
the communication routines via GPI.


Communicator Factory
^^^^^^^^^^^^^^^^^^^^

.. code-block:: c++

   CommunicatorPtr gpiComm = Communicator::getCommunicatorPtr( Communicator::GPI );

