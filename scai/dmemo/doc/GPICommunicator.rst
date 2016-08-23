.. _GPICommunicator:

GPICommunicator
===============

GPICommunicator is a derived Communicator class that implements
the communication routines via GPI.

Communicator Factory
^^^^^^^^^^^^^^^^^^^^

A GPI communicator is available via the factory as follows:

.. code-block:: c++

   CommunicatorPtr gpiComm = Communicator::getCommunicatorPtr( Communicator::GPI );


.. code-block:: c++

   CommunicatorPtr comm = Communicator::getCommunicatorPtr();

