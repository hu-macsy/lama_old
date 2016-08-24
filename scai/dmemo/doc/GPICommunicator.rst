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

For portability, we recommend using in your applications the default communicator:

.. code-block:: c++

   CommunicatorPtr comm = Communicator::getCommunicatorPtr();

The actual communicator can be specified at runtime.

.. code-block:: c++

    gaspi_run -m machines $PWD/lamaApplication.exe --SCAI_COMMUNICATOR=GPI


Some pitfalls regarding gaspi_run:

* The full path must be specified for the executable
* In a script, you should add a ``sleep 1`` command between two immediate gaspi_run commands

