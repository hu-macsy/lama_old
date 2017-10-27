.. _Communicator:

Communicator
============

A communicator provides a set of routines for exchanging non-local data between processors.

Communication with Heterogeneous Arrays
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here is an example of a high-level communication routine that shifts a heterogeneous array
between neighbored processors.

.. code-block:: c++

   hmemo::HArray<ValueType> recvArray, sendArray;

   sendArray = ...
   int direction = 1;  // shift to right, -1 for shift to left
   comm->shiftArray( recvArray, sendArray, direction );

* The method is a template method for different value types of the array.
* The method must be called by all processors (SPMD mode).

For a complete list of available routines we refer to the Doxygen system documentation.

Derived Classes
^^^^^^^^^^^^^^^

* NoCommunicator is a dummy class for a single processor.
* MPICommunicator uses MPI for the implementation of the communication routines.

Communicator Factory
^^^^^^^^^^^^^^^^^^^^

.. code-block:: c++

   CommunicatorPtr mpiComm = Communicator::getCommunicatorPtr( Communicator::MPI );

