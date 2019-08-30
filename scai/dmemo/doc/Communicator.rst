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
* The method must be called by all processors of the communicator (SPMD mode) 

.. figure:: _images/shift.*
    :width: 500px
    :align: center
    :alt: ShiftArray

    Shifting a heterogeneous array within a communicator.

The following example is a more complicated one where arbitrary data might
be exchanged between processors. For each value of an array the new owner is 
specified.

.. code-block:: c++

   HArray<ValueType> sendArray( ... );  // values to be exchanged
   HArray<PartitionId> onwers( ... );   // tell for each send element where to go

   hmemo::HArray<ValueType> recvArray;  // colllects the values from other processors

   dmemo::globalExchange( recvArray, sendArray, owners, comm );

.. figure:: _images/global_exchange.*
    :width: 700px
    :align: center
    :alt: Global exchange.

    Exchanging arbitrary values between processors.

Per default, communication takes only place between host memory, i.e. read access
for the send array and write access for the received array use both the Host context.
If the send array has only a valid copy in CUDA memory, the data is copied from the
GPU device to the host memory implicitly by the communication routine.
If the communicator supports CUDA memory (e.g. for CUDA-aware MPI) or any other
memory, the read and write access might also use the corresponding data directly.

For a complete list of available routines we refer to the Doxygen system documentation.

Derived Classes
^^^^^^^^^^^^^^^

* NoCommunicator is a dummy class for a single processor.
* MPICommunicator uses MPI for the implementation of the communication routines.

The Communicator provides a factory that returns for each supported communicator type
a corresponding object.

.. code-block:: c++

   CommunicatorPtr mpiComm = Communicator::getCommunicatorPtr( CommunicatorType::MPI );
   CommunicatorPtr noComm = Communicator::getCommunicatorPtr( CommunicatorType::NO );

In many methods the communicator object is an optional argument. If it is omitted,
the current communicator is taken.

.. code-block:: c++

   CommunicatorPtr comm = Communicator::getCommunicatorPtr();

The  communicator is either the default communicator or it might be 
any other communicator that has been set as current communicator.

Communicator Splitting
^^^^^^^^^^^^^^^^^^^^^^

Most applications might work well with the default communicator that uses
all available processors. In some situations it might be useful to work only on
a subset of processors.
Therefore it is possible to split up a communicator into sub-communicators by a
color argument.

.. code-block:: c++

   CommunicatorPtr commWorld = Communicator::getCommunicatorPtr();
   PartitionId color = commWorld->getRank() % 2; 
   CommunicatorPtr commTask = commWorld->split( color );

All processors with the same color will build a new communicator. 

.. figure:: _images/splitting.*
    :width: 500px
    :align: center
    :alt: CommunicatorSplit

    Splitting a communicator into sub-communicators.

After the processors are divided up into subgroups with a new communicator,
the macro ``SCAI_DMEMO_TASK`` should be used to set this communicator as
the current communicator within the actual scope.

.. code-block:: c++

    CommunicatorPtr commWorld = Communicator::getCommunicatorPtr();

    {
        // commTask becomes new current communicator

        SCAI_DMEMO_TASK( commTask )

        // define a block distribution, uses current communicator commTask

        auto dist = std::make_shared<BlockDistribution>( N ); 
        ...
    }

    // this distribution now uses the original communicator commWorld

    auto dist = std::make_shared<BlockDistribution>( N );
