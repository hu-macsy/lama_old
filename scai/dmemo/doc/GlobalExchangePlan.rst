.. _GlobalExchangePlan:

GlobalExchangePlan
==================

A global exchange plan is a data structure specifies an irregular communication pattern
for the elements of an array. Typically such a plan is built
by specifying an target array that tells where each element has to
go.

.. code-block:: c++

    HArray<ValueType> data = ...;   
    HArray<PartitionId> owners = ...;   // specifiy for each element of data where it goes
    GlobalExchangePlay plan( owners, comm );
    plan.exchange( otherData, data, comm );

A global exchange plan has exactly three entries:

 - a permuatation that is used to pack the data into contiguous sections for each processor
 - the communication schedule for sending the data
 - the corresponding communication schedule for receiving the data (transpose)

.. figure:: _images/global_exchange_plan.*
    :width: 700px
    :align: center
    :alt: Global exchange.

    Exchanging arbitrary values between processors.

The global exchange can also be used for communication in the other direction. 
The permutation is used to scatter the received data in the final array.

.. code-block:: c++

    otherData = ... ;                           // prepare data to send back
    plan.exchangeBack( data, otherData, comm );

GlobalAddressingPlan
====================

A GlobalAddressingPlan is an extension of a GlobalExchangePlan with an additional
entry that tells how to match the data from other processors with the local data, 
i.e. it tells for each entry how it matches with local entries.

A typically use of a global exchange plan is for gathering data from a distributed array.

.. code-block:: c++

     HArray<double> sourceArray 
     HArray<IndexType> indexes = ...
     HArray<ValueType> targetArray;
 
     // targetArray[i] = sourceArray[ indexes[i] ]  for 0 <= i < indexes.size()
     HArrayUtils::gather( targetArray, sourceArray, indexes, ... ); 

In the following we assume that the above arrays are just the local parts of 
corresponding distributed arrays.
The serial solution does not work for the distributed arrays because the indexes are global 
and might result in out-of-range addressing for the source array that contains only the local part.

.. code-block:: c++

    const Distribution& dist =  ...                      // distribution of source array
    const Communicator& comm =  dist.getCommunicator();

    HArray<PartitionId> owners;
    sourceDistribution.computeOwners( owners, indexes );
    GlobalExchangePlan plan( owners, comm );
    HArray<IndexType> myIndexes;
    plan.exchange( myIndexes, indexes, comm );   // I am owner of myIndexes
    sourceDistribution.global2LocalV( myIndexes, myIndexes );
    HArray<ValueType> sendValues;  // values to send from my source values
    HArrayUtils::gather( sendValues, sourceArray, myIndexes, BinaryOp::COPY );
    plan.exchangeBack(  targetArray, sendValues, distribution.getCommunicator() );

