.. _GlobalAddressingPlan:

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

    HArray<PartitionId> owners;
    distribution.computeOwners( owners, indexes );
    GlobalExchangePlan plan( owners, distribution->getCommunicatorPtr() );
    HArray<IndexType> myIndexes;
    plan.exchange( myIndexes, indexes, comm );   // I am owner of myIndexes
    sourceDistribution.global2LocalV( myIndexes, myIndexes );

The array myIndexes becomes now part of the global addressing plan.

.. code-block:: c++

    GlobalExchangePlan plan( indexes, distribution );

.. code-block:: c++

    plan.gather( targetArray, sourceArray, BinaryOp::COPY );

    HArray<ValueType> sendValues;  // values to send from my source values
    HArrayUtils::gather( sendValues, sourceArray, plan.myIndexes, BinaryOp::COPY );
    plan.exchangeBack(  targetArray, sendValues );

.. figure:: _images/global_gather.*
    :width: 700px
    :align: center
    :alt: Global gather.

    Gathering from distributed data.

.. code-block:: c++

    plan.scatter( targetArray, sourceArray, BinaryOp::COPY );

    HArray<ValueType> recvValues;  // values to send from my source values
    plan.exchange(  recvValues, sourceArray );
    HArrayUtils::scatter( targetArray, plan.myIndexes, recvValues, BinaryOp::COPY );

.. figure:: _images/global_scatter.*
    :width: 700px
    :align: center
    :alt: Global gather.

    Scattering into distributed data.
