.. _GlobalAddressingPlan:

GlobalAddressingPlan
====================

A GlobalAddressingPlan is built when a distributed array is globally accessed,
i.e. each processor has a set of global indexes for which it wants to read
values (gather) or to write values (scatter) of the distributed array.

In a first step a GlobalExchangePlan is built by the owners of the global
indexes. Via a global exchange each processor sends the owning processors
the global indexes for which it wants to read or write elements. Each 
processor translates the received global indexes to local indexes.
The global exchange plan together with the local indexes builds the
global addressing plan.

.. code-block:: c++

    GlobalAddressingPlan globalAddressingPlan( 
      const Distribution& dist, 
      const HArray<IndexType>& globalIndexes )
    {
       auto owners = dist.owners( globalIndexes );
       auto plan = globalExchangePlan plan( owners, dist>getCommunicatorPtr() );
       HArray<IndexType> myIndexes;
       plan.exchange( myIndexes, globalIndexes );   // I am owner of myIndexes
       dist.global2LocalV( myIndexes, myIndexes );
       return GlobalAddressingPlan( plan, myIndexes );
    }

The following figure shows a example how the global addressing plan is 
built by a distribution and a set of global indexes for each processor.

.. figure:: _images/global_addressing_plan.*
    :width: 700px
    :align: center
    :alt: Global gather.

    GlobalAddressingPlan = GlobalExchangePlan + localIndexes

Therefore the class GlobalAddressingPlan is derived from GlobalExchangePlan and 
contains an additional member variable for the localized global indexes.

Gather
^^^^^^

A gather operation reads values from an array at arbitrary positions. For an 
array of indexes it returns an array of the values at the corresponding index
positions. The module ``utilskernel`` provides a corresponding function for
heterogeneous arrays.

.. code-block:: c++

   // targetArray[i] = sourceArray[indexes[i]], 0 <= i < indexes.size()
   HArrayUtils::gather( targetArray, sourceArray, indexes, BinaryOp::COPY );

If now the source array is a distributed array (dist), the gather method can be implemented
as follows:

.. code-block:: c++

    auto plan = globalAddressingPlan( dist, indexes );
    plan.gather( targetArray, sourceArray, BinaryOp::COPY );

The gather method of a global addressing plan gathers the values from the local part of the
array via the local indexes. This array is exchanged via the communication plans. The
received values are now scattered into the target array via the same permutation that
sorted the global indexes for the the owning processors.

.. code-block:: c++
  
    // implementation of GlobalExchangePlan::gather( targetArray, sourceArray, op )
    HArray<ValueType> sendArray, recvArray;   // temporay data
    sendArray = sourceArray[ plan.myIndexes ];
    plan.comm->exchangeByPlan( recvArray, plan.sendPlan, sendArray, plan.recvPlan );
    targetArray[ plan.packPerm ] = recvValues;

The following example shows how

.. figure:: _images/global_gather.*
    :width: 700px
    :align: center
    :alt: Global gather.

    Gathering from distributed array.

Scatter
^^^^^^^

Scatter is writing values at certain positions in an array.

.. code-block:: c++

   // targetArray[indexes[i]] = sourceArray[i]
   HArrayUtils::scatter( targetArray, indexes, unique, sourceArray, BinaryOp::COPY );

If indexes are global indexes for a distributed target array, global scatter is implemented
via a global addressing plan as follows:

.. code-block:: c++

    auto plan = globalAddressingPlan( dist, indexes, unique );
    plan.scatter( targetArray, sourceArray, BinaryOp::COPY );

.. code-block:: c++

    // implementation of GlobalExchangePlan::scatter( targetArray, sourceArray, op )
    HArray<ValueType> sendArray, recvArray;   // temporay data
    sendArray = sourceArray[ plan.packPerm ];
    plan.comm->exchangeByPlan( recvArray, plan.recvPlan, sendArray, plan.sendPlan );
    targetArray[plan.localIndexes] = recvValues;

Here is an example of scattering values into a distributed array. The send permutation
is used to build the send array so communication data is contiguous. The received 
data now is scattered into the local part of the distributed array via the local indexes.

.. figure:: _images/global_scatter.*
    :width: 700px
    :align: center
    :alt: Global gather.

    Scattering into distributed array.

The ``unique`` flag is helpful for scatter operations. It indicates that no (global)
index appears twice and therefore no array element is updated more than once. Only
if elements might be updated more than once, additional (atomic) synchronization is required
for parallel updates. For global scattering on distributed data, the flag can be set when
building the plan.

Redistribution
^^^^^^^^^^^^^^

A redistribution of data is a very common operation for distributed memory programming,
i.e. a distributed array with a source distribution becomes assigned to a distributed
array with a target distribution.

.. code-block:: c++

    HArray<ValueType> sourceArray = ...;   // local part of distributed array via sourceDist
    HArray<ValueType> targetArray;         // becomes local part of distributed array via targetDist
 
    globalAssign( targetArray, targetDist, sourceArray, sourceDist );


Actually it is straightforward to implement this operation by using a global addressing plan.
The first possibility is to gather the new local target data from the 'distributed' source array.

.. code-block:: c++

    // targetArray[i] = sourceArray[targetDist.local2Global(i)]
    auto plan = globalAddressingPlan( sourceDist, targetDist.ownedGlobalIndexes() );
    plan.gather( targetArray, sourceArray );

But it can be implemented also via a global scatter operation, where each local source data is
scattered into the 'distributed' target array.

.. code-block:: c++

    // targetArray[sourceDist.local2Global(i)] = sourceArray[i]
    auto plan = GlobalAddressingPlan( targetDist, sourceDist.ownedGlobalIndexes() );
    plan.scatter( targetArray, sourceArray );

Actually the both plans are just the reverse of each other which is obvious as the
plans represent global permutations that are exactly inverse to each other.
The following figure shows an example of a GlobalAddressingPlan used for a permutation.


.. figure:: _images/global_redistribute.*
    :width: 700px
    :align: center
    :alt: Global redistribute.

    Redistribution with a Global Addressing Plan

Nevertheless, there will be an own plan data structure for redistribution that is optimized 
for local data movement, i.e. that avoids self communication.
