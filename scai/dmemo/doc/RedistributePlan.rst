.. _RedistributePlan:

RedistributePlan
================

A redistribution of data is a very common operation for distributed memory programming,
i.e. a distributed array with a source distribution becomes assigned to a distributed
array with a target distribution. 

.. code-block:: c++

    HArray<ValueType> sourceArray = ...;   // local part of distributed array via sourceDist
    HArray<ValueType> targetArray;         // becomes local part of distributed array via targetDist
 
    globalAssign( targetArray, targetDist, sourceArray, sourceDist );

Redistribution via GlobalAddressingPlan
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Actually it is straightforward to implement this operation by using a global addressing plan.
The first possibility is to gather the new local target data from the 'distributed' source array.

.. code-block:: c++

    // targetArray[i] = sourceArray[targetDist.local2Global(i)]
    auto plan = GlobalAddressingPlan( sourceDist, targetDist.ownedIndexes() );
    plan.gather( targetArray, sourceArray );

But it can be implemented also via a global scatter operation, where each local source data is
scattered into the 'distributed' target array.

.. code-block:: c++

    // targetArray[sourceDist.local2Global(i)] = sourceArray[i]
    auto plan = GlobalAddressingPlan( targetDist, sourceDist.ownedIndexes() );
    plan.scatter( targetArray, sourceArray );

Actually the both plans are just the reverse of each other which is obvious as the
plans represent global permutations that are exactly inverse to each other.

Member Variables of a Redistribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In contrary to a GlobalAddressingPlan we split up the local communication, i.e. the
data that remains local.

Construction of a RedistributePlan
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Usually, the new target distribution is given by
the new owner for each local index that a processor owns in the current source distribution.

.. figure:: _images/redistribution.*
    :width: 700px
    :align: center
    :alt: Redistribution

    Definition of a new distribution by new owners for the local indexes of the current one.

A redistribute plan is now an object that describes the global permutation/communication
pattern that is used for redistribution of distributed data from the source to the target 
distribution.

A redistribute plan is very close to a global exchange plan but has the following 
differences:

 *  it does not contain only a pack permutation but also an unpack permutation to resort
    data incomming from other processors
 *  it is optimized for packing and unpacking local data

.. code-block:: c++

    HArray<PartitionId> newOwners = ...;   // specifiy for each element of data where it goes
    RedistributePlan plan( newOwners, dist );

The array with the new owners is used like in the global exchange plan to set up a global
exchange plan. This global exchange plan is used to send each processor its new global indexes
that it will own for the new distribution. The new global indexes coming from the different
processors form the new target distribution and will be sorted in ascending order. The
corresponding permutation can be used to unpack the redistributed data.

.. figure:: _images/redistribute_plan.*
    :width: 700px
    :align: center
    :alt: RedistributePlan

    Setting up a redistribution plan.

Redistribution
^^^^^^^^^^^^^^

The redistribute plan can be applied to redistribute multiple vectors and/or matrices that have 
the corresponding source distribution and have afterwards the new target distribution.

.. code-block:: c++

    plan.redistribute( targetArray, sourceArray );
    plan.redistributeN( targetArray, sourceArray );
    plan.redistributeV( targetArray, targetSizes, sourceArray, sourceSizes );
