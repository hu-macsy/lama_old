.. _RedistributePlan:

RedistributePlan
================

A redistribute plan is very close to a global exchange plan but has the following 
differnces:

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

The redistribute plan can be applied to redistribute vectors or matrices that have 
the corresponding source distribution and have afterwards the new target distribution.
