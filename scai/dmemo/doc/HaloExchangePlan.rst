.. _HaloExchangePlan:

HaloExchangePlan
================

Halo
^^^^

For distributed arrays it is often the case that other processors need
non-local values for their computations. A halo is a data structure that
keeps these non-local values from other processor where each entry in the 
halo stands for a corresponding entry from other processors.

.. figure:: _images/halo.* 
    :width: 600px
    :align: center
    :alt: HaloUpdate

    Halo as array of elements for non-local entries required from other processors.

Halo Exchange Plan
^^^^^^^^^^^^^^^^^^

The halo plan is a data structure that specifies the structure of the halo as well
as how to exchange the corresponding values between the processors.
The halo plan is built by the required indexes for an array with a given distribution.

.. code-block:: c++

     auto plan = haloExchangePlan( distribution, requiredIndexes );

The halo exchange plan contains the following individual information on each processor:

 * An array for the global (non-local) indexes that specify the positions of the halo entries.
   The entries in the halo are sorted by the owners so that communication can be done with
   contiguous sections.
 * A halo communication plan for communication in the halo array
 * A local communication plan for communication in the local parts of the distributed array
 * An array with local indexes for the data that has halo incarnations. It is used to pack
   or unpack data needed for the halo on other processors.

The following figure shows an example of a halo exchange plan for a distributed
array, here 40 elements block distributed onto 4 processors.

.. figure:: _images/halo_plan.* 
    :width: 600px
    :align: center
    :alt: HaloExchangePlan

    Halo exchange plan with structure of the halo and communication plans for data exchange.

The HaloExchangePlan is very similiar to a GlobalAddressingPlan, so here are the differences:

- The halo is organized in contiguous sections for non-local data from each other processor
  so updata of the halo does not require any unpacking of data from other processors.
- Therefore the pack permutation is replaced by the array halo2GlobalIndexes
- It contains a hash map to get for a global index its corresponding halo index.

.. code-block:: c++

    halo2GlobalIndexes[ packPerm ] = requiredIndexes;

Furthermore, a halo is usually built by only non-local indexes. Even if a halo might also
contain copies of local data, this causes overhead and is not the intention of the halo.

Halo Update
^^^^^^^^^^^

Such a halo exchange plan might be used to update the halo with the actual values
from the local parts of other processors, i.e. each halo entry contains the actual value
of the 'global' array.

The corresponding code is as follows:

.. code-block:: c++

     HArray<double> localArray = ...

     HArray<double> haloArray;
     haloExchangePlan.update( haloArray, localArray, comm );

The following figure shows the halo arrays with updated values.

.. figure:: _images/halo_update.* 
    :width: 600px
    :align: center
    :alt: HaloUpdate

    Halo update with actual values from the local parts of a distributed array.

For the halo update each processor has to gather the local data that is needed
on other processors as given by the local indexes array of the halo exchange plan.
Sending is done via the local communication plan, receiving via the global communication 
plan. 

.. code-block:: c++

     HArray<double> sendArray;
     gather( sendArray, localArray, haloPlan.getLocalIndexes() );
     comm.exchangeByPlan( haloArray, haloPlan.getHaloCommunicationPlan(),
                          sendArray, haloPlan.getLocalCommunicationPlan() );

Halo Reduce
^^^^^^^^^^^

A halo plan can also be used to write back halo values to the corresponding positions
of the local counterparts. The communication is exactly in the reverse order: sending
data with the halo communication plan, receiving with the local communication plan
and scattering the received values in the local array. As a local entry might have
a halo counterpart on multiple processors, a reduction operation has to be specified
how to reduce these mutliple values.

.. code-block:: c++

     HArray<double> haloArray = ...
     haloExchangePlan.updateByHalo( localArray, haloArray, common::BinaryOp::ADD, comm );


.. figure:: _images/halo_reduce.* 
    :width: 600px
    :align: center
    :alt: HaloReduce

    Update of local array by halo entries, multiple entries are added.





