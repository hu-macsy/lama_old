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

The following figure shows how the exchange plan is used to communication values
of an array.

.. figure:: _images/global_exchange.*
    :width: 700px
    :align: center
    :alt: Global exchange.

    Global exchange of data by using a global exchange plan.

The global exchange can also be used for communication in the other direction. 
The permutation is used to scatter the received data in the final array.

.. code-block:: c++

    otherData = ... ;                           // prepare data to send back
    plan.exchangeBack( data, otherData, comm );

