.. _GlobalExchangePlan:

GlobalExchangePlan
==================

A global exchange plan is used for unstructured communication, i.e. for
each element (usually of an array) an individual target processor is
specified.

Typically such a plan is built by specifying an target array that tells 
where each element has to go.

.. code-block:: c++

    HArray<PartitionId> owners = ...;   // specifiy for each element of data where it goes
    auto plan = globalExchangePlay( owners );

A global exchange plain contains these entries:

 - a permuatation that is used to pack the data into contiguous sections for each processor
 - the communication plan for sending the (contiguous) data
 - the corresponding communication schedule for receiving the data (transpose)
 - the communicator that specifies the involved processors

.. figure:: _images/global_exchange_plan.*
    :width: 700px
    :align: center
    :alt: Global exchange.

    Exchanging arbitrary values between processors.

After such a global exchange plan is built, it can be used to exchange data.
The number of entries in the send data must match the size of the owners.
Of course this plan can be used for multiple data exchange.

.. code-block:: c++

    HArray<ValueType> source = ...;   // same size as owners
    HArray<ValueType> target;
    SCAI_ASSERT_EQ_ERROR( source.size(), plan.sendSize() );
    plan.exchange( target, source );
    SCAI_ASSERT_EQ_ERROR( target.size(), plan.recvSize() );

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

    HArray<ValueType> source = ...;   // same size as the receive size of the plan
    HArray<ValueType> target;
    SCAI_ASSERT_EQ_ERROR( source.size(), plan.recvSize() );
    plan.exchangeBack( data, otherData, comm );
    SCAI_ASSERT_EQ_ERROR( target.size(), plan.sendSize() );

The GlobalExchangePlan is used in some way for all kind of unstructured communication,
especially for the other plan classes like GlobalAddressingPlan (gather, scatter
of distributed arrays), RedistributePlan, and HaloExchangePlan. A typical use
of a global exchange is the assembly of matrix data where each processor might 
collect or compute matrix data and only when the final coordinate matrix is built
the coordinate data is sent to the owning processors.

Setting up a an unstructured communication is usually rather expensive. By using
a plan it is possible to reuse the communication pattern and only apply it just
for other data.
