.. _CommunicationPlan:

CommunicationPlan
=================

The CommunicationPlan is the data structure used in the dmemo module
for communication of contiguous sections between arbitrary processors.
It specifies for each pair (p, q) of processors how many data is sent
from processor p to processor q. It is built by a call on
each processor with the number of quantities this processor wants to
exchange with other processors.

.. code-block:: c++

    std::vector<IndexType> quantities;
    for ( PartitionId p = 0; p < size; ++p )
    {
        quantities.push( ... );  // number of entries to exchange with processor p 
    }
    CommunicationPlan sendPlan( quantities );

Therefore a communication plan is just a dense matrix of size np x np where
np is the number of processors of the corresponding communicator. Each processor
stores only one row of this matrix.

A communication plan might be either used for sending or for receiving data. 
The matching receive plan to a send plan is just the transpose of the corresponding
matrix and can be computed as follows:

.. code-block:: c++

    auto comm = Communicator::getCommunicatorPtr();
    auto sendPlan = comm->transpose( recvPlan );

The following figure shows an example of a send and receive plan for four processors.

.. figure:: _images/communication_plan.* 
    :width: 700px
    :align: center
    :alt: CommunicationPlan

A communicator provides a method to exchange data between processors by using the 
communication plans. The send and receive plan must match in the above sense that one
it the transpose of the other one. 

.. code-block:: c++

    comm.exchangeByPlan( recvArray, recvPlan, sendArray, sendPlan );

.. figure:: _images/exchange_by_plan.* 
    :width: 700px
    :align: center
    :alt: ExchangeByPlan 

The communication plans are used in more complex plans for arbitrary communication,
e.g. GlobalExchangePlan, GlobalAddressingPlan, HaloPlan and RedistributePlan.

