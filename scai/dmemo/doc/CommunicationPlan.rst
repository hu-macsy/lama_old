.. _CommunicationPlan:

CommunicationPlan
=================

The CommunicationPlan is one central data structure used in the dmemo module
for irregular communication. It specifies for each other processor the number of entries
that will be either sent to or received from it.

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
    :width: 800px
    :align: center
    :alt: ExchangeByPlan 


