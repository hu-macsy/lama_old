.. _lama_SetDistribution:

Setting a Distribution
======================

For writing algorithms you can use the generic data structures ``Vector`` and ``Matrix`` without the decision how to divide their data amoung mulitple processes, but when implementing dedicated applications you need to choose a dedicated ``Distribution``. Otherwise the default is a ``NoDistribution``, which results in replicated data on all processes when executing parallely. If you write single-node applications, you do not have to do anything.

You can choose between these types of distributions. For detailed informations about them, please refer to :ref:`scaidmemo:main-page_dmemo`.

* BlockDistribution: builds one continued block for each process with (nearly) equal block size for each process
* CyclicDistribution: builds slices of a specified size (stride), assigning the slices to the processes
* GenBlockDistribution: builds continued blocks for the processes with given block sizes
* GeneralDistribution: assigning each row a specific process
* MetisDistribution: assigment by graph partitioning of a given matrix resulting in a GeneralDistribution

For creating a ``Distribution`` you always also need a ``Communicator``, that later on will take care of the communication between the node. Available communicators beneath the ``NoCommunicator`` are **MPI** and **GPI**. You can get one from the factory by calling ``getCommunicatorPtr`` with the choosen ``Communicator``, or without a specification to get the default communicator (defined by environment variable or following a predifined order):

.. code-block:: c++

	CommunicatorPtr noComm  = Communicator::getCommunicatorPtr( Communicator::NO );
	CommunicatorPtr mpiComm = Communicator::getCommunicatorPtr( Communicator::MPI );
	CommunicatorPtr gpiComm = Communicator::getCommunicatorPtr( Communicator::GPI );

	CommunicatorPtr comm = Communicator::getCommunicatorPtr(); // returning the default communicator

Now you can create a ``Distribution`` by calling the constructor or you can get one from the factory by calling ``getDistributionPtr``.

.. code-block:: c++

	IndexType size = 71;
	float weight = 1.0;
	DistributionPtr dist ( Distribution::getDistribution( "CYCLIC", comm, size, weight ) );

You only can get an instance of a context by the factory by calling ``getContextPtr`` with a context name:



.. code-block:: c++

    z.redistribute( dist );
