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

For creating a ``Distribution`` you always also need a ``Communicator``, that later on will take care of the communication between the node. In addition to ``NoCommunicator``, a communicator for **MPI** is also available. You can get one from the factory by calling ``getCommunicatorPtr`` with the choosen ``Communicator``, or without a specification to get the default communicator (defined by environment variable or following a predifined order):

.. code-block:: c++

	dmemo::CommunicatorPtr noComm  = dmemo::Communicator::getCommunicatorPtr( dmemo::CommunicatorType::NO );
	dmemo::CommunicatorPtr mpiComm = dmemo::Communicator::getCommunicatorPtr( dmemo::CommunicatorType::MPI );

	dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr(); // returning the default communicator

Now you can create a ``Distribution`` by calling the constructor or you can get one from the factory by calling ``getDistributionPtr``.

.. code-block:: c++

	IndexType size = 71;
	float weight = 1.0;
	
	// calling the constructor (resulting in a distribution variable)
	dmemo::BlockDistribution block ( size, comm );

	// calling the constructor to initialize a distribution pointer
	dmemo::DistribtionPtr blockPtr ( new dmemo::BlockDistribution( size, comm ) );

	// calling the factory (resulting in distribution pointer)
	dmemo::DistributionPtr cyclicDist ( dmemo::Distribution::getDistribution( "CYCLIC", comm, size, weight ) );


You can set a ``Distribution`` to a ``Vector`` or ``Matrix`` at initialization time by calling the constructor. For the different constructors please refer to :ref:`lama_Vector` or :ref:`lama_Matrix`. To set a distribution after creation, e.g. after initialization, or just to change it later in your application, you can call ``redistribute`` on both with one (``Vector``) or two (``Matrix``) ``Distributions``.

.. code-block:: c++

	DenseVector<float> x( ... );
	CSRSparseMatrix<float> m( ... );

    x.redistribute( blockPtr );
    m.redistribute( cyclicDist, blockPtr );

Note: what you have to keep in mind, distributions in parallel applications are just like matrix and vector sizes in serial applications and have to match for the usage in expressions, e.g. in a matrix-vector-multiplication with **m** and **x** x must have the same ``Distribution`` as m has for its column ``Distribution``. Furthermore the resulting vector has the same ``Distribution`` as m has for its row ``Distribution``.
