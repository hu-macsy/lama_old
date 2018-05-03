:orphan:

Task 4: Let the CG Iterative Solver run MPI parallel
----------------------------------------------------

In this task we will exploit the data parallelism of LAMA applications. Therefore
the data structures like matrices and vector are distributed among the different
processes/processors and each processor executes the operations on its data.
Therefore the next exercise is to run our CG-Solver MPI parallel. 

LAMA provides the concept of Communicator and
Distribution to handle the complexity of distributed memory parallelism 
Therefore we need to create an object of type Communicator to get access to 
a parallel environment. We can obtain a Communicator from the CommunicatorFactory.
Communicator objects are always managed by shared pointers as this guarantees
that the parallel environment will not be released before all objects using this
communicator have been destructed.

.. code-block:: c++

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

For mapping data to the processors of the parallel environment, we need a Distribution.
In this tutorial we want to work with the class BlockDistribution that maps elements
in equally sized chunks to the processors. There are many more distributions, e.g. CyclicDistribution,
GenBlockDistribution and NoDistribution. Distribution objects are also always managed
by shared pointers, mainly due to the fact that one distribution might be used for multiple
vectors and/or matrices.

.. code-block:: c++

    IndexType size = matrix.getNumRows ( );
    DistributionPtr dist( new BlockDistribution( size, comm ) );

For more explanation we refer to :ref:`scaidmemo:dmemo-distributions`.

An object like matrix or vector can be distributed in one of the following ways:

* Instead of a size argument (size for vectors, number of rows and/or columns for matrices) you take
  the distribution (pointer) as argument.

.. code-block:: c++

    DenseVector<double> v( size, 1.0 );       DenseVector<double> v( dist, 1.0 )
    CSRSparseMatrix<double> m( size, size );  CSRSparseMatrix<double> m( dist, dist );

* If the object has already been constructed, you can redistribute it.

.. code-block:: c++

    v.redistribute( dist );
    m.redistribute( dist, dist );

Operations using vectors and matrices have not to be rewritten at all because
LAMA exploits data parallelism on distributed vectors and distributed matrices by default.

A LAMA program runs parallel with distributed data when 
you start it with mpirun and multiple processes.

.. code-block:: c++

   # disable OpenMP parallelism when using multiple processes on the same node
   export OMP_NUM_THREADS=1          
   mpirun -np 4 task4 gr_30_30.mtx

**Exercise**: Takee on of the applications (our proposed solution takes the original program from task 0)
and rewrite it in such a way that data parallelism is exploited.

.. csv-table:: 
   :header: "previous", "Solution", "next"
   :widths: 330, 340, 330

   ":doc:`task_3`", ":doc:`solution_task_4`", ":doc:`task_5`"

