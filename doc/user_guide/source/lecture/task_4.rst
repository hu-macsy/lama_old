:orphan:

Task 4: Let the CG Iterative Solver run MPI parallel
====================================================

LAMA offers different solvers for equation systems, e.g. Jacobi, AMG, or SOR.
Using one of these LAMA-Solvers is an alternative to integrating a
self-provided implementation. To demonstrate the use of a LAMA-Solver, you can
see an example in the appropriated solution of task 2. You can either use the
LAMA provided CG-Solver or the one you have created in Task 2 to solve the
equation system in parallel. LAMA provides the concept of Communicators and
Distributions to handle the complexity of distributed memory parallelism 
through the classes Communicator and Distribution.

Our next exercise is to run a CG-Solver MPI parallel. To run a program parallel,
we need to create an object of type CommunicatorPtr to get access to a parallel
environment. We can obtain a CommunicatorPtr from the CommunicatorFactory.
Besides this object, we need a DistributionPtr. A DistributionPtr is a shared
pointer to a Distribution. You can create one by passing the pointer to a
Distribution to its constructor. In this tutorial we want to work with the class 
BlockDistribution. BlockDistribution is one possible solution in LAMA to handle
distributed data types (there are many more: CyclicDistribution,
GenBlockDistribution and NoDistribution).

You have to redistribute the CSRSparseMatrix and the DenseVectors solution and
rhs from task 1 and 2 by calling their redistribute()-method. After this the
program runs implicitly parallel when you start it with mpirun and multiple
processes.

.. csv-table:: 
   :header: "previous", "Solution", "next"
   :widths: 330, 340, 330

   ":doc:`task_3`", ":doc:`solution_task_4`", ":doc:`task_5`"