:orphan:

Solution Task 4
---------------

Here is the solution of task 4. The code demonstrate a CG-Solver running with MPI. 

.. literalinclude:: ../../../../solver/examples/lecture/task4.cpp 
   :language: c++
   :lines: 37-101
   :emphasize-lines: 34,35,36,47

(1) Creating a CommunicationPointer to get access to a parallel environment.
(2) Creating a DistributionPointer of a BlockDistribution.
(3) Redistributing the matrix with row and column distribution. In this case it is a square matrix, so row and column distributions are equal. In a matrix vector multiplication the column distribution must match the distribution of the vector.
(4) Redistributing the rhs Vector.

To execute task4 in parallel use mpirun

.. code-block::bash

   mpirun -np <num-procs> ./task4 <input-file>

:download:`Download complete solution Task 4 <../../../../solver/examples/lecture/task4.cpp>`

**Excursion 2:**

Be careful by printing out values of a parallel run. Use
getValue(IndexType globalindex) to get a correct access with a global index.

.. code-block:: c++

   for ( IndexType i = 0; i < solution.size(); ++i )
   {
       std::cout << solution.getValue( i ) << " ";
   }

   
.. csv-table::
   :header: "back to this Task", "Index", "next Task"
   :widths: 330, 340, 330

   ":doc:`task_4`", ":doc:`../lecture`", ":doc:`task_5`"
