:orphan:

Solution Task 4
---------------

Here is the solution of task 4. The code is the original CG-Solver from task 0 but now running with MPI. 

.. literalinclude:: ../../../solver/examples/lecture/task4.cpp 
   :language: c++
   :lines: 34-107
   :emphasize-lines: 30-31,37-39,40,44,48

(1) Creating a Communicator (pointer) to get access to a parallel environment.
(2) Creating a Distribution (pointer), here a block distribution.
(3) Redistributing the matrix with row and column distribution. In this case it is a square matrix, so row and column distributions are equal. In a matrix vector multiplication the column distribution must match the distribution of the vector.
(4) The solution vector is created directly with the block distribution.
(5) The rhs vector inherits its distribution from the matrix-vector multiplication, i.e. the row distribution of the matrix.

Remark: The column distribution of a matrix specifies the distribution of the source vector when this matrix is used as a linear
mapping. The row distributon specifies the distribution of the target vector. Internally, only the row distribution is used to map the
matrix coordinates to the processors where the column distribution is used to determine a communication schedule for exchanging
non-local data required for matrix-vector mutliplication.

To execute task4 in parallel use mpirun

.. code-block::bash

   mpirun -np <num-procs> ./task4 <input-file>

:download:`Download complete solution Task 4 <../../../solver/examples/lecture/task4.cpp>`

**Excursion**

When running the executable on multiple processors, all processors write their output to console.
This is not very helpful especially when they are printing all the same message.
Here are some hints how to avoid multiple output.

The most convenient way is to guard output in such a way that it is only printed by
one processor. Usually you take processor 0 so your application remains independent on the
number of processors on which it will run.

.. code-block:: c++

   if ( comm->getRank() == 0 )
   {
       std::cout << "This message is only printed once." << std::endl;
   } 

But you have to be careful. When you run a data parallel operation only with a subset of processors
your application will probably stay in a deadlock. Here is a very typical example that results 
in such a problem.

.. code-block:: c++

   // ATTENTION: ERRONEOUS code
   if ( comm->getRank() == 0 )
   {
       std::cout << "ERROR: never do this: norm of vector = " << l2Norm( v ) << std::endl;
   } 

This is also the reason that you should not set a logger only for a subset of processors.
To avoid output of multiple processors you can set a disabled flag for a logger.

.. code-block:: c++

   bool isDisabled = comm->getRank() > 0;
   LoggerPtr logger( new CommonLogger( "CGLogger: ", LogLevel::convergenceHistory, LoggerWriteBehaviour::toConsoleOnly, isDisabled ) );
   
In some situations, it is indeed helpful that each processor writes its own output. This
is the case if it prints information about its local, processor-specific data.
Good practice is here to print the communicator at first so it is easy to identify which
output comes from which processor.

.. code-block:: c++

   std::cout << *comm << ": my part of vector rhs = " << rhs << std::endl;

:download:`Download complete solution Task 4a <../../../solver/examples/lecture/task4a.cpp>`

.. csv-table::
   :header: "back to this Task", "Index", "next Task"
   :widths: 330, 340, 330

   ":doc:`task_4`", ":doc:`../lecture`", ":doc:`task_5`"
