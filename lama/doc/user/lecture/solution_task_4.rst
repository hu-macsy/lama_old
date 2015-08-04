:orphan:

Solution Task 4
===============

Here is the solution of task 4. The code demonstrate a CG-Solver running with MPI. 

.. code-block:: c++
   :emphasize-lines: 26,27,28,40

   #include <lama.hpp>

   #include <lama/storage/SparseAssemblyStorage.hpp>
   #include <lama/matrix/CSRSparseMatrix.hpp>
   #include <lama/DenseVector.hpp>
   #include <tracing.hpp>

   #include <lama/solver/CG.hpp>
   #include <lama/solver/criteria/IterationCount.hpp>

   #include <lama/CommunicatorFactory.hpp>
   #include <lama/distribution/BlockDistribution.hpp>

   #include <iostream>

   using namespace lama;

   int main( int argc, char* argv[] )
   {
      //TASK 1:
      CSRSparseMatrix<double> m( argv[1] );
      IndexType size = m.getNumRows();
   
      //TASK 4:
      CommunicatorPtr comm( CommunicatorFactory::get( "MPI" ) ); /*(1)*/
      DistributionPtr dist( new BlockDistribution( size, comm ) ); /*(2)*/
      m.redistribute( dist, dist ); /*(3)*/

      //Creation of Vector (task 2):
      DenseVector<double> rhs( size, 0.0 );
      HostWriteAccess<double> hwarhs( rhs.getLocalValues() );

      for ( IndexType i = 0; i < size; ++i )
      {
         hwarhs[i] = double( i + 1 );
      }
      hwarhs.release();

      rhs.redistribute( dist ); /*(4)*/

      DenseVector<double> solution( 0.0, dist );

      //TASK 2:
      //Here is the self-provided implementation of task 2
      {...}

      for ( IndexType i = 0; i < solution.size(); ++i ) 
      {
         std::cout << solution.getValue( i ) << " ";
      }
      return 0;
    }

(1) Creating a CommunicationPointer to get access to a parallel environment.
(2) Creating a DistributionPointer of a BlockDistribution.
(3) Redistributing the matrix with row and column distribution. In this case it is a square matrix, so row and column distributions are equal. In a matrix vector multiplication the column distribution must match the distribution of the vector.
(4) Redistributing the rhs Vector.

To execute task4 in parallel use mpirun

.. code-block::bash

   mpirun -np <num-procs> ./task4 <input-file>

:download:`Download complete solution Task 4 <../../../examples/lecture/task4.cpp>`

**Excursion 2:**

Be careful by printing out values of a parallel run. Use
getValue(IndexType globalindex) to get a correct access with a global index.

::

   for ( IndexType i = 0; i < solution.size(); ++i )
   {
       std::cout << solution.getValue( i ) << " ";
   }

   
.. csv-table::
   :header: "back to this Task", "Index", "next Task"
   :widths: 330, 340, 330

   ":doc:`task_4`", ":doc:`../lecture`", ":doc:`task_5`"
   
