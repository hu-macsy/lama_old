:orphan:

Solution Task 2
===============

Here is the solution of task 2. The code demonstrate a self-provided CG-Solver. 

.. code-block:: c++
   :emphasize-lines: 37,43,48

   #include <lama.hpp>

   #include <lama/storage/SparseAssemblyStorage.hpp>
   #include <lama/matrix/CSRSparseMatrix.hpp>
   #include <lama/DenseVector.hpp>
   #include <lama/tracing.hpp>

   #include <lama/expression/MatrixVectorExpressions.hpp>
   #include <lama/expression/VectorExpressions.hpp>
   #include <lama/norm/L2Norm.hpp>

   #include <iostream>

   using namespace lama;

   int main( int argc, char* argv[] )
   {
      if ( argc < 2 )
      {
          std::cerr << "No input file specified" << std::endl;
          exit( -1 );
      }

      CSRSparseMatrix<double> m(argv[1]);
      IndexType size=m.getNumRows();

      DenseVector<double> rhs( size, 0.0 );
      HostWriteAccess<double> hwarhs( rhs.getLocalValues() );

      for ( IndexType i = 0; i < size; ++i )
      {
         hwarhs[i] = double( i + 1 );
      }
      hwarhs.release();

      DenseVector<double> solution( size, 0.0 );

      //TASK 2
      /*(1)*/     Scalar rNew;
      Scalar rOld;
      Scalar alpha;
      Scalar beta;
      Scalar eps = 0.00001;

      /*(2)*/     DenseVector<double> r ( size, 0.0 );
      DenseVector<double> help ( size, 0.0 );
      DenseVector<double> d ( size, 0.0 );
      DenseVector<double> Ad( size, 0.0 );

      /*(3)*/     //d = r = solution - m * rhs
      help = m * solution;
      r = rhs - help;
      d = r;

      rOld = r * r;

      L2Norm norm;

      for ( int k = 1; k < 10; k++ )
      {
         Scalar rnorm = norm( r );

         std::cout << "Iter k = " << k << " : norm( r ) = " << rnorm.getValue<double>() << std::endl;

         if ( norm( r ) < eps )
         {
            break;
         }
         else
         {
            Ad = m * d;
            alpha = rOld / ( d * Ad );
   
            solution = solution + ( alpha * d );
            r = r - alpha * Ad;

            rNew = r * r;
            beta = rNew / rOld;

            d = r + beta * d;
            rOld = rNew;
         }
      }
      return 0;
    }

(1) Creation of Scalars and epsilon.
(2) Creation of Vectors, that are needed by your implementation.
(3) Here is the self-provided algorithm of a CG-Solver.

The solution of task 2 can be downloaded `here`__.

__ http://libama.sourceforge.net/tutorial/solutions/task2.cpp

A alternative solution for task 2 is the lama++-CG-class, which contains the same calculation.

::

   Scalar eps = 0.00001;
   L2Norm norm;

   CG cgSolver( "CGTestSolver" );
   const ResidualThreshold rt( norm, eps, ResidualThreshold::Absolute );
   cgSolver.setStoppingCriterion( rt );
   cgSolver.initialize( m );
   cgSolver.solve( solution, rhs );

The solution of task 2 alternative can be downloaded `here`__.

__ http://libama.sourceforge.net/tutorial/solutions/task2a.cpp

**Excursion:**

To see results, it is possible to print out the values of vectors. By the way, an explicit function to print out results does not exist.
But it is possile to use a HostReadAccess to get reading access to the values.

::

    std::cout << "The solution is: ";
    HostReadAccess<double> hra( solution.getLocalValues() );
    for ( int i = 0; i < solution.size(); i++ )
    {
        std::cout << hra[i] << " ";
    }
    std::cout << std::endl;


.. csv-table::
   :header: "back to this Task", "Index", "next Task"
   :widths: 330, 340, 330

   ":doc:`task_2`", ":doc:`index`", ":doc:`task_3`"
