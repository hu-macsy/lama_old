:orphan:

Solution Task 3
===============

.. code-block:: c++
   :emphasize-lines: 46,47,50

   #include <lama.hpp>

   #include <lama/matrix/CSRSparseMatrix.hpp>
   #include <lama/DenseVector.hpp>
   #include <lama/expression/MatrixVectorExpressions.hpp>
   #include <lama/expression/VectorExpressions.hpp>

   #include <lama/solver/CG.hpp>
   #include <lama/solver/criteria/ResidualThreshold.hpp>
   #include <lama/norm/L2Norm.hpp>

   #include <lama/solver/logger/CommonLogger.hpp>
   #include <lama/solver/logger/Timer.hpp>

   #include <iostream>

   using namespace lama;

   int main( int argc, char* argv[] )
   {
       if ( argc < 2 )
       {
           std::cerr << "No input file specified" << std::endl;
           exit( -1 );
       }
   
       //Read a sparse matrix from the passed input file
       CSRSparseMatrix<double> m( argv[1] );
       std::cout << "Read matrix m : " << m << std::endl;
       IndexType size = m.getNumRows();
   
       //Create rhs vector
       DenseVector<double> rhs( size , 0.0 );
       std::cout << "Vector rhs : " << rhs << std::endl;
   
       //Create solution vector
       DenseVector<double> solution( size , 1.0 );
       std::cout << "Vector solution : " << solution << std::endl;
   
       //Compute the rhs that fits our solution to be able to calculate the error later
       rhs = m * solution;
   
       //Reset soution to zero so that there is something to solve
       solution = 0.0;
   
       /* 1 */      std::auto_ptr<Timer> timer( new Timer() );
       /* 2 */      std::auto_ptr<Logger> logger( new CommonLogger("CGLogger", LogLevel::convergenceHistory, LoggerWriteBehaviour::toConsoleOnly, timer) );

       //Create a CG solver using the Constructor with logger as an argument
       /* 3 */      CG cgSolver( "CGTestSolver", logger );
   
       //Create a stopping criterion for the iterative solver cgSolver
       NormPtr norm = NormPtr( new L2Norm() );
       CriterionPtr criterion( new ResidualThreshold( norm, 1E-8, ResidualThreshold::Absolute ) );
       cgSolver.setStoppingCriterion( criterion );
   
       //Initialize the solver
       cgSolver.initialize( m );
       //Solve m * solution = rhs
       cgSolver.solve( solution, rhs );
   
       //calculate the error and its L2-Norm
       DenseVector<double> error( size, 1.0 );
   
       error = error - solution;

       std::cout<<"L2-Norm of error is "<<l2Norm( error )<<std::endl;

       return 0;
   }

(1) Creating a Timer
(2) Creating a CommonLogger
(3) Using Constructor with logger as an argument

The solution of task 3 can be downloaded `here`__.

__ http://libama.sourceforge.net/tutorial/solutions/task3.cpp


.. csv-table::
   :header: "back to this Task", "Index", "next Task"
   :widths: 330, 340, 330

   ":doc:`task_3`", ":doc:`index`", ":doc:`task_4`"

