Solver Example
==============

The following C++ program shows how to use a linear equation solver like the Conjugate Gradient 
Method.

::
   
   typedef double ValueType;

   DenseVector<ValueType> solution( 10000, 2.0 );
   const DenseVector<ValueType> rhs( solution.getDistributionPtr(), 1.0 );

   CSRSparseMatrix<ValueType> matrix;
   MatrixCreator<ValueType>::buildPoisson2D( matrix, 9, 100, 100 );

   LoggerPtr logger(
      new CommonLogger(
         "<CG>: ",
         LogLevel::solverInformation,
         LoggerWriteBehaviour::toConsoleOnly,
         std::auto_ptr<Timer>( new Timer() ) ) );

   CriterionPtr criterion( new IterationCount( 10 ) );

   CG solver( "mySolverName", logger );

   solver.setStoppingCriterion( criterion );
   solver.initialize( matrix );
   solver.solve( solution, rhs );  

:download:`Download source file <../../../lama/examples/tutorial/solver.cpp>`

First we need the required vectors and matrices for the linear equation :math:`Ax=b` where A is the 
matrix, x is the solution vector and b is the right-hand-side (rhs).

Then we define the :ref:`solver-logging`. This is optional. By default there is no logging. The example 
creates a logger that shows all solver information like current iteration, residual, etc.

Every iterative solver needs a stopping criterion to define how exact the solution should be or how 
many iterations the solver should run through. These criteria are connectable by using logical 
connectives. For more information see :ref:`stopping-criteria`. 

Afterwards we are ready to create a solver, by defining an identifier and handing the logger over to 
the CGs constructor. Now we can initialize the solver with the matrix. This is a required step before a 
solve method call is possible. The solve method requires the two vectors solution and rhs of the  
equation system. Internally the CG calls its iterate method until the defined stopping criterion 
is satisfied. 
