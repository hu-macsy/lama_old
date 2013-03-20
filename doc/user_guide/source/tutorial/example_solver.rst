.. _example:

Solver Example
==============

The following C++ program shows how to use a linear equation solver like the Conjugate Gradient 
Method.

::
   
   typedef double ValueType;
   
   DenseVector<ValueType> solution( ... );
   const DenseVector<ValueType> rhs( solution );
   CSRSparseMatrix<ValueType> matrix( ... );
   
   LoggerPtr logger( 
      new CommonLogger( 
         "<CG>: ",        
         LogLevel::solverInformation,
         LoggerWriteBehaviour::toConsoleOnly,
         std::auto_ptr<Timer>( new OpenMPTimer() ) ) );
   
   BooleanConditionPtr criterion( new IterationCount( 10 ) );
   
   CG solver( "mySolverName", logger );
   
   solver.setStoppingCriterion( criteria );
   solver.intitialize( matrix ); 
   solver.solve( solution, rhs );   


First we need the required vectors and matrices for the lineare equation :math:`Ax=b` where A is the 
matrix, x is the solution vector and b is the right-hand-side (rhs).

Then we define the :ref:`solver-logging`. This is optional. By default there is no logging. The example 
creates a logger that shows all solver information like current iteration, residual, etc.

Every iterative solver needs a stopping criterion to define how exact the solution should be or how 
many iterations the solver should run through. These criterions are connectable by using logical 
connectives. For more information see :ref:`stopping-criteria`. 

Afterwards we are ready to create a solver, by defining a identifier and handing the logger over to 
the CGs constructor. Now we can intialize the solver with the matrix. This a required step before a 
solve method call is possible. The solve method requires the two vectors solution and rhs of the  
equation system. Internally the CG calls its iterate method until the defined stopping criterion 
is satisfied. 