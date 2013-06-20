Solver
======

LAMA provides different linear equation solvers. These are available in the ``src/lama/solver/`` directory.

See the example for a short introduction how the solvers can be used. For writing your own solver see
:doc:`here <solver/writingSolver>`.

Available Solvers
-----------------

- Direct Solver

.. - LUSolver
 
 - InverseSolver

- Iterative Solver

  - Jacobi 
 
    - DefaultJacobi
   
    - SpecializedJacobi
   
  - CG
 
  - GMRES
 
  - AMG
 
  - SOR
 
  - SimpleAMG

- Generative Solver

  - :doc:`MetaSolver <solver/metasolver>`

Usage
-----

For using one of the predifined solver you create an instance of it, initialize it with the matrix, set a stopping
criteria and start to solve it.

::

    CG cgSolver( "CGSolver" );
    cgSolver.intialize( matrix );
    // set stopping critia (see below)
    cgSolver.solve( solutionVector, rhsVector );

.. _stopping-criteria:

Stopping Criteria
-----------------

LAMA has three different stopping criteria that can be combined by conjunction (&&) and disjunction (||).
The solver stops if the criteria can be evaluated to true after an iteration.

Possible stopping criteria:

- IterationCount: the iteration count starts with zero with the last solve call and is incremented in every iteration

- ResidualStagnation: the residual stagnation is the absolute difference between the residual of the last two iterations

- ResidualThreshold: the residual threshold is the absolute residual

::

    Scalar eps = 0.00001;
    NormPtr norm = NormPtr( new L2Norm() );
    CriterionPtr rt( new ResidualThreshold( norm, eps, ResidualThreshold::Absolute ) );

    mySolver.setStoppingCriterion( rt );

.. _solver-preconditioning:

Preconditioning
---------------

Every iterative solver can be preconditioned by another solver. Create a solver as preconditioner and set it for your
new solver.

::

    SolverPtr precon( new CG( "CG preconditioner" ) );
    GMRES preconGMRES( "preconGMRES" );
    preconGMRES.setPreconditioner( precon );


.. _solver-logging:

Solver Logging
--------------

The solver logging is an additional logggin to the general LAMA logging. They can be used both together or standalone.
The solver logging gives information about the state of solver, e.g. actual iteration count or residual.

Create a *CommonLogger* and pass it to the solver.

::

   LoggerPtr logger( new CommonLogger( "<LoggerName>: ", LogLevel::solverInformation,
         LoggerWriteBehaviour::toConsoleOnly, std::auto_ptr<Timer>( new Timer() ) ) );
         
   CG cgWithLogger( "CGSolver", logger );

You can choose between five log level:

- noLogging: no solver logging (same as used passing no LoggerPtr)

- convergenceHistory: logs the convergence history of the solver - iterations and residuals

- solverInformation: more information about the solver will be logged

- advancedInformation: advanced solver information like residual requests and stopping criteria checks will be logged

- completeInformation: logs every log message of the solver

There you three write behaviour:

- toConsoleOnly: log messages will be written to standard out only

- toFileOnly: log messages will be written to the log file only

- toFileAndConsole: log messages will be written to the console and the logfile

The Timer will time dedicated calls, that will be logged in with advanced information. 
