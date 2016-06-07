.. _solver-logging:

Solver Logging
--------------

The solver logging is an additional logging to the general LAMA logging. They can be used both together or standalone.
The solver logging gives information about the state of the solver, e.g. actual iteration count or residual.

Create a *CommonLogger* and pass it to the solver.

.. code-block:: c++

   LoggerPtr logger( new CommonLogger( "<LoggerName>: ", LogLevel::solverInformation,
         LoggerWriteBehaviour::toConsoleOnly, std::auto_ptr<Timer>( new Timer() ) ) );
         
   CG cgWithLogger( "CGSolver", logger );

You can choose between five log levels:

- noLogging: no solver logging (same as used passing no LoggerPtr)

- convergenceHistory: logs the convergence history of the solver - iterations and residuals

- solverInformation: more information about the solver will be logged

- advancedInformation: advanced solver information like residual requests and stopping criteria checks will be logged

- completeInformation: logs every log message of the solver

Additonally you have to select one of the three logger behaviours:

- toConsoleOnly: log messages will be written to standard out only

- toFileOnly: log messages will be written to the log file only

- toFileAndConsole: log messages will be written to the console and the logfile

The Timer will time dedicated calls, that will be logged in with advanced information. 