.. _solver-logging:

Solver Logging
--------------

The solver logging gives information about the state of the solver, e.g. actual iteration count or residual.
By default a solver has no Logger and nothing will be logged. If you want to see additional information about the status of a solver you can create a ``CommonLogger`` and pass it to the solver with the constructor of the solver or afterwards by ``setLogger``.

NOTE: The solver logging is an additional logging to the general :ref:`scailogging:main-page_logging`. They can be used both together or standalone. 

.. code-block:: c++

   LoggerPtr cgLogger( new CommonLogger( "<CGLogger>: ", LogLevel::convergenceHistory, LoggerWriteBehaviour::toConsoleOnly ) );
         
   CG cgWithLogger( "CGSolver", cgLogger ); // setting logger with constructor

   LoggerPtr bicgLogger( new CommonLogger( "<BiCGLogger>: ", LogLevel::completeInformation, LoggerWriteBehaviour::toFileOnly ) );

   BiCG biCG ( "BiCGSolver" );
   biCG.setLogger( logger );

A solver logger has five ``LogLevel`` and three ``LoggerWriteBehaviour`` you can select.

Available log levels:

* noLogging: no solver logging (same as used passing no LoggerPtr)

* convergenceHistory: logs the convergence history of the solver - iterations and residuals

* solverInformation: more information about the solver will be logged

* advancedInformation: advanced solver information like residual requests and stopping criteria checks will be logged

* completeInformation: logs every log message of the solver

Available write behaviours:

* toConsoleOnly: log messages will be written to standard out only

* toFileOnly: log messages will be written to the log file only

* toFileAndConsole: log messages will be written to the console and the logfile

.. The Timer will time dedicated calls, that will be logged in with advanced information. 
