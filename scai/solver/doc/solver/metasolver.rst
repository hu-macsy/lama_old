:orphan:

MetaSolver
----------

The *MetaSolver* is a solver that creates solver instances and connects them, by a user defined configuration file.
Your benefit is, that you can write one program using the metasolver, compile it once and switch between the solver and
solver configurations in every run by just editing the configuration file.

Writing a configuration
^^^^^^^^^^^^^^^^^^^^^^^

To create a configuration you define a text file with the content or hand it over to the configure method as a std::string.

::

   CommonLogger loggerDef( <CGpre>, solverInformation, toConsoleOnly, Timer );

   SOR sorSolver {
      logger = loggerDef; 
      omega = 1.5;
   }
   
   CG root {
      logger = loggerDef;
      preconditioner = sorSolver;
      stoppingCriteria = 
         (IterationCount(10) OR ResidualThreshold(L2Norm, 0.001, Absolute));
   }

Run a MetaSolver
^^^^^^^^^^^^^^^^

::

   void userFunction( 
      const std::string& confgiuration, 
      const Matrix& matrix, 
      const Vector& rhs,
      Vector& solution )
   {
      MetaSolver mt("<MetaSolverName>");
   
      mt.configure( configuration );
      mt.intialize( matrix );
      mt.solve( solution, rhs );
   }
 