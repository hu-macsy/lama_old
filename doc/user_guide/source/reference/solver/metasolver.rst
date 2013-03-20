MetaSolver
----------

The MetaSolver is a solver that creates solver instances and connects them, by a user defined 
configuration file.

Writing a configuration
.......................

To create a configuration you define a text file with the content or hand it over to the 
configure method as a std::string.

code block::

   CommonLogger loggerDef( <CGpre>, solverInformation, toConsoleOnly, OpenMPTimer );

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
................

code block::

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
 