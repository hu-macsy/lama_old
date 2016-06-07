.. _solver-usage:

Solver Usage
============

SCAI solver provides different linear equation solvers. Using one of the predefined solvers is really easy: you just create an instance of it, initialize it with the equation matrix and can start to solve it. By default the solver just makes one iteration. For running till a predefined criteria set a (combination of) stopping criteria (see :ref:`here <stopping-criteria>`). You can call solve for different right hand sides and initialize the solver with different equation matrices as you want.

The general call to a CG solver is:

.. code-block:: c++

    CG cgSolver( "myCGSolver" );
    cgSolver.intialize( matrix );
    cgSolver.solve( solutionVector, rhsVector );

You may also create the solver instance by calling the factory:

.. code-block:: c++

    SolverPtr solver( Solver::create( "CG", "myCGsolver" ) );
    solver->intialize( matrix );
    solver->solve( solutionVector, rhsVector );
