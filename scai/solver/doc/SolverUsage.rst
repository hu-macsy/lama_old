.. _solver-usage:

Solver
======

LAMA provides different linear equation solvers. These are available in the ``scai/lama/solver/`` directory.

See the example for a short introduction how the solvers can be used. For writing your own solver see
:doc:`here <solver/writingSolver>`.

Usage
-----

For using one of the predefined solver you create an instance of it, initialize it with the matrix, set a stopping
criteria and start to solve it.

.. code-block:: c++

    CG cgSolver( "CGSolver" );
    cgSolver.intialize( matrix );
    // set stopping critia (see below)
    cgSolver.solve( solutionVector, rhsVector );
