Solver
======

LAMA provides different linear equation solvers. These are available in the ``src/lama/solver/`` directory.
See the example for a short introduction how the solvers can be used.
For writing your own solver see :doc:`here <writingSolver>`.

Available Direct Solvers
------------------------

- Direct Solver

 - LUSolver
 
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

 - :doc:`MetaSolver <metasolver>`

Solver Interface
----------------

::

    Solver::intialize( Matrix& matrix );

::

    Solver::solve( Vector& solution, Vector& rhs );

.. _stopping-criteria-label:

Stopping Criteria
-----------------


.. _solver-logging-label:

Solver Logging
--------------


