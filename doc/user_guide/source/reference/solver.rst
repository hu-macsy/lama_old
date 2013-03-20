Solver
======

LAMA provides different linear euquation solvers. These are available in the ``src/lama/solver/`` 
directory. See the example for a short introduction how the solvers can be used.s

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

.. toctree::
   :titlesonly:
   :maxdepth: 2

   solver/metasolver

Solver Interface
----------------

::

    Solver::intialize( Matrix& matrix );

::

    Solver::solve( Vector& solution, Vector& rhs );

.. _stopping-criteria:

Stopping Criteria
-----------------


.. _solver-logging:

Solver Logging
--------------


