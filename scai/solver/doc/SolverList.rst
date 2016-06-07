.. _solver_solver:

***********
Solver List
***********

Direct solvers
--------------

* InverseSolver (NOTE: only single-node implementation; passing solve to LAPACK)

Iterative solver
----------------

Splitting methods
^^^^^^^^^^^^^^^^^

* Jacobi 
* Richardson

Krylow subspace methods
^^^^^^^^^^^^^^^^^^^^^^^

* BiCG
* BiCGstab
* CG
* CGNE
* CGNR
* CGS
* GMRES
* MINRES
* QMR
* TFQMR

Multigrid methods
^^^^^^^^^^^^^^^^^

* SimpleAMG

NOTE: SCAI solver does not support a full AMG. We prepare an interface to |SAMG| - another (commercial) Fraunhofer SCAI library. Please contact us via lama[at]scai.fraunhofer.de if you are interested in using SCAI solver with SAMG.

.. |SAMG| raw:: html

   <a href="https://www.scai.fraunhofer.de/de/geschaeftsfelder/schnelle-loeser/produkte/samg.html" target="_blank"> SAMG </a>
