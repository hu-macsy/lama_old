.. _solver-preconditioning:

Preconditioning
---------------

Every iterative solver can be preconditioned by another solver. Create a solver as preconditioner and set it for your
new solver.

.. code-block:: c++

    SolverPtr precon( new CG( "CG preconditioner" ) );
    GMRES preconGMRES( "preconGMRES" );
    preconGMRES.setPreconditioner( precon );
    