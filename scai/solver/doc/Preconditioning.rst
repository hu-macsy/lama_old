.. _solver-preconditioning:

Preconditioning
---------------

SCAI solver prepares a simple preconditioning concept: every iterative solver can be preconditioned by another solver.
You just create a solver and can do everything you know about a solver (e.g. setting a solver logger, initalize with a equation matrix). Afterwards you set your precondition solver with ``setPreconditioner`` for your new solver. The preconditioner will be called before every iteration of your new solver. 

.. code-block:: c++

    SolverPtr prec( new CG( "CG preconditioner" ) );
    GMRES preGMRES( "preconGMRES" );
    preconGMRES.setPreconditioner( pre );
