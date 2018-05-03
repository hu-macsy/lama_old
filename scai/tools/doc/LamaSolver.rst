.. _LamaSolver:

**********
lamaSolver
**********

This executable can be used to try the different solvers for given matrices.


.. code-block:: c++

    Usage: <solver> <matrix_filename> [ rhs ] [start_solution] [ final_solution_filename ] [ options ] 

The command line arguments rhs, start_solution can be either a value or a filename.

.. code-block:: c++

         --SCAI_DISTRBUTION=BLOCK|<dist_filename>
         --SCAI_SOLVER=BiCG|BiCGstab|CG|CGNE|CGNR|CGS|DecompositionSolver|GMRES|InverseSolver|Jacobi|Kaczmarz|MINRES|QMR|Richardson|SimpleAMG|TFQMR|TrivialPreconditioner
         --SCAI_SOLVER_LOG=noLogging|convergenceHistory|solverInformation|advancedInformation|completeInformation
         --SCAI_MAX_ITER=<int_val>
         --SCAI_NORM=L1|L2|Max
         --SCAI_REL_TOL=<val>
         --SCAI_ABS_TOL=<val>
         --SCAI_DIV_TOL=<val>
         --SCAI_FORMAT=[CSR|ELL|JDS|DIA|COO]
         --SCAI_TYPE=[float|double|LongDouble|ComplexFloat|ComplexDouble|ComplexLongDouble]
         --SCAI_NUM_THREADS=...
         --SCAI_USE_METIS=<flag>
         --SCAI_ASYNCHRONOUS=<flag>
         or general options:
         --SCAI_COMMUNICATOR=[MPI|NO]
         --SCAI_CONTEXT=[Host|CUDA]
         --SCAI_DEVICE=[0|1|...]
         --SCAI_CUDA_USE_TEXTURE=[0|1]
         --SCAI_CUDA_USE_SHARED_MEM=[0|1]
         --SCAI_CUDA_BLOCK_SIZE=[64|128|...]

